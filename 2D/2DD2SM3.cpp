// 2DD2SM3.cpp

#include <Eigen/Dense>
#include <algorithm>  
#include <cmath>       
#include <fstream>     
#include <iomanip>    
#include <iostream>   
#include <random>      
#include <string>     
#include <vector>      
#include <locale>
#include <omp.h>

using namespace std;
using namespace Eigen;

// Use Eigen's 2D vector instead of State struct
using State = Vector2d;

// ========================================
// Function definitions
// ========================================

// sgn 関数の定義
constexpr int sgn(double x) {
    if (isnan(x)) return 0;
    return (x > 0) - (x < 0);
}

// Compute sum involving delta value
inline double compute_sum_state(double delta_val,double delta_val_sq) {
    return delta_val - 0.5 * delta_val_sq;  
}

// Local time approximation kernel function (for scalar x)
inline double phi_n(double z_const, double alpha, double x, double h){
    if (h <= 0.0) throw invalid_argument("h must be positive!");
    const double t = pow(h, alpha);
    const double PI = 3.141592653589793238462643383279502884;
    const double norm = 1.0 / sqrt(2.0 * PI * t);
    const double z_const_sq = (x - z_const) * (x - z_const);
    return norm * exp(-z_const_sq / (2.0 * t));
}

// D.2 Test Functional: f(x) = arctan(x)
inline double f(double L) {
    return atan(L);
}

// s(x) = 2 + sin(x)
inline double s_func(double x) { return 2.0 + sin(x); }
// s'(x) = cos(x)
inline double ds_func(double x) { return cos(x); }
// s''(x) = -sin(x)
inline double dds_func(double x) { return -sin(x); }

// ========================================
// Hermite Polynomials Helper Functions
// ========================================

inline double H11(double t, double r1, double S2) {
    // H^{11} = r1^2 / (t^2 * S2^2) - 1 / (t * S2)
    // Note: S2 here represents s(x2)^2
    return (r1 * r1) / (t * t * S2 * S2) - 1.0 / (t * S2);
}

inline double H22(double t, double r2, double S1) {
    // H^{22} = r2^2 / (t^2 * S1^2) - 1 / (t * S1)
    // Note: S1 here represents s(x1)^2
    return (r2 * r2) / (t * t * S1 * S1) - 1.0 / (t * S1);
}

inline double H12(double t, double r1, double r2, double S1, double S2) {
    // H^{12} = (r1 * r2) / (t^2 * S2 * S1)
    return (r1 * r2) / (t * t * S2 * S1);
}

inline double H112(double t, double r1, double r2, double S1, double S2) {
    // H^{112} = r1^2*r2 / (t^3 * S2^2 * S1) - r2 / (t^2 * S2 * S1)
    double term1 = (r1 * r1 * r2) / (t * t * t * S2 * S2 * S1);
    double term2 = r2 / (t * t * S2 * S1);
    return term1 - term2;
}

inline double H221(double t, double r1, double r2, double S1, double S2) {
    // H^{221} = r2^2*r1 / (t^3 * S1^2 * S2) - r1 / (t^2 * S1 * S2)
    double term1 = (r2 * r2 * r1) / (t * t * t * S1 * S1 * S2);
    double term2 = r1 / (t * t * S1 * S2);
    return term1 - term2;
}

inline double H1111(double t, double r1, double S2) {
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t2 * t2;
    double S2_2 = S2 * S2;
    double S2_3 = S2_2 * S2;
    double S2_4 = S2_2 * S2_2;
    
    // H^{1111}
    return pow(r1, 4) / (t4 * S2_4) - 6.0 * pow(r1, 2) / (t3 * S2_3) + 3.0 / (t2 * S2_2);
}

inline double H2222(double t, double r2, double S1) {
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t2 * t2;
    double S1_2 = S1 * S1;
    double S1_3 = S1_2 * S1;
    double S1_4 = S1_2 * S1_2;

    // H^{2222}
    return pow(r2, 4) / (t4 * S1_4) - 6.0 * pow(r2, 2) / (t3 * S1_3) + 3.0 / (t2 * S1_2);
}

inline double H1112(double t, double r1, double r2, double S1, double S2) {
    double t3 = t * t * t;
    double t4 = t3 * t;
    // H^{1112}
    return (pow(r1, 3) * r2) / (t4 * pow(S2, 3) * S1) - (3.0 * r1 * r2) / (t3 * pow(S2, 2) * S1);
}

inline double H1222(double t, double r1, double r2, double S1, double S2) {
    double t3 = t * t * t;
    double t4 = t3 * t;
    // H^{1222}
    return (r1 * pow(r2, 3)) / (t4 * S2 * pow(S1, 3)) - (3.0 * r1 * r2) / (t3 * S2 * pow(S1, 2));
}

inline double H1122(double t, double r1, double r2, double S1, double S2) {
    // H^{1122} = H^{11} * H^{22}
    return H11(t, r1, S2) * H22(t, r2, S1);
}

// ========================================
// Delta Calculations (Delta^1 & Delta^2)
// ========================================
// Delta^1 Formula
inline double Delta1(const State& x, const State& y, double t) {
    // 1. Prepare variables
    double s_x1 = s_func(x(0));
    double s_x2 = s_func(x(1));
    double ds_x1 = ds_func(x(0));
    double ds_x2 = ds_func(x(1));

    // a(x) = [x2, -x1]
    double a1 = x(1);
    double a2 = -x(0);

    // r = y - x - a(x)t
    double r1 = y(0) - x(0) - a1 * t;
    double r2 = y(1) - x(1) - a2 * t;

    // S1 = s(x1)^2, S2 = s(x2)^2
    double S1 = s_x1 * s_x1;
    double S2 = s_x2 * s_x2;

    // 2. Compute Hermite terms
    double h112 = H112(t, r1, r2, S1, S2);
    double h221 = H221(t, r1, r2, S1, S2);

    // 3. Compute Delta^1
    // Term 1: 0.5 * t^2 * s'(x2) * s(x2) * s(x1)^2 * H^{112}
    double term1 = 0.5 * t * t * ds_x2 * s_x2 * S1 * h112;

    // Term 2: 0.5 * t^2 * s'(x1) * s(x1) * s(x2)^2 * H^{221}
    double term2 = 0.5 * t * t * ds_x1 * s_x1 * S2 * h221;

    return term1 + term2;
}

// Delta^2 Formula
inline double Delta2(const State& x, const State& y, double t) {
    // 1. Prepare variables
    double s_x1 = s_func(x(0));
    double s_x2 = s_func(x(1));
    double ds_x1 = ds_func(x(0));
    double ds_x2 = ds_func(x(1));
    double dds_x1 = dds_func(x(0));
    double dds_x2 = dds_func(x(1));

    // a(x) = [x2, -x1]
    double a1 = x(1);
    double a2 = -x(0);
    
    // Derivatives of drift: 
    // ∂1 a1 = 0, ∂2 a1 = 1
    // ∂1 a2 = -1, ∂2 a2 = 0
    double da1_dx1 = 0.0; double da1_dx2 = 1.0;
    double da2_dx1 = -1.0; double da2_dx2 = 0.0;

    // r = y - x - a(x)t
    double r1 = y(0) - x(0) - a1 * t;
    double r2 = y(1) - x(1) - a2 * t;

    // S1 = s(x1)^2, S2 = s(x2)^2
    double S1 = s_x1 * s_x1;
    double S2 = s_x2 * s_x2;
    // S'(x) = 2*s(x)*s'(x)
    double S_prime_x1 = 2.0 * s_x1 * ds_x1;
    double S_prime_x2 = 2.0 * s_x2 * ds_x2;

    // 2. Compute Hermite terms
    double h11 = H11(t, r1, S2);
    double h22 = H22(t, r2, S1);
    double h12 = H12(t, r1, r2, S1, S2);
    double h1111 = H1111(t, r1, S2);
    double h2222 = H2222(t, r2, S1);
    double h1112 = H1112(t, r1, r2, S1, S2);
    double h1222 = H1222(t, r1, r2, S1, S2);
    double h1122 = H1122(t, r1, r2, S1, S2);

    double t2 = t * t;
    double t3 = t * t * t;

    // 3. Compute Delta^2 lines (matching the formula)
    
    // Line 1: 0.5*t^2 [ ∂1 a1 * S2 * H11 + ∂2 a1 * S1 * H12 ]
    double L1 = 0.5 * t2 * (da1_dx1 * S2 * h11 + da1_dx2 * S1 * h12);

    // Line 2: 0.5*t^2 [ ∂1 a2 * S2 * H12 + ∂2 a2 * S1 * H22 ]
    double L2 = 0.5 * t2 * (da2_dx1 * S2 * h12 + da2_dx2 * S1 * h22);

    // Line 3: 0.25*t^2 [ S'(x2)*a2 * H11 + S'(x1)*a1 * H22 ]
    double L3 = 0.25 * t2 * (S_prime_x2 * a2 * h11 + S_prime_x1 * a1 * h22);

    // Line 4: 0.25*t^2 [ s''(x2)s(x2)S1 * H11 + s''(x1)s(x1)S2 * H22 ]
    double L4 = 0.25 * t2 * (dds_x2 * s_x2 * S1 * h11 + dds_x1 * s_x1 * S2 * h22);

    // Line 5: 0.125*t^2 [ (s'(x2)s(x1))^2 * H11 + (s'(x1)s(x2))^2 * H22 ] - 0.25*t^2 [ s'(x1)s(x1)s'(x2)s(x2) * H12 ]
    double L5 = 0.125 * t2 * (pow(ds_x2 * s_x1, 2) * h11 + pow(ds_x1 * s_x2, 2) * h22)
              - 0.25 * t2 * (ds_x1 * s_x1 * ds_x2 * s_x2 * h12);

    // Line 6: t^3/24 [ (s'(x2))^2 S1 S2 * H1111 + (s'(x1))^2 S1 S2 * H2222 ]
    double L6 = (t3 / 24.0) * (pow(ds_x2, 2) * S1 * S2 * h1111 + pow(ds_x1, 2) * S1 * S2 * h2222);

    // Line 7: t^3/12 [ s'(x1)s'(x2) S1 S2^1.5 * H1112 + s'(x1)s'(x2) S1^1.5 S2 * H1222 ]
    // Note: S1*S2^1.5 = s(x1)^2 * s(x2)^3, S1^1.5*S2 = s(x1)^3 * s(x2)^2
    double L7 = (t3 / 12.0) * (ds_x1 * ds_x2 * S1 * pow(s_x2, 3) * h1112 + ds_x1 * ds_x2 * pow(s_x1, 3) * S2 * h1222);

    // Line 8: Big bracket term * H1122
    // Bracket = t^3/6 s''(x1)s(x1)S2^2 + t^3/6 s''(x2)s(x2)S1^2 + t^3/24 (s'(x2))^2 S1^2 + t^3/24 (s'(x1))^2 S2^2
    double term8_1 = (t3 / 6.0) * dds_x1 * s_x1 * pow(S2, 2);
    double term8_2 = (t3 / 6.0) * dds_x2 * s_x2 * pow(S1, 2);
    double term8_3 = (t3 / 24.0) * pow(ds_x2, 2) * pow(S1, 2);
    double term8_4 = (t3 / 24.0) * pow(ds_x1, 2) * pow(S2, 2);
    double L8 = (term8_1 + term8_2 + term8_3 + term8_4) * h1122;

    return L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8;
}

// ========================================
// Simulation Schemes from D.1
// ========================================

// Euler-Maruyama (Order 0.5)
inline State A0(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_func(curr(1)) * dW(0), s_func(curr(0)) * dW(1));
    return curr + drift * dt + diffusion; 
}

// Milstein (Order 1.0)
inline State A1(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    
    double s_x0 = s_func(curr(0));
    double s_x1 = s_func(curr(1));
    double ds_x0 = ds_func(curr(0));
    double ds_x1 = ds_func(curr(1));

    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_x1 * dW(0), s_x0 * dW(1));
    State base = curr + drift * dt + diffusion;

    Vector2d milstein_corr(
        0.5 * ds_x1 * s_x0 * dW(0) * dW(1),
        0.5 * ds_x0 * s_x1 * dW(0) * dW(1)
    );
    return base + milstein_corr;
}

// Strong Order 1.5 Scheme
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double w1 = dW(0), w2 = dW(1);
    double w1_sq = w1 * w1;
    double w2_sq = w2 * w2;

    double s0 = s_func(curr(0));
    double s1 = s_func(curr(1));
    double ds0 = ds_func(curr(0));
    double ds1 = ds_func(curr(1));
    double dds0 = dds_func(curr(0));
    double dds1 = dds_func(curr(1));
    
    Vector2d a(curr(1), -curr(0));

    double w111 = (w1 * w1_sq) - 3.0 * dt * w1;
    double w222 = (w2 * w2_sq) - 3.0 * dt * w2;
    double w112 = (w1_sq - dt) * w2;
    double w122 = (w2_sq - dt) * w1;

    State next = curr;

    // --- X1 Component ---
    next(0) += a(0) * dt + s1 * w1;
    next(0) += 0.5 * ds1 * s0 * w1 * w2;
    next(0) += 0.5 * dt * (s0 * w2 + ds1 * a(1) * w1);
    
    double term_t8_1 = 2.0 * dds1 * s0 * s0 * w1;
    double term_t8_2 = (ds1 * s0 * ds1 * s0 / s1) * w1;
    double term_t8_3 = -ds0 * ds1 * s1 * w2;
    next(0) += 0.125 * dt * (term_t8_1 + term_t8_2 + term_t8_3);

    double term_h1 = (1.0/24.0) * (ds1*ds1 / s1) * s0*s0 * w111;
    double term_h2 = (1.0/12.0) * ds0 * ds1 * s1 * w112;
    double term_h3 = (1.0/6.0)  * dds1 * s0*s0 * w122;
    double term_h4 = (1.0/24.0) * (ds1*ds1 / s1) * s0*s0 * w122;
    next(0) += term_h1 + term_h2 + term_h3 + term_h4;

    // --- X2 Component ---
    next(1) += a(1) * dt + s0 * w2;
    next(1) += 0.5 * ds0 * s1 * w1 * w2;
    next(1) += 0.5 * dt * (-s1 * w1 + ds0 * a(0) * w2);

    double term2_t8_1 = 2.0 * dds0 * s1 * s1 * w1;
    double term2_t8_2 = ds0 * s1 * ds0 * s1 * w1/ s0;
    double term2_t8_3 = -ds1 * ds0 * s0 * w1;
    next(1) += 0.125 * dt * (term2_t8_1 + term2_t8_2 + term2_t8_3);

    double term2_h1 = (1.0/24.0) * (ds0*ds0 / s0) * s1*s1 * w222;
    double term2_h2 = (1.0/12.0) * ds0 * ds1 * s0 * w122;
    double term2_h3 = (1.0/6.0)  * dds0 * s1*s1 * w112;
    double term2_h4 = (1.0/24.0) * (ds0*ds0 / s0) * s1*s1 * w112;
    next(1) += term2_h1 + term2_h2 + term2_h3 + term2_h4;

    return next;
}


// ========================================
// Main Program
// ========================================

int main() {  

    constexpr double t_start = 0.0;
    constexpr double t_end = 1.0;
    constexpr double mu = 0.0;
    constexpr double sigma = 1.0;
    
    const State x0_state = Vector2d(1.0, 1.0);
    constexpr int max_n = 9;

    vector<double> A(max_n + 1, 0.0);
    vector<double> Am(max_n + 1, 0.0);
    vector<double> A_1_5(max_n + 1, 0.0);
    vector<double> E(max_n + 1, 0.0);
    vector<double> Em(max_n + 1, 0.0);
    vector<double> E_1_5(max_n + 1, 0.0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str()); 
    const string csv_path = dir_path + "/2DD2SM3_100_1000_data.csv"; 
    ofstream ofs(csv_path, ios::out | ios::trunc);
    
    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }
    
    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1.5,A,Am,A_1.5\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n; 
        const int paths = 10 * points * points;
        
        const double dt = (t_end - t_start) / (points - 1);
        const double sqrt_dt = sqrt(dt);

        double S = 0.0, Sm = 0.0, S_1_5 = 0.0;
        double B = 0.0, Bm = 0.0, B_1_5 = 0.0;


        #pragma omp parallel reduction(+:S,Sm,S_1_5,B,Bm,B_1_5) 
        {
            mt19937 rng(42); 
            mt19937 rng1(30);
            normal_distribution<double> dist(mu, sigma);
         
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < paths; ++p) {
                State st_em  = x0_state;
                State st_mil = x0_state;
                State st_15  = x0_state;
                
                double D_A0  = 0.0, D_A1  = 0.0, D_A2  = 0.0;
                double sum_A0 = 0.0, sum_A1 = 0.0, sum_A2 = 0.0;
                
                for (int idx = 1; idx < points; ++idx) {
                    double Z1 = dist(rng);
                    double Z2 = dist(rng1);
                    
                    State next_em  = A0(st_em, dt, Z1, Z2);
                    State next_mil = A1(st_mil, dt, Z1, Z2);
                    State next_15  = A2(st_15, dt, Z1, Z2);

                    D_A0  += Delta2(st_em, next_em, dt);
                    D_A1 += Delta2(st_mil, next_mil, dt);
                    D_A2  += Delta2(st_15, next_15, dt);

                    st_em  = next_em;
                    st_mil = next_mil;
                    st_15  = next_15;
                    // cout << "Path " << p << " Step(idx) " << idx << "\n";
                    // cout << "next_em = [" << next_em(0) << ", " << next_em(1) << "]\n";
                    // cout << "next_mil = [" << next_mil(0) << ", " << next_mil(1) << "]\n";
                    // cout << "next_15 = [" << next_15(0) << ", " << next_15(1) << "]\n";
                    // cout << "D_A0 = " << D_A0 << ", D_A1 = " << D_A1 << ", D_A2 = " << D_A2 << "\n";
                }
                
                sum_A0 = D_A0;
                sum_A1 = D_A1;
                sum_A2 = D_A2;

                S += sgn(sum_A0) / sqrt(dt);
                Sm += sgn(sum_A1) / sqrt(dt);
                S_1_5 += sgn(sum_A2) / sqrt(dt);
                B += sgn(sum_A0) * sgn(sum_A0) / dt;
                Bm += sgn(sum_A1) * sgn(sum_A1) /dt;
                B_1_5 += sgn(sum_A2) * sgn(sum_A2) /dt;

                // cout << "Path " << p << " Step(paths) " << paths << "\n";
                // cout << "sum_A0 = " << sum_A0 << ", sum_A1 = " << sum_A1 << ", sum_A2 = " << sum_A2 << "\n";
                // cout << "S = " << S << ", Sm = " << Sm << ", S_1_5 = " << S_1_5 << "\n";
            }
       } // End of parallel region

        const double inv_paths = 1.0 / paths;
        A[n] = S * inv_paths;
        Am[n] = Sm * inv_paths;
        A_1_5[n] = S_1_5 * inv_paths;
        
        E[n] = B * inv_paths - A[n] * A[n];
        Em[n] = Bm * inv_paths - Am[n] * Am[n];
        E_1_5[n] = B_1_5 * inv_paths - A_1_5[n] * A_1_5[n];

        cout << "-------------------------------------------------" << n << "\n";      
        cout << setprecision(10) << "points = " << points << "\n";       
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(10) << "Var EM      = " << E[n] << "\n";         
        cout << setprecision(10) << "Var Milstein = " << Em[n] << "\n";         
        cout << setprecision(10) << "Var 1.5     = " << E_1_5[n] << "\n";      
        cout << setprecision(10) << "Mean EM     = " << A[n] << "\n";          
        cout << setprecision(10) << "Mean Milstein = " << Am[n] << "\n";         
        cout << setprecision(10) << "Mean 1.5    = " << A_1_5[n] << "\n";      

        ofs << n << "," << points << ","  
            << fixed << setprecision(18) 
            << E[n] << "," << Em[n] << "," << E_1_5[n] << "," 
            << A[n] << "," << Am[n] << "," << A_1_5[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}
