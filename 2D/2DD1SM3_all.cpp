// 2DD1SM3_all.cpp
// Error simulation for 2DD2SM3 scheme with Limit Simulation I_T^0
// Option: s(x) = 2 + sin(x) or s(x) = 2 + min(max(-1,x),4)

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

inline double f_sgn(double x, double min_val = -100.0, double max_val = 0.0) {
    return max(min_val, min(x, max_val));
}


// Compute sum involving delta value
inline double compute_sum_state(double delta_val,double delta_val_sq) {
    return delta_val - 0.5 * delta_val_sq;  
}

// s(x) = 2 + sin(x)
inline double s_func(double x) { return 2.0 + sin(x); }

// s'(x) = cos(x)
inline double ds_func(double x) { return cos(x); }

// s''(x) = -sin(x) (Approximation kept from original code)
inline double dds_func(double x) { return -sin(x); }

// // Clipping function
// inline double f_clip(double x) {
//     double min_val = -1.0, max_val = 4.0;
//     return min(max_val, max(x, min_val));
// }

// // s(x) = 2 + min(max(x,-1),4)
// inline double s_func(double x) { return 2.0 + f_clip(x); }
// // s'(x) = 1 if -1 <= x <= 4 else 0
// inline double ds_func(double x) {
//     if (x >= -1.0 && x <= 4.0) {
//         return 1.0;
//     } else {
//         return 0.0;
//     }
// }

// // s''(x) = 0
// inline double dds_func(double x) { return 0.0; }


// ========================================
// Hermite Polynomials Helper Functions
// ========================================

inline double H11(double t, double r1, double S2) {
    return (r1 * r1) / (t * t * S2 * S2) - 1.0 / (t * S2);
}

inline double H22(double t, double r2, double S1) {
    return (r2 * r2) / (t * t * S1 * S1) - 1.0 / (t * S1);
}

inline double H12(double t, double r1, double r2, double S1, double S2) {
    return (r1 * r2) / (t * t * S2 * S1);
}

inline double H112(double t, double r1, double r2, double S1, double S2) {
    double term1 = (r1 * r1 * r2) / (t * t * t * S2 * S2 * S1);
    double term2 = r2 / (t * t * S2 * S1);
    return term1 - term2;
}

inline double H221(double t, double r1, double r2, double S1, double S2) {
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
    return pow(r1, 4) / (t4 * S2_4) - 6.0 * pow(r1, 2) / (t3 * S2_3) + 3.0 / (t2 * S2_2);
}

inline double H2222(double t, double r2, double S1) {
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t2 * t2;
    double S1_2 = S1 * S1;
    double S1_3 = S1_2 * S1;
    double S1_4 = S1_2 * S1_2;
    return pow(r2, 4) / (t4 * S1_4) - 6.0 * pow(r2, 2) / (t3 * S1_3) + 3.0 / (t2 * S1_2);
}

inline double H1112(double t, double r1, double r2, double S1, double S2) {
    double t3 = t * t * t;
    double t4 = t3 * t;
    return (pow(r1, 3) * r2) / (t4 * pow(S2, 3) * S1) - (3.0 * r1 * r2) / (t3 * pow(S2, 2) * S1);
}

inline double H1222(double t, double r1, double r2, double S1, double S2) {
    double t3 = t * t * t;
    double t4 = t3 * t;
    return (r1 * pow(r2, 3)) / (t4 * S2 * pow(S1, 3)) - (3.0 * r1 * r2) / (t3 * S2 * pow(S1, 2));
}

inline double H1122(double t, double r1, double r2, double S1, double S2) {
    return H11(t, r1, S2) * H22(t, r2, S1);
}

// ========================================
// Delta Calculations (Delta^1 & Delta^2)
// ========================================
// Delta^1 Formula
inline double Delta1(const State& x, const State& y, double t) {
    double s_x1 = s_func(x(0));
    double s_x2 = s_func(x(1));
    double ds_x1 = ds_func(x(0));
    double ds_x2 = ds_func(x(1));

    double a1 = x(0);
    double a2 = x(1);

    double r1 = y(0) - x(0) - a1 * t;
    double r2 = y(1) - x(1) - a2 * t;

    double S1 = s_x1 * s_x1;
    double S2 = s_x2 * s_x2;

    double h112 = H112(t, r1, r2, S1, S2);
    double h221 = H221(t, r1, r2, S1, S2);

    double term1 = 0.5 * t * t * ds_x2 * s_x2 * S1 * h112;
    double term2 = 0.5 * t * t * ds_x1 * s_x1 * S2 * h221;

    return term1 + term2;
}

// ========================================
// Simulation Schemes
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

    double term2_t8_1 = 2.0 * dds0 * s1 * s1 * w2;
    double term2_t8_2 = ds0 * s1 * ds0 * s1 * w2 / s0;
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
    
    const State x0_state = Vector2d(1, 1); //X_start
    constexpr int max_n = 3;

    vector<double> A(max_n + 1, 0.0);
    vector<double> Am(max_n + 1, 0.0);
    vector<double> A_1_5(max_n + 1, 0.0);
    vector<double> E(max_n + 1, 0.0);
    vector<double> Em(max_n + 1, 0.0);
    vector<double> E_1_5(max_n + 1, 0.0);
    
    // Stats for Limit Simulation
    vector<double> A_lim(max_n + 1, 0.0);
    vector<double> E_lim(max_n + 1, 0.0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str()); 
    const string csv_path = dir_path + "/2DD1SM3_all_s2sin_100_1000_data.csv"; // or "_min_max_100_1000_data.csv" for min-max option
    ofstream ofs(csv_path, ios::out | ios::trunc);
    
    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }
    
    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1.5,E_lim,A,Am,A_1.5,A_lim\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n; 
        const int paths = 10 * points * points;
        
        const double dt = (t_end - t_start) / (points - 1);
        const double dtm = dt / (points - 1); 
        const double sqrt_dtm = sqrt(dtm);
        const double sqrt_dt = sqrt(dt);

        double S = 0.0, Sm = 0.0, S_1_5 = 0.0,S_lim = 0.0;
        double B = 0.0, Bm = 0.0, B_1_5 = 0.0,B_lim = 0.0;
        
        // #pragma omp parallel reduction(+:S,Sm,S_1_5,S_lim,B,Bm,B_1_5,B_lim) 
        // {
            mt19937 rng_nm(40);
            mt19937 rng1_nm(50);
            
            // Random generator for the independent scalar W_tilde
            mt19937 rng_lim(1000);

            normal_distribution<double> dist(mu, sigma);
         
            // #pragma omp for schedule(static) nowait
            for (int p = 0; p < paths; ++p) {

                State st_em  = x0_state;
                State st_mil = x0_state;
                State st_15  = x0_state;
                State st_nm  = x0_state;

                double D_A0  = 0.0, D_A1  = 0.0, D_A2  = 0.0;
                double D_A0_sq  = 0.0, D_A1_sq  = 0.0, D_A2_sq  = 0.0;
                double D_nm = 0.0, D_nm_sq = 0.0;
                
                // Variable for the Limit Simulation Z
                double I_T = 0.0, I_T_inner = 0.0;
                double sum_A0 = 0.0, sum_A1 = 0.0, sum_A2 = 0.0, sum_nm = 0.0;
                
                for (int idx = 1; idx < points; ++idx) {
                    double Z1 = 0.0;
                    double Z2 = 0.0;
                    
                    for (int m = 0; m < points; ++m){
                        
                        // Generate noise for Benchmark X
                        double Z1_nm = dist(rng_nm);
                        double Z2_nm = dist(rng1_nm);
                        
                        // Generate INDEPENDENT noise for Limit Simulation (scalar W_tilde)
                        double Z_tilde = dist(rng_lim);

                        // Compute Limit Integral Term
                        // Formula: sqrt(3/2) * Integral( sqrt( (s'(x1)/s(x1) * s(x2))^2 + (s'(x2)/s(x2) * s(x1))^2 ) dW_tilde )
                        double x1 = st_nm(0);
                        double x2 = st_nm(1);
                        double s1 = s_func(x1);
                        double s2 = s_func(x2);
                        double ds1 = ds_func(x1);
                        double ds2 = ds_func(x2);
                        
                        // Term 1: (s'(x1)/s(x1)) * s(x2)
                        double term1 = (ds1 / s1) * s2;
                        
                        // Term 2: (s'(x2)/s(x2)) * s(x1)
                        double term2 = (ds2 / s2) * s1;
                        
                        double integrand = sqrt(term1 * term1 + term2 * term2);
                        
                        // Update I^0_T = I^0_T + sqrt(3/2) * integrand * Z_tilde * sqrt(dtm)
                        I_T += sqrt(1.5) * integrand * Z_tilde * sqrt_dtm;
                        I_T_inner += 1.5 * integrand * integrand * dtm;

                        // Update Benchmark State
                        State nm_benchmark = A1(st_nm, dtm, Z1_nm, Z2_nm);
                        
                        //compute Delta_1 of Benchmark
                        double delta_nm = Delta1(st_nm, nm_benchmark, dtm);
                        D_nm += delta_nm;
                        D_nm_sq += delta_nm * delta_nm;
                        st_nm = nm_benchmark;
                        Z1 += Z1_nm / sqrt(points);
                        Z2 += Z2_nm / sqrt(points);
                    }
                    
                    State next_em  = A0(st_em, dt, Z1, Z2);
                    State next_mil = A1(st_mil, dt, Z1, Z2);
                    State next_15  = A2(st_15, dt, Z1, Z2);

                    D_A0  += Delta1(st_em, next_em, dt);
                    D_A1 += Delta1(st_mil, next_mil, dt);
                    D_A2  += Delta1(st_15, next_15, dt);

                    D_A0_sq  += Delta1(st_em, next_em, dt) * Delta1(st_em, next_em, dt);
                    D_A1_sq += Delta1(st_mil, next_mil, dt) * Delta1(st_mil, next_mil, dt);
                    D_A2_sq  += Delta1(st_15, next_15, dt) * Delta1(st_15, next_15, dt);

                    st_em  = next_em;
                    st_mil = next_mil;
                    st_15  = next_15;
                }
                sum_A0 = compute_sum_state(D_A0, D_A0_sq);
                sum_A1 = compute_sum_state(D_A1, D_A1_sq);
                sum_A2 = compute_sum_state(D_A2, D_A2_sq);
                sum_nm = compute_sum_state(D_nm, D_nm_sq);

                // Update States for Limition iteration
                double computed_I_T = exp(-I_T - 0.5 * I_T_inner) - 1.0; 

                //absolute value for the computed_I_T
                double abs_computed_I_T = fabs(computed_I_T);

                S += sgn(sum_nm) - sgn(sum_A0) ;
                Sm += sgn(sum_nm) - sgn(sum_A1) ;
                S_1_5 += sgn(sum_nm) - sgn(sum_A2) ;
                S_lim += abs_computed_I_T;
                B += (sgn(sum_A0) - sgn(sum_nm)) * (sgn(sum_A0) - sgn(sum_nm));
                Bm += (sgn(sum_A1) - sgn(sum_nm)) * (sgn(sum_A1) - sgn(sum_nm));
                B_1_5 += (sgn(sum_A2) - sgn(sum_nm)) * (sgn(sum_A2) - sgn(sum_nm));
                B_lim += abs_computed_I_T * abs_computed_I_T;
                
                
                

            }
       //} // End of parallel region

       // Compute Means and Variances
       //calculate inverse of paths 
        const double inv_paths = 1.0 / paths;

        // Compute Means
        A[n] = S * inv_paths;
        Am[n] = Sm * inv_paths;
        A_1_5[n] = S_1_5 * inv_paths;

        //Compute Variances
        E[n] = B * inv_paths - A[n] * A[n];
        Em[n] = Bm * inv_paths - Am[n] * Am[n];
        E_1_5[n] = B_1_5 * inv_paths - A_1_5[n] * A_1_5[n];
        
        // Compute Limit Stats
        A_lim[n] = S_lim * inv_paths;
        E_lim[n] = B_lim * inv_paths - A_lim[n] * A_lim[n];

        cout << "-------------------------------------------------" << n << "\n";      
        cout << setprecision(10) << "points = " << points << "\n";       
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(15) << "Var EM      = " << E[n] << "\n";         
        cout << setprecision(15) << "Var Milstein = " << Em[n] << "\n";         
        cout << setprecision(15) << "Var 1.5     = " << E_1_5[n] << "\n";      
        cout << setprecision(15) << "Var Lim   = " << E_lim[n] << "\n"; 
        cout << setprecision(15) << "Mean EM     = " << A[n] << "\n";          
        cout << setprecision(15) << "Mean Milstein = " << Am[n] << "\n";         
        cout << setprecision(15) << "Mean 1.5    = " << A_1_5[n] << "\n";      
        cout << setprecision(15) << "Mean Lim   = " << A_lim[n] << "\n"; 
        
        ofs << n << "," << points << ","  
            << fixed << setprecision(15) 
            << E[n] << "," << Em[n] << "," << E_1_5[n] << "," << E_lim[n] << ","
            << A[n] << "," << Am[n] << "," << A_1_5[n] << ","<< A_lim[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}