//2DD4_check

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

// Local time approximation kernel function (for scalar x)
inline double phi_n(double z_const, double alpha, double x, double h){
    if (h <= 0.0) throw invalid_argument("h must be positive!");
    const double t = pow(h, alpha);
    const double PI = 3.141592653589793238462643383279502884;
    const double norm = 1.0 / sqrt(2.0 * PI * t);
    const double z_const_sq = (x - z_const) * (x - z_const);
    return norm * exp(-z_const_sq / (2.0 * t));
}

// inline double phi_n_2d(double z_const, double alpha, const Vector2d& x_vec, double h){
//     if (h <= 0.0) throw invalid_argument("h must be positive!");
//     const double t = pow(h, alpha);
//     const double PI = 3.141592653589793238462643383279502884;
    
//     Vector2d diff = x_vec - Vector2d::Constant(z_const);
//     double dist_sq = diff.squaredNorm();  // ||x - z||^2
    
//     const double norm = 1.0 / (2.0 * PI * t);
//     return norm * exp(-dist_sq / (2.0 * t));
// }

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
// 模拟格式 (Simulation Schemes from D.1)
// ========================================

// Euler-Maruyama (Order 0.5)
inline State A0(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    // Noise increments
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);

    // Drift term: a(X) = [x2, -x1]^T
    Vector2d drift(curr(1), -curr(0));
    
    // Diffusion term: σ(X) * dW = [s(x2)*dW1, s(x1)*dW2]^T
    Vector2d diffusion(s_func(curr(1)) * dW(0), s_func(curr(0)) * dW(1));

    return curr + drift * dt + diffusion; 
}

// Milstein (Order 1.0)
inline State A1(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    // Noise increments
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    
    double s_x0 = s_func(curr(0));  // s(x1)
    double s_x1 = s_func(curr(1));  // s(x2)
    double ds_x0 = ds_func(curr(0)); // s'(x1)
    double ds_x1 = ds_func(curr(1)); // s'(x2)

    // Euler base
    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_x1 * dW(0), s_x0 * dW(1));
    State base = curr + drift * dt + diffusion;

    // Milstein Correction: 1/2 * s'(x_j)s(x_i) * dW1*dW2
    Vector2d milstein_corr(
        0.5 * ds_x1 * s_x0 * dW(0) * dW(1),  // Row 1
        0.5 * ds_x0 * s_x1 * dW(0) * dW(1)   // Row 2
    );

    return base + milstein_corr;
}

// Strong Order 1.5 Scheme
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    // 1. Basic definitions
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double w1 = dW(0), w2 = dW(1);
    double w1_sq = w1 * w1;
    double w2_sq = w2 * w2;

    double s0 = s_func(curr(0));     // s(x1)
    double s1 = s_func(curr(1));     // s(x2)
    double ds0 = ds_func(curr(0));   // s'(x1)
    double ds1 = ds_func(curr(1));   // s'(x2)
    double dds0 = dds_func(curr(0)); // s''(x1)
    double dds1 = dds_func(curr(1)); // s''(x2)
    
    Vector2d a(curr(1), -curr(0));   // a(x) = [x2, -x1]^T

    // 2. High order noise terms
    double w111 = (w1 * w1_sq) - 3.0 * dt * w1;
    double w222 = (w2 * w2_sq) - 3.0 * dt * w2;
    double w112 = (w1_sq - dt) * w2;
    double w122 = (w2_sq - dt) * w1;

    // 3. Initialize result
    State next = curr;

    // --- X1 Component  ---
    
    // Base Euler
    next(0) += a(0) * dt + s1 * w1;
    
    // Milstein (Order 1.0)
    next(0) += 0.5 * ds1 * s0 * w1 * w2; // = Base Euler + 0.5 * ds1 * s0 * w1 * w2
    
    // Order 1.5 Drift Correction
    next(0) += 0.5 * dt * (s0 * w2 + ds1 * a(1) * w1);
    
    // Order 1.5 Sigma Correction
    double term_t8_1 = 2.0 * dds1 * s0 * s0 * w1;
    double term_t8_2 = (ds1 * s0 * ds1 * s0 / s1) * w1;
    double term_t8_3 = -ds0 * ds1 * s1 * w2;
    next(0) += 0.125 * dt * (term_t8_1 + term_t8_2 + term_t8_3);

    // Order 1.5 Higher Order Noise
    double term_h1 = (1.0/24.0) * (ds1*ds1 / s1) * s0*s0 * w111;
    double term_h2 = (1.0/12.0) * ds0 * ds1 * s1 * w112;
    double term_h3 = (1.0/6.0)  * dds1 * s0*s0 * w122;
    double term_h4 = (1.0/24.0) * (ds1*ds1 / s1) * s0*s0 * w122;
    next(0) += term_h1 + term_h2 + term_h3 + term_h4;

    // --- X2 Component  ---
    
    // Base Euler
    next(1) += a(1) * dt + s0 * w2;

    // Milstein (Order 1.0)
    next(1) += 0.5 * ds0 * s1 * w1 * w2;

    // Order 1.5 Drift Correction
    next(1) += 0.5 * dt * (-s1 * w1 + ds0 * a(0) * w2);

    // Order 1.5 Sigma Correction
    double term2_t8_1 = 2.0 * dds0 * s1 * s1 * w1;
    double term2_t8_2 = ds0 * s1 * ds0 * s1 * w1/ s0;
    double term2_t8_3 = -ds1 * ds0 * s0 * w1;
    next(1) += 0.125 * dt * (term2_t8_1 + term2_t8_2 + term2_t8_3);

    // Order 1.5 Higher Order Noise
    double term2_h1 = (1.0/24.0) * (ds0*ds0 / s0) * s1*s1 * w222;
    double term2_h2 = (1.0/12.0) * ds0 * ds1 * s0 * w122;
    double term2_h3 = (1.0/6.0)  * dds0 * s1*s1 * w112;
    double term2_h4 = (1.0/24.0) * (ds0*ds0 / s0) * s1*s1 * w112;
    next(1) += term2_h1 + term2_h2 + term2_h3 + term2_h4;

    return next;
}

// ========================================
// メインプログラム
// ========================================

int main() {  

    double z_const = 1.38;//z_const = 1.38 \approx 1.4
    constexpr double alpha = 1.0;
    constexpr double t_start = 0.0;
    constexpr double t_end = 1.0;
    constexpr double mu = 0.0;
    constexpr double sigma = 1.0;
    
    // Initial State x_0 = (1, 1) from D.2
    const State x0_state = Vector2d(1.0, 1.0);
    //const State x0_nm = Vector2d(0.0, 0.0);
    constexpr int max_n = 9;
    


    // 配列の初期化
    vector<double> A(max_n + 1, 0.0);
    vector<double> Am(max_n + 1, 0.0);
    vector<double> A_1_5(max_n + 1, 0.0);
    vector<double> Ab(max_n + 1, 0.0);
    vector<double> E(max_n + 1, 0.0);
    vector<double> Em(max_n + 1, 0.0);
    vector<double> E_1_5(max_n + 1, 0.0);
    vector<double> Eb(max_n + 1, 0.0);

    // CSV ファイル名の設定
    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str()); 
    const string csv_path = dir_path + "/2DD4_check_100_1000_data.csv"; 
    ofstream ofs(csv_path, ios::out | ios::trunc);
    
    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }
    
    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1.5,E_b,A,Am,A_1.5,A_b\n";

    // 時間ステップ数のループ
    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n; 
        const int paths = 10 * points * points;
        //const double nm = n;
        const double dt = (t_end - t_start) / (points - 1);
        const double dtm = dt / (points - 1);
        const double sqrt_dt = sqrt(dt);

        // パラメータ初期化
        double S     = 0.0;
        double Sm    = 0.0; 
        double S_1_5 = 0.0;
        double Sb    = 0.0;
        double Bb    = 0.0;
        double B     = 0.0;
        double Bm    = 0.0;
        double B_1_5 = 0.0;

        // OpenMP threadの並列化
        #pragma omp parallel reduction(+:S,Sm,S_1_5,Sb,B,Bm,B_1_5,Bb) 
        {
            // mt19937 rng(42); 
            // mt19937 rng1(30);
            mt19937 rng_nm(30); 
            mt19937 rng1_nm(42);

           
            normal_distribution<double> dist(mu, sigma);
         
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < paths; ++p) {
                // Initialize States
                State st_em  = x0_state;
                State st_mil = x0_state;
                State st_15  = x0_state;
                State st_nm  = x0_state;
                
                // Accumulators for Local Time (L) of the *first component* X1
                double L_em  = 0.0;
                double L_mil = 0.0;
                double L_15  = 0.0;
                double L_nm  = 0.0;
                
                for (int idx = 1; idx < points; ++idx) {
                    // Random numbers
                    // double Z1 = dist(rng);
                    // double Z2 = dist(rng1);
                    double Z1 = 0.0;
                    double Z2 = 0.0;
                    double Z1_nm = dist(rng_nm);
                    double Z2_nm = dist(rng1_nm);

                    //slide benchmark fomula of D4 (Slide Milstein)
                    for (int m = 0; m < points; ++m){
                        State nm_benchmark = A1(st_nm, dtm, Z1_nm, Z2_nm);
                        st_nm = nm_benchmark;
                        Z1 += Z1_nm;
                        Z2 += Z2_nm;
                        
                    }
                    
                    // 1. Update States
                    State next_em  = A0(st_em, dt, Z1, Z2);
                    State next_mil = A1(st_mil, dt, Z1, Z2);
                    State next_15  = A2(st_15, dt, Z1, Z2);

                    // 2. Accumulate Local Time for X1 component (index 0)
                    L_em  += dt * phi_n(z_const, alpha, next_em(0), dt);
                    L_mil += dt * phi_n(z_const, alpha, next_mil(0), dt);
                    L_15  += dt * phi_n(z_const, alpha, next_15(0), dt);
                    L_nm += dt * phi_n(z_const, alpha, st_nm(0), dt);

                    // 3. Move forward
                    st_em  = next_em;
                    st_mil = next_mil;
                    st_15  = next_15;
                }

                // Apply Test Functional f(L) = arctan(L)
                double val_em  = f(L_em);
                double val_mil = f(L_mil);
                double val_15  = f(L_15);
                double val_nm  = f(L_nm);
                

                // Expectation Accumulation
                S += val_em ;
                Sm += val_mil ;
                S_1_5 += val_15 ;
                Sb += val_nm;

                // Variance Accumulation
                B += (val_em ) * (val_em );
                Bm += (val_mil ) * (val_mil );
                B_1_5 += (val_15 ) * (val_15 );
                Bb += (val_nm ) * (val_nm );
            }
        }

        // 期待値の計算
        const double inv_paths = 1.0 / paths;
        A[n] = S * inv_paths;
        Am[n] = Sm * inv_paths;
        A_1_5[n] = S_1_5 * inv_paths;
        Ab[n] = Sb * inv_paths;
        // 分散の計算
        E[n] = B * inv_paths - A[n] * A[n];
        Em[n] = Bm * inv_paths - Am[n] * Am[n];
        E_1_5[n] = B_1_5 * inv_paths - A_1_5[n] * A_1_5[n];
        Eb[n] = Bb * inv_paths - Ab[n] * Ab[n];
        // 出力
        cout << "-------------------------------------------------" << n << "\n";      
        cout << setprecision(10) << "points = " << points << "\n";       
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(15) << "Var EM      = " << E[n] << "\n";         
        cout << setprecision(15) << "Var Milstein= " << Em[n] << "\n";         
        cout << setprecision(15) << "Var 1.5     = " << E_1_5[n] << "\n";  
        cout << setprecision(15) << "Var benchmark    = " << Eb[n] << "\n";    
        cout << setprecision(15) << "Mean EM     = " << A[n] << "\n";          
        cout << setprecision(15) << "Mean Mil    = " << Am[n] << "\n";         
        cout << setprecision(15) << "Mean 1.5    = " << A_1_5[n] << "\n";      
        cout << setprecision(15) << "Mean benchmark    = " << Ab[n] << "\n";      
        // CSVファイルに書き込み
        ofs << n << "," << points << ","  
            << fixed << setprecision(15)
            << E[n] << "," << Em[n] << "," << E_1_5[n] << "," << Eb[n] << "," 
            << A[n] << "," << Am[n] << "," << A_1_5[n] << "," << Ab[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}