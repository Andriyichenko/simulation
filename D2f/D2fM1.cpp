// D2fM1 with x_0=1, T=1, a=0.5, b=0.5

#include <algorithm>  
#include <cmath>       
#include <fstream>     
#include <iomanip>    
#include <iostream>   
#include <random>      
#include <string>     
#include <vector>      
#include <locale>
#include <omp.h>       // 事前にOpenMPをインストールしてください %macOS: brew install libomp

using namespace std;

// ========================================
// 関数の定義
// ========================================

// sgn 関数の定義
constexpr int sgn(double x) {
    if (isnan(x)) return 0;
    return (x > 0) - (x < 0);
}

// delta_1(t,x,y) の定義　順番直し
inline double delta_1(double sigma_prime, double sigma, double dt, 
                      double x, double y, double mu) {
    const double diff = y - x - mu * dt;
    const double sigma_inv = 1.0 / sigma;
    const double diff_sigma_inv = diff * sigma_inv;
    const double diff_sq = diff * diff;
    const double dt_inv = 1.0 / dt;


    const double A = diff_sq * diff * sigma_inv * sigma_inv * sigma_inv * dt_inv;
    const double B = diff_sigma_inv;
    
    return 0.5 * sigma_prime * (A - 3.0 * B);
}

// delta_2(t,x,y) の定義 
inline double delta_2(double a_x, double a_x_prime, double sigma_prime, double sigma, 
                      double sigma_prime_2, double dt, double x, double y) {
    double diff = y - x - a_x * dt;
    double diff_sq = diff * diff;
    double sigma_sq = sigma * sigma;
    double sigma_inv = 1.0 / sigma;
    double sigma_inv_sq = sigma_inv * sigma_inv;
    double sigma_inv_cube = sigma_inv_sq * sigma_inv;
    
    double A = 0.5 * a_x_prime * sigma + 0.5 * a_x * sigma_prime + 
               0.25 * sigma_prime_2 * sigma_sq;
    double B = diff_sq * sigma_inv_cube;
    double C = dt * sigma_inv;
    double D = sigma_prime_2 * sigma + sigma_prime * sigma_prime;
    double E = diff_sq * diff_sq * sigma_inv_sq * sigma_inv_sq / dt;    
    double F = diff_sq * sigma_inv_sq;
    double G = 3.0 * dt;
    return A * (B - C) + D * (E - 6.0 * F + G) / 6.0;
}

// ========================================
// パラメータの設定
// ========================================

struct StateCoeff {
    double sqrt_W_sq_plus_1, sqrt_X_b_sq_plus_1;                              //sqrt(W^2 + 1), sqrt(X_b^2 + 1)
    double drift, drift_deriv, drift_X_b, drift_X_b_deriv;                       // a(W), a'(W), a_m(W), a_m'(W)
    double sigma, sigma_deriv,sigma_deriv2, sigma_inv;                       // sigma(W), sigma'(W), sigma_m(W), sigma_m'(W)
    double sigma_X_b, sigma_X_b_deriv;                                       // sigma(X_b), sigma''(X_b)
    double sigma_sq, sigma_cube,sigma_sq_deriv, sigma_sq_deriv2;        // sigma(W)

// 係数の計算
inline void compute( double a, double b, double W_state) {
        const double w_sq = W_state * W_state;
        const double W_sq_plus_1 = w_sq + 1.0;
        const double b_sq = b * b, b_quad = b * b * b, a_b = a * b;
    
    
        sqrt_W_sq_plus_1 = sqrt(W_sq_plus_1);

        //a(x)の計算部分
        drift = 0.5 * b_sq * W_state + 0.5 * a_b * sqrt_W_sq_plus_1; //a_x
        drift_deriv = 0.5 * b_sq + (0.5 * a_b * W_state / sqrt_W_sq_plus_1); //a_x'
        
        // sigma(x)
        if (fabs(b) < 1e-12) {
            sigma = 0.0;
            sigma_inv = 0.0;
            sigma_sq = 0.0;
            sigma_deriv = 0.0;
            sigma_deriv2 = 0.0;
            sigma_sq_deriv = 0.0;
            sigma_sq_deriv2 = 0.0;
        } else {
            sigma = b * sqrt_W_sq_plus_1; // sigma = b_const *sqrt(W^2 + 1)
            sigma_inv = 1.0 / sigma; //1/sigma
            sigma_sq = sigma * sigma; //b = sigma^2
            sigma_deriv = b * b * W_state * sigma_inv; //sigma'
            sigma_deriv2 = b * b * b * b * sigma_inv * sigma_inv * sigma_inv; //sigma''
            sigma_sq_deriv = 2.0 * b * b * W_state; // b'
            sigma_sq_deriv2 = 2.0 * b * b; // b''
        }


    }
};




// ========================================
// 近似の更新関数の定義
// ========================================

// Euler-Maruyama
inline double A0(double W_state, const StateCoeff& coef,
                                         double dt, double Z) {
    const double sqrt_dt = sqrt(dt);

    return W_state + coef.drift * dt + coef.sigma * sqrt_dt * Z;
}

// Milstein
inline double A1(double W_state, const StateCoeff& coef,
                                   double dt, double Z) {
    const double sqrt_dt = sqrt(dt);
    const double dt_sq = dt * dt;
    const double Z_sq = Z * Z;
    return W_state + coef.drift * dt + coef.sigma * sqrt_dt * Z + 
           0.5 * coef.sigma_deriv * coef.sigma * (Z_sq * dt - dt);
           
}

// 1.5 
inline double A2(double W_state, const StateCoeff& coef,
                            double dt, double Z) {
    const double dt_sq = dt * dt;
    const double Z_sq = Z * Z;
    const double sqrt_dt = sqrt(dt);
    const double dW = Z * sqrt_dt;
    const double dW_cube = dW * dW * dW;
    const double base = W_state + coef.drift * dt + coef.sigma * dW;

    const double milstein_term = 0.5 * coef.sigma_deriv * coef.sigma * (Z_sq * dt - dt);

    const double mixed_deriv = 0.5 * coef.drift_deriv * coef.sigma + 
                               0.5 * coef.sigma_deriv * coef.drift + 
                               0.25 * coef.sigma_deriv2 * coef.sigma_sq;

    const double term3 = mixed_deriv * dW * dt;

    const double triple_deriv = coef.sigma_deriv2 * coef.sigma_sq + 
                                coef.sigma_deriv * coef.sigma_deriv * coef.sigma;

    const double term4 = triple_deriv * (dW_cube - 3 * dW * dt) / 6.0;
    
    return base + milstein_term + term3 + term4;
}

// M1のbenchmark関数の定義
inline double benchmark(double X_b, double dt, double Z, double b, double a) {
    const double asinh_Xb = asinh(X_b);
    const double sqrt_dt = sqrt(dt);
    const double a_b = a * b;
   

    return sinh((asinh_Xb) + 0.5 * a_b * dt + b * sqrt_dt * Z);

}

inline double compute_sum_state(double delta_val) {
    return delta_val - 0.5 * delta_val * delta_val;  
}

inline double f(double x, double min_val = -100.0, double max_val = 0.0) {
    return max(min_val, min(x, max_val)); 
}

inline double c_4_sq(const StateCoeff& coef) { 
    return coef.sigma_sq_deriv2 * coef.sigma_sq * coef.sigma_sq / 12.0; // b'' * b^2 /12.0
}

inline double c_2_sq(const StateCoeff& coef) {
    return coef.drift_deriv * coef.sigma_sq * 0.5 + coef.sigma_sq_deriv * coef.drift * 0.25
           + coef.sigma_sq_deriv2 * coef.sigma_sq * 0.125 - coef.sigma_sq_deriv * coef.sigma_sq_deriv * 0.0625;
           // 0.5 * a' * b + 0.25 * b' * a + 0.125 * b'' * b - 0.0625 * (b')^2
}

// ========================================
// メインプログラム
// ========================================

int main() {  

    // 定数の定義
    constexpr double t_start = 0.0;
    constexpr double t_end = 1.0;
    constexpr double mu = 0.0;
    constexpr double sigma = 1.0;
    constexpr double b = 0.5;
    constexpr double a = 0.5;
    constexpr double b_sq = b * b;          
    constexpr double b_quad = b_sq * b_sq;
    constexpr double x_0 = 1.0;
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
    system(("mkdir -p " + dir_path).c_str()); //フォルダーの確認 
    const string csv_path = dir_path + "/D2fM1_100_1000_data.csv"; //data sourceのファイル名指定
    ofstream ofs(csv_path, ios::out | ios::trunc);
    
    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }
    
    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1_5,Eb,A,Am,A_1_5,Ab\n";

    // 時間ステップ数のループ
    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n; //(10-50-100-200-400-600-800-1000)
        const int paths = 10 * points * points;
        
        const double dt = (t_end - t_start) / (points - 1);
        const double sqrt_dt = sqrt(dt);
        const double dt_sqrt_dt = dt * sqrt_dt;


        // パラメータ初期化
        double S = 0.0, Sm = 0.0, S_1_5 = 0.0, Sb = 0.0;
        double B = 0.0, Bm = 0.0, B_1_5 = 0.0, Bb = 0.0;

        // OpenMP threadの並列化
       #pragma omp parallel reduction(+:S, Sm, S_1_5, Sb, B, Bm, B_1_5, Bb)
        {
            // 各threadは独自の乱数生成器を持つ
            mt19937 rng(42);
            mt19937 rng1(30);
            mt19937 rng2(50);
            normal_distribution<double> dist(mu, sigma);
            
            #pragma omp for schedule(static) nowait//threadごとに均等に計算を分配
            for (int p = 0; p < paths; ++p) {
                // 変数の初期化
                double sum_W = 0.0, sum_W1 = 0.0, sum_W2 = 0.0, sum_Xb = 0.0;
                double I_W_stateb_1 = 0.0, I_W_stateb_2 = 0.0;
                double delta_W = 0.0, delta_W1 = 0.0, delta_W2 = 0.0, delta_Xb = 0.0;
                double X_b_Y,W_state_Y,W_state1_Y,W_state2_Y;
                double W_state = x_0, W_state1 = x_0, W_state2 = x_0,X_b = x_0;
                double c2_sq = 0.0, c4_sq = 0.0;
                
                for (int idx = 1; idx < points; ++idx) {
                    // ランダム数の生成
                    const double Z = dist(rng);
                    const double Z1 = dist(rng1);
                    const double Z2 = dist(rng2);
                    const double Z1_sqrt_dt = Z1 * sqrt_dt, Z2_sqrt_dt = Z2 * sqrt_dt;
                    const double dW = sqrt_dt * Z;
                    const double Z_sq = Z * Z;
                    const double Z_sq_minus_1 = Z_sq - 1.0;
                    const double Z_cube_minus_3Z = Z * (Z_sq - 3.0);
                    
                    // 係数の計算
                    StateCoeff coef_em, coef_m, coef_1_5, coef_X_b;
                    coef_em.compute(a, b, W_state);
                    coef_m.compute(a, b, W_state1);
                    coef_1_5.compute(a, b, W_state2);
                    coef_X_b.compute(a, b, X_b);

                    
                    // 状態の更新
                    W_state_Y = A0(W_state, coef_em, dt, Z);
                    W_state1_Y = A1(W_state1, coef_m, dt, Z);
                    W_state2_Y = A2(W_state2, coef_1_5, dt, Z);
                    double sp_W_stateb = 1.0 / coef_X_b.sigma_sq;
                    X_b_Y = benchmark(X_b, dt, Z, b, a);
                    
                    delta_W  = delta_2(coef_em.drift, coef_em.drift_deriv, coef_em.sigma_deriv, coef_em.sigma, 
                                                   coef_em.sigma_deriv2, dt, W_state, W_state_Y);
                    delta_W1 = delta_2(coef_m.drift, coef_m.drift_deriv, coef_m.sigma_deriv, coef_m.sigma, 
                                                   coef_m.sigma_deriv2, dt, W_state1, W_state1_Y);
                    delta_W2 = delta_2(coef_1_5.drift, coef_1_5.drift_deriv, coef_1_5.sigma_deriv, coef_1_5.sigma, 
                                                   coef_1_5.sigma_deriv2, dt, W_state2, W_state2_Y);
                    delta_Xb = delta_2(coef_X_b.drift, coef_X_b.drift_deriv, coef_X_b.sigma_deriv, coef_X_b.sigma, 
                                                   coef_X_b.sigma_deriv2, dt, X_b, X_b_Y);

                    c2_sq = fabs(c_2_sq(coef_X_b));
                    c4_sq = fabs(c_4_sq(coef_X_b));

                    sum_W += delta_W ;
                    sum_W1 += delta_W1 ;
                    sum_W2 += delta_W2 ;
                    sum_Xb += delta_Xb ;
              
                    I_W_stateb_1 += c2_sq * fabs(sp_W_stateb) * Z1_sqrt_dt;
                    I_W_stateb_2 += c4_sq * sp_W_stateb * sp_W_stateb * Z2_sqrt_dt;
                    
                    //更新過程
                    W_state = W_state_Y;
                    W_state1 = W_state1_Y;
                    W_state2 = W_state2_Y;
                    X_b = X_b_Y;


                }
 

                    // 指数項の計算
                    double I_T = sqrt(2) * I_W_stateb_1 + 2 * sqrt(6) * I_W_stateb_2;
                    double inner = f(I_T);
                    double limit = inner * I_T;

                    //期待値の計算
                    S  += (f(sum_Xb / sqrt(dt)) - f(sum_W / sqrt(dt))) / sqrt(dt);       
                    Sm += (f(sum_Xb / sqrt(dt)) - f(sum_W1 / sqrt(dt))) / sqrt(dt);      
                    S_1_5 += (f(sum_Xb / sqrt(dt)) - f(sum_W2 / sqrt(dt))) / sqrt(dt);
                    Sb += limit; 

                    //分散の計算
                    B += (f(sum_Xb / sqrt(dt)) - f(sum_W / sqrt(dt))) * (f(sum_Xb / sqrt(dt)) - f(sum_W / sqrt(dt))) / dt;        
                    Bm += (f(sum_Xb / sqrt(dt)) - f(sum_W1 / sqrt(dt))) * (f(sum_Xb / sqrt(dt)) - f(sum_W1 / sqrt(dt))) / dt;     
                    B_1_5 += (f(sum_Xb / sqrt(dt)) - f(sum_W2 / sqrt(dt))) * (f(sum_Xb / sqrt(dt)) - f(sum_W2 / sqrt(dt))) / dt;  
                    Bb += limit * limit;

            }
        }// end of parallel region

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
        Eb[n] = Bb * inv_paths - (Ab[n] * Ab[n]);


        // 出力用
        cout << "-------------------------------------------------" << n << "\n";      
        cout << setprecision(10) << "points = " << points << "\n";       
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(10) << "Var E      = " << E[n] << "\n";         
        cout << setprecision(10) << "Var E_m    = " << Em[n] << "\n";        
        cout << setprecision(10) << "Var E_1.5  = " << E_1_5[n] << "\n";
        cout << setprecision(10) << "Var limit_E_b    = " << Eb[n] << "\n";
        cout << setprecision(10) << "Mean A      = " << A[n] << "\n";         
        cout << setprecision(10) << "Mean A_m    = " << Am[n] << "\n";        
        cout << setprecision(10) << "Mean A_1.5  = " << A_1_5[n] << "\n"; 
        cout << setprecision(10) << "Mean limit_A_b    = " << Ab[n] << "\n";

        // CSVファイルに書き込み
        ofs << n << "," << points << ","  
            << fixed << setprecision(10) 
            << E[n] << "," << Em[n] << "," << E_1_5[n] << "," << Eb[n] << ","
            << A[n] << "," << Am[n] << "," << A_1_5[n] << "," << Ab[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}
