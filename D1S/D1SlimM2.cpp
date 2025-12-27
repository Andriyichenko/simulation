// D1SlimM2 with x_0=1, T=1, a=0.5, b=0.5
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
// パラメータの設定
// ========================================

struct StateCoeff {
    double sqrt_W_sq_plus_1, sqrt_X_b_sq_plus_1;                              //sqrt(W^2 + 1), sqrt(X_b^2 + 1)
    double drift, drift_deriv, drift_X_b, drift_X_b_deriv;                       // a(W), a'(W), a_m(W), a_m'(W)
    double sigma, sigma_deriv,sigma_deriv2, sigma_inv;                       // sigma(W), sigma'(W), sigma_m(W), sigma_m'(W)
    double sigma_X_b, sigma_X_b_deriv;                                       // sigma(X_b), sigma''(X_b)
    double sigma_sq, sigma_cube; 

    // 係数の計算
inline void compute(double a, double b, double W_state) {
    const double w_sq = W_state * W_state;
    const double W_sq_plus_1 = w_sq + 1.0;
    const double b_sq = b * b;
    
    sqrt_W_sq_plus_1 = sqrt(W_sq_plus_1);

    // a(x)
    drift = 0.5 * b_sq * W_state + a * sqrt_W_sq_plus_1 * asinh(W_state);
    drift_deriv = 0.5 * b_sq + a + (a * W_state * asinh(W_state) / sqrt_W_sq_plus_1);
    
    // sigma(x) = b*sqrt(x^2+1)
    if (fabs(b) < 1e-12) {
        sigma = 0.0;
        sigma_inv = 0.0;
        sigma_sq = 0.0;
        sigma_deriv = 0.0;
        sigma_deriv2 = 0.0;
    } else {
        sigma = b * sqrt_W_sq_plus_1;
        sigma_inv = 1.0 / sigma;
        sigma_sq = sigma * sigma;
        
        // sigma'(x) = b*x / sqrt(x^2+1)
        sigma_deriv = b * W_state / sqrt_W_sq_plus_1;
        
        // sigma''(x) = b / (x^2+1)^(3/2)
        const double W_sq_plus_1_pow_1_5 = W_sq_plus_1 * sqrt_W_sq_plus_1;
        sigma_deriv2 = b / W_sq_plus_1_pow_1_5;
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

// M2のbenchmark関数の定義
inline double benchmark(double X_b, double t, double dW, double dW1, double b, double a) {
    const double Y_0 = asinh(X_b);
    const double a_t = a * t, t_a_inv = 1.0/a_t,t_2a_inv = 1.0/(2*a_t);
    const double exp_at = exp(a * t), exp_2at = exp(2 * a * t);
    const double alpha_t = (exp_at - 1) * t_a_inv; 
    const double beta_first = (exp_2at - 1) * t_2a_inv;
    const double beta_t = sqrt(beta_first - (alpha_t * alpha_t));

    return sinh(exp_at * Y_0 + b * (alpha_t * dW + beta_t * dW1));

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
    const string csv_path = dir_path + "/D1SlimM2_100_1000_data.csv"; //data sourceのファイル名指定
    ofstream ofs(csv_path, ios::out | ios::trunc);
    
    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }
    
    ofs.imbue(locale::classic());
    ofs << "n,points,Eb,Ab\n";

    // 時間ステップ数のループ
    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n; //(10-50-100-200-400-600-800-1000)
        const int paths = 100 * points * points;
        
        const double dt = (t_end - t_start) / (points - 1);
        const double sqrt_dt = sqrt(dt);
        const double dt_sqrt_dt = dt * sqrt_dt;


        // パラメータ初期化
        double Sb = 0.0, Bb = 0.0;


    


        // OpenMP threadの並列化
        #pragma omp parallel reduction(+:Sb, Bb)
        {
            // 各threadは独自の乱数生成器を持つ
            mt19937 rng(42);
            mt19937 rng1(30);
            mt19937 rng2(56);
            normal_distribution<double> dist(mu, sigma);
            
            #pragma omp for schedule(dynamic, 64) nowait
            for (int p = 0; p < paths; ++p) {
                // 変数の初期化
                double X_b = x_0,X_b_Y=x_0;
                // 積分変数の初期化
                double I_W_stateb = 0.0;
                double I_quad_W_stateb = 0.0;
                
                for (int idx = 1; idx < points; ++idx) {
                    // ランダム数の生成
                    const double Z = dist(rng);
                    const double Z1 = dist(rng1);
                    const double Z2 = dist(rng2);
                    double Z2_sqrt_dt = Z2 * sqrt_dt;
                    double dW = sqrt_dt * Z;
                    double dW1 = sqrt_dt * Z1;
                    double dW2 = sqrt_dt * Z2;
             
                    
                    // 係数の計算
                    StateCoeff coefb;
                    coefb.compute(a, b, X_b_Y);

                    // diffusion係数の導数の絶対値
                    double sp_W_stateb = fabs(coefb.sigma_deriv);

                    // 状態の更新
                    X_b_Y = benchmark(X_b, dt, dW, dW1, b, a);
                    
                    // 積分項の更新
                    I_W_stateb += sqrt(1.5) * sp_W_stateb * Z2_sqrt_dt;
                    I_quad_W_stateb += 1.5 * sp_W_stateb * sp_W_stateb * dt;
                    X_b = X_b_Y;


                }

                // 指数項の計算
                double inner_b = exp(-I_W_stateb - 0.5 * I_quad_W_stateb) - 1.0;

                // 絶対値の計算
                const double abs_inner_b = fabs(inner_b);

                // 誤差の累計
                Sb += abs_inner_b;
                Bb += abs_inner_b * abs_inner_b;

            }
           
        } // end of parallel region

        // 期待値の計算
        const double inv_paths = 1.0 / paths;
        Ab[n] = Sb * inv_paths;
        
        // 分散の計算
        Eb[n] = Bb * inv_paths - Ab[n] * Ab[n];

        // 出力用
        cout << "-------------------------------------------------" << n << "\n";      
        cout << setprecision(10) << "points = " << points << "\n";       
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(10) << "Eb     = " << Eb[n] << "\n";     
        cout << setprecision(10) << "Ab     = " << Ab[n] << "\n";    

        // CSVファイルに書き込み
        ofs << n << "," << points << ","  
            << fixed << setprecision(10) 
            << Eb[n] << "," << Ab[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}