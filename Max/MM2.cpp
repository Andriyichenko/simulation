// MM2 with x_0=1, T=1, a=0.5, b=0.5

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
    double sigma_sq, sigma_cube;                                             // sigma(W)^2, sigma(W)^3
    
    
    // 係数の計算
inline void compute( double a, double b, double W_state) {
        const double w_sq = W_state * W_state;
        const double W_sq_plus_1 = w_sq + 1.0;
        const double b_sq = b * b, b_quad = b * b * b, a_b = a * b;
    
    
        sqrt_W_sq_plus_1 = sqrt(W_sq_plus_1);

        //a(x)の計算部分
        drift = 0.5 * b_sq * W_state + a * sqrt_W_sq_plus_1 * asinh(W_state); //a_x
        drift_deriv = 0.5 * b_sq + a + (a * W_state * asinh(W_state) / sqrt_W_sq_plus_1); //a_x'
        
        // sigma(x)
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
            sigma_deriv = b * b * W_state * sigma_inv;
            sigma_deriv2 = b * b * b * sigma_inv * sigma_inv * sigma_inv;
        }

       
    }
};

// ========================================
// 近似の更新関数の定義
// ========================================

// Euler-Maruyama
inline double A0(double W_state, const StateCoeff& coef,
                                         double dt, double sqrt_dt, double Z) {
    return W_state + coef.drift * dt + coef.sigma * sqrt_dt * Z;
}

// Milstein
inline double A1(double W_state, const StateCoeff& coef,
                                   double dt, double Z, double sqrt_dt) {

    const double Z_sq = Z * Z;
    return W_state + coef.drift * dt + coef.sigma * sqrt_dt * Z + 
           0.5 * coef.sigma_deriv * coef.sigma * (Z_sq * dt - dt);
}

// 1.5 
inline double A2(double W_state, const StateCoeff& coef,
                            double dt, double Z, double sqrt_dt) {
    const double dt_sq = dt * dt;
    const double Z_sq = Z * Z;
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
inline double benchmark_opt(double x_0, double t, double dW, double dW_prime, double b, double a) {
    const double Y_0 = asinh(x_0);
    const double a_t = a * t, t_a_inv = 1.0/a_t,t_2a_inv = 1.0/(2*a_t);
    const double exp_at = exp(a * t), exp_2at = exp(2 * a * t);
    const double alpha_t = (exp_at - 1) * t_a_inv; 
    const double beta_first = (exp_2at - 1) * t_2a_inv;
    const double beta_t = sqrt(beta_first - (alpha_t * alpha_t));

    return sinh(exp_at * Y_0 + b * (alpha_t * dW + beta_t * dW_prime));

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
    const string dir_path = "data_source";
    system(("mkdir -p " + dir_path).c_str()); //フォルダーの確認 
    const string csv_path = dir_path + "/MM2_100_1000_error_data.csv"; //data sourceのファイル名指定
    ofstream ofs(csv_path, ios::out | ios::trunc);
    
    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }
    
    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1_5,A,Am,A_1_5\n";

    // 時間ステップ数のループ
    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + (n * 100); //(10-50-100-200-400-600-800-1000)
        const int paths = 8 * points * points;
        
        const double dt = (t_end - t_start) / (points - 1);
        const double sqrt_dt = sqrt(dt);
        const double dt_sqrt_dt = dt * sqrt_dt;


        // パラメータ初期化
        double S = 0.0, Sm = 0.0, S_1_5 = 0.0, Sb = 0.0;
        double B = 0.0, Bm = 0.0, B_1_5 = 0.0, Bb = 0.0;

        // OpenMP threadの並列化
        #pragma omp parallel reduction(+:S,Sm,S_1_5,B,Bm,B_1_5) //OpenMPのreduction句で指定の変数の和(+:)を並列計算 i.e. reduction(operator : variable_list)
        {
            mt19937 rng(42);
            mt19937 rng1(27);
            normal_distribution<double> dist(mu, sigma);
            #pragma omp for schedule(static) nowait//threadごとに均等に計算を分配
            for (int p = 0; p < paths; ++p) {
                // 変数の初期化
                double W_state = x_0, W_state1 = x_0, W_state2 = x_0;
                double X_b = x_0,X_b_Y=x_0,W_state_Y=x_0,W_state1_Y=x_0,W_state2_Y=x_0;
                double dX0 = 0.0, dX1 = 0.0, dX2 = 0.0;

                for (int idx = 1; idx < points; ++idx) {
                    // ランダム数の生成
                    const double Z = dist(rng), Z1 = dist(rng1);
                    double dW = sqrt_dt * Z;
                    double dWp = sqrt_dt * Z1;
   
                    // 係数の計算
                    StateCoeff coef_em, coef_m, coef_1_5;
                    coef_em.compute(a, b, W_state);
                    coef_m.compute(a, b, W_state1);
                    coef_1_5.compute(a, b, W_state2);
                    
                    // 状態の更新
                    W_state_Y = A0(W_state, coef_em, dt, sqrt_dt, Z);
                    W_state1_Y = A1(W_state1, coef_m, dt, Z, sqrt_dt);
                    W_state2_Y = A2(W_state2, coef_1_5, dt,Z, sqrt_dt);

                    //benchmarkの更新
                    X_b_Y = benchmark_opt(X_b, dt, dW, dWp, b, a);

                    //X_bとWの誤差の計算
                    double diff0  = X_b_Y - W_state_Y;    
                    double diff1 = X_b_Y - W_state1_Y;   
                    double diff2 = X_b_Y - W_state2_Y; 

                    // max二乗誤差の更新
                    dX0  = max(dX0,  diff0  * diff0);   
                    dX1 = max(dX1, diff1 * diff1);   
                    dX2 = max(dX2, diff2 * diff2);  

                    //wとx_bの更新
                    W_state1 = W_state1_Y;
                    W_state2 = W_state2_Y;
                    W_state = W_state_Y;
                    X_b = X_b_Y;
                    

                }

                // 期待値の累計
                S += dX0;
                Sm += dX1;
                S_1_5 += dX2;

                // 分散用の累計
                B += dX0 * dX0;
                Bm += dX1 * dX1;
                B_1_5 += dX2 * dX2;
            }
            
           
        } // end of parallel region

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

        // 出力用
        cout << "-------------------------------------------------" << n << "\n";      
        cout << setprecision(10) << "points = " << points << "\n";       
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(10) << "E0      = " << E[n] << "\n";        
        cout << setprecision(10) << "E1      = " << Em[n] << "\n";       
        cout << setprecision(10) << "E2      = " << E_1_5[n] << "\n";    
        cout << setprecision(10) << "A0      = " << A[n] << "\n";        
        cout << setprecision(10) << "A1      = " << Am[n] << "\n";       
        cout << setprecision(10) << "A2      = " << A_1_5[n] << "\n";    

        // CSVファイルに書き込み
        ofs << n << "," << points << ","  
            << fixed << setprecision(10) 
            << E[n] << "," << Em[n] << "," << E_1_5[n] << ","
            << A[n] << "," << Am[n] << "," << A_1_5[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}
