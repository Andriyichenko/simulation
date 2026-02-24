// LTM2_EFF with x_0=1, T=1, z=0.5, a=0.5, b=0.5,alpha=1
//E[F(\bar{X}) - F(\bar{X}^\alpha)]

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



inline double phi_n(double z_const,double alpha,double x,double h){
    if (h <=0.0) throw invalid_argument("h must be positive!");
    const double t =pow(h,alpha);
    const double PI = 3.141592653589793238462643383279502884;
    const double norm = 1.0 / sqrt(2.0 * PI * t);
    const double z_const_sq = (x - z_const) * (x - z_const);
    return norm * exp(-z_const_sq / (2.0 * t));

}

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

// テスト関数 f(L) = arctan(L)
inline double f(double x) {
    return atan(sqrt(x));
}






// ========================================
// メインプログラム
// ========================================

int main() {  

    constexpr double z_const=2.4;//z_constは小さいほど精度が良い
    constexpr double alpha=1.0;//alpha \in (0,2) 1 と　0.1（分散が小さいが、数値漏れる） と1.9（分散が大きい）
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
    const string csv_path = dir_path + "/LTM2_EFF_100_1000_data.csv"; //data sourceのファイル名指定
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
        const int paths = 10 * points * points;
        
        const double dt = (t_end - t_start) / (points - 1);
        const double sqrt_dt = sqrt(dt);
        const double dt_sqrt_dt = dt * sqrt_dt;


        // パラメータ初期化
        double S = 0.0, Sm = 0.0, S_1_5 = 0.0, Sb = 0.0;
        double B = 0.0, Bm = 0.0, B_1_5 = 0.0, Bb = 0.0;

        // OpenMP threadの並列化
        #pragma omp parallel reduction(+:S,Sm,S_1_5,B,Bm,B_1_5) //OpenMPのreduction句で指定の変数の和(+:)を並列計算 i.e. reduction(operator : variable_list)
        {
            normal_distribution<double> dist(mu, sigma);
            #pragma omp for schedule(static) nowait//threadごとに均等に計算を分配
            for (int p = 0; p < paths; ++p) {
                seed_seq ss0{42u, 0u, (uint32_t)p};
                seed_seq ss1{27u, 1u, (uint32_t)p};
                mt19937 rng(ss0);
                mt19937 rng1(ss1);
                // 変数の初期化
                double W_state = x_0, W_state1 = x_0, W_state2 = x_0;
                double X_b = x_0,X_b_Y=x_0,W_state_Y=x_0,W_state1_Y=x_0,W_state2_Y=x_0;
                double L_x0 = 0.0, L_xb = 0.0, L_x1 = 0.0, L_x2 = 0.0;

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
                    W_state_Y= A0(W_state, coef_em, dt, Z);
                    W_state1_Y = A1(W_state1, coef_m, dt, Z);
                    W_state2_Y = A2(W_state2, coef_1_5, dt, Z);

                    //benchmarkの更新
                    X_b_Y = benchmark(X_b, dt, dW, dWp, b, a);

                    // 誤差の累計
                    L_xb += dt * phi_n(z_const, alpha, X_b_Y, dt);
                    L_x0 += dt * phi_n(z_const, alpha, W_state_Y, dt);
                    L_x1 += dt * phi_n(z_const, alpha, W_state1_Y, dt);
                    L_x2 += dt * phi_n(z_const, alpha, W_state2_Y, dt);

                    //wとx_bの更新
                    W_state1 = W_state1_Y;
                    W_state2 = W_state2_Y;
                    W_state = W_state_Y;
                    X_b = X_b_Y;
                    

                }

                // LTのF関数の期待値の累計
                S += f(L_xb) - f(L_x0);
                Sm += f(L_xb) - f(L_x1);
                S_1_5 += f(L_xb) - f(L_x2);
                // LTのF関数の分散用の累計
                B += (f(L_xb) - f(L_x0)) * (f(L_xb) - f(L_x0));
                Bm += (f(L_xb) - f(L_x1)) * (f(L_xb) - f(L_x1));
                B_1_5 += (f(L_xb) - f(L_x2)) * (f(L_xb) - f(L_x2));
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
        cout << setprecision(10) << "points = " << points << "\n";       // ステップ数 points をコンソール表示
        cout << "-------------------------------------------------" <<  "\n";      
        cout << setprecision(10) << "E0      = " << E[n] << "\n";         // EM 法の推定値 A[n] をコンソール表示
        cout << setprecision(10) << "E1      = " << Em[n] << "\n";         // Milstein 法の推定値 Am[n] をコンソール表示
        cout << setprecision(10) << "E2      = " << E_1_5[n] << "\n";      // 1.5 次法の推定値 A_1_5[n] をコンソール表示
        cout << setprecision(10) << "A0      = " << A[n] << "\n";         // EM 法の推定値 A[n] をコンソール表示
        cout << setprecision(10) << "A1      = " << Am[n] << "\n";         // Milstein 法の推定値 Am[n] をコンソール表示
        cout << setprecision(10) << "A2      = " << A_1_5[n] << "\n";      // 1.5 次法の推定値 A_1_5[n] をコンソール表示

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
