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

// ========================================
// φ_n 関数
// ========================================
inline double phi_n(double z_const, double alpha, double x, double h) {
    if (h <= 0.0) throw invalid_argument("h must be positive!");
    const double t = pow(h, alpha);
    const double PI = 3.141592653589793238462643383279502884;
    const double norm = 1.0 / sqrt(2.0 * PI * t);
    const double z_const_sq = (x - z_const) * (x - z_const);
    return norm * exp(-z_const_sq / (2.0 * t));
}

// ========================================
// StateCoeff 
// ========================================
struct StateCoeff {
    double sqrt_W_sq_plus_1;
    double drift, drift_deriv;
    double sigma, sigma_deriv, sigma_deriv2;
    double sigma_sq, sigma_cube;

    inline void compute(double a, double b, double W_state) {
        const double w_sq = W_state * W_state;
        const double W_sq_plus_1 = w_sq + 1.0;
        const double b_sq = b * b, b_quad = b_sq * b_sq;

        sqrt_W_sq_plus_1 = sqrt(W_sq_plus_1);
        drift = 0.5 * b_sq * W_state + a * sqrt_W_sq_plus_1 * asinh(W_state);
        drift_deriv = 0.5 * b_sq + a + (a * W_state * asinh(W_state) / sqrt_W_sq_plus_1);
        sigma = b * sqrt_W_sq_plus_1;
        sigma_sq = sigma * sigma;
        sigma_deriv = b_sq * W_state / sigma;
        sigma_deriv2 = b_quad / (sigma * sigma * sigma);
    }
};

// Euler-Maruyama
inline double A0(double W_state, const StateCoeff &coef, double dt, double sqrt_dt, double Z) {
    return W_state + coef.drift * dt + coef.sigma * sqrt_dt * Z;
}

// Milstein
inline double A1(double W_state, const StateCoeff &coef, double dt, double Z, double sqrt_dt) {
    const double Z_sq = Z * Z;
    return W_state + coef.drift * dt + coef.sigma * sqrt_dt * Z +
           0.5 * coef.sigma_deriv * coef.sigma * (Z_sq * dt - dt);
}

// 1.5 阶随机泰勒展开
inline double A2(double W_state, const StateCoeff &coef, double dt, double Z, double sqrt_dt) {
    const double Z_sq = Z * Z;
    const double dW = Z * sqrt_dt;
    const double dW_cube = dW * dW * dW;
    const double base = W_state + coef.drift * dt + coef.sigma * dW;
    const double milstein_term = 0.5 * coef.sigma_deriv * coef.sigma * (Z_sq * dt - dt);

    const double mixed_deriv = 0.5 * coef.drift_deriv * coef.sigma +
                               0.5 * coef.sigma_deriv * coef.drift +
                               0.25 * coef.sigma_deriv2 * coef.sigma_sq;

    const double term3 = mixed_deriv * dW * dt;
    const double triple_deriv =
        coef.sigma_deriv2 * coef.sigma_sq + coef.sigma_deriv * coef.sigma_deriv * coef.sigma;
    const double term4 = triple_deriv * (dW_cube - 3 * dW * dt) / 6.0;
    return base + milstein_term + term3 + term4;
}

// Y_t = e^{a t} asinh(x0) + b (α_t W_t + β_t W_t')
inline double benchmark(double X_b, double t, double W, double Wp, double b, double a) {
    const double asinh_Xb = asinh(X_b);
    const double exp_at = exp(a * t), exp_2at = exp(2 * a * t), t_a_inv = 1.0 / (a * t), t_2a_inv = 1.0 / (2 * a * t);
    const double alpha_t =(exp_at - 1) * t_a_inv ;   
    const double beta_first = (exp_2at - 1) * t_2a_inv;                  // α_t = (e^{a t} - 1) / (a t)
    const double beta_t = sqrt( beta_first - alpha_t * alpha_t);
    return exp_at * asinh_Xb + b * (alpha_t * W + beta_t * Wp);
}

// ========================================
// メインプログラム
// ========================================
int main() {

    constexpr double z_const = 1.0;
    constexpr double alpha = 1.0;
    constexpr double t_start = 0.0;
    constexpr double t_end = 1.0;
    constexpr double mu = 0.0;
    constexpr double sigma = 1.0;
    constexpr double b = 0.5;
    constexpr double a = 0.5;
    constexpr double x_0 = 1.0;
    constexpr int max_n = 10;

    vector<double> A(max_n + 1, 0.0), Am(max_n + 1, 0.0), A_1_5(max_n + 1, 0.0), Ab(max_n + 1, 0.0);
    vector<double> E(max_n + 1, 0.0), Em(max_n + 1, 0.0), E_1_5(max_n + 1, 0.0), Eb(max_n + 1, 0.0);

    const string dir_path = "data_source";
    system(("mkdir -p " + dir_path).c_str());
    const string csv_path = dir_path + "/LTM2_test2_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);

    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }

    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1_5,Eb,A,Am,A_1_5,Ab\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + (n * 10);
        const int paths = 8 * points * points;
        const double dt = (t_end - t_start) / (points - 1);
        const double sqrt_dt = sqrt(dt);

        double S = 0.0, Sm = 0.0, S_1_5 = 0.0, Sb = 0.0;
        double B = 0.0, Bm = 0.0, B_1_5 = 0.0, Bb = 0.0;

#pragma omp parallel reduction(+ : S, Sm, S_1_5, Sb, B, Bm, B_1_5, Bb)
        {
            mt19937 rng(42);
            mt19937 rng1(30);
            normal_distribution<double> dist(mu, sigma);

            #pragma omp for schedule(static)
            for (int p = 0; p < paths; ++p) {
                double W_state = x_0, W_state1 = x_0, W_state2 = x_0;
                double X_b_Y = x_0, X_b = x_0;
                double dX0 = 0.0, dX1 = 0.0, dX2 = 0.0, dX_b = 0.0;
                double W_sum = 0.0, Wp_sum = 0.0; // 独立したブラウン運動の和

                for (int idx = 1; idx < points; ++idx) {
                    double Z = dist(rng);
                    double Z1 = dist(rng1);
                    double dW = sqrt_dt * Z;
                    double dWp = sqrt_dt * Z1;
                    W_sum += dW;
                    Wp_sum += dWp;
                    double t = idx * dt;

                    StateCoeff coef_em, coef_m, coef_1_5;
                    coef_em.compute(a, b, W_state);
                    coef_m.compute(a, b, W_state1);
                    coef_1_5.compute(a, b, W_state2);

                    double W_state_Y = A0(W_state, coef_em, dt, sqrt_dt, Z);
                    double W_state1_Y = A1(W_state1, coef_m, dt, Z, sqrt_dt);
                    double W_state2_Y = A2(W_state2, coef_1_5, dt, Z, sqrt_dt);
                    double X_b_Y = benchmark(X_b, t, W_sum, Wp_sum, b, a);

                    dX0 += dt * phi_n(z_const, alpha, W_state_Y, dt);
                    dX1 += dt * phi_n(z_const, alpha, W_state1_Y, dt);
                    dX2 += dt * phi_n(z_const, alpha, W_state2_Y, dt);
                    dX_b += dt * phi_n(z_const, alpha, X_b_Y, dt);

                    W_state = W_state_Y;
                    W_state1 = W_state1_Y;
                    W_state2 = W_state2_Y;
                }

                S += dX0;
                Sm += dX1;
                S_1_5 += dX2;
                Sb += dX_b;
                B += dX0 * dX0;
                Bm += dX1 * dX1;
                B_1_5 += dX2 * dX2;
                Bb += dX_b * dX_b;
            }
        }

        const double inv_paths = 1.0 / paths;
        A[n] = S * inv_paths;
        Am[n] = Sm * inv_paths;
        A_1_5[n] = S_1_5 * inv_paths;
        Ab[n] = Sb * inv_paths;

        E[n] = B * inv_paths - A[n] * A[n];
        Em[n] = Bm * inv_paths - Am[n] * Am[n];
        E_1_5[n] = B_1_5 * inv_paths - A_1_5[n] * A_1_5[n];
        Eb[n] = Bb * inv_paths - Ab[n] * Ab[n];

        cout << "-------------------------------------------------" << n << "\n";
        cout << "points = " << points << "\n";
        cout << "E0  = " << E[n] << "\nE1  = " << Em[n] << "\nE2  = " << E_1_5[n] << "\nEb  = " << Eb[n] << "\n";
        cout << "A0  = " << A[n] << "\nA1  = " << Am[n] << "\nA2  = " << A_1_5[n] << "\nA_b= " << Ab[n] << "\n";

        ofs << n << "," << points << "," << fixed << setprecision(10)
            << E[n] << "," << Em[n] << "," << E_1_5[n] << "," << Eb[n] << ","
            << A[n] << "," << Am[n] << "," << A_1_5[n] << "," << Ab[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}