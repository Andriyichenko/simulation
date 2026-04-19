// 2DD0F0M3_all_with_NV.cpp
// Error simulation for Model 3 with EM, Milstein, Strong 1.5, and Ninomiya-Victoir
// The code also keeps the limit simulation part from the base program.

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <random>
#include <string>
#include <vector>
// #include <omp.h>

using namespace std;
using namespace Eigen;

// Use Eigen's 2D vector as the state
using State = Vector2d;

// ========================================
// Utility functions
// ========================================

// Sign function
constexpr int sgn(double x) {
    if (std::isnan(x)) return 0;
    return (x > 0.0) - (x < 0.0);
}

// Clipped helper (kept from the base style, not used below)
inline double f_sgn(double x, double min_val = -100.0, double max_val = 0.0) {
    return std::max(min_val, std::min(x, max_val));
}

// Compute the sum used in the error functional
inline double compute_sum_state(double delta_val, double delta_val_sq) {
    return delta_val - 0.5 * delta_val_sq;
}

// s(x) = 2 + sin(x)
inline double s_func(double x) { return 2.0 + std::sin(x); }

// s'(x) = cos(x)
inline double ds_func(double x) { return std::cos(x); }

// s''(x) = -sin(x)
inline double dds_func(double x) { return -std::sin(x); }

// ========================================
// Hermite polynomial helper functions
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
    return std::pow(r1, 4) / (t4 * S2_4)
         - 6.0 * std::pow(r1, 2) / (t3 * S2_3)
         + 3.0 / (t2 * S2_2);
}

inline double H2222(double t, double r2, double S1) {
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t2 * t2;
    double S1_2 = S1 * S1;
    double S1_3 = S1_2 * S1;
    double S1_4 = S1_2 * S1_2;
    return std::pow(r2, 4) / (t4 * S1_4)
         - 6.0 * std::pow(r2, 2) / (t3 * S1_3)
         + 3.0 / (t2 * S1_2);
}

inline double H1112(double t, double r1, double r2, double S1, double S2) {
    double t3 = t * t * t;
    double t4 = t3 * t;
    return (std::pow(r1, 3) * r2) / (t4 * std::pow(S2, 3) * S1)
         - (3.0 * r1 * r2) / (t3 * std::pow(S2, 2) * S1);
}

inline double H1222(double t, double r1, double r2, double S1, double S2) {
    double t3 = t * t * t;
    double t4 = t3 * t;
    return (r1 * std::pow(r2, 3)) / (t4 * S2 * std::pow(S1, 3))
         - (3.0 * r1 * r2) / (t3 * S2 * std::pow(S1, 2));
}

inline double H1122(double t, double r1, double r2, double S1, double S2) {
    return H11(t, r1, S2) * H22(t, r2, S1);
}

// ========================================
// Delta calculation
// ========================================

// Delta^1 formula
inline double Delta1(const State& x, const State& y, double t) {
    double s_x1 = s_func(x(0));
    double s_x2 = s_func(x(1));
    double ds_x1 = ds_func(x(0));
    double ds_x2 = ds_func(x(1));

    // Model 3 drift: a(x) = (x2, -x1)
    double a1 = x(1);
    double a2 = -x(0);

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
// Simulation schemes
// ========================================

// Euler-Maruyama (order 0.5)
inline State A0(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = std::sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_func(curr(1)) * dW(0), s_func(curr(0)) * dW(1));
    return curr + drift * dt + diffusion;
}

// Milstein (order 1.0)
inline State A1(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = std::sqrt(dt);
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

// Strong order 1.5 scheme
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = std::sqrt(dt);
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

    // X1 component
    next(0) += a(0) * dt + s1 * w1;
    next(0) += 0.5 * ds1 * s0 * w1 * w2;
    next(0) += 0.5 * dt * (s0 * w2 + ds1 * a(1) * w1);

    double term_t8_1 = 2.0 * dds1 * s0 * s0 * w1;
    double term_t8_2 = (ds1 * s0 * ds1 * s0 / s1) * w1;
    double term_t8_3 = -ds0 * ds1 * s1 * w2;
    next(0) += 0.125 * dt * (term_t8_1 + term_t8_2 + term_t8_3);

    double term_h1 = (1.0 / 24.0) * (ds1 * ds1 / s1) * s0 * s0 * w111;
    double term_h2 = (1.0 / 12.0) * ds0 * ds1 * s1 * w112;
    double term_h3 = (1.0 / 6.0)  * dds1 * s0 * s0 * w122;
    double term_h4 = (1.0 / 24.0) * (ds1 * ds1 / s1) * s0 * s0 * w122;
    next(0) += term_h1 + term_h2 + term_h3 + term_h4;

    // X2 component
    next(1) += a(1) * dt + s0 * w2;
    next(1) += 0.5 * ds0 * s1 * w1 * w2;
    next(1) += 0.5 * dt * (-s1 * w1 + ds0 * a(0) * w2);

    double term2_t8_1 = 2.0 * dds0 * s1 * s1 * w2;
    double term2_t8_2 = ds0 * s1 * ds0 * s1 * w2 / s0;
    double term2_t8_3 = -ds1 * ds0 * s0 * w1;
    next(1) += 0.125 * dt * (term2_t8_1 + term2_t8_2 + term2_t8_3);

    double term2_h1 = (1.0 / 24.0) * (ds0 * ds0 / s0) * s1 * s1 * w222;
    double term2_h2 = (1.0 / 12.0) * ds0 * ds1 * s0 * w122;
    double term2_h3 = (1.0 / 6.0)  * dds0 * s1 * s1 * w112;
    double term2_h4 = (1.0 / 24.0) * (ds0 * ds0 / s0) * s1 * s1 * w112;
    next(1) += term2_h1 + term2_h2 + term2_h3 + term2_h4;

    return next;
}

// Ninomiya-Victoir scheme
inline State A3(const State& curr, double dt, double Z1, double Z2, int xi) {
    const double sqrt_dt = std::sqrt(dt);
    double dW1 = sqrt_dt * Z1;
    double dW2 = sqrt_dt * Z2;

    // Exact drift flow: phi_0
    auto flow_drift = [](double t, const State& x) -> State {
        double c = std::cos(t);
        double s = std::sin(t);
        return State(c * x(0) + s * x(1), -s * x(0) + c * x(1));
    };

    // Exact sigma_1 flow
    auto flow_sigma1 = [](double z, const State& x) -> State {
        return State(x(0) + z * s_func(x(1)), x(1));
    };

    // Exact sigma_2 flow
    auto flow_sigma2 = [](double z, const State& x) -> State {
        return State(x(0), x(1) + z * s_func(x(0)));
    };

    State y = flow_drift(0.5 * dt, curr);

    if (xi >= 0) {
        y = flow_sigma1(dW1, y);
        y = flow_sigma2(dW2, y);
    } else {
        y = flow_sigma2(dW2, y);
        y = flow_sigma1(dW1, y);
    }

    y = flow_drift(0.5 * dt, y);
    return y;
}

// ========================================
// Main program
// ========================================

int main() {
    constexpr double t_start = 0.0;
    constexpr double t_end   = 1.0;
    constexpr double mu      = 0.0;
    constexpr double sigma   = 1.0;

    const State x0_state(1.0, 1.0);
    constexpr int max_n = 9;

    vector<double> A(max_n + 1, 0.0);
    vector<double> Am(max_n + 1, 0.0);
    vector<double> A_1_5(max_n + 1, 0.0);
    vector<double> A_nv(max_n + 1, 0.0);

    vector<double> E(max_n + 1, 0.0);
    vector<double> Em(max_n + 1, 0.0);
    vector<double> E_1_5(max_n + 1, 0.0);
    vector<double> E_nv(max_n + 1, 0.0);

    // Statistics for the limit simulation
    vector<double> A_lim(max_n + 1, 0.0);
    vector<double> E_lim(max_n + 1, 0.0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str());

    const string csv_path = dir_path + "/TEST_2DD0F0M3_all_with_NV_100_1000_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);

    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }

    ofs.imbue(locale::classic());
    ofs << "n,points,E,Em,E_1_5,E_nv,E_lim,A,Am,A_1_5,A_nv,A_lim\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n;
        const int paths  = 10 * points * points;

        const double dt = (t_end - t_start) / (points - 1);

        // Fine benchmark subdivision count
        const int M = points - 1;
        const double dtm = dt / M;
        const double sqrt_dtm = std::sqrt(dtm);

        long double S      = 0.0L;
        long double Sm     = 0.0L;
        long double S_1_5  = 0.0L;
        long double S_nv   = 0.0L;
        long double S_lim  = 0.0L;

        long double B      = 0.0L;
        long double Bm     = 0.0L;
        long double B_1_5  = 0.0L;
        long double B_nv   = 0.0L;
        long double B_lim  = 0.0L;

        #pragma omp parallel reduction(+:S,Sm,S_1_5,S_nv,S_lim,B,Bm,B_1_5,B_nv,B_lim)
        {
            normal_distribution<double> dist(mu, sigma);
            bernoulli_distribution dist_xi(0.5);

            #pragma omp for schedule(static) 
            for (int p = 0; p < paths; ++p) {
                dist.reset();

                seed_seq ss0{40u, 0u, static_cast<uint32_t>(p)};
                seed_seq ss1{50u, 1u, static_cast<uint32_t>(p)};
                seed_seq ss2{1000u, 2u, static_cast<uint32_t>(p)};
                seed_seq ss3{2000u, 3u, static_cast<uint32_t>(p)};

                mt19937 rng_nm(ss0);
                mt19937 rng1_nm(ss1);
                mt19937 rng_lim(ss2);
                mt19937 rng_xi(ss3);

                State st_em  = x0_state;
                State st_mil = x0_state;
                State st_15  = x0_state;
                State st_nv  = x0_state;
                State st_nm  = x0_state;

                double D_A0    = 0.0;
                double D_A1    = 0.0;
                double D_A2    = 0.0;
                double D_A3    = 0.0;
                double D_nm    = 0.0;

                double D_A0_sq = 0.0;
                double D_A1_sq = 0.0;
                double D_A2_sq = 0.0;
                double D_A3_sq = 0.0;
                double D_nm_sq = 0.0;

                double I_T = 0.0;
                double I_T_inner = 0.0;

                for (int idx = 1; idx < points; ++idx) {
                    double Z1 = 0.0;
                    double Z2 = 0.0;

                    // Fine-grid benchmark and limit simulation
                    for (int m = 0; m < M; ++m) {
                        double Z1_nm = dist(rng_nm);
                        double Z2_nm = dist(rng1_nm);
                        double Z_tilde = dist(rng_lim);

                        double x1 = st_nm(0);
                        double x2 = st_nm(1);
                        double s1 = s_func(x1);
                        double s2 = s_func(x2);
                        double ds1 = ds_func(x1);
                        double ds2 = ds_func(x2);

                        double term1 = (ds1 / s1) * s2;
                        double term2 = (ds2 / s2) * s1;
                        double integrand = std::sqrt(term1 * term1 + term2 * term2);

                        I_T += std::sqrt(0.5) * integrand * Z_tilde * sqrt_dtm;
                        I_T_inner += 0.5 * integrand * integrand * dtm;

                        State nm_benchmark = A1(st_nm, dtm, Z1_nm, Z2_nm);

                        double delta_nm = Delta1(st_nm, nm_benchmark, dtm);
                        D_nm += delta_nm;
                        D_nm_sq += delta_nm * delta_nm;

                        st_nm = nm_benchmark;

                        Z1 += Z1_nm / std::sqrt(static_cast<double>(M));
                        Z2 += Z2_nm / std::sqrt(static_cast<double>(M));
                    }

                    int xi = dist_xi(rng_xi) ? 1 : -1;

                    State next_em  = A0(st_em, dt, Z1, Z2);
                    State next_mil = A1(st_mil, dt, Z1, Z2);
                    State next_15  = A2(st_15, dt, Z1, Z2);
                    State next_nv  = A3(st_nv, dt, Z1, Z2, xi);

                    double delta_em  = Delta1(st_em,  next_em,  dt);
                    double delta_mil = Delta1(st_mil, next_mil, dt);
                    double delta_15  = Delta1(st_15,  next_15,  dt);
                    double delta_nv  = Delta1(st_nv,  next_nv,  dt);

                    D_A0 += delta_em;
                    D_A1 += delta_mil;
                    D_A2 += delta_15;
                    D_A3 += delta_nv;

                    D_A0_sq += delta_em  * delta_em;
                    D_A1_sq += delta_mil * delta_mil;
                    D_A2_sq += delta_15  * delta_15;
                    D_A3_sq += delta_nv  * delta_nv;

                    st_em  = next_em;
                    st_mil = next_mil;
                    st_15  = next_15;
                    st_nv  = next_nv;
                }

                double sum_A0 = compute_sum_state(D_A0, D_A0_sq);
                double sum_A1 = compute_sum_state(D_A1, D_A1_sq);
                double sum_A2 = compute_sum_state(D_A2, D_A2_sq);
                double sum_A3 = compute_sum_state(D_A3, D_A3_sq);
                double sum_nm = compute_sum_state(D_nm, D_nm_sq);

                double computed_I_T = std::exp(-I_T - 0.5 * I_T_inner) - 1.0;
                double abs_computed_I_T = std::fabs(computed_I_T);

                S     += sgn(sum_nm) - sgn(sum_A0);
                Sm    += sgn(sum_nm) - sgn(sum_A1);
                S_1_5 += sgn(sum_nm) - sgn(sum_A2);
                S_nv  += sgn(sum_nm) - sgn(sum_A3);
                S_lim += abs_computed_I_T;

                B     += (sgn(sum_A0) - sgn(sum_nm)) * (sgn(sum_A0) - sgn(sum_nm));
                Bm    += (sgn(sum_A1) - sgn(sum_nm)) * (sgn(sum_A1) - sgn(sum_nm));
                B_1_5 += (sgn(sum_A2) - sgn(sum_nm)) * (sgn(sum_A2) - sgn(sum_nm));
                B_nv  += (sgn(sum_A3) - sgn(sum_nm)) * (sgn(sum_A3) - sgn(sum_nm));
                B_lim += abs_computed_I_T * abs_computed_I_T;
            }
        }// End of parallel region

        const double inv_paths = 1.0 / static_cast<double>(paths);

        A[n]     = static_cast<double>(S * inv_paths);
        Am[n]    = static_cast<double>(Sm * inv_paths);
        A_1_5[n] = static_cast<double>(S_1_5 * inv_paths);
        A_nv[n]  = static_cast<double>(S_nv * inv_paths);
        A_lim[n] = static_cast<double>(S_lim * inv_paths);

        E[n]     = std::max(0.0, static_cast<double>(B * inv_paths - A[n] * A[n]));
        Em[n]    = std::max(0.0, static_cast<double>(Bm * inv_paths - Am[n] * Am[n]));
        E_1_5[n] = std::max(0.0, static_cast<double>(B_1_5 * inv_paths - A_1_5[n] * A_1_5[n]));
        E_nv[n]  = std::max(0.0, static_cast<double>(B_nv * inv_paths - A_nv[n] * A_nv[n]));
        E_lim[n] = std::max(0.0, static_cast<double>(B_lim * inv_paths - A_lim[n] * A_lim[n]));

        cout << "-------------------------------------------------" << n << "\n";
        cout << setprecision(10) << "points = " << points << "\n";
        cout << "-------------------------------------------------\n";
        cout << setprecision(15) << "Var EM        = " << E[n] << "\n";
        cout << setprecision(15) << "Var Milstein  = " << Em[n] << "\n";
        cout << setprecision(15) << "Var 1.5       = " << E_1_5[n] << "\n";
        cout << setprecision(15) << "Var NV        = " << E_nv[n] << "\n";
        cout << setprecision(15) << "Var Limit     = " << E_lim[n] << "\n";
        cout << setprecision(15) << "Mean EM       = " << A[n] << "\n";
        cout << setprecision(15) << "Mean Milstein = " << Am[n] << "\n";
        cout << setprecision(15) << "Mean 1.5      = " << A_1_5[n] << "\n";
        cout << setprecision(15) << "Mean NV       = " << A_nv[n] << "\n";
        cout << setprecision(15) << "Mean Limit    = " << A_lim[n] << "\n";

        ofs << n << "," << points << ","
            << fixed << setprecision(15)
            << E[n] << "," << Em[n] << "," << E_1_5[n] << "," << E_nv[n] << "," << E_lim[n] << ","
            << A[n] << "," << Am[n] << "," << A_1_5[n] << "," << A_nv[n] << "," << A_lim[n] << "\n";
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}