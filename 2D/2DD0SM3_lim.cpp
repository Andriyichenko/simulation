// 2DD1SM3_lim.cpp
// Limit Simulation only (I_T^0)
// s(x) = 2 + sin(x)
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <omp.h>
#include <random>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

using State = Vector2d;

// s(x) = 2 + sin(x)
inline double s_func(double x) { return 2.0 + sin(x); }
// s'(x) = cos(x)
inline double ds_func(double x) { return cos(x); }

// Milstein (Order 1.0) used as benchmark state evolution for the limit integrand
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

int main() {
    constexpr double t_start = 0.0;
    constexpr double t_end   = 1.0;
    constexpr double mu      = 0.0;
    constexpr double sigma   = 1.0;

    const State x0_state = Vector2d(1, 1); // X_start
    constexpr int max_n = 9;

    // Limit Simulation stats only
    vector<double> A_lim(max_n + 1, 0.0);
    vector<double> E_lim(max_n + 1, 0.0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str());
    const string csv_path = dir_path + "/2DD1SM3_lim_100_1000_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);

    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }

    ofs.imbue(locale::classic());
    ofs << "n,points,E_lim,A_lim\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n;
        const int paths  = 10 * points * points;

        const double dt   = (t_end - t_start) / (points - 1);
        const double dtm  = dt / (points - 1);
        const double sqrt_dtm = sqrt(dtm);

        double S_lim = 0.0;
        double B_lim = 0.0;

        #pragma omp parallel reduction(+:S_lim,B_lim)
        {
            normal_distribution<double> dist(mu, sigma);

            #pragma omp for schedule(static) nowait
            for (int p = 0; p < paths; ++p) {
                seed_seq ss0{40u, 0u, (uint32_t)p};
                seed_seq ss1{50u, 1u, (uint32_t)p};
                seed_seq ss2{1000u, 2u, (uint32_t)p};
                mt19937 rng_nm(ss0);
                mt19937 rng1_nm(ss1);
                mt19937 rng_lim(ss2);
                State st_nm = x0_state;

                // Limit Simulation variable
                double I_T = 0.0;

                for (int idx = 1; idx < points; ++idx) {
                    for (int m = 0; m < points; ++m) {
                        // Noise for benchmark state evolution
                        double Z1_nm = dist(rng_nm);
                        double Z2_nm = dist(rng1_nm);

                        // Independent noise for Limit Simulation (scalar W_tilde)
                        double Z_tilde = dist(rng_lim);

                        // Compute integrand based on current benchmark state st_nm
                        double x1 = st_nm(0);
                        double x2 = st_nm(1);
                        double s1 = s_func(x1);
                        double s2 = s_func(x2);
                        double ds1 = ds_func(x1);
                        double ds2 = ds_func(x2);

                        double term1 = (ds1 / s1) * s2;
                        double term2 = (ds2 / s2) * s1;
                        double integrand = sqrt(term1 * term1 + term2 * term2);

                        // Update I_T
                        I_T += sqrt(1.5) * integrand * Z_tilde * sqrt_dtm;

                        // Update benchmark state
                        st_nm = A1(st_nm, dtm, Z1_nm, Z2_nm);
                    }
                }

                S_lim += I_T;
                B_lim += I_T * I_T;
            }
        } // End parallel region

        const double inv_paths = 1.0 / paths;
        A_lim[n] = S_lim * inv_paths;
        E_lim[n] = B_lim * inv_paths - A_lim[n] * A_lim[n];

        cout << "-------------------------------------------------" << n << "\n";
        cout << setprecision(10) << "points = " << points << "\n";
        cout << "-------------------------------------------------\n";
        cout << setprecision(15) << "Var Lim  = " << E_lim[n] << "\n";
        cout << setprecision(15) << "Mean Lim = " << A_lim[n] << "\n";

        ofs << n << "," << points << ","
            << fixed << setprecision(15)
            << E_lim[n] << "," << A_lim[n] << "\n";
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}