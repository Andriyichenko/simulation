// 2DD2SM3_lim.cpp
// Limit Simulation I_T^1 only (OLD one)

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
inline double ds_func(double x) { return cos(x); }
inline double dds_func(double x) { return -sin(x); }

// Milstein (Order 1.0) — used to generate the benchmark path for the limit integrand
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
    constexpr double t_end = 1.0;
    constexpr double mu = 0.0;
    constexpr double sigma = 1.0;

    const State x0_state = Vector2d(1, 1);
    constexpr int max_n = 9;

    vector<double> A_lim(max_n + 1, 0.0);
    vector<double> E_lim(max_n + 1, 0.0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str());
    const string csv_path = dir_path + "/2DD2SM3_lim2_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);

    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }

    ofs.imbue(locale::classic());
    ofs << "n,points,E_lim,A_lim\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n;
        const int paths = 10 * points * points;

        const double dt = (t_end - t_start) / (points - 1);
        const double dtm = dt / (points - 1);
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
                State st_nm_y = x0_state;

                double Z_limit = 0.0; // I_T^1 (accumulator)

                for (int idx = 1; idx < points; ++idx) {
                    st_nm = st_nm_y;

                    for (int m = 0; m < points; ++m) {
                        double Z1_nm = dist(rng_nm);
                        double Z2_nm = dist(rng1_nm);

                        // independent noise dW~
                        double Z_tilde = dist(rng_lim);

                        // coefficients for integrand
                        double x1 = st_nm(0);
                        double x2 = st_nm(1);

                        double s1 = s_func(x1);
                        double s2 = s_func(x2);
                        double ds1 = ds_func(x1);
                        double ds2 = ds_func(x2);
                        double dds1 = dds_func(x1);
                        double dds2 = dds_func(x2);

                        // drift a(x) = (x2, -x1)
                        double a1_val = x2;
                        double a2_val = -x1;

                        // m=2 coefficients
                        double c2_12 = 0.5 * (1.0 / s2) + 0.5 * (-1.0 / s1) - 0.25 * ds1 * ds2;

                        double term11_2 = 0.25 * (ds2 * a2_val) / (s2 * s2);
                        double term11_3 = 0.25 * (dds2 * s1 * s1) / s2;
                        double term11_4 = 0.125 * (ds2 * s1) * (ds2 * s1) / (s2 * s2);
                        double c2_11 = term11_2 + term11_3 + term11_4;

                        double term22_2 = 0.25 * (ds1 * a1_val) / (s1 * s1);
                        double term22_3 = 0.25 * (dds1 * s2 * s2) / s1;
                        double term22_4 = 0.125 * (ds1 * s2) * (ds1 * s2) / (s1 * s1);
                        double c2_22 = term22_2 + term22_3 + term22_4;

                        //double Sum_m2 = 2.0 * (pow(c2_11, 2) + pow(c2_22, 2) + 2.0 * pow(c2_12, 2));
                        double Sum_m2 = 2.0 * (pow(c2_11, 2) + pow(c2_22, 2) + 0.5 * pow(c2_12, 2));

                        // m=4 coefficients
                        double c2_1111 = (1.0 / 24.0) * (pow(ds2, 2) * pow(s1, 2)) / pow(s2, 2);
                        double c2_2222 = (1.0 / 24.0) * (pow(ds1, 2) * pow(s2, 2)) / pow(s1, 2);

                        double c2_1112 = (1.0 / 12.0) * ds1 * ds2;
                        double c2_1222 = c2_1112;

                        double term1122_1 = (1.0 / 6.0) * dds1 * (pow(s2, 2) / s1);
                        double term1122_2 = (1.0 / 6.0) * dds2 * (pow(s1, 2) / s2);
                        double term1122_3 = (1.0 / 24.0) * pow(ds2, 2) * (pow(s1, 2) / pow(s2, 2));
                        double term1122_4 = (1.0 / 24.0) * pow(ds1, 2) * (pow(s2, 2) / pow(s1, 2));
                        double c2_1122 = term1122_1 + term1122_2 + term1122_3 + term1122_4;

                        // double sum_c4_sq = pow(c2_1111, 2) + pow(c2_2222, 2)
                        //                  + 4.0 * pow(c2_1112, 2)
                        //                  + 4.0 * pow(c2_1222, 2)
                        //                  + 6.0 * pow(c2_1122, 2);

                        double sum_c4_sq = 24*pow(c2_1111, 2) + 24*pow(c2_2222, 2)
                                         + 6.0 * pow(c2_1112, 2)
                                         + 6.0 * pow(c2_1222, 2)
                                         + 4.0 * pow(c2_1122, 2);

                        double Sum_m4 = 24.0 * sum_c4_sq;

                        double integrand = sqrt(Sum_m2 + Sum_m4);

                        // accumulate limit integral
                        Z_limit += integrand * Z_tilde * sqrt_dtm;

                        // update benchmark state
                        st_nm = A1(st_nm, dtm, Z1_nm, Z2_nm);
                    }

                    st_nm_y = st_nm;
                }

                S_lim += fabs(Z_limit);
                B_lim += Z_limit * Z_limit;
            }
        }

        const double inv_paths = 1.0 / paths;
        A_lim[n] = S_lim * inv_paths;
        E_lim[n] = B_lim * inv_paths - A_lim[n] * A_lim[n];

        cout << "-------------------------------------------------" << n << "\n";
        cout << setprecision(10) << "points = " << points << "\n";
        cout << "-------------------------------------------------\n";
        cout << setprecision(15) << "Var Limit  = " << E_lim[n] << "\n";
        cout << setprecision(15) << "Mean Limit = " << A_lim[n] << "\n";

        ofs << n << "," << points << ","
            << fixed << setprecision(15)
            << E_lim[n] << "," << A_lim[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}