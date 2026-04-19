// 2DD1SNVM3_all.cpp
// Integrated with Euler-Maruyama (A0), Milstein (A1), Order 1.5 (A2), and Ninomiya-Victoir (A3)

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

using State = Vector2d; 

// ========================================================================
// Math & Precision Helpers
// ========================================================================

// Safely compute variance to avoid catastrophic cancellation and negative values
static inline double var_nonneg_from_sums(long double sum, long double sumsq, long double N) {
    const long double mean  = sum   / N;   
    const long double mean2 = sumsq / N;   
    long double var = std::fma(-mean, mean, mean2); // Equivalently: mean2 - mean*mean
    if (var < 0.0L) var = 0.0L;                     // Force non-negative constraint
    return static_cast<double>(var);
}

// ========================================================================
// CONFIGURATION s(x), a(x) Analytical Functions
// ========================================================================

inline double s_func(double x)   { return 2.0 + sin(x); }
inline double ds_func(double x)  { return cos(x); }
inline double dds_func(double x) { return -sin(x); }

// ========================================================================
// Simulation Schemes
// ========================================================================

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
    double s0 = s_func(curr(0)), s1 = s_func(curr(1));
    double ds0 = ds_func(curr(0)), ds1 = ds_func(curr(1));

    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s1 * dW(0), s0 * dW(1));
    State base = curr + drift * dt + diffusion;

    Vector2d milstein_corr(
        0.5 * ds1 * s0 * dW(0) * dW(1),
        0.5 * ds0 * s1 * dW(0) * dW(1)
    );
    return base + milstein_corr;
}

// Strong Order 1.5 Scheme
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double w1 = dW(0), w2 = dW(1);

    double s0 = s_func(curr(0)), s1 = s_func(curr(1));
    double ds0 = ds_func(curr(0)), ds1 = ds_func(curr(1));
    double dds0 = dds_func(curr(0)), dds1 = dds_func(curr(1));
    Vector2d a(curr(1), -curr(0));

    double w111 = (w1 * w1 * w1) - 3.0 * dt * w1;
    double w222 = (w2 * w2 * w2) - 3.0 * dt * w2;
    double w112 = (w1 * w1 - dt) * w2;
    double w122 = (w2 * w2 - dt) * w1;

    State next = curr;
    // X1 component
    next(0) += a(0) * dt + s1 * w1 + 0.5 * ds1 * s0 * w1 * w2;
    next(0) += 0.5 * dt * (s0 * w2 + ds1 * a(1) * w1);
    next(0) += 0.125 * dt * (2.0 * dds1 * s0 * s0 * w1 + (ds1 * s0 * ds1 * s0 / s1) * w1 - ds0 * ds1 * s1 * w2);
    next(0) += (1.0/24.0)*(ds1*ds1/s1)*s0*s0*w111 + (1.0/12.0)*ds0*ds1*s1*w112 + (1.0/6.0)*dds1*s0*s0*w122 + (1.0/24.0)*(ds1*ds1/s1)*s0*s0*w122;

    // X2 component
    next(1) += a(1) * dt + s0 * w2 + 0.5 * ds0 * s1 * w1 * w2;
    next(1) += 0.5 * dt * (-s1 * w1 + ds0 * a(0) * w2);
    next(1) += 0.125 * dt * (2.0 * dds0 * s1 * s1 * w2 + (ds0 * s1 * ds0 * s1 / s0) * w2 - ds1 * ds0 * s0 * w1);
    next(1) += (1.0/24.0)*(ds0*ds0/s0)*s1*s1*w222 + (1.0/12.0)*ds0*ds1*s0*w122 + (1.0/6.0)*dds0*s1*s1*w112 + (1.0/24.0)*(ds0*ds0/s0)*s1*s1*w112;

    return next;
}

// Ninomiya-Victoir Scheme (A3)
inline State A3(const State& curr, double dt, double Z1, double Z2, int xi) {
    const double sqrt_dt = sqrt(dt);
    double dW1 = sqrt_dt * Z1;
    double dW2 = sqrt_dt * Z2;

    // V0: Exact drift flow (Rotation)
    auto flow0 = [](double t, const State& x) -> State {
        double c = cos(t); double s = sin(t);
        return State(x(0) * c + x(1) * s, -x(0) * s + x(1) * c);
    };
    // V1: Exact flow for sigma column 1
    auto flow1 = [](double z, const State& x) -> State {
        return State(x(0) + z * s_func(x(1)), x(1));
    };
    // V2: Exact flow for sigma column 2
    auto flow2 = [](double z, const State& x) -> State {
        return State(x(0), x(1) + z * s_func(x(0)));
    };

    State y = flow0(0.5 * dt, curr);
    if (xi > 0) {
        y = flow1(dW1, y);
        y = flow2(dW2, y);
    } else {
        y = flow2(dW2, y);
        y = flow1(dW1, y);
    }
    return flow0(0.5 * dt, y);
}

// ========================================
// Main Simulation
// ========================================

int main() {
    const double z_const = 1.38; 
    const double t_end = 1.0;
    const State x0_state(1.0, 1.0);
    const int max_n = 9;

    vector<double> A(max_n + 1), Am(max_n + 1), A_1_5(max_n + 1), A_nv(max_n + 1);
    vector<double> E(max_n + 1), Em(max_n + 1), E_1_5(max_n + 1), E_nv(max_n + 1);

    const string csv_path = "../data_source/2DD1SM3_all_data.csv";
    system("mkdir -p ../data_source");
    ofstream ofs(csv_path, ios::out | ios::trunc);
    ofs << "n,points,E,Em,E_1.5,E_nv,A,Am,A_1.5,A_nv\n";

    auto phi_n = [](double z, double x, double h) {
        const double t = h; // alpha = 1.0
        const double PI = 3.14159265358979323846;
        return (1.0 / sqrt(2.0 * PI * t)) * exp(-(x - z) * (x - z) / (2.0 * t));
    };

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n;
        const long double N = 10.0L * points * points;
        const double dt = t_end / (points - 1);
        const double dtm = dt / (points - 1);

        long double S = 0, Sm = 0, S_1_5 = 0, S_nv = 0;
        long double B = 0, Bm = 0, B_1_5 = 0, B_nv = 0;

        #pragma omp parallel reduction(+:S,Sm,S_1_5,S_nv,B,Bm,B_1_5,B_nv)
        {
            mt19937 rng_z1, rng_z2, rng_xi;
            normal_distribution<double> dist_z1(0, 1), dist_z2(0, 1);
            bernoulli_distribution dist_xi(0.5);

            #pragma omp for schedule(static)
            for (int p = 0; p < (int)N; ++p) {
                rng_z1.seed(30u + p); rng_z2.seed(42u + p); rng_xi.seed(54u + p);
                
                State st_em = x0_state, st_mil = x0_state, st_15 = x0_state, st_nv = x0_state, st_nm = x0_state;
                double L_em = 0, L_mil = 0, L_15 = 0, L_nv = 0, L_nm = 0;

                for (int idx = 1; idx < points; ++idx) {
                    double Z1 = 0, Z2 = 0;
                    for (int m = 0; m < points; ++m) {
                        double z1_sub = dist_z1(rng_z1), z2_sub = dist_z2(rng_z2);
                        st_nm = A1(st_nm, dtm, z1_sub, z2_sub);
                        Z1 += z1_sub / sqrt(points);
                        Z2 += z2_sub / sqrt(points);
                    }
                    
                    int xi = dist_xi(rng_xi) ? 1 : -1;

                    st_em  = A0(st_em, dt, Z1, Z2);
                    st_mil = A1(st_mil, dt, Z1, Z2);
                    st_15  = A2(st_15, dt, Z1, Z2);
                    st_nv  = A3(st_nv, dt, Z1, Z2, xi);

                    L_em  += dt * phi_n(z_const, st_em(0), dt);
                    L_mil += dt * phi_n(z_const, st_mil(0), dt);
                    L_15  += dt * phi_n(z_const, st_15(0), dt);
                    L_nv  += dt * phi_n(z_const, st_nv(0), dt);
                    L_nm  += dt * phi_n(z_const, st_nm(0), dt);
                }

                auto f = [](double L) { return atan(L); };
                double v_nm = f(L_nm), v_em = f(L_em), v_mil = f(L_mil), v_15 = f(L_15), v_nv = f(L_nv);

                S += (v_nm - v_em);    B += (long double)(v_nm - v_em) * (v_nm - v_em);
                Sm += (v_nm - v_mil);  Bm += (long double)(v_nm - v_mil) * (v_nm - v_mil);
                S_1_5 += (v_nm - v_15); B_1_5 += (long double)(v_nm - v_15) * (v_nm - v_15);
                S_nv += (v_nm - v_nv); B_nv += (long double)(v_nm - v_nv) * (v_nm - v_nv);
            }
        }

        A[n] = (double)(S / N); Am[n] = (double)(Sm / N); A_1_5[n] = (double)(S_1_5 / N); A_nv[n] = (double)(S_nv / N);
        E[n] = var_nonneg_from_sums(S, B, N);
        Em[n] = var_nonneg_from_sums(Sm, Bm, N);
        E_1_5[n] = var_nonneg_from_sums(S_1_5, B_1_5, N);
        E_nv[n] = var_nonneg_from_sums(S_nv, B_nv, N);

        cout << "n=" << n << " pts=" << points << " VarNV=" << scientific << E_nv[n] << " MeanNV=" << A_nv[n] << endl;
        ofs << n << "," << points << "," << fixed << setprecision(10) << E[n] << "," << Em[n] << "," << E_1_5[n] << "," << E_nv[n] 
            << "," << A[n] << "," << Am[n] << "," << A_1_5[n] << "," << A_nv[n] << endl;
    }
    return 0;
}