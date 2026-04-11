// 2DD2SSM3_all_test.cpp
// Scheme: Projected 1.5 Ito-Taylor (X^{2,*}), Delta^{2,*}, Sigma_{2,*}

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

static inline double var_nonneg_from_sums(long double sum, long double sumsq, long double N) {
    const long double mean  = sum   / N;
    const long double mean2 = sumsq / N;
    long double var = std::fma(-mean, mean, mean2);
    if (var < 0.0L) var = 0.0L;
    return static_cast<double>(var);
}

constexpr int sgn(double x) {
    if (isnan(x)) return 0;
    return (x > 0) - (x < 0);
}

// ========================================================================
// Model: s(x), a(x) and their analytical derivatives
// ========================================================================
namespace Model {

    inline double s(double x)   { return 1.0 + 0.8 * sin(x); }
    inline double ds(double x)  { return 0.8 * cos(x); }
    inline double dds(double x) { return -0.8 * sin(x); }

    inline State a(const State& x) {
        return State(2.0 + sin(x(1)), 1.0 + cos(x(0)));
    }

    inline double da1_dx1(const State& x) { (void)x; return 0.0; }
    inline double da1_dx2(const State& x) { return  cos(x(1)); }
    inline double da2_dx1(const State& x) { return -sin(x(0)); }
    inline double da2_dx2(const State& x) { (void)x; return 0.0; }


}

// ========================================================================
// Sigma_{2,*}(x)  — limit coefficient for kappa_{2,*}
//
// Sigma_{2,*}(x) = (29/288) * ( (s(x2)*s'(x1))^2/s(x1)^2
//                                + (s(x1)*s'(x2))^2/s(x2)^2 )^2
// ========================================================================
inline double Sigma2star(const State& x) {
    const double s1  = Model::s(x(0));
    const double s2  = Model::s(x(1));
    const double ds1 = Model::ds(x(0));
    const double ds2 = Model::ds(x(1));

    const double A = (s2 * ds1) * (s2 * ds1) / (s1 * s1);   // (s(x2)*s'(x1))^2 / s(x1)^2
    const double B = (s1 * ds2) * (s1 * ds2) / (s2 * s2);   // (s(x1)*s'(x2))^2 / s(x2)^2

    return (29.0 / 288.0) * (A + B) * (A + B);
}

// ========================================================================
// Hermite Polynomial Formulas
// ========================================================================
inline double H11(double t, double r1, double S2)
    { return r1*r1 / (t*t*S2*S2) - 1.0/(t*S2); }

inline double H22(double t, double r2, double S1)
    { return r2*r2 / (t*t*S1*S1) - 1.0/(t*S1); }

inline double H12(double t, double r1, double r2, double S1, double S2)
    { return r1*r2 / (t*t*S2*S1); }

inline double H1111(double t, double r1, double S2) {
    double t2=t*t, t3=t2*t, t4=t2*t2;
    double S2_2=S2*S2, S2_3=S2_2*S2, S2_4=S2_2*S2_2;
    return pow(r1,4)/(t4*S2_4) - 6.0*r1*r1/(t3*S2_3) + 3.0/(t2*S2_2);
}

inline double H2222(double t, double r2, double S1) {
    double t2=t*t, t3=t2*t, t4=t2*t2;
    double S1_2=S1*S1, S1_3=S1_2*S1, S1_4=S1_2*S1_2;
    return pow(r2,4)/(t4*S1_4) - 6.0*r2*r2/(t3*S1_3) + 3.0/(t2*S1_2);
}

inline double H1112(double t, double r1, double r2, double S1, double S2) {
    double t3=t*t*t, t4=t3*t;
    return pow(r1,3)*r2 / (t4*pow(S2,3)*S1) - 3.0*r1*r2 / (t3*S2*S2*S1);
}

inline double H1222(double t, double r1, double r2, double S1, double S2) {
    double t3=t*t*t, t4=t3*t;
    return r1*pow(r2,3) / (t4*S2*pow(S1,3)) - 3.0*r1*r2 / (t3*S2*S1*S1);
}

inline double H1122(double t, double r1, double r2, double S1, double S2) {
    return H11(t, r1, S2) * H22(t, r2, S1);
}

// ========================================================================
// Delta^{2,*}_t(x, y) =
//     + (t^2/8) * (s(x1)*s'(x2))^2 * H^{11}
//     + (t^2/8) * (s(x2)*s'(x1))^2 * H^{22}
//     - (t^2/4) * s(x1)*s(x2)*s'(x1)*s'(x2)* H^{12}
//     + (t^3/24)*(s(x1)*s'(x2))^2 *s(x2)^2 * H^{1111}
//     + (t^3/24)*(s(x2)*s'(x1))^2 * s(x1)^2 * H^{2222}
//     - (t^3/12)* s(x1)*s(x2)^3 *s'(x1)*s'(x2) * H^{1112}
//     - (t^3/12)* s(x1)^3 *s(x2)*s'(x1)*s'(x2) * H^{1222}
//     + (t^3/12)*((s(x2)^4 * s'(x1))^2 + s(x1)^4 * s'(x2)^2) * H^{1122}
// ========================================================================
inline double Delta2star(const State& x, const State& y, double t) {
    const double s1  = Model::s(x(0));
    const double s2  = Model::s(x(1));
    const double ds1 = Model::ds(x(0));
    const double ds2 = Model::ds(x(1));

    const State a_val = Model::a(x);

    // Residuals after removing drift (needed for Hermite polynomials)
    const double r1 = y(0) - x(0) - a_val(0)*t;
    const double r2 = y(1) - x(1) - a_val(1)*t;

    // Variance factors S1 = s(x1)^2, S2 = s(x2)^2
    const double S1 = s1 * s1;
    const double S2 = s2 * s2;

    // Coefficient combinations
    const double C11 = (s1*ds2)*(s1*ds2) ;   // (s(x1)*s'(x2))^2 
    const double C22 = (s2*ds1)*(s2*ds1);   // (s(x2)*s'(x1))^2 
    const double Cds = ds1 * ds2;
    const double C12 = s1 * s2 * Cds;              // s(x1)*s(x2)*s'(x1)*s'(x2)

    const double t2 = t * t;
    const double t3 = t2 * t;

    double result = 0.0;

    // t^2 terms
    result += (t2/8.0) * C11 * H11(t, r1, S2);
    result += (t2/8.0) * C22 * H22(t, r2, S1);
    result -= (t2/4.0) * C12 * H12(t, r1, r2, S1, S2);

    // t^3 terms
    result += (t3/24.0) * C11 * S2 * H1111(t, r1, S2);
    result += (t3/24.0) * C22 * S1 * H2222(t, r2, S1);
    result -= (t3/12.0) * C12 * S2 * H1112(t, r1, r2, S1, S2);
    result -= (t3/12.0) * C12 * S1 * H1222(t, r1, r2, S1, S2);
    result += (t3/12.0) * (C22 * S2 + C11 * S1) * H1122(t, r1, r2, S1, S2);

    return result;
}

// ========================================================================
// Simulation Schemes
// ========================================================================

// Euler-Maruyama (Order 0.5)
inline State A0(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    Vector2d drift = Model::a(curr);
    Vector2d diffusion(Model::s(curr(1)) * dW(0), Model::s(curr(0)) * dW(1));
    return curr + drift * dt + diffusion;
}

// Milstein (Order 1.0)
inline State A1(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);

    double s_x0  = Model::s(curr(0));
    double s_x1  = Model::s(curr(1));
    double ds_x0 = Model::ds(curr(0));
    double ds_x1 = Model::ds(curr(1));

    Vector2d drift     = Model::a(curr);
    Vector2d diffusion(s_x1 * dW(0), s_x0 * dW(1));
    State base = curr + drift * dt + diffusion;

    // Milstein cross-correction: (1/2) s'(x2)*s(x1)*w1*w2 and (1/2) s'(x1)*s(x2)*w1*w2
    Vector2d milstein_corr(
        0.5 * ds_x1 * s_x0 * dW(0) * dW(1),
        0.5 * ds_x0 * s_x1 * dW(0) * dW(1)
    );
    return base + milstein_corr;
}

// Strong Order 1.5 Scheme
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double w1 = dW(0), w2 = dW(1);
    double w1_sq = w1 * w1;
    double w2_sq = w2 * w2;

    double s0 = Model::s(curr(0));
    double s1 = Model::s(curr(1));
    double ds0 = Model::ds(curr(0));
    double ds1 = Model::ds(curr(1));
    double dds0 = Model::dds(curr(0));
    double dds1 = Model::dds(curr(1));
    
    Vector2d a_val = Model::a(curr);
    
    // Decoupled drift geometry mapping via analytical derivatives
    double da1_dx1 = Model::da1_dx1(curr);
    double da1_dx2 = Model::da1_dx2(curr);
    double da2_dx1 = Model::da2_dx1(curr);
    double da2_dx2 = Model::da2_dx2(curr);

    double w111 = (w1 * w1_sq) - 3.0 * dt * w1;
    double w222 = (w2 * w2_sq) - 3.0 * dt * w2;
    double w112 = (w1_sq - dt) * w2;
    double w122 = (w2_sq - dt) * w1;

    State next = curr;

    // --- X1 Component ---
    next(0) += a_val(0) * dt + s1 * w1; 
    next(0) += 0.5 * ds1 * s0 * w1 * w2; 
    next(0) += 0.5 * dt * (s1 * da1_dx1 * w1 + s0 * da1_dx2 * w2 + ds1 * a_val(1) * w1); 
    
    double term_t8_1 = 2.0 * dds1 * s0 * s0 * w1;
    double term_t8_2 = (ds1 * s0 * ds1 * s0 / s1) * w1;
    double term_t8_3 = -ds0 * ds1 * s1 * w2;
    next(0) += 0.125 * dt * (term_t8_1 + term_t8_2 + term_t8_3);

    double term_h1 = (1.0/24.0) * (ds1*ds1 / s1) * s0*s0 * w111;
    double term_h2 = (1.0/12.0) * ds0 * ds1 * s1 * w112;
    double term_h3 = (1.0/6.0)  * dds1 * s0*s0 * w122;
    double term_h4 = (1.0/24.0) * (ds1*ds1 / s1) * s0*s0 * w122;
    next(0) += term_h1 + term_h2 + term_h3 + term_h4;

    // --- X2 Component ---
    next(1) += a_val(1) * dt + s0 * w2; 
    next(1) += 0.5 * ds0 * s1 * w1 * w2; 
    next(1) += 0.5 * dt * (s1 * da2_dx1 * w1 + s0 * da2_dx2 * w2 + ds0 * a_val(0) * w2); 

    double term2_t8_1 = 2.0 * dds0 * s1 * s1 * w2;
    double term2_t8_2 = ds0 * s1 * ds0 * s1 * w2 / s0;
    double term2_t8_3 = -ds1 * ds0 * s0 * w1;
    next(1) += 0.125 * dt * (term2_t8_1 + term2_t8_2 + term2_t8_3);

    double term2_h1 = (1.0/24.0) * (ds0*ds0 / s0) * s1*s1 * w222;
    double term2_h2 = (1.0/12.0) * ds0 * ds1 * s0 * w122;
    double term2_h3 = (1.0/6.0)  * dds0 * s1*s1 * w112;
    double term2_h4 = (1.0/24.0) * (ds0*ds0 / s0) * s1*s1 * w112;
    next(1) += term2_h1 + term2_h2 + term2_h3 + term2_h4;

    return next;
}

// ========================================================================
// Projected 1.5 Ito-Taylor scheme  X^{2,*}  
// ========================================================================
inline State A2star(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    const double w1 = sqrt_dt * Z1;
    const double w2 = sqrt_dt * Z2;

    const double s1   = Model::s(curr(0));
    const double s2   = Model::s(curr(1));
    const double ds1  = Model::ds(curr(0));
    const double ds2  = Model::ds(curr(1));
    const double dds1 = Model::dds(curr(0));
    const double dds2 = Model::dds(curr(1));

    const State  a_val  = Model::a(curr);
    const double a1     = a_val(0);
    const double a2     = a_val(1);

    const double d1a1 = Model::da1_dx1(curr);
    const double d2a1 = Model::da1_dx2(curr);
    const double d1a2 = Model::da2_dx1(curr);
    const double d2a2 = Model::da2_dx2(curr);

    // Third-order Hermite monomials
    const double w112 = (w1*w1 - dt) * w2;   // (w1^2 - t)*w2
    const double w122 = (w2*w2 - dt) * w1;   // (w2^2 - t)*w1

    State next = curr;

    // --- Component 1 (X1) ---
    // Euler
    next(0) += a1 * dt + s2 * w1;
    // Milstein cross
    next(0) += 0.5 * ds2 * s1 * w1 * w2;
    // t/2 drift-diffusion interaction
    next(0) += (dt/2.0) * (a2*ds2*w1 + d1a1*s2*w1 + d2a1*s1*w2);
    // t/4 second-derivative term
    next(0) += (dt/4.0) * (s1*s1 * dds2 * w1);
    // 1/6 cubic Hermite terms
    next(0) += (1.0/6.0) * (s2*ds1*ds2*w112 + s1*s1*dds2*w122);

    // --- Component 2 (X2) ---
    // Euler
    next(1) += a2 * dt + s1 * w2;
    // Milstein cross
    next(1) += 0.5 * ds1 * s2 * w1 * w2;
    // t/2 drift-diffusion interaction
    next(1) += (dt/2.0) * (a1*ds1*w2 + d1a2*s2*w1 + d2a2*s1*w2);
    // t/4 second-derivative term
    next(1) += (dt/4.0) * (s2*s2 * dds1 * w2);
    // 1/6 cubic Hermite terms
    next(1) += (1.0/6.0) * (dds1*s2*s2*w112 + s1*ds1*ds2*w122);

    return next;
}

// ========================================================================
// Main Program
// ========================================================================
int main() {

    constexpr double t_start = 0.0;
    constexpr double t_end   = 0.5;
    constexpr double mu      = 0.0;
    constexpr double sigma   = 1.0;

    const State x0_state = Vector2d(0.5, 0.5);
    constexpr int max_n = 9;

    vector<double> A(max_n+1,0.0), Am(max_n+1,0.0), A_2s(max_n+1,0.0), A_lim(max_n+1,0.0), A_1_5(max_n+1,0.0);
    vector<double> E(max_n+1,0.0), Em(max_n+1,0.0), E_2s(max_n+1,0.0), E_lim(max_n+1,0.0), E_1_5(max_n+1,0.0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str());
    const string csv_path = dir_path + "/2DD2SSM3_all_100_1000_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);

    if (!ofs) {
        cerr << "Failed to open CSV file: " << csv_path << endl;
        return 1;
    }

    ofs.imbue(locale::classic());
    ofs << "n,points,E_EM,E_Mil,E_2star,E_lim,E_1_5,A_EM,A_Mil,A_2star,A_lim,A_1_5\n";

    for (int n = 0; n <= max_n; ++n) {
        const int points = 100 + 100 * n;
        const int paths  = 10 * points * points;

        const double dt  = (t_end - t_start) / (points - 1);
        const double dtm = dt / points;

        long double S=0,Sm=0,S_2s=0,S_lim=0,S_1_5=0;
        long double B=0,Bm=0,B_2s=0,B_lim=0,B_1_5=0;

        #pragma omp parallel reduction(+:S,Sm,S_2s,S_lim,S_1_5,B,Bm,B_2s,B_lim,B_1_5)
        {
            normal_distribution<double> dist(mu, sigma);

            #pragma omp for schedule(static) nowait
            for (int p = 0; p < paths; ++p) {
                seed_seq ss0{40u, 0u, (uint32_t)p};
                seed_seq ss1{50u, 1u, (uint32_t)p};
                mt19937 rng_nm(ss0);
                mt19937 rng1_nm(ss1);

                State st_em  = x0_state;
                State st_mil = x0_state;
                State st_1_5 = x0_state;
                State st_2s  = x0_state;   // X^{2,*} scheme
                State st_nm  = x0_state;
                State st_nm_y = x0_state;

                double D_A0=0, D_A1=0, D_A2=0, D_A2s=0, D_bench=0;
                double Integrated_Sigma2s = 0.0;   // integrates Sigma_{2,*}

                for (int idx = 1; idx < points; ++idx) {
                    double Z1 = 0.0, Z2 = 0.0;
                    st_nm = st_nm_y;

                    for (int m = 0; m < points; ++m) {
                        double Z1_nm = dist(rng_nm);
                        double Z2_nm = dist(rng1_nm);

                        // Accumulate Sigma_{2,*} integral
                        Integrated_Sigma2s += Sigma2star(st_nm) * dtm;

                        //iteration of X_t
                        st_nm = A1(st_nm, dtm, Z1_nm, Z2_nm);

                        Z1 += Z1_nm / sqrt(points);
                        Z2 += Z2_nm / sqrt(points);
                    }
                    State nm_bench_y = st_nm;

                    State next_em  = A0(st_em,  dt, Z1, Z2);
                    State next_mil = A1(st_mil, dt, Z1, Z2);
                    State next_1_5 = A2(st_1_5, dt, Z1, Z2);
                    State next_2s  = A2star(st_2s, dt, Z1, Z2);

                    // Accumulate Delta^{2,*} test statistics
                    D_A0   += Delta2star(st_em,   next_em,   dt);
                    D_A1   += Delta2star(st_mil,  next_mil,  dt);
                    D_A2   += Delta2star(st_1_5,  next_1_5,  dt);
                    D_A2s  += Delta2star(st_2s,   next_2s,   dt);
                    D_bench+= Delta2star(st_nm_y, nm_bench_y,dt);

                    st_em   = next_em;
                    st_mil  = next_mil;
                    st_1_5  = next_1_5;
                    st_2s   = next_2s;
                    st_nm_y = nm_bench_y;
                }

                // F^{2,*}_n statistics: h^{-1/2} * sgn(sum Delta^{2,*})
                S    += (long double)(sgn(D_bench) - sgn(D_A0))  / sqrt(dt);
                Sm   += (long double)(sgn(D_bench) - sgn(D_A1))  / sqrt(dt);
                S_1_5 += (long double)(sgn(D_bench) - sgn(D_A2))  / sqrt(dt);
                S_2s += (long double)(sgn(D_bench) - sgn(D_A2s)) / sqrt(dt);


                B    += (long double)(sgn(D_A0)  - sgn(D_bench)) * (sgn(D_A0)  - sgn(D_bench)) / dt;
                Bm   += (long double)(sgn(D_A1)  - sgn(D_bench)) * (sgn(D_A1)  - sgn(D_bench)) / dt;
                B_1_5 += (long double)(sgn(D_A2)  - sgn(D_bench)) * (sgn(D_A2)  - sgn(D_bench)) / dt;
                B_2s += (long double)(sgn(D_A2s) - sgn(D_bench)) * (sgn(D_A2s) - sgn(D_bench)) / dt;

                // kappa_{2,*} sample: E[ sqrt(2/pi * integral Sigma_{2,*} dt) ]
                const double PI = 3.14159265358979323846;
                double inside = 2.0 * Integrated_Sigma2s / PI;
                if (inside < 0.0) inside = 0.0;
                double kappa_sample = sqrt(inside);
                S_lim += (long double)kappa_sample;
                B_lim += (long double)kappa_sample * (long double)kappa_sample;
            }
        } // end parallel

        const long double N = (long double)paths;

        A[n]    = static_cast<double>(S    / N);
        Am[n]   = static_cast<double>(Sm   / N);
        A_1_5[n] = static_cast<double>(S_1_5 / N);
        A_2s[n] = static_cast<double>(S_2s / N);
        A_lim[n]= static_cast<double>(S_lim/ N);

        E[n]    = var_nonneg_from_sums(S,    B,    N);
        Em[n]   = var_nonneg_from_sums(Sm,   Bm,   N);
        E_1_5[n] = var_nonneg_from_sums(S_1_5, B_1_5, N);
        E_2s[n] = var_nonneg_from_sums(S_2s, B_2s, N);
        E_lim[n]= var_nonneg_from_sums(S_lim,B_lim,N);

        cout << "------------------------------------------------- n=" << n << "\n";
        cout << setprecision(10) << "points = " << points << "\n";
        cout << "-------------------------------------------------\n";
        cout << setprecision(15) << "Var EM        = " << E[n]    << "\n";
        cout << setprecision(15) << "Var Milstein  = " << Em[n]   << "\n";
        cout << setprecision(15) << "Var A2        = " << E_1_5[n] << "\n";
        cout << setprecision(15) << "Var A2_star   = " << E_2s[n] << "\n";
        cout << setprecision(15) << "Var Limit     = " << E_lim[n]<< "\n";
        cout << setprecision(15) << "Mean EM       = " << A[n]    << "\n";
        cout << setprecision(15) << "Mean Milstein = " << Am[n]   << "\n";
        cout << setprecision(15) << "Mean A2       = " << A_1_5[n] << "\n";
        cout << setprecision(15) << "Mean A2_star  = " << A_2s[n] << "\n";
        cout << setprecision(15) << "Mean Limit    = " << A_lim[n]<< "\n";
        

        ofs << n << "," << points << ","
            << fixed << setprecision(15)
            << E[n]    << "," << Em[n]   << "," << E_2s[n] << "," << E_lim[n] << "," << E_1_5[n] << ","
            << A[n]    << "," << Am[n]   << "," << A_2s[n] << "," << A_lim[n] << "," << A_1_5[n] << endl;
    }

    ofs.close();
    cout << "CSV written to: " << csv_path << endl;
    return 0;
}