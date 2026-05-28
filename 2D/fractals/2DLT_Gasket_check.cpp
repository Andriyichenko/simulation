//2DLTM3_gasket
// Sierpinski Gasket local time simulation
//   d=2, q=3, p=2, d_F = log_2(3) ≈ 1.5850
//   Leb(F^(k)) = (sqrt(3)/4) * (3/4)^k
//   k_n = floor( (1+log2_3) / (log2_3*(2-log2_3)+2) * log2(n) )
//         ≈ 0.972 * log2(n)
//   sigma ≈ 0.596
//
// Approximation:
//   L_t^{mu_n,n} = 1/(n*Leb(F^(k_n))) * sum_{k=1}^{floor(nt)} 1_{F^(k_n)}(X_{k/n})


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <string>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <locale>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

using State = Vector2d;

// ================================================================
// 拡散係数とその導関数
// s(x) = 2 + sin(x), s'(x) = cos(x), s''(x) = -sin(x)
// ================================================================
inline double s_func(double x)   { return 2.0 + sin(x); }
inline double ds_func(double x)  { return cos(x); }
inline double dds_func(double x) { return -sin(x); }

// ================================================================
// テスト関数 f(L) = arctan(L)
// ================================================================
inline double f(double L) { return atan(L); }

// ================================================================
// 数値格式
// ================================================================

// A0: Euler-Maruyama (次数 0.5)
inline State A0(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_func(curr(1)) * dW(0), s_func(curr(0)) * dW(1));
    return curr + drift * dt + diffusion;
}

// A1: Milstein (次数 1.0)
inline State A1(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double s_x0 = s_func(curr(0)), s_x1 = s_func(curr(1));
    double ds_x0 = ds_func(curr(0)), ds_x1 = ds_func(curr(1));
    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_x1 * dW(0), s_x0 * dW(1));
    State base = curr + drift * dt + diffusion;
    Vector2d milstein_corr(
        0.5 * ds_x1 * s_x0 * dW(0) * dW(1),
        0.5 * ds_x0 * s_x1 * dW(0) * dW(1)
    );
    return base + milstein_corr;
}

// A2: 1.5次
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double w1 = dW(0), w2 = dW(1);
    double w1_sq = w1*w1, w2_sq = w2*w2;
    double s0 = s_func(curr(0)), s1 = s_func(curr(1));
    double ds0 = ds_func(curr(0)), ds1 = ds_func(curr(1));
    double dds0 = dds_func(curr(0)), dds1 = dds_func(curr(1));
    Vector2d a(curr(1), -curr(0));
    double w111 = (w1*w1_sq) - 3.0*dt*w1;
    double w222 = (w2*w2_sq) - 3.0*dt*w2;
    double w112 = (w1_sq - dt)*w2;
    double w122 = (w2_sq - dt)*w1;
    State next = curr;
    next(0) += a(0)*dt + s1*w1;
    next(0) += 0.5*ds1*s0*w1*w2;
    next(0) += 0.5*dt*(s0*w2 + ds1*a(1)*w1);
    next(0) += 0.125*dt*(2.0*dds1*s0*s0*w1 + (ds1*s0*ds1*s0/s1)*w1 - ds0*ds1*s1*w2);
    next(0) += (1.0/24.0)*(ds1*ds1/s1)*s0*s0*w111
             + (1.0/12.0)*ds0*ds1*s1*w112
             + (1.0/6.0)*dds1*s0*s0*w122
             + (1.0/24.0)*(ds1*ds1/s1)*s0*s0*w122;
    next(1) += a(1)*dt + s0*w2;
    next(1) += 0.5*ds0*s1*w1*w2;
    next(1) += 0.5*dt*(-s1*w1 + ds0*a(0)*w2);
    next(1) += 0.125*dt*(2.0*dds0*s1*s1*w2 + ds0*s1*ds0*s1*w2/s0 - ds1*ds0*s0*w1);
    next(1) += (1.0/24.0)*(ds0*ds0/s0)*s1*s1*w222
             + (1.0/12.0)*ds0*ds1*s0*w122
             + (1.0/6.0)*dds0*s1*s1*w112
             + (1.0/24.0)*(ds0*ds0/s0)*s1*s1*w112;
    return next;
}

// A3: Ninomiya-Victoir
inline State A3(const State& curr, double dt, double Z1, double Z2, int xi) {
    const double sqrt_dt = sqrt(dt);
    double dW1 = sqrt_dt * Z1;
    double dW2 = sqrt_dt * Z2;
    auto flow_drift = [](double t, const State& x) -> State {
        double c = cos(t), s = sin(t);
        return State(c*x(0)+s*x(1), -s*x(0)+c*x(1));
    };
    auto flow_sigma1 = [](double z, const State& x) -> State {
        return State(x(0) + z*s_func(x(1)), x(1));
    };
    auto flow_sigma2 = [](double z, const State& x) -> State {
        return State(x(0), x(1) + z*s_func(x(0)));
    };
    State y = flow_drift(0.5*dt, curr);
    if (xi >= 0) { y = flow_sigma1(dW1, y); y = flow_sigma2(dW2, y); }
    else          { y = flow_sigma2(dW2, y); y = flow_sigma1(dW1, y); }
    y = flow_drift(0.5*dt, y);
    return y;
}

// ================================================================
// Sierpinski Gasket 構築と指示関数
// ================================================================

// 点 x から線分 [a,b] への距離
inline double dist_point_to_segment(const Vector2d& x,
                                    const Vector2d& a,
                                    const Vector2d& b) {
    Vector2d v = b - a;
    double len2 = v.squaredNorm();
    if (len2 < 1e-18) return (x - a).norm();
    double theta = (x - a).dot(v) / len2;
    theta = min(1.0, max(0.0, theta));
    return (x - (a + theta * v)).norm();
}

inline bool point_in_triangle(const Vector2d& x,
                               const Vector2d& v0,
                               const Vector2d& v1,
                               const Vector2d& v2) {
    double d1 = (x - v1).x()*(v0 - v1).y() - (v0 - v1).x()*(x - v1).y();
    double d2 = (x - v2).x()*(v1 - v2).y() - (v1 - v2).x()*(x - v2).y();
    double d3 = (x - v0).x()*(v2 - v0).y() - (v2 - v0).x()*(x - v0).y();
    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(has_neg && has_pos);
}

// 点 x から三角形 [v0,v1,v2] への距離 (内部なら 0)
inline double dist_point_to_triangle(const Vector2d& x,
                                     const Vector2d& v0,
                                     const Vector2d& v1,
                                     const Vector2d& v2) {
    if (point_in_triangle(x, v0, v1, v2)) return 0.0;
    return min({ dist_point_to_segment(x, v0, v1),
                 dist_point_to_segment(x, v1, v2),
                 dist_point_to_segment(x, v2, v0) });
}

// レベル k の Sierpinski Gasket F^(k) を構築
// Leb(F^(k)) = (sqrt(3)/4)*(3/4)^k
using Triangle = tuple<Vector2d, Vector2d, Vector2d>;

vector<Triangle> build_gasket(int k) {
    const Vector2d v0(0.0, 0.0);
    const Vector2d v1(1.0, 0.0);
    const Vector2d v2(0.5, sqrt(3.0) / 2.0);
    vector<Triangle> tris;
    tris.push_back({v0, v1, v2});
    for (int iter = 0; iter < k; ++iter) {
        vector<Triangle> nxt;
        nxt.reserve(tris.size() * 3);
        for (auto& [a, b, c] : tris) {
            Vector2d mab = 0.5 * (a + b);
            Vector2d mbc = 0.5 * (b + c);
            Vector2d mac = 0.5 * (a + c);
            nxt.push_back({a,   mab, mac});
            nxt.push_back({mab, b,   mbc});
            nxt.push_back({mac, mbc, c  });
        }
        tris = move(nxt);
    }
    return tris;
}

// 1_{F^(k)}(x):
inline bool in_gasket(const Vector2d& x, const vector<Triangle>& tris) {
    for (auto& [a, b, c] : tris)
        if (point_in_triangle(x, a, b, c)) return true;
    return false;
}

// ================================================================
//   k_n = floor( (1+log2(3)) / (log2(3)*(2-log2(3))+2) * log2(n) )
//       ≈ 0.972 * log2(n)
//   上限 8 (3^8=6561 三角形)
// ================================================================
inline int choose_level_gasket(int n) {
    const double log2_3 = log(3.0) / log(2.0);
    const double coeff  = (1.0 + log2_3) / (log2_3 * (2.0 - log2_3) + 2.0);
    int k = (int)floor(coeff * log2((double)n));
    return max(1, min(k, 8));
}

// ================================================================
// Leb(F^(k)) = (sqrt(3)/4) * (3/4)^k   
// ================================================================
inline double leb_gasket(int k) {
    return (sqrt(3.0) / 4.0) * pow(3.0 / 4.0, (double)k);
}

// ================================================================
// Main
// ================================================================
int main() {

    constexpr double alpha   = 1.0;
    constexpr double t_start = 0.0;
    constexpr double t_end   = 1.0;
    constexpr double mu      = 0.0;
    constexpr double sigma_v = 1.0;

    const double dF_gasket = log(3.0) / log(2.0);  // ≈ 1.58496

    const State x0_state = Vector2d(0.1, 0.9);
    constexpr int max_n = 9;

    vector<double> A(max_n+1,0), Am(max_n+1,0);
    vector<double> A_1_5(max_n+1,0), A_nv(max_n+1,0), A_nm(max_n+1,0);
    vector<double> E(max_n+1,0), Em(max_n+1,0);
    vector<double> E_1_5(max_n+1,0), E_nv(max_n+1,0), E_nm(max_n+1,0);

    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str());
    const string csv_path = dir_path + "/2DLTM3_gasket_100_1000_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);
    if (!ofs) { cerr << "CSV を開けません: " << csv_path << endl; return 1; }
    ofs.imbue(locale::classic());
    ofs << "n,points,k_gasket,leb_gasket,"
        << "Var_em,Var_mil,Var_1_5,Var_nv,Var_nm,"
        << "Mean_em,Mean_mil,Mean_1_5,Mean_nv,Mean_nm\n";

    for (int n = 0; n <= max_n; ++n) {

        const int    points     = 100 + 100 * n;
        const int    paths      = 10 * points * points;
        const double dt         = (t_end - t_start) / (points - 1);
        const double dtm        = dt / (points - 1);

        // k_n
        const int    k_n        = choose_level_gasket(points);
        // Leb(F^(k_n))
        const double leb_kn     = leb_gasket(k_n);
        // norm_coeff = 1 / (n * Leb(F^(k_n)))
        const double norm_coeff = 1.0 / ((double)points * leb_kn);

        // F^(k_n) を構築
        auto tris_gasket = build_gasket(k_n);

        long double S=0, Sm=0, S_1_5_acc=0, S_nv_acc=0, S_nm_acc=0;
        long double B_em=0, Bm=0, B_1_5=0, B_nv_acc=0, B_nm_acc=0;

#pragma omp parallel \
    reduction(+: S, Sm, S_1_5_acc, S_nv_acc, S_nm_acc) \
    reduction(+: B_em, Bm, B_1_5, B_nv_acc, B_nm_acc)
        {
        normal_distribution<double> dist_nm0(mu, sigma_v);
        normal_distribution<double> dist_nm1(mu, sigma_v);
        bernoulli_distribution dist_xi(0.5);

#pragma omp for schedule(static) nowait
        for (int p = 0; p < paths; ++p) {

            seed_seq ss0{30u, 0u, (uint32_t)p};
            seed_seq ss1{42u, 1u, (uint32_t)p};
            seed_seq ss2{54u, 2u, (uint32_t)p};
            mt19937 rng_nm(ss0), rng1_nm(ss1), rng_xi(ss2);

            State st_em=x0_state, st_mil=x0_state;
            State st_15=x0_state, st_nv =x0_state;
            State st_nm=x0_state;

            // 指示関数の累積カウント
            double cnt_em=0, cnt_mil=0, cnt_15=0, cnt_nv=0, cnt_nm=0;

            for (int idx = 1; idx < points; ++idx) {

                double Z1=0.0, Z2=0.0;

                // benchmark: Milstein
                for (int m = 0; m < points; ++m) {
                    double Z1_nm = dist_nm0(rng_nm);
                    double Z2_nm = dist_nm1(rng1_nm);
                    st_nm = A1(st_nm, dtm, Z1_nm, Z2_nm);
                    Z1 += Z1_nm / sqrt((double)points);
                    Z2 += Z2_nm / sqrt((double)points);
                }

                // 1_{F^(k_n)}(X_benchmark)
                if (in_gasket(st_nm, tris_gasket)) cnt_nm += 1.0;

                int xi = dist_xi(rng_xi) ? 1 : -1;

                // 格式を更新
                State ne   = A0(st_em,  dt, Z1, Z2);
                State nm_  = A1(st_mil, dt, Z1, Z2);
                State n15  = A2(st_15,  dt, Z1, Z2);
                State nnv  = A3(st_nv,  dt, Z1, Z2, xi);

                // 1_{F^(k_n)}(X_scheme)
                if (in_gasket(ne,  tris_gasket)) cnt_em  += 1.0;
                if (in_gasket(nm_, tris_gasket)) cnt_mil += 1.0;
                if (in_gasket(n15, tris_gasket)) cnt_15  += 1.0;
                if (in_gasket(nnv, tris_gasket)) cnt_nv  += 1.0;

                st_em=ne; st_mil=nm_; st_15=n15; st_nv=nnv;

            } // end idx

            //  L = cnt / (n * Leb(F^(k_n)))
            double L_em  = cnt_em  * norm_coeff;
            double L_mil = cnt_mil * norm_coeff;
            double L_15  = cnt_15  * norm_coeff;
            double L_nv  = cnt_nv  * norm_coeff;
            double L_nm  = cnt_nm  * norm_coeff;

            // テスト関数 f(L) = arctan(L)
            double val_em  = f(L_em);
            double val_mil = f(L_mil);
            double val_15  = f(L_15);
            double val_nv  = f(L_nv);
            double val_nm  = f(L_nm);

            S         += val_em;
            Sm        += val_mil;
            S_1_5_acc += val_15;
            S_nv_acc  += val_nv;
            S_nm_acc  += val_nm;

            B_em     += (val_em  * val_em );
            Bm       += (val_mil * val_mil);
            B_1_5    += (val_15  * val_15 );
            B_nv_acc += (val_nv  * val_nv );
            B_nm_acc += (val_nm  * val_nm );

        } // end path
        } // end omp parallel

        const double inv = 1.0 / paths;
        A[n]     = (double)(S         * inv);
        Am[n]    = (double)(Sm        * inv);
        A_1_5[n] = (double)(S_1_5_acc * inv);
        A_nv[n]  = (double)(S_nv_acc  * inv);
        A_nm[n]  = (double)(S_nm_acc  * inv);
        E[n]     = (double)(B_em      * inv) - A[n]     * A[n];
        Em[n]    = (double)(Bm       * inv) - Am[n]    * Am[n];
        E_1_5[n] = (double)(B_1_5    * inv) - A_1_5[n] * A_1_5[n];
        E_nv[n]  = (double)(B_nv_acc * inv) - A_nv[n]  * A_nv[n];
        E_nm[n]  = (double)(B_nm_acc * inv) - A_nm[n]  * A_nm[n];
        cout << "-------------------------------------------------" << n << "\n";
        cout << setprecision(10) << "points = " << points << "\n";
        cout << "-------------------------------------------------" << "\n";
        cout << setprecision(6)  << "Gasket k_n       = " << k_n    << "\n";
        cout << setprecision(10) << "Leb(F^(k_n))     = " << leb_kn << "\n";
        cout << setprecision(10) << "norm_coeff       = " << norm_coeff << "\n";
        cout << setprecision(15) << "Gasket Var EM        = " << E[n]     << "\n";
        cout << setprecision(15) << "Gasket Var Milstein  = " << Em[n]    << "\n";
        cout << setprecision(15) << "Gasket Var 1.5       = " << E_1_5[n] << "\n";
        cout << setprecision(15) << "Gasket Var NV        = " << E_nv[n]  << "\n";
        cout << setprecision(15) << "Gasket Var benchmark        = " << E_nm[n]  << "\n";
        cout << setprecision(15) << "Gasket Mean EM       = " << A[n]     << "\n";
        cout << setprecision(15) << "Gasket Mean Milstein = " << Am[n]    << "\n";
        cout << setprecision(15) << "Gasket Mean 1.5      = " << A_1_5[n] << "\n";
        cout << setprecision(15) << "Gasket Mean NV       = " << A_nv[n]  << "\n";
        cout << setprecision(15) << "Gasket Mean benchmark       = " << A_nm[n]  << "\n";
        ofs << n << "," << points << "," << k_n << ","
            << fixed << setprecision(10) << leb_kn << ","
            << E[n]     << "," << Em[n]    << ","
            << E_1_5[n] << "," << E_nv[n]  << ","
            << E_nm[n]  << ","
            << A[n]     << "," << Am[n]    << ","
            << A_1_5[n] << "," << A_nv[n]  << "," << A_nm[n] << endl;

    } // end n

    ofs.close();
    cout << "\nCSV を書き込みました: " << csv_path << endl;
    return 0;
}