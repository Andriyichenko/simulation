//2DLTM3_koch

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <locale>
#include <algorithm>
#include <Eigen/Dense>
#include <omp.h>
#include <cstdlib>
#include <cstdint>

using namespace std;
using namespace Eigen;

using State = Vector2d;

// ================================================================
// 拡散係数とその導関数 
//   s(x)   = 2 + sin(x)
//   s'(x)  = cos(x)
//   s''(x) = -sin(x)
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

// A0: Euler-Maruyama 法 (強収束次数 0.5)
inline State A0(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    Vector2d drift(curr(1), -curr(0));
    Vector2d diffusion(s_func(curr(1)) * dW(0), s_func(curr(0)) * dW(1));
    return curr + drift * dt + diffusion;
}

// A1: Milstein 法 (強収束次数 1.0)
inline State A1(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double s_x0  = s_func(curr(0)),  s_x1  = s_func(curr(1));
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

// A2: 1.5 次元
inline State A2(const State& curr, double dt, double Z1, double Z2) {
    const double sqrt_dt = sqrt(dt);
    Vector2d dW(sqrt_dt * Z1, sqrt_dt * Z2);
    double w1 = dW(0), w2 = dW(1);
    double w1_sq = w1*w1, w2_sq = w2*w2;
    double s0   = s_func(curr(0)),   s1   = s_func(curr(1));
    double ds0  = ds_func(curr(0)),  ds1  = ds_func(curr(1));
    double dds0 = dds_func(curr(0)), dds1 = dds_func(curr(1));
    Vector2d a(curr(1), -curr(0));
    double w111 = (w1*w1_sq) - 3.0*dt*w1;
    double w222 = (w2*w2_sq) - 3.0*dt*w2;
    double w112 = (w1_sq - dt)*w2;
    double w122 = (w2_sq - dt)*w1;
    State next = curr;
    // X1 成分
    next(0) += a(0)*dt + s1*w1;
    next(0) += 0.5*ds1*s0*w1*w2;
    next(0) += 0.5*dt*(s0*w2 + ds1*a(1)*w1);
    next(0) += 0.125*dt*(2.0*dds1*s0*s0*w1 + (ds1*s0*ds1*s0/s1)*w1 - ds0*ds1*s1*w2);
    next(0) += (1.0/24.0)*(ds1*ds1/s1)*s0*s0*w111
             + (1.0/12.0)*ds0*ds1*s1*w112
             + (1.0/6.0)*dds1*s0*s0*w122
             + (1.0/24.0)*(ds1*ds1/s1)*s0*s0*w122;
    // X2 成分
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

// A3: Ninomiya-Victoir 法 
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
    else         { y = flow_sigma2(dW2, y); y = flow_sigma1(dW1, y); }
    y = flow_drift(0.5*dt, y);
    return y;
}

// ================================================================
// Koch 曲線 
// ================================================================

// 点 x から線分 [a,b] までの距離
inline double dist_point_to_segment(const Vector2d& x,
                                    const Vector2d& a,
                                    const Vector2d& b) {
    Vector2d v    = b - a;
    double   len2 = v.squaredNorm();
    if (len2 < 1e-18) return (x - a).norm();
    double theta = (x - a).dot(v) / len2;
    theta = min(1.0, max(0.0, theta));
    return (x - (a + theta * v)).norm();
}

//  レベル m の Koch 折線 F_m を構築
// 初期: a_0=(0,0), b_0=(1,0)
// m 回反復 => |Ω_m| = 4^m 本
vector<pair<Vector2d, Vector2d>> build_koch(int m) {
    const double c60 = 0.5;
    const double s60 = sqrt(3.0) / 2.0;
    vector<pair<Vector2d, Vector2d>> segs;
    segs.push_back({ Vector2d(0.0, 0.0), Vector2d(1.0, 0.0) });
    for (int k = 0; k < m; ++k) {
        vector<pair<Vector2d, Vector2d>> nxt;
        nxt.reserve(segs.size() * 4);
        for (auto& [a, b] : segs) {
            Vector2d v  = b - a;
            Vector2d Rv(c60*v(0) - s60*v(1), s60*v(0) + c60*v(1));
            Vector2d p2 = a + v / 3.0;
            Vector2d p3 = p2 + Rv / 3.0;
            Vector2d p4 = a + 2.0 * v / 3.0;
            nxt.push_back({ a,  p2 });
            nxt.push_back({ p2, p3 });
            nxt.push_back({ p3, p4 });
            nxt.push_back({ p4, b  });
        }
        segs = move(nxt);
    }
    return segs;
}

// d(x, F_m) = min_{ω} d(x, [a_ω, b_ω])
inline double dist_koch(const Vector2d& x,
                        const vector<pair<Vector2d, Vector2d>>& segs) {
    double dmin = 1e18;
    for (auto& [a, b] : segs)
        dmin = min(dmin, dist_point_to_segment(x, a, b));
    return dmin;
}

// ================================================================
//   1 ステップの寄与 = ε_n^{-(2-d_F)} · exp(-d²/(2ε_n²))
//   ε_n² = h^α,  ε_n^{2-d_F} = h^{α(2-d_F)/2}
//   Koch: d_F = log4/log3 ≈ 1.2619
// ================================================================
inline double fractal_kernel(double d, double dF, double alpha, double h) {
    double eps2     = pow(h, alpha);
    double eps_norm = pow(h, alpha * (2.0 - dF) / 2.0);
    return (1.0 / eps_norm) * exp(-d * d / (2.0 * eps2));
}

// ================================================================
//   Koch: r=1/3, m=⌈log(1/(c·ε_n))/log3⌉, 上限 6
// ================================================================
inline int choose_level_koch(double alpha, double h, double c = 2.0) {
    const double r     = 1.0 / 3.0;
    const double eps_n = pow(h, alpha / 2.0);
    if (eps_n <= 0.0 || c <= 0.0) return 1;
    double val = log(1.0 / (c * eps_n)) / log(1.0 / r);
    int m = (int)ceil(val);
    return max(1, min(m, 6));
}

// ================================================================
// メインプログラム
// ================================================================
int main() {

    constexpr double alpha   = 1.0;
    constexpr double t_start = 0.0;
    constexpr double t_end   = 100.0;
    constexpr double mu      = 0.0;
    constexpr double sigma_v = 1.0;

    // Koch 曲線 Hausdorff 次元: d_F = log4/log3
    const double dF_koch = log(4.0) / log(3.0);

    // 初期状態 X_0 = (1, 1)
    const State   x0_state = Vector2d(1.0, 1.0);
    constexpr int max_n    = 9;

    // ---- 結果格納配列 ----
    vector<double> A(max_n+1,0),     Am(max_n+1,0);
    vector<double> A_1_5(max_n+1,0), A_nv(max_n+1,0);
    vector<double> E(max_n+1,0),     Em(max_n+1,0);
    vector<double> E_1_5(max_n+1,0), E_nv(max_n+1,0);

    // ---- CSV 設定 ----
    const string dir_path = "../data_source";
    system(("mkdir -p " + dir_path).c_str());
    const string csv_path = dir_path + "/2DLTM3_koch_100_1000_data.csv";
    ofstream ofs(csv_path, ios::out | ios::trunc);
    if (!ofs) { cerr << "CSV を開けません: " << csv_path << endl; return 1; }
    ofs.imbue(locale::classic());
    ofs << "n,points,m_koch,"
        << "Var_em,Var_mil,Var_1_5,Var_nv,"
        << "Mean_em,Mean_mil,Mean_1_5,Mean_nv\n";

    // ================================================================
    // 外側ループ n = 0..max_n
    // ================================================================
    for (int n = 0; n <= max_n; ++n) {

        const int    points = 100 + 100 * n;
        const int    paths  = 10 * points * points;
        const double dt     = (t_end - t_start) / (points - 1);  // 粗いステップ h_n
        const double dtm    = dt / (points - 1);                  // ベンチマーク細かいステップ

        // int  m_koch    = choose_level_koch(alpha, dt);
        int  m_koch    = 10;
        auto segs_koch = build_koch(m_koch);

        // reduction 
        long double S=0,    Sm=0,    S_1_5_acc=0,  S_nv_acc=0;
        long double B_=0,   Bm=0,    B_1_5=0,      B_nv_acc=0;

#pragma omp parallel \
    reduction(+: S, Sm, S_1_5_acc, S_nv_acc) \
    reduction(+: B_, Bm, B_1_5, B_nv_acc)
        {
            normal_distribution<double> dist_nm0(mu, sigma_v);
            normal_distribution<double> dist_nm1(mu, sigma_v);
            bernoulli_distribution      dist_xi(0.5);

#pragma omp for schedule(static) nowait
            for (int p = 0; p < paths; ++p) {

                // 再現性のあるシード 
                seed_seq ss0{30u, 0u, (uint32_t)p};
                seed_seq ss1{42u, 1u, (uint32_t)p};
                seed_seq ss2{54u, 2u, (uint32_t)p};
                mt19937 rng_nm(ss0), rng1_nm(ss1), rng_xi(ss2);

                State st_em=x0_state, st_mil=x0_state;
                State st_15=x0_state, st_nv =x0_state;
                State st_nm=x0_state;

                double L_em=0, L_mil=0, L_15=0, L_nv=0, L_nm=0;

                for (int idx = 1; idx < points; ++idx) {

                    double Z1=0.0, Z2=0.0;

                    // ================================================
                    //　benchmark
                    // ================================================
                    for (int m = 0; m < points; ++m) {
                        double Z1_nm = dist_nm0(rng_nm);
                        double Z2_nm = dist_nm1(rng1_nm);
                        State nm_benchmark = A1(st_nm, dtm, Z1_nm, Z2_nm);
                        st_nm = nm_benchmark;
                        Z1 += Z1_nm / sqrt((double)points);
                        Z2 += Z2_nm / sqrt((double)points);
                    }

                    // L_nm += dt * phi_n(..., st_nm(0), dt)
                    L_nm += dt * fractal_kernel(
                        dist_koch(st_nm, segs_koch), dF_koch, alpha, dt);

                    int xi = dist_xi(rng_xi) ? 1 : -1;

                    // ================================================
                    // 状態更新 
                    // ================================================
                    State ne  = A0(st_em,  dt, Z1, Z2);
                    State nm_ = A1(st_mil, dt, Z1, Z2);
                    State n15 = A2(st_15,  dt, Z1, Z2);
                    State nnv = A3(st_nv,  dt, Z1, Z2, xi);

                    // ================================================
                    //   L_Koch += h · K_ε(d(X_{t_k}, F_m))
                    // ================================================
                    L_em  += dt * fractal_kernel(dist_koch(ne,  segs_koch), dF_koch, alpha, dt);
                    L_mil += dt * fractal_kernel(dist_koch(nm_, segs_koch), dF_koch, alpha, dt);
                    L_15  += dt * fractal_kernel(dist_koch(n15, segs_koch), dF_koch, alpha, dt);
                    L_nv  += dt * fractal_kernel(dist_koch(nnv, segs_koch), dF_koch, alpha, dt);

                    st_em=ne; st_mil=nm_; st_15=n15; st_nv=nnv;

                } // end idx

                // テスト関数 f(L) = arctan(L) を適用
                double val_em  = f(L_em);
                double val_mil = f(L_mil);
                double val_15  = f(L_15);
                double val_nv  = f(L_nv);
                double val_nm  = f(L_nm);

                // 平均誤差の累積: S += f(L_nm) - f(L_scheme)
                S        += val_nm - val_em;
                Sm       += val_nm - val_mil;
                S_1_5_acc += val_nm - val_15;
                S_nv_acc  += val_nm - val_nv;

                // 分散計算: B += (f(L_scheme) - f(L_nm))²
                B_        += (val_em  - val_nm) * (val_em  - val_nm);
                Bm        += (val_mil - val_nm) * (val_mil - val_nm);
                B_1_5     += (val_15  - val_nm) * (val_15  - val_nm);
                B_nv_acc  += (val_nv  - val_nm) * (val_nv  - val_nm);

            } // end path
        } // end omp parallel

        // ================================================================
        // 統計量の計算
        //   A[n]  = (1/N) Σ( f(L_nm) - f(L_scheme) )
        //   E[n]  = (1/N) Σ(f(L_scheme)-f(L_nm))² - A[n]²
        // ================================================================
        const double inv = 1.0 / paths;

        A[n]     = (double)(S         * inv);
        Am[n]    = (double)(Sm        * inv);
        A_1_5[n] = (double)(S_1_5_acc * inv);
        A_nv[n]  = (double)(S_nv_acc  * inv);

        E[n]     = (double)(B_        * inv) - A[n]     * A[n];
        Em[n]    = (double)(Bm        * inv) - Am[n]    * Am[n];
        E_1_5[n] = (double)(B_1_5     * inv) - A_1_5[n] * A_1_5[n];
        E_nv[n]  = (double)(B_nv_acc  * inv) - A_nv[n]  * A_nv[n];

        // ================================================================
        // コンソール出力 (指定フォーマットに完全準拠)
        // ================================================================
        cout << "-------------------------------------------------" << n << "\n";
        cout << setprecision(10) << "points = " << points << "\n";
        cout << "-------------------------------------------------" << "\n";
        cout << setprecision(6) <<  "Koch m       = " << m_koch     << "\n";
        cout << setprecision(15) << "Koch Var EM       = " << E[n]     << "\n";
        cout << setprecision(15) << "Koch Var Milstein = " << Em[n]    << "\n";
        cout << setprecision(15) << "Koch Var 1.5      = " << E_1_5[n] << "\n";
        cout << setprecision(15) << "Koch Var NV       = " << E_nv[n]  << "\n";
        cout << setprecision(15) << "Koch Mean EM      = " << A[n]     << "\n";
        cout << setprecision(15) << "Koch Mean Milstein= " << Am[n]    << "\n";
        cout << setprecision(15) << "Koch Mean 1.5     = " << A_1_5[n] << "\n";
        cout << setprecision(15) << "Koch Mean NV      = " << A_nv[n]  << "\n";

        // ================================================================
        // CSV 出力
        // ================================================================
        ofs << n << "," << points << "," << m_koch << ","
            << fixed << setprecision(10)
            << E[n]     << "," << Em[n]    << ","
            << E_1_5[n] << "," << E_nv[n]  << ","
            << A[n]     << "," << Am[n]    << ","
            << A_1_5[n] << "," << A_nv[n]  << endl;

    } // end n

    ofs.close();
    cout << "\nCSV を書き込みました: " << csv_path << endl;
    return 0;
}