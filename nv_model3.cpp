#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

struct Vec2 {
    double x1{0.0};
    double x2{0.0};

    Vec2 operator+(const Vec2& other) const { return {x1 + other.x1, x2 + other.x2}; }
    Vec2 operator-(const Vec2& other) const { return {x1 - other.x1, x2 - other.x2}; }
    Vec2 operator*(double c) const { return {c * x1, c * x2}; }
};

inline Vec2 operator*(double c, const Vec2& v) { return v * c; }

struct Model3Config {
    // Drift a(x) = (a1(x), a2(x)). Replace this lambda to test other drifts.
    std::function<Vec2(const Vec2&)> drift = [](const Vec2& x) {
        return Vec2{x.x2, -x.x1};
    };

    // Volatility profile s(x). Replace this lambda to test other s.
    std::function<double(double)> s = [](double x) {
        return 2.0 + std::sin(x);
    };

    // Optional exact flow for the drift ODE x' = a(x).
    // If empty, the code falls back to RK4.
    std::function<Vec2(double, const Vec2&)> exact_drift_flow = [](double t, const Vec2& x) {
        const double c = std::cos(t);
        const double s = std::sin(t);
        return Vec2{c * x.x1 + s * x.x2, -s * x.x1 + c * x.x2};
    };

    // Number of RK4 substeps used if exact_drift_flow is not supplied.
    int rk4_substeps = 8;
};

class NinomiyaVictoir2D {
public:
    explicit NinomiyaVictoir2D(Model3Config cfg) : cfg_(std::move(cfg)) {}

    Vec2 flow_drift(double dt, const Vec2& x) const {
        if (cfg_.exact_drift_flow) {
            return cfg_.exact_drift_flow(dt, x);
        }
        return rk4_flow(dt, x, cfg_.drift, cfg_.rk4_substeps);
    }

    // sigma_1(x) = (s(x2), 0)
    Vec2 flow_sigma1(double z, const Vec2& x) const {
        return Vec2{x.x1 + z * cfg_.s(x.x2), x.x2};
    }

    // sigma_2(x) = (0, s(x1))
    Vec2 flow_sigma2(double z, const Vec2& x) const {
        return Vec2{x.x1, x.x2 + z * cfg_.s(x.x1)};
    }

    // One Ninomiya-Victoir step.
    // If xi = +1, use order sigma_1 then sigma_2.
    // If xi = -1, use order sigma_2 then sigma_1.
    Vec2 step(double h, double dW1, double dW2, int xi, const Vec2& x) const {
        Vec2 y = flow_drift(0.5 * h, x);
        if (xi >= 0) {
            y = flow_sigma1(dW1, y);
            y = flow_sigma2(dW2, y);
        } else {
            y = flow_sigma2(dW2, y);
            y = flow_sigma1(dW1, y);
        }
        y = flow_drift(0.5 * h, y);
        return y;
    }

    std::vector<Vec2> simulate_path(double T, int n_steps, const Vec2& x0, std::uint64_t seed = 123456789ULL) const {
        if (n_steps <= 0) {
            throw std::invalid_argument("n_steps must be positive");
        }
        const double h = T / static_cast<double>(n_steps);
        const double sqrt_h = std::sqrt(h);

        std::mt19937_64 rng(seed);
        std::normal_distribution<double> normal(0.0, 1.0);
        std::uniform_int_distribution<int> rademacher(0, 1);

        std::vector<Vec2> path;
        path.reserve(static_cast<std::size_t>(n_steps) + 1);
        Vec2 x = x0;
        path.push_back(x);

        for (int k = 0; k < n_steps; ++k) {
            const double dW1 = sqrt_h * normal(rng);
            const double dW2 = sqrt_h * normal(rng);
            const int xi = (rademacher(rng) == 0 ? -1 : 1);
            x = step(h, dW1, dW2, xi, x);
            path.push_back(x);
        }
        return path;
    }

private:
    Model3Config cfg_;

    static Vec2 rk4_flow(double dt, const Vec2& x0,
                         const std::function<Vec2(const Vec2&)>& f,
                         int substeps) {
        if (substeps <= 0) {
            throw std::invalid_argument("rk4_substeps must be positive");
        }
        const double h = dt / static_cast<double>(substeps);
        Vec2 x = x0;
        for (int i = 0; i < substeps; ++i) {
            const Vec2 k1 = f(x);
            const Vec2 k2 = f(x + 0.5 * h * k1);
            const Vec2 k3 = f(x + 0.5 * h * k2);
            const Vec2 k4 = f(x + h * k3);
            x = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        return x;
    }
};

int main() {
    Model3Config cfg;

    // Default model from the note:
    //   s(x) = 2 + sin(x)
    //   a(x) = (x2, -x1)
    // To test another model, replace cfg.s and/or cfg.drift.
    // If the new drift has no exact flow, set cfg.exact_drift_flow = {};
    // and the code will use RK4 for the drift half-steps.

    NinomiyaVictoir2D nv(cfg);

    const Vec2 x0{1.0, 1.0};
    const double T = 1.0;
    const int n_steps = 1000;
    const std::uint64_t seed = 20260405ULL;

    const std::vector<Vec2> path = nv.simulate_path(T, n_steps, x0, seed);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "# k x1 x2\n";
    for (std::size_t k = 0; k < path.size(); ++k) {
        std::cout << k << ' ' << path[k].x1 << ' ' << path[k].x2 << '\n';
    }

    return 0;
}
