#include <iostream>
#include <random>
#include <vector>
#include <omp.h>
#include <iomanip>

int main() {
    const int paths = 4;   
    const int threads = 8; // スレッド数を指定

    omp_set_num_threads(threads);

    // =========================================================
    // A: p番号ベースのシード
    // =========================================================
    std::cout << "========================================\n";
    std::cout << "【A】p番号ベースのシード\n";
    std::cout << "========================================\n";

    struct ResultA {
        int p, tid;
        double val0, val1, val2;
    };
    std::vector<ResultA> resultsA(paths);

    #pragma omp parallel for schedule(static) num_threads(threads)
    for (int p = 0; p < paths; ++p) {
        int tid = omp_get_thread_num();

        std::seed_seq ss0{42u, 0u, (uint32_t)p};
        std::seed_seq ss1{42u, 1u, (uint32_t)p};
        std::seed_seq ss2{42u, 2u, (uint32_t)p};
        std::mt19937 rng(ss0);
        std::mt19937 rng1(ss1);
        std::mt19937 rng2(ss2);

        std::normal_distribution<double> dist(0.0, 1.0);
        double v0 = dist(rng);
        double v1 = dist(rng1);
        double v2 = dist(rng2);

        resultsA[p] = {p, tid, v0, v1, v2};
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "p=0とp=1のrng値が違う → パスごとに独立していることを確認\n\n";
    for (auto& r : resultsA) {
        std::cout << "  path=" << r.p
                  << "  担当thread=" << r.tid
                  << "  rng=" << std::setw(10) << r.val0
                  << "  rng1=" << std::setw(10) << r.val1
                  << "  rng2=" << std::setw(10) << r.val2
                  << "\n";
    }

    // =========================================================
    // Aは同じpなら何スレッドでも同じ値になるか？
    // =========================================================
    std::cout << "\n----------------------------------------\n";
    std::cout << "【A 重要確認】同じp=0を8スレッドで独立計算 → 全部同じ値になるはず\n";
    std::cout << "----------------------------------------\n";

    std::vector<double> checkA(threads);
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (int t = 0; t < threads; ++t) {
        int p = 0; // 全員同じpで計算
        std::seed_seq ss0{42u, 0u, (uint32_t)p};
        std::mt19937 rng(ss0);
        std::normal_distribution<double> dist(0.0, 1.0);
        checkA[t] = dist(rng);
    }

    for (int t = 0; t < threads; ++t) {
        std::cout << "  thread=" << t << "  p=0のrng値=" << checkA[t];
        if (t == 0) std::cout << "  ← 基準値";
        else {
            if (checkA[t] == checkA[0])
                std::cout << "  ✓ 基準値と一致";
            else
                std::cout << "  ✗ 違う！";
        }
        std::cout << "\n";
    }

    // =========================================================
    // B: スレッドベースのシード（非推奨方式）
    // =========================================================
    std::cout << "\n========================================\n";
    std::cout << "【B】スレッドベースのシード（非推奨方式）\n";
    std::cout << "========================================\n";

    struct ResultB {
        int p, tid;
        double val0, val1, val2;
    };
    std::vector<ResultB> resultsB(paths);

    #pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        std::mt19937 rng(42);
        std::mt19937 rng1(30);
        std::mt19937 rng2(50);
        std::normal_distribution<double> dist(0.0, 1.0);

        #pragma omp for schedule(static)
        for (int p = 0; p < paths; ++p) {
            double v0 = dist(rng);
            double v1 = dist(rng1);
            double v2 = dist(rng2);
            resultsB[p] = {p, tid, v0, v1, v2};
        }
    }

    std::cout << "全スレッドが同じシード(42/30/50)で初期化 → スレッドが変わると値が変わる\n\n";
    for (auto& r : resultsB) {
        std::cout << "  path=" << r.p
                  << "  担当thread=" << r.tid
                  << "  rng=" << std::setw(10) << r.val0
                  << "  rng1=" << std::setw(10) << r.val1
                  << "  rng2=" << std::setw(10) << r.val2
                  << "\n";
    }

    // =========================================================
    // 【重要テスト】Bは同じpでもスレッドが変わると値が変わるか？
    // =========================================================
    std::cout << "\n----------------------------------------\n";
    std::cout << "【B 重要確認】p=0を8スレッドで独立計算 → スレッドは全員同じシードで初期化されるので全部同じ値\n";
    std::cout << "（ただし実際のシミュレーションではスレッドごとに担当パスが違うため\n";
    std::cout << " 同じパスを再計算すると担当スレッドが変わり結果が変わる恐れがある）\n";
    std::cout << "----------------------------------------\n";

    std::vector<double> checkB(threads);
    #pragma omp parallel for schedule(static) num_threads(threads)
    for (int t = 0; t < threads; ++t) {
        std::mt19937 rng(42);
        std::normal_distribution<double> dist(0.0, 1.0);
        checkB[t] = dist(rng);
    }

    for (int t = 0; t < threads; ++t) {
        std::cout << "  thread=" << t << "  rng値=" << checkB[t];
        if (t == 0) std::cout << "  ← 基準値";
        else {
            if (checkB[t] == checkB[0])
                std::cout << "  ✓ 同じ（シードが同じなので当然）";
            else
                std::cout << "  ✗ 違う！";
        }
        std::cout << "\n";
    }

    // =========================================================
    // まとめ
    // =========================================================
    std::cout << "\n========================================\n";
    std::cout << "【まとめ】\n";
    std::cout << "========================================\n";
    std::cout << "A方式: 同じpなら何スレッドでも同じ値 → 再現性あり ✓\n";
    std::cout << "B方式: 全スレッドが同じシードで初期化される\n";
    std::cout << "       → スレッド数や担当パスが変わると再現性が壊れる ✗\n";

    return 0;
}