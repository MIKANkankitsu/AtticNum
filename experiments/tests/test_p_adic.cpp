#include <gtest/gtest.h>
#include <random>
#include <vector>
#include <functional>
#include <gmpxx.h>
#include "p_adic.hpp"

using namespace attic;

class PIntStressTest : public ::testing::Test {
protected:
    std::mt19937_64 rng{1337};
    const std::vector<uint64_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23};

    uint64_t getRandomPrime() {
        std::uniform_int_distribution<size_t> dist(0, primes.size() - 1);
        return primes[dist(rng)];
    }

    uint64_t getRandomInt(uint64_t min_v, uint64_t max_v) {
        std::uniform_int_distribution<uint64_t> dist(min_v, max_v);
        return dist(rng);
    }

    mpq_class getRandompIntegralRational(uint64_t p) {
        while (true) {
            uint64_t num = getRandomInt(0, 7);
            uint64_t den = getRandomInt(1, 8);
            mpq_class q(num, den);
            q.canonicalize();
            if (q.get_den() % p != 0) return q;
        }
    }
};

TEST_F(PIntStressTest, RunAllAssertionFormulas10000) {
    using PInt = attic::p_int;
    
    // p_int 同士の加減算
    auto check_add_sub_identity = [](const PInt& A, const PInt& B, const PInt&, int64_t, uint64_t p, uint64_t prec) {
        EXPECT_EQ((A + B) - A, B);
        EXPECT_EQ((A - B) + B, A);
        PInt res = A;
        res += B;
        EXPECT_EQ(A + B, res);
        res = B;
        res -= A;
        EXPECT_EQ(B - A, res);
        EXPECT_EQ(A - A, PInt(0, p, prec));
        EXPECT_TRUE((A - A).is_zero());
        EXPECT_TRUE((A +(-A)).is_zero());
    };

    // integer_like との加算
    auto check_integer_like_add = [](const PInt& A, const PInt&, const PInt&, int64_t val, uint64_t, uint64_t) {
        PInt res1 = A + val;
        PInt res2 = val + A;
        PInt res3 = A;
        res3 += val;

        EXPECT_EQ(res1, res2);
        EXPECT_EQ(res1, res3);
        EXPECT_EQ(res1 - val, A);
    };

    // integer_like との減算
    auto check_integer_like_sub = [](const PInt& A, const PInt&, const PInt&, int64_t val, uint64_t, uint64_t) {
        PInt res1 = A - val;
        PInt res2 = val - A;
        PInt res3 = A;
        res3 -= val;

        EXPECT_EQ(res1, -res2);
        EXPECT_EQ(res1, res3);
        EXPECT_EQ(res1 + val, A);
    };

    // p_int 同士の乗除算
    auto check_mul_div_identity = [](const PInt& A, const PInt& B, const PInt&, int64_t, uint64_t p, uint64_t prec) {
        if (!A.is_unit() || !B.is_unit()) return;

        EXPECT_EQ((A * B) / A, B);
        EXPECT_EQ((A / B) * B, A);
        PInt res = A;
        res *= B;
        EXPECT_EQ(A * B, res);
        res = B;
        res /= A;
        EXPECT_EQ(B / A, res);
        EXPECT_EQ(A / A, PInt(1, p, prec));
        EXPECT_TRUE((A / A).is_unit());
    };

    // integer_like との乗算
    auto check_integer_like_mul = [](const PInt& A, const PInt&, const PInt&, int64_t val, uint64_t p, uint64_t) {
        PInt res1 = A * val;
        PInt res2 = val * A;
        PInt res3 = A;
        res3 *= val;

        EXPECT_EQ(res1, res2);
        EXPECT_EQ(res1, res3);
        if (val % p != 0) {
            EXPECT_EQ(res1 / val, A);
        }  
    };

    // integer_like との除算
    auto check_integer_like_div = [](const PInt& A, const PInt&, const PInt&, int64_t val, uint64_t p, uint64_t) {
        if (val % p == 0 || !A.is_unit()) return;
        PInt res1 = A / val;
        PInt res2 = val / A;
        PInt res3 = A;
        res3 /= val;

        EXPECT_EQ(res1, res2.inverse());
        EXPECT_EQ(res1, res3);
        EXPECT_EQ(res1 * val, A);
    };

    // 分配律 (やりたいだけ)
    auto check_distributive = [](const PInt& A, const PInt& B, const PInt& C, int64_t, uint64_t, uint64_t) {
        EXPECT_EQ(A * (B+C), A*B + A*C);
        EXPECT_EQ(A * (B-C), A*B - A*C);
        EXPECT_EQ((A + B) * (A - B), A*A - B*B);
        EXPECT_EQ(A*A*A + B*B*B + C*C*C - 3*A*B*C, (A+B+C) * (A*A + B*B + C*C - A*B - B*C - C*A));
        if (C.is_unit()) {
            EXPECT_EQ((A + B/C) * C, A*C + B);
        }
    };

    // 乗除算の単則と可換性 (可逆元の場合)
    auto check_mul_div_inverse = [](const PInt& A, const PInt& B, const PInt&, int64_t, uint64_t p, uint64_t prec) {
        if (A.is_unit() && B.is_unit()) {
            PInt prod = A * B;
            prod /= A;
            prod /= B;
            EXPECT_EQ(prod, PInt(1, p, prec));
            EXPECT_EQ(A * B, B * A);
        } else {
            if (!A.is_unit()) {
                EXPECT_THROW(A.inverse(), std::domain_error);
                EXPECT_THROW(PInt(1,p,prec)/(p * A), std::domain_error);
                EXPECT_THROW(PInt(1,p,prec)/(0 * A), std::domain_error);
            }
        }
    };

    using TestFormula = std::function<void(const PInt&, const PInt&, const PInt&, int64_t, uint64_t, uint64_t)>;
    const std::vector<TestFormula> test_formulas = {
        check_add_sub_identity,
        check_integer_like_add,
        check_integer_like_sub,
        check_mul_div_identity,
        check_integer_like_mul,
        check_integer_like_div,
        check_distributive,
        check_mul_div_inverse
    };

    // -------------------------------------------------------------------------
    // メインの 10,000 回ループ（ランダム生成と各式の全適用）
    // -------------------------------------------------------------------------
    constexpr int NUM_TRIALS = 10000;

    for (int i = 0; i < NUM_TRIALS; ++i) {
        const uint64_t p = getRandomPrime();
        const uint64_t prec = getRandomInt(3, 5);

        const PInt A(getRandompIntegralRational(p), p, prec);
        const PInt B(getRandompIntegralRational(p), p, prec);
        const PInt C(getRandompIntegralRational(p), p, prec);
        const int64_t int_val = getRandomInt(0, 216);

        // 登録されたすべての検証式を順番に実行
        for (const auto& formula : test_formulas) {
            formula(A, B, C, int_val, p, prec);
        }
    }
}

// --- 異精度・Exact/Finite 混合演算のテスト ---
TEST(PIntMixedTest, PrecisionCombining) {
    p_int exact_a(3, 5, 0);     // Exact
    p_int finite_b(4, 5, 2);    // mod 25

    // Exact + Finite -> 低い方の精度 (precision = 2) に合わせられる
    p_int res = exact_a + finite_b;
    EXPECT_FALSE(res.is_exact());
    EXPECT_EQ(res.precision(), 2);
    EXPECT_EQ(res.rep(), mpz_class(7));
}

// --- 例外処理 (@test_throws に相当) のテスト ---
TEST(PIntExceptionTest, DomainAndInvalidArgument) {
    // 素数でない値の指定
    EXPECT_THROW(p_int(10, 4, 0), std::invalid_argument);

    // Exact モードで分母が p で割れる有理数は不可
    // 1/5 は p=5 で domain_error
    EXPECT_THROW(p_int(mpq_class(1, 5), 5, 0), std::domain_error);

    // 0 での割り算
    p_int a(5, 5, 2);
    p_int zero(0, 5, 2);
    EXPECT_THROW(a / zero, std::domain_error);

    // p 進単元でない値 (5 の倍数) での有限精度割り算
    p_int non_unit(5, 5, 2);
    EXPECT_THROW(a / non_unit, std::domain_error);
}