#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <gmp.h>
#include <gmpxx.h>
#include <stdexcept>
#include <type_traits>
#include <utility>

// sizeof(unsigned long) が 8 バイト未満ならコンパイルエラーにしてビルドを中断
static_assert(sizeof(unsigned long) >= 8,
    "unsigned long int must be at least 64-bit (8 bytes). "
    "Windows (LLP64) is not directly supported without 64-bit truncation helpers!");

namespace attic {

struct unchecked_prime_t {
    explicit unchecked_prime_t() = default;
};

inline constexpr unchecked_prime_t unchecked_prime{};

template <class T>
concept integer_like =
    std::integral<std::remove_cvref_t<T>> ||
    std::same_as<std::remove_cvref_t<T>, mpz_class>;

class p_int {
    private:
        uint64_t prime_;
        uint64_t precision_;

        mpz_class modulus_;
        mpz_class int_p_;

        mpq_class exact_q_;

        //素数判定
        static bool is_prime_u64(uint64_t p) {
            if (p < 2) return false;

            mpz_class z = p;
            return mpz_probab_prime_p(z.get_mpz_t(), 25) != 0;
        }

        static uint64_t checked_prime(uint64_t prime) {
            if (!is_prime_u64(prime)) {
                throw std::invalid_argument("p_int prime must be prime");
            }

            return prime;
        }

        //コンストラクタ関連
        static mpz_class make_modulus(uint64_t prime, uint64_t precision) {
            mpz_class mod;
            mpz_ui_pow_ui(mod.get_mpz_t(), prime, precision);
            return mod;
        }

        //計算に関すること
        static uint64_t combine_precision(uint64_t a, uint64_t b) {
            if (a == 0) return b;
            if (b == 0) return a;
            return std::min(a, b);
        }

        static bool divisible_by_prime(const mpz_class& value, uint64_t prime) {
            return mpz_divisible_ui_p(value.get_mpz_t(), prime) != 0;
        }

        static void require_p_integral(const mpq_class& value, uint64_t prime) {
            if (divisible_by_prime(value.get_den(), prime)) {
                throw std::domain_error("p_int exact rational denominator is divisible by prime");
            }
        }

        static void normalize_mod(mpz_class& value, const mpz_class& modulus) {
            mpz_mod(value.get_mpz_t(), value.get_mpz_t(), modulus.get_mpz_t());
        }

        static void normalize_mod(mpz_class& r, const mpz_class& a, const mpz_class& m) {
            mpz_mod(r.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t());
        }

        //Hensel持ち上げを利用した高速な有限精度除算
        static mpz_class inverse_mod_prime_power(
                const mpz_class& value,
                uint64_t prime,
                uint64_t precision,
                const mpz_class& modulus
                ) {
            if (precision == 0) {
                throw std::logic_error("inverse modulo p^precision requires finite precision");
            }

            mpz_class p = prime;
            mpz_class value_mod_p;
            normalize_mod(value_mod_p, value, p);

            mpz_class inv;
            if (mpz_invert(inv.get_mpz_t(), value_mod_p.get_mpz_t(), p.get_mpz_t()) == 0) {
                throw std::domain_error("p_int value is not invertible modulo prime");
            }

            uint64_t lifted_precision = 1;
            mpz_class lifted_modulus = p;
            mpz_class value_mod, temp;
            while (lifted_precision < precision) {
                if (lifted_precision < (precision >> 1)) {
                    lifted_precision <<= 1;
                    lifted_modulus *= lifted_modulus;
                } else {
                    lifted_precision = precision;
                    lifted_modulus = modulus;
                }
                
                normalize_mod(value_mod, value, lifted_modulus);

                // temp = 2 - value_mod * inv
                mpz_set_ui(temp.get_mpz_t(), 2);
                mpz_submul(temp.get_mpz_t(), value_mod.get_mpz_t(), inv.get_mpz_t());

                inv *= temp;
                normalize_mod(inv, lifted_modulus);
            }

            return inv;
        }

        //exactから有限精度
        static mpz_class finite_rep_from_q(const mpq_class& value, uint64_t prime, const mpz_class& modulus, uint64_t precision) {
            mpz_class rep = inverse_mod_prime_power(value.get_den(), prime, precision, modulus);
            rep *= value.get_num();
            normalize_mod(rep, modulus);
            return rep;
        }

        void normalize() {
            if (is_exact()) return;

            normalize_mod(int_p_, modulus_);
        }

        void require_same_prime(const p_int& rhs) const {
            if (prime_ != rhs.prime_) {
                throw std::invalid_argument("p_int operation requires matching primes");
            }
        }

    public:
        //コンストラクタ

        // uncheck mod z
        p_int(unchecked_prime_t, mpz_class value, uint64_t prime, uint64_t precision, mpz_class modulus)
            : prime_(prime),
              precision_(precision),
              modulus_(precision == 0 ? mpz_class(0) : std::move(modulus)),
              int_p_(precision == 0 ? mpz_class(0) : std::move(value)),
              exact_q_(precision == 0 ? mpq_class(value) : mpq_class(0)) {
            
            normalize();
        }

        // uncheck mod q
        p_int(unchecked_prime_t, mpq_class value, uint64_t prime, uint64_t precision, mpz_class modulus)
            : prime_(prime),
              precision_(precision),
              modulus_(precision == 0 ? mpz_class(0) : std::move(modulus)),
              int_p_(0),
              exact_q_(0) {
            
            value.canonicalize();
            require_p_integral(value, prime_);

            if (is_exact()) {
                exact_q_ = std::move(value);
            } else {
                int_p_ = finite_rep_from_q(value, prime_, modulus_, precision_);
            }
        }

        // uncheck nomod z  -->  uncheck mod z
        p_int(unchecked_prime_t tag, mpz_class value, uint64_t prime, uint64_t precision)
            : p_int(tag, std::move(value), prime, precision, precision == 0 ? mpz_class(0) : make_modulus(prime, precision)) {}

        // uncheck nomod q  -->  uncheck mod q
        p_int(unchecked_prime_t tag, mpq_class value, uint64_t prime, uint64_t precision)
            : p_int(tag, std::move(value), prime, precision, precision == 0 ? mpz_class(0) : make_modulus(prime, precision)) {}

        // check mod z  -->  uncheck mod z
        p_int(mpz_class value, uint64_t prime, uint64_t precision, mpz_class modulus)
            : p_int(unchecked_prime, std::move(value), checked_prime(prime), precision, std::move(modulus)) {}

        // check mod q  -->  uncheck mod q
        p_int(mpq_class value, uint64_t prime, uint64_t precision, mpz_class modulus)
            : p_int(unchecked_prime, std::move(value), checked_prime(prime), precision, std::move(modulus)) {}

        // check nomod z  -->  uncheck nomod z
        p_int(mpz_class value, uint64_t prime, uint64_t precision)
            : p_int(unchecked_prime, std::move(value), checked_prime(prime), precision) {}

        // check nomod q  -->  uncheck nomod q
        p_int(mpq_class value, uint64_t prime, uint64_t precision)
            : p_int(unchecked_prime, std::move(value), checked_prime(prime), precision) {}

        // uncheck nomod T  -->  uncheck nomod z
        template <std::integral T>
        p_int(unchecked_prime_t tag, T value, uint64_t prime, uint64_t precision)
            : p_int(tag, mpz_class(value), prime, precision) {}

        // uncheck mod T  -->  uncheck mod z
        template <std::integral T>
        p_int(unchecked_prime_t tag, T value, uint64_t prime, uint64_t precision, mpz_class modulus)
            : p_int(tag, mpz_class(value), prime, precision, std::move(modulus)) {}

        // check nomod T  -->  check nomod z
        template <std::integral T>
        p_int(T value, uint64_t prime, uint64_t precision)
            : p_int(mpz_class(value), prime, precision) {}

        // check mod T -->  check mod z
        template <std::integral T>
        p_int(T value, uint64_t prime, uint64_t precision, mpz_class modulus)
            : p_int(mpz_class(value), prime, precision, std::move(modulus)) {}

        //各要素取得
        bool is_exact() const { return precision_ == 0; }
        uint64_t prime() const { return prime_; }
        uint64_t precision() const { return precision_; }

        //representative: 代表元
        const mpz_class& rep() const {
            if (is_exact()) {
                throw std::logic_error("p_int::rep() is only valid for finite values");
            }

            return int_p_;
        }

        const mpz_class& modulus() const {
            if (is_exact()) {
                throw std::logic_error("p_int::modulus() is only valid for finite values");
            }

            return modulus_;
        }

        const mpq_class& exact_value() const {
            if (!is_exact()) {
                throw std::logic_error("p_int::exact_value() is only valid for exact values");
            }

            return exact_q_;
        }

        //精度変更、より精度が低い方に変更する
        p_int to_precision(uint64_t target_precision) const {
            if (precision_ == target_precision) {
                return *this;
            }

            if (target_precision == 0) {
                throw std::domain_error("cannot recover an exact p_int from a finite value");
            }

            if (is_exact()) {
                return p_int(unchecked_prime, exact_q_, prime_, target_precision);
            }

            if (precision_ < target_precision) {
                throw std::domain_error("cannot increase finite p_int precision");
            }

            return p_int(unchecked_prime, int_p_, prime_, target_precision);
        }

        //足し算
        p_int& operator+=(const p_int& rhs) {
            require_same_prime(rhs);

            if (precision_ == rhs.precision_) {
                if (is_exact()) {
                    exact_q_ += rhs.exact_q_;
                } else {
                    int_p_ += rhs.int_p_;
                    if (int_p_ >= modulus_) {
                        int_p_ -= modulus_;
                    }
                }
                return *this;
            }

            if (precision_ < rhs.precision_) {
                if (is_exact()) {
                    // exact += finite 右の設定を引きつぐ
                    int_p_ = finite_rep_from_q(exact_q_, prime_, rhs.modulus_, rhs.precision_) + rhs.int_p_;
                    normalize_mod(int_p_, rhs.modulus_);
                    exact_q_ = 0;
                    precision_ = rhs.precision_;
                    modulus_ = rhs.modulus_;
                } else {
                    // low += high
                    int_p_ += rhs.int_p_;
                    normalize_mod(int_p_, modulus_);
                }
                return *this;
            } else { // precision_ > rhs.precision_
                if (rhs.is_exact()) {
                    // finite += exact
                    int_p_ += finite_rep_from_q(rhs.exact_q_, prime_, modulus_, precision_);
                    normalize_mod(int_p_, modulus_);
                } else {
                    // high += low
                    int_p_ += rhs.int_p_;
                    normalize_mod(int_p_, rhs.modulus_);
                    precision_ = rhs.precision_;
                    modulus_ = rhs.modulus_;
                }
                return *this;
            }
        }

        template <integer_like T>
        p_int& operator+=(T&& rhs) {
            if (is_exact()) {
                exact_q_ += rhs;
            } else {
                int_p_ += rhs;
                normalize_mod(int_p_, modulus_);
            }
            return *this;
        }

        friend p_int operator+(p_int lhs, const p_int& rhs) {
            lhs += rhs;
            return lhs;
        }

        template <integer_like T>
        friend p_int operator+(p_int lhs, T&& rhs) {
            lhs += std::forward<T>(rhs);
            return lhs;
        }

        template <integer_like T>
        friend p_int operator+(T&& lhs, p_int rhs) {
            rhs += std::forward<T>(lhs);
            return rhs;
        }

        //引き算
        p_int& operator-=(const p_int& rhs) {
            require_same_prime(rhs);

            if (is_exact() && rhs.is_exact()) {
                exact_q_ -= rhs.exact_q_;
                exact_q_.canonicalize();
                return *this;
            }

            uint64_t precision = combine_precision(precision_, rhs.precision_);
            p_int l = to_precision(precision);
            p_int r = rhs.to_precision(precision);
            l.int_p_ -= r.int_p_;
            l.normalize();
            *this = std::move(l);
            return *this;
        }

        template <integer_like T>
        p_int& operator-=(T&& rhs) {
            return *this -= p_int(
                unchecked_prime,
                std::forward<T>(rhs),
                prime_,
                0
            );
        }

        friend p_int operator-(p_int lhs, const p_int& rhs) {
            lhs -= rhs;
            return lhs;
        }

        template <integer_like T>
        friend p_int operator-(p_int lhs, T&& rhs) {
            lhs -= std::forward<T>(rhs);
            return lhs;
        }

        template <integer_like T>
        friend p_int operator-(T&& lhs, const p_int& rhs) {
            p_int result(
                unchecked_prime,
                std::forward<T>(lhs),
                rhs.prime_,
                0
            );

            result -= rhs;
            return result;
        }

        p_int operator-() const {
            p_int res(*this);
            if (res.is_exact()) {
                res.exact_q_ = -res.exact_q_;
                res.exact_q_.canonicalize();
            } else {
                res.int_p_ = -res.int_p_;
                res.normalize();
            }
            return res;
        }

        //掛け算
        p_int& operator*=(const p_int& rhs) {
            require_same_prime(rhs);

            if (is_exact() && rhs.is_exact()) {
                exact_q_ *= rhs.exact_q_;
                exact_q_.canonicalize();
                return *this;
            }

            uint64_t precision = combine_precision(precision_, rhs.precision_);
            p_int l = to_precision(precision);
            p_int r = rhs.to_precision(precision);
            l.int_p_ *= r.int_p_;
            l.normalize();
            *this = std::move(l);
            return *this;
        }

        template <integer_like T>
        p_int& operator*=(T&& rhs) {
            return *this *= p_int(
                unchecked_prime,
                std::forward<T>(rhs),
                prime_,
                0
            );
        }

        friend p_int operator*(p_int lhs, const p_int& rhs) {
            lhs *= rhs;
            return lhs;
        }

        template <integer_like T>
        friend p_int operator*(p_int lhs, T&& rhs) {
            lhs *= std::forward<T>(rhs);
            return lhs;
        }

        template <integer_like T>
        friend p_int operator*(T&& lhs, const p_int& rhs) {
            p_int result(
                unchecked_prime,
                std::forward<T>(lhs),
                rhs.prime_,
                0
            );

            result *= rhs;
            return result;
        }

        //割り算
        bool is_zero() const {
            return is_exact() ? exact_q_ == 0 : int_p_ == 0;
            // return (exact_q_ == 0) && (int_p_ == 0);
        }

        bool is_unit() const {
            if (is_exact()) {
                if (exact_q_ == 0) return false;
                return !divisible_by_prime(exact_q_.get_num(), prime_);
            }

            return !divisible_by_prime(int_p_, prime_);
        }

        p_int inverse() const {
            if (!is_unit()) {
                throw std::domain_error("p_int inverse is only defined for units");
            }

            if (is_exact()) {
                mpq_class inv = 1 / exact_q_;
                // もうunitと確定しているので、逆数の分母はpで割り切れない
                return p_int(unchecked_prime, std::move(inv), prime_, 0);
            }

            mpz_class inv = inverse_mod_prime_power(int_p_, prime_, precision_, modulus_);
            return p_int(unchecked_prime, std::move(inv), prime_, precision_);
        }

        p_int& operator/=(const p_int& rhs) {
            require_same_prime(rhs);
            if (rhs.is_zero()) {
                throw std::domain_error("cannot divide p_int by zero");
            }

            if (is_exact() && rhs.is_exact()) {
                mpq_class result = exact_q_ / rhs.exact_q_;
                result.canonicalize();
                require_p_integral(result, prime_);
                *this = p_int(unchecked_prime, std::move(result), prime_, 0);
                return *this;
            }

            uint64_t precision = combine_precision(precision_, rhs.precision_);
            p_int l = to_precision(precision);
            p_int r = rhs.to_precision(precision);
            if (!r.is_unit()) {
                throw std::domain_error("finite p_int division currently requires a unit divisor");
            }

            l *= r.inverse();
            *this = std::move(l);
            return *this;
        }

        template <integer_like T>
        p_int& operator/=(T&& rhs) {
            return *this /= p_int(
                unchecked_prime,
                std::forward<T>(rhs),
                prime_,
                0
            );
        }

        friend p_int operator/(p_int lhs, const p_int& rhs) {
            lhs /= rhs;
            return lhs;
        }

        template <integer_like T>
        friend p_int operator/(p_int lhs, T&& rhs) {
            lhs /= std::forward<T>(rhs);
            return lhs;
        }

        template <integer_like T>
        friend p_int operator/(T&& lhs, const p_int& rhs) {
            p_int result(
                unchecked_prime,
                std::forward<T>(lhs),
                rhs.prime_,
                0
            );

            result /= rhs;
            return result;
        }

        //テスト用整数回累乗
        p_int pow_uint(uint64_t n) const {
            p_int base(*this);
            p_int result = is_exact()
                ? p_int(unchecked_prime, mpz_class(1), prime_, 0)
                : p_int(unchecked_prime, mpz_class(1), prime_, precision_);

            while (n > 0) {
                if (n & 1) result *= base;
                base *= base;
                n >>= 1;
            }
            return result;
        }
};

}//namespace
