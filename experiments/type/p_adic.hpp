#include <concepts>
#include <cstdint>
#include <gmp.h>
#include <gmpxx.h>
#include <stdexcept>
#include <sys/types.h>
#include <type_traits>
#include <utility>

namespace attic{

struct unchecked_prime_t {
    explicit unchecked_prime_t() = default;
};

inline constexpr unchecked_prime_t unchecked_prime{};

class p_int {
    private:
        uint64_t prime_;
        uint64_t precision_;
        mpz_class modulus_;
        mpz_class int_p_;
        // 気持ちとしてはint_p mod p^k

        static bool is_prime_u64(uint64_t p) {
            if (p<2) return false;

            mpz_class z;
            mpz_import(z.get_mpz_t(), 1, 1, sizeof(p), 0, 0, &p);

            int r = mpz_probab_prime_p(z.get_mpz_t(), 25);
            return r!=0;
        }
        static mpz_class from_uint64(uint64_t x) {
            mpz_class z;
            mpz_import(z.get_mpz_t(), 1, 1, sizeof(x), 0, 0, &x);
            return z;
        } 

        static uint64_t checked_prime(uint64_t prime) {
            if (!is_prime_u64(prime)) {
                throw std::invalid_argument("p_int prime must be prime");
            }

            return prime;
        }

        static mpz_class make_modulus(uint64_t prime, uint64_t precision) {
            mpz_class mod;
            mpz_ui_pow_ui(mod.get_mpz_t(), prime, precision);
            return mod;
        }

        static mpz_class from_int64(int64_t x) {
            if (x>= 0) {
                return from_uint64(static_cast<uint64_t>(x));
            }

            //int64_minとかでもウゴクようにするため
            uint64_t mag = static_cast<uint64_t>(-(x+1)) + 1;
            return -from_uint64(mag);
        }

        template <std::integral T>
        static mpz_class to_mpz(T x) {
            if constexpr (std::is_signed_v<T>) {
                return from_int64(static_cast<int64_t>(x));
            } else {
                return from_uint64(static_cast<uint64_t>(x));
            }
        }

        void normalize() {
            int_p_ %= modulus_;
            if (int_p_<0) int_p_ += modulus_;
        }

        void reduce_precision_for_binary_op(const p_int& rhs) {
            if (rhs.precision_ >= precision_) return;

            precision_ = rhs.precision_;
            modulus_ = rhs.modulus_;
        }
    public:
        //変数宣言

        p_int(unchecked_prime_t, mpz_class int_p, uint64_t prime, uint64_t precision)
            : prime_(prime),
              precision_(precision),
              modulus_(make_modulus(prime_, precision_)),
              int_p_(std::move(int_p)) {
            normalize();
        }

        p_int(mpz_class int_p, uint64_t prime, uint64_t precision)
            : p_int(unchecked_prime, std::move(int_p), checked_prime(prime), precision) {}

        template <std::integral T>
        p_int(T int_p, uint64_t prime, uint64_t precision) : p_int(to_mpz(int_p), prime, precision) {}

        template <std::integral T>
        p_int(unchecked_prime_t tag, T int_p, uint64_t prime, uint64_t precision)
            : p_int(tag, to_mpz(int_p), prime, precision) {}
        //足し算
        p_int& operator+=(const p_int& rhs) {
            if (prime_ != rhs.prime_) {
                throw std::invalid_argument("cannot add p_int with different primes");
            }

            reduce_precision_for_binary_op(rhs);
            int_p_ += rhs.int_p_;
            normalize();

            return *this;
        }

        friend p_int operator+(p_int lhs, const p_int& rhs) {
            lhs += rhs;
            return lhs;
        }

        //引き算
        p_int& operator-=(const p_int& rhs) {
            if (prime_ != rhs.prime_) {
                throw std::invalid_argument("cannt subtract p_int with different primes");
            }

            reduce_precision_for_binary_op(rhs);
            int_p_ -= rhs.int_p_;
            normalize();

            return *this;
        }

        friend p_int operator-(p_int lhs, const p_int& rhs) {
            lhs -= rhs;
            return lhs;
        }

        p_int operator-() const {
            p_int res(*this);
            res.int_p_ = -res.int_p_;
            res.normalize();
            return res;
        }
        //掛け算
        p_int& operator*=(const p_int& rhs) {
            if (prime_ != rhs.prime_) {
                throw std::invalid_argument("cannot multiply p_int with different primes");
            }

            reduce_precision_for_binary_op(rhs);
            int_p_ *= rhs.int_p_;
            normalize();

            return *this;
        }

        friend p_int operator*(p_int lhs, const p_int& rhs) {
            lhs *= rhs;
            return lhs;
        }

        //確認用整数累乗

        p_int pow_uint(uint64_t n) const {
            p_int base(*this);
            p_int result(unchecked_prime, 1, prime_, precision_);

            while (n>0) {
                if (n&1) result *= base;
                base *= base;
                n>>=1;
            }
            return result;
        }

};
}
