#pragma once

#include <boost/multiprecision/gmp.hpp>
#include <cstdint>
#include <stdexcept>
#include <utility>

//共通の整数仕様
namespace attic {
    namespace mp = boost::multiprecision;

   using bigint = mp::mpz_int;
}

//実数
namespace attic::real {
    class Real {
        private:
            attic::bigint A_;
            attic::bigint P_;

            static attic::bigint pow10(std::int64_t n) {
                if (n<0) throw std::invalid_argument("pow10: negative exponent");

                attic::bigint r = 1;
                attic::bigint b = 10;

                while (n>0) {
                    if (n&1) r*=b;
                    b *= b;
                    n>>=1;
                }

                return r;
            }

            static attic::bigint rescale_down_trunc(
                    const attic::bigint& A,
                    std::int64_t oldP,
                    std::int64_t newP
                    ) {
                if (newP > oldP) return A * pow10(newP - oldP);

                return A / pow10(oldP - newP);
            }
    };
}
