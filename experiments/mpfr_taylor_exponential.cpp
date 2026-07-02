#include <iostream>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <iomanip>

//これはtaylorの定理を用いて愚直に指数関数の計算をするものです。

namespace mp = boost::multiprecision;

mp::mpfr_float make_error(int digits, int extra = 5) {
    return 1 / pow(mp::mpfr_float(10), digits + extra);
}

mp::mpfr_float exp_taylor_positive(const mp::mpfr_float& x, int digits) {
    mp::mpfr_float error = make_error(digits);

    int n = 0;

    mp::mpfr_float term = 1; // x^0 / 0!
    mp::mpfr_float sum = 0;

    mp::mpfr_float M = pow(mp::mpfr_float(3), ceil(x));
    mp::mpfr_float rem = M;

    while (rem > error) {
        sum += term;

        ++n;
        term *= x / n;
        rem *= x / n;
    }

    return sum;
}

mp::mpfr_float exp_taylor(const mp::mpfr_float& x, int digits) {
    if (x < 0) {
        return 1 / exp_taylor_positive(-x, digits);
    }

    return exp_taylor_positive(x, digits);
}

int main(void) {
    int digits;
    std::cin >> digits;
    
    int bits = static_cast<int>(digits * 3.32192809488736234787) + 64;
    mp::mpfr_float::default_precision(bits);

    mp::mpfr_float x;
    std::cin >> x;

    mp::mpfr_float y = exp_taylor(x, digits);

    std::cout << std::setprecision(digits) << y << '\n';

    return 0;
}
