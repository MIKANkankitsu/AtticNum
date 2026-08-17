#include <iostream>
#include <cassert>
#include <gmpxx.h>

// テスト対象の関数
static void normalize_mod(mpz_class& value, const mpz_class& modulus) {
    mpz_mod(value.get_mpz_t(), value.get_mpz_t(), modulus.get_mpz_t());
    // value %= modulus;
    // if (value < 0) value += modulus;
}

int main() {
    const mpz_class mod = 5;

    // パターン 1: 正の数（そのまま余りが取れるか）
    {
        mpz_class val = 7;
        normalize_mod(val, mod);
        assert(val == 2); // 7 % 5 = 2
    }

    // パターン 2: 負の数（正の余りに補正されるか）★ここが一番重要
    {
        mpz_class val = -7;
        normalize_mod(val, mod);
        assert(val == 3); // -7 ≡ 3 (mod 5)
    }

    // パターン 3: ちょうど割り切れる負の数（-0 にならず 0 になるか）
    {
        mpz_class val = -10;
        normalize_mod(val, mod);
        assert(val == 0); // -10 % 5 = 0, if (0 < 0) は走らない
    }

    // パターン 4: モジュラスより小さい負の数
    {
        mpz_class val = -2;
        normalize_mod(val, mod);
        assert(val == 3); // -2 ≡ 3 (mod 5)
    }

    std::cout << "All normalize_mod tests passed successfully!" << std::endl;
    std::cout << (12345 >> 1) << std::endl;
    return 0;
}