// #include <iostream>
// #include <cassert>
// #include <gmpxx.h>

// // テスト対象の関数
// static void normalize_mod(mpz_class& value, const mpz_class& modulus) {
//     mpz_mod(value.get_mpz_t(), value.get_mpz_t(), modulus.get_mpz_t());
//     // value %= modulus;
//     // if (value < 0) value += modulus;
// }

// int main() {
//     const mpz_class mod = 5;

//     // パターン 1: 正の数（そのまま余りが取れるか）
//     {
//         mpz_class val = 7;
//         normalize_mod(val, mod);
//         assert(val == 2); // 7 % 5 = 2
//     }

//     // パターン 2: 負の数（正の余りに補正されるか）★ここが一番重要
//     {
//         mpz_class val = -7;
//         normalize_mod(val, mod);
//         assert(val == 3); // -7 ≡ 3 (mod 5)
//     }

//     // パターン 3: ちょうど割り切れる負の数（-0 にならず 0 になるか）
//     {
//         mpz_class val = -10;
//         normalize_mod(val, mod);
//         assert(val == 0); // -10 % 5 = 0, if (0 < 0) は走らない
//     }

//     // パターン 4: モジュラスより小さい負の数
//     {
//         mpz_class val = -2;
//         normalize_mod(val, mod);
//         assert(val == 3); // -2 ≡ 3 (mod 5)
//     }

//     std::cout << "All normalize_mod tests passed successfully!" << std::endl;
//     std::cout << (12345 >> 1) << std::endl;
//     return 0;
// }

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <gmp.h>
#include <gmpxx.h>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <iostream>

struct Tracker {
    Tracker() { std::cout << "Default\n"; }
    Tracker(const Tracker&) { std::cout << "Copy\n"; }
    Tracker(Tracker&&) noexcept { std::cout << "Move\n"; }
};

void foo(Tracker t) {
    Tracker b = std::move(t);
}

// int main() {
//     Tracker a;
//     foo(a); // 何回 Copy / Move が出るか観察する！
// }

class Foo {
private:
    long secret_data_; // private なメンバ変数

public:
    explicit Foo(int val) : secret_data_(val) {}

    // 別インスタンス (rhs) の private メンバ secret_data_ に直接アクセスできるか？
    Foo& operator+=(const Foo& rhs) {
        std::cout << "my secret_data_: " << this->secret_data_ << "\n";
        std::cout << "your's (rhs) secret_data_: " << rhs.secret_data_ << "\n";

        // ★ ここがコンパイル通るか・アクセスできるかの実験！
        this->secret_data_ += rhs.secret_data_;

        return *this;
    }

    int get_val() const { return secret_data_; }
};

// int main() {
//     Foo a(10);
//     Foo b(20);

//     a += b; // 実行！

//     std::cout << "after added a: " << a.get_val() << "\n";
//     return 0;
// }

#include <cstdint>

// sizeof(unsigned long) が 8 バイト未満ならコンパイルエラーにしてビルドを中断
static_assert(sizeof(unsigned long) >= 8, 
    "unsigned long int must be at least 64-bit (8 bytes). "
    "Windows (LLP64) is not directly supported without 64-bit truncation helpers!");

// int main() {
//     std::cout << sizeof(unsigned long int) << " " << sizeof(uint64_t) << std::endl;
// }

int main() {
    mpz_class a = 12346;

    std::cout << mpz_divisible_ui_p(a.get_mpz_t(),(uint64_t)5) << std::endl;
    // mpz_divisible_ui_pは、割り切れるなら0以外,そうでないなら0を返す

    return 0;
}