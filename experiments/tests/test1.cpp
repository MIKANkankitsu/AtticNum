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

int main() {
    std::cout << sizeof(unsigned long int) << " " << sizeof(uint64_t) << std::endl;
}