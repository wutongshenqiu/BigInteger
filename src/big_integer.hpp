//
// Created by qiufeng on 2020/12/7.
//

#ifndef SRC_BIG_INTEGER_HPP
#define SRC_BIG_INTEGER_HPP

#endif //SRC_BIG_INTEGER_HPP


// 最快的大整数乘法是一个没有解决的问题
//      1. Complex multiplication algorithm
//      2. Karatsuba multiplication
//      3. Toom–Cook
//      4. Fourier transform methods

// 可以参考的文章和代码
// 1. https://gist.github.com/ar-pa/957297fb3f88996ead11
// 2. https://www.codementor.io/@arpitbhayani/how-python-implements-super-long-integers-12icwon5vk
// 3. https://math.stackexchange.com/questions/111150/changing-a-number-between-arbitrary-bases


#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include <sstream>


class BigInteger {
private:
    constexpr static uint32_t KShift = 30;
    constexpr static uint32_t KBase = 1 << KShift;
    constexpr static uint32_t KMask = KBase - 1;
    // use for output
    constexpr static uint32_t KDecimalBase = 1e9;
    constexpr static uint32_t KDecimalShift = 9;
    // use for multiplication
    constexpr static uint32_t KMulShift = 20;
    constexpr static uint32_t KMulBase = 1 << KMulShift;
    constexpr static uint32_t KMulMask = KMulBase - 1;
    // use for miller rabin test
    constexpr static std::array<uint32_t, 9> KMillerRabinTestBase = {
        2, 3, 5, 7, 11, 13, 17, 19, 23,
    };

private:
    // 符号标识
    // 1 代表正数
    // -1 代表负数
    int sign;
    std::vector<uint32_t> digits;

public:
    // constructors:
    BigInteger();
    // TODO
    // usage of explicit?
    explicit BigInteger(int64_t v);
    BigInteger(const BigInteger &v);
    explicit BigInteger(const std::string &s);
    BigInteger(const std::vector<uint32_t> &d, uint32_t base, int sign);
    // move constructors
    BigInteger(BigInteger &&v) noexcept;

    // Assignment operators:
    BigInteger &operator=(int64_t v);
    BigInteger &operator=(const std::string &s);
    // TODO
    // default usage?
    BigInteger &operator=(const BigInteger &v);
    // move assign
    BigInteger &operator=(BigInteger &&v) noexcept;

    // Unary arithmetic operators:
    BigInteger operator+() const;
    BigInteger operator-() const;

    // Binary arithmetic operators:
    BigInteger operator+(const BigInteger &v) const;
    BigInteger operator-(const BigInteger &v) const;
    BigInteger operator*(const BigInteger &v) const;
    BigInteger operator/(const BigInteger &v) const;
    BigInteger operator%(const BigInteger &v) const;
    BigInteger operator+(const int64_t &v) const;
    BigInteger operator-(const int64_t &v) const;
    BigInteger operator*(const int64_t &v) const;
    // the operator `/` here is different from python's operator `//`
    BigInteger operator/(const int32_t &v) const;
    int32_t operator%(const int32_t &v) const;
    BigInteger operator+(const std::string &v) const;
    BigInteger operator-(const std::string &v) const;
    BigInteger operator*(const std::string &v) const;
    BigInteger operator/(const std::string &v) const;
    BigInteger operator%(const std::string &v) const;

    // Arithmetic-assignment operators:
    BigInteger &operator+=(const BigInteger &v);
    BigInteger &operator-=(const BigInteger &v);
    BigInteger &operator*=(const BigInteger &v);
    BigInteger &operator/=(const BigInteger &v);
    BigInteger &operator%=(const BigInteger &v);
    BigInteger &operator+=(const int64_t &v);
    BigInteger &operator-=(const int64_t &v);
    BigInteger &operator*=(const int64_t &v);
    BigInteger &operator/=(const int32_t &v);
    BigInteger &operator+=(const std::string &v);
    BigInteger &operator-=(const std::string &v);
    BigInteger &operator*=(const std::string &v);
    BigInteger &operator/=(const std::string &v);
    BigInteger &operator%=(const std::string &v);

    // Relational operators:
    bool operator>(const BigInteger &v) const;
    bool operator>=(const BigInteger &v) const;
    bool operator<(const BigInteger &v) const;
    bool operator<=(const BigInteger &v) const;
    bool operator==(const BigInteger &v) const;
    bool operator!=(const BigInteger &v) const;
    bool operator>(const int64_t &v) const;
    bool operator>=(const int64_t &v) const;
    bool operator<(const int64_t &v) const;
    bool operator<=(const int64_t &v) const;
    bool operator==(const int64_t &v) const;
    bool operator!=(const int64_t &v) const;
    bool operator>(const std::string &v) const;
    bool operator>=(const std::string &v) const;
    bool operator<(const std::string &v) const;
    bool operator<=(const std::string &v) const;
    bool operator==(const std::string &v) const;
    bool operator!=(const std::string &v) const;

    // I/O stream operators:
    friend std::ostream &operator<<(std::ostream &os, const BigInteger &v);

public:
    // get absolute value
    [[nodiscard]] BigInteger Abs() const;

    // convert integer to string with base 10
    [[nodiscard]] std::string ToString() const;

    // convert integer to string with base 10
    [[nodiscard]] std::string ToHexString() const;

    // load from string with base10
    void FromString(const std::string &s);

    [[nodiscard]] inline bool IsZero() const {
        return digits.empty() || (digits.size() == 1 && !digits[0]);
    }

    [[nodiscard]] inline bool IsEven() const {
        return digits.empty() || (!(digits[0] & 1));
    }

public:
    // ascii char to digit
    inline static uint32_t CharToDigit(char c) {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;

        return 0;
    }

    // digit to hex char
    inline static char DigitToHex(uint32_t t) {
        if (t >= 0 && t <= 9) return static_cast<char>(t + '0');
        else if (t <= 15) return static_cast<char>(t - 10 + 'A');
        else return 0;
    }
    // add two vectors
    template<class T>
    static std::vector<T> QuickAdd(const std::vector<T> &a, const std::vector<T> &b, uint32_t shift) {
        if (a.size() < b.size()) return QuickAdd(b, a, shift);

        std::vector<T> res(a.size() + 1);
        T mask = (1 << shift) - 1;

        T carry = 0;
        for (auto i = 0; i < b.size(); i++) {
            carry += a[i] + b[i];
            res[i] = carry & mask;
            carry >>= shift;
        }

        for (auto i = b.size(); i < a.size(); i++) {
            carry += a[i];
            res[i] = carry & mask;
            carry >>= shift;
        }
        if (carry) res[res.size()-1] = carry;
        else TrimZero(res);

        return res;
    }
    // Convert `old_digits` with `old_base` to `new_base`
    //
    // Args:
    //      old_digits: a vector that stores digit in descending order
    //
    // Returns:
    //      a vector that stores digit with base `new_base` in descending order
    static std::vector<uint32_t> ConvertBase(const std::vector<uint32_t> &old_digits,
                                             uint32_t old_base,
                                             uint32_t new_base);

    // `old_base` must equal to 2^old_shift
    static std::vector<uint32_t> QuickConvertBase(const std::vector<uint32_t> &old_digits,
                                                  uint32_t old_shift,
                                                  uint32_t new_base);

    // pop zero in the end of vector
    template<class T>
    inline static void TrimZero(std::vector<T> &a) {
        while (!a.empty() && !a.back()) a.pop_back();
    }

    // below are multiplication algorithms
    // see https://en.wikipedia.org/wiki/Multiplication_algorithm

    // note:
    //      1. before use the result, call TrimZero
    inline static std::vector<uint64_t> GridMultiplication(const std::vector<uint64_t> &a,
                                                           const std::vector<uint64_t> &b) {
        std::vector<uint64_t> res(a.size() + b.size(), 0);

        for (auto i = 0; i < a.size(); i++) {
            for (auto j = 0; j < b.size(); j++) {
                res[i + j] += a[i] * b[j];
            }
        }

        return res;
    }

    inline static std::vector<uint32_t> GridMultiplication(const std::vector<uint32_t> &a,
                                                           const std::vector<uint32_t> &b) {
        std::vector<uint32_t> res(a.size() + b.size(), 0);

        for (auto i = 0; i < a.size(); i++) {
            uint64_t carry = 0;
            auto pres = res.begin() + i;
            auto pb = b.begin();
            auto pbend = b.end();

            while (pb != pbend) {
                carry += *pres + *pb++ * static_cast<uint64_t>(a[i]);
                *pres++ = static_cast<uint32_t>(carry & KMask);
                carry >>= KShift;
            }
            if (carry) {
                *pres += static_cast<uint32_t>(carry & KMask);
            }
        }

        return res;
    }



    // note:
    //      1. we assume vector `a` and `b` have the same size
    //      2. size of `a` better to be power of 2
    //      3. before use the result, call TrimZero
    // for more information, see https://zhuanlan.zhihu.com/p/144813558
    static std::vector<uint64_t> KaratsubaMultiplication(const std::vector<uint64_t> &a,
                                                         const std::vector<uint64_t> &b);

    // (quotient, remainder)
    static std::pair<BigInteger, BigInteger> Divmod(const BigInteger &a1, const BigInteger &b1);

    // greatest common divisor
    static BigInteger Gcd(const BigInteger &a, const BigInteger &b);

    // extended euclidean algorithm
    // give `a`, `b`;
    // output `s`, `t`, `g` s.t. sa + tb = (a, b) = g
    static std::tuple<BigInteger, BigInteger, BigInteger> ExtendedEuclidean(const BigInteger &a,
                                                                            const BigInteger &b);

    // calculate `k` that satisfies b^e = k (mod m)
    static BigInteger ModularExponentiation(const BigInteger &b, const BigInteger &e, const BigInteger &m);

    // miller-rabin primality test
    static bool MillerRabinTest(const BigInteger &v);

//    static bool IsPrime(const BigInteger &v);
};