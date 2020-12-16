//
// Created by qiufeng on 2020/12/7.
//

#include "big_integer.hpp"


BigInteger::BigInteger() : sign(1) {}

BigInteger::BigInteger(int64_t v)  {
    *this = v;
}

BigInteger::BigInteger(const BigInteger &v)  {
    *this = v;
}

BigInteger::BigInteger(const std::string &s) {
    FromString(s);
}

BigInteger::BigInteger(const std::vector<uint32_t> &d, uint32_t base, int sign) {
    this->sign = sign > 0 ? 1 : -1;
    digits = ConvertBase(d, base, KBase);
}

BigInteger & BigInteger::operator=(int64_t v)  {
    digits.clear();
    sign = 1;

    if (v < 0) {
        sign = -1, v = -v;
    }
    for (; v > 0; v /= KBase) {
        digits.push_back(v & KMask);
    }

    return *this;
}

BigInteger & BigInteger::operator=(const std::string &s) {
    FromString(s);

    return *this;
}

BigInteger & BigInteger::operator=(const BigInteger &v) = default;

std::string BigInteger::ToString() const {
    std::stringstream ss;
    ss << *this;

    std::string s;
    ss >> s;

    return s;
}

std::ostream& operator<<(std::ostream& os, const BigInteger& v) {
    if (v.digits.empty()) {
        os << '0';
        return os;
    }
    if (v.sign == -1) os << '-';

    auto digits_base_decimal = BigInteger::QuickConvertBase(v.digits, BigInteger::KShift, BigInteger::KDecimalBase);
    os << digits_base_decimal.back();

    for (int i = static_cast<int>(digits_base_decimal.size()) - 2; i >= 0; --i) {
        os << std::setw(BigInteger::KDecimalShift) << std::setfill('0') << digits_base_decimal[i];
    }

    return os;
}

void BigInteger::FromString(const std::string &s) {
    digits.clear();
    sign = 1;

    size_t pos = 0;
    // get sign of the string
    if (!s.empty()) {
        if (s[pos] == '-') {
            sign = -sign;
            ++pos;
        } else if (s[pos] == '+') {
            ++pos;
        }
    }
    uint32_t base = 10;
    // get base of the string
    if (pos < s.size() - 1) {
        if (s[pos] == '0' && s[pos+1] == 'x') {
            base = 16;
            pos += 2;
        }
    }

    std::vector<uint32_t> str_digits(s.size() - pos, 0);
    std::transform(s.rbegin(), s.rend() - pos, str_digits.begin(), CharToDigit);

    digits = ConvertBase(str_digits, base, KBase);
}

std::vector<uint32_t> BigInteger::ConvertBase(const std::vector<uint32_t> &old_digits, uint32_t old_base, uint32_t new_base)  {
    // TODO
    // this algorithm are copied from `long_to_decimal_string_internal()` in `Objects/longobject.c`
    // I haven't fully understood it yet

    // TODO
    // preallocate for efficiency
    std::vector<uint32_t> new_digits;
    for (int32_t i = old_digits.size(); --i >= 0;) {
        uint32_t hi = old_digits[i];
        for (uint32_t & new_digit : new_digits) {
            uint64_t z = static_cast<uint64_t>(new_digit) * old_base + hi;
            hi = static_cast<uint32_t>(z / new_base);
            new_digit = static_cast<uint32_t>(z - static_cast<uint64_t>(hi) * new_base);
        }
        while (hi) {
            new_digits.push_back(hi % new_base);
            hi /= new_base;
        }
    }

    return new_digits;
}

std::vector<uint32_t> BigInteger::QuickConvertBase(const std::vector<uint32_t> &old_digits,
                                                   uint32_t old_shift,
                                                   uint32_t new_base) {
    // TODO
    // preallocate for efficiency
    std::vector<uint32_t> new_digits;
    for (int32_t i = old_digits.size(); --i >= 0;) {
        uint32_t hi = old_digits[i];
        for (uint32_t & new_digit : new_digits) {
            uint64_t z = static_cast<uint64_t>(new_digit) << old_shift | hi;
            hi = static_cast<uint32_t>(z / new_base);
            new_digit = static_cast<uint32_t>(z - static_cast<uint64_t>(hi) * new_base);
        }
        while (hi) {
            new_digits.push_back(hi % new_base);
            hi /= new_base;
        }
    }

    return new_digits;
}

BigInteger BigInteger::operator+() const {
    BigInteger res = *this;

    return res;
}

BigInteger BigInteger::operator-() const {
    BigInteger res = *this;
    res.sign = -sign;

    return res;
}

BigInteger BigInteger::operator+(const BigInteger &v) const {
    if (sign == v.sign) {
        // choose longer one may avoid allocation in vector
        BigInteger res;
        res.sign = sign;
        res.digits = QuickAdd(this->digits, v.digits, KShift);

        return res;
    }

    return *this - (-v);
}

BigInteger BigInteger::operator-(const BigInteger &v) const {
    if (sign == v.sign) {
        // TODO
        // more efficiency
        if (this->Abs() < v.Abs()) return -(v - *this);

        uint32_t borrow = 0;
        BigInteger res = *this;

        for (auto i = 0; i < v.digits.size(); i++) {
            borrow = digits[i] - v.digits[i] - borrow;
            res.digits[i] = borrow & KMask;
            // if borrow is negative(signed int), than borrow >> 30 will always be 3
            borrow >>= KShift;
            borrow &= 1;
        }

        for (auto i = v.digits.size(); i < digits.size(); i++) {
            borrow = digits[i] - borrow;
            res.digits[i] = borrow & KMask;
            borrow >>= KShift;
            borrow &= 1;
        }
        TrimZero(res.digits);

        return res;
    }

    return *this + (-v);
}

BigInteger BigInteger::Abs() const {
    BigInteger res = *this;
    res.sign = std::abs(res.sign);

    return res;
}

bool BigInteger::operator<(const BigInteger &v) const {
    if (sign != v.sign) return sign < v.sign;
    if (digits.size() != v.digits.size()) return sign * digits.size() < v.sign * v.digits.size();

    for (int32_t i = digits.size() - 1; i >= 0; i--) {
        if (digits[i] != v.digits[i]) return sign * digits[i] < v.sign * v.digits[i];
    }

    return false;
}

bool BigInteger::operator>(const BigInteger &v) const {
    return v < *this;
}

bool BigInteger::operator==(const BigInteger &v) const {
    if (!(digits.size() | v.digits.size())) return true;
    if (sign != v.sign) return false;
    if (digits.size() != v.digits.size()) return false;

    for (auto i = 0; i < digits.size(); i++) {
        if (digits[i] != v.digits[i]) return false;
    }

    return true;
}

bool BigInteger::operator!=(const BigInteger &v) const {
    return !(*this == v);
}

bool BigInteger::operator<=(const BigInteger &v) const {
    return !(*this > v);
}

bool BigInteger::operator>=(const BigInteger &v) const {
    return !(*this < v);
}

BigInteger BigInteger::operator*(const int64_t &v) const {
    BigInteger res = *this;

    if (v < 0) {
        res.sign = -res.sign;
        return res * -v;
    }
    if (v > KBase) {
        // TODO
        // * to >>
        return res * (v >> KShift) * KBase + res * (v & KMask);
    }
    uint64_t carry = 0;
    for (uint32_t & digit : res.digits) {
        carry = digit * v + carry;
        digit = carry & KMask;
        carry >>= KShift;
    }
    if (carry) res.digits.push_back(carry);

    TrimZero(res.digits);

    return res;
}

BigInteger & BigInteger::operator+=(const BigInteger &v) {
    *this = *this + v;

    return *this;
}

BigInteger & BigInteger::operator-=(const BigInteger &v) {
    *this = *this - v;

    return *this;
}

BigInteger & BigInteger::operator*=(const int64_t &v) {
    *this = *this * v;

    return *this;
}

std::vector<uint64_t> BigInteger::KaratsubaMultiplication(const std::vector<uint64_t> &a,
                                                          const std::vector<uint64_t> &b) {
    auto n = a.size();
    std::vector<uint64_t> res(n << 1);
    if (n <= 32) {
        return GridMultiplication(a, b);
    }

    auto k = n >> 1;
    std::vector<uint64_t> a1(a.begin(), a.begin() + k);
    std::vector<uint64_t> a2(a.begin() + k, a.end());
    std::vector<uint64_t> b1(b.begin(), b.begin() + k);
    std::vector<uint64_t> b2(b.begin() + k, b.end());

    auto a1b1 = KaratsubaMultiplication(a1, b1);
    auto a2b2 = KaratsubaMultiplication(a2, b2);

    for (auto i = 0; i < k; i++) {
        a2[i] += a1[i];
        b2[i] += b1[i];
    }

    std::vector<uint64_t> r = KaratsubaMultiplication(a2, b2);
    for (auto i = 0; i < a1b1.size(); i++) r[i] -= a1b1[i];
    for (auto i = 0; i < a2b2.size(); i++) r[i] -= a2b2[i];

    for (auto i = 0; i < r.size(); i++) res[i + k] += r[i];
    for (auto i = 0; i < a1b1.size(); i++) res[i] += a1b1[i];
    for (auto i = 0; i < a2b2.size(); i++) res[i + n] += a2b2[i];

    return res;
}

BigInteger BigInteger::operator*(const BigInteger &v) const {
    // change base to 2^20 to prevent overflow
    std::vector<uint32_t> a20 = QuickConvertBase(digits, KShift, KMulBase);
    std::vector<uint32_t> b20 = QuickConvertBase(v.digits, KShift, KMulBase);

    std::vector<uint64_t> a(a20.begin(), a20.end());
    std::vector<uint64_t> b(b20.begin(), b20.end());

    // make a have the same size as b
    while (a.size() < b.size())
        a.push_back(0);
    while (b.size() < a.size())
        b.push_back(0);
    // make size of a is power of 2
    while (a.size() & (a.size() - 1))
        a.push_back(0), b.push_back(0);

    auto c = KaratsubaMultiplication(a, b);
    std::vector<uint32_t> res_digits(c.size(), 0);
    uint64_t carry = 0;
    for (auto i = 0; i < c.size(); i++) {
        carry = carry + c[i];
        res_digits[i] = static_cast<uint32_t>(carry & KMulMask);
        carry >>= KMulShift;
    }
    TrimZero(res_digits);

    BigInteger res;
    res.sign = sign * v.sign;
    res.digits = QuickConvertBase(res_digits, KMulShift, KBase);

    return res;
}

BigInteger & BigInteger::operator*=(const BigInteger &v) {
    *this = *this * v;

    return *this;
}

BigInteger BigInteger::operator/(const int32_t &v) const {
    BigInteger res = *this;
    if (v < 0) {
        res.sign = -res.sign;
        return res / -v;
    }

    uint64_t rem = 0;
    for (int i = static_cast<int>(res.digits.size() - 1); i >= 0; i--) {
        rem = res.digits[i] + (rem << KShift);
        res.digits[i] = static_cast<uint32_t>(rem / v);
        rem = rem % v;
    }
    TrimZero(res.digits);

    return res;
}

BigInteger & BigInteger::operator/=(const int32_t &v) {
    *this = *this / v;

    return *this;
}

int32_t BigInteger::operator%(const int32_t &v) const {
    if (v < 0) return *this % -v;

    int64_t m = 0;
    for (int i = static_cast<int>(digits.size() - 1); i >= 0; i--) {
        m = (digits[i] + (m << KShift)) % v;
    }

    return static_cast<int32_t>(m) * sign;
}

BigInteger BigInteger::operator+(const int64_t &v) const {
    return *this + BigInteger(v);
}

BigInteger & BigInteger::operator+=(const int64_t &v) {
    *this = *this + v;

    return *this;
}

BigInteger BigInteger::operator-(const int64_t &v) const {
    return *this - BigInteger(v);
}

BigInteger & BigInteger::operator-=(const int64_t &v) {
    *this = *this - v;

    return *this;
}

bool BigInteger::operator<=(const int64_t &v) const {
    return *this <= BigInteger(v);
}

bool BigInteger::operator>=(const int64_t &v) const {
    return *this >= BigInteger(v);
}

bool BigInteger::operator==(const int64_t &v) const {
    return *this == BigInteger(v);
}

bool BigInteger::operator>(const int64_t &v) const {
    return *this > BigInteger(v);
}

bool BigInteger::operator<(const int64_t &v) const {
    return *this < BigInteger(v);
}

bool BigInteger::operator!=(const int64_t &v) const {
    return *this != BigInteger(v);
}

// TODO
// I haven't fully understood the algorithm
// information can be found in [The Art of Computer Programming, Vol. 2 (3rd
// edn.), section 4.3.1, Algorithm D]
std::pair<BigInteger, BigInteger> BigInteger::Divmod(const BigInteger &a1, const BigInteger &b1) {
    int32_t norm = KBase / (b1.digits.back() + 1);
    auto a = a1.Abs() * norm;
    auto b = b1.Abs() * norm;
    BigInteger q, r;
    q.digits.resize(a.digits.size());

    for (int i = static_cast<int>(a.digits.size()) - 1; i >= 0; i--) {
        r *= KBase;
        r += a.digits[i];

        uint32_t s1 = r.digits.size() <= b.digits.size() ? 0 : r.digits[b.digits.size()];
        uint32_t s2 = r.digits.size() <= b.digits.size() - 1 ? 0 : r.digits[b.digits.size() - 1];

        auto d = static_cast<uint32_t>(((static_cast<uint64_t>(s1) << KShift) + s2) / b.digits.back());
        r -= b * d;

        while (r < 0) {
            r += b;
            --d;
        }
        q.digits[i] = d;
    }

    q.sign = a1.sign * b1.sign;
    r.sign = a1.sign;
    TrimZero(q.digits);
    TrimZero(r.digits);

    return std::make_pair(q, r / norm);
}

BigInteger BigInteger::operator/(const BigInteger &v) const {
    return Divmod(*this, v).first;
}

BigInteger & BigInteger::operator/=(const BigInteger &v) {
    *this = *this / v;

    return *this;
}

BigInteger BigInteger::operator%(const BigInteger &v) const {
    return Divmod(*this, v).second;
}

BigInteger & BigInteger::operator%=(const BigInteger &v) {
    *this = (*this % v);

    return *this;
}

BigInteger BigInteger::operator+(const std::string &v) const {
    return *this + BigInteger(v);
}

BigInteger BigInteger::operator-(const std::string &v) const {
    return *this - BigInteger(v);
}

BigInteger BigInteger::operator*(const std::string &v) const {
    return *this * BigInteger(v);
}

BigInteger BigInteger::operator/(const std::string &v) const {
    return *this / BigInteger(v);
}

BigInteger BigInteger::operator%(const std::string &v) const {
    return *this % BigInteger(v);
}

BigInteger & BigInteger::operator+=(const std::string &v) {
    return *this += BigInteger(v);
}

BigInteger & BigInteger::operator-=(const std::string &v) {
    return *this -= BigInteger(v);
}

BigInteger & BigInteger::operator*=(const std::string &v) {
    return *this *= BigInteger(v);
}

BigInteger & BigInteger::operator/=(const std::string &v) {
    return *this /= BigInteger(v);
}

BigInteger & BigInteger::operator%=(const std::string &v) {
    return *this %= BigInteger(v);
}

bool BigInteger::operator>(const std::string &v) const {
    return *this > BigInteger(v);
}

bool BigInteger::operator>=(const std::string &v) const {
    return *this >= BigInteger(v);
}

bool BigInteger::operator<(const std::string &v) const {
    return *this < BigInteger(v);
}

bool BigInteger::operator<=(const std::string &v) const {
    return *this <= BigInteger(v);
}

bool BigInteger::operator==(const std::string &v) const {
    return *this == BigInteger(v);
}

bool BigInteger::operator!=(const std::string &v) const {
    return *this != BigInteger(v);
}

BigInteger BigInteger::Gcd(const BigInteger &a, const BigInteger &b) {
    return b.IsZero() ? a : Gcd(b, a % b);
}

std::tuple<BigInteger, BigInteger, BigInteger> BigInteger::ExtendedEuclidean(const BigInteger &a,
                                                                             const BigInteger &b) {
    BigInteger a0 = a;
    BigInteger b0 = b;
    std::vector<std::pair<BigInteger, BigInteger>> st;

    while (!b0.IsZero()) {
        st.emplace_back(a0, b0);
        auto tmp = b0;
        b0 = a0 % b0;
        a0 = tmp;
    }
    // gcd(a, b)
    BigInteger g = st.back().second;

    BigInteger s(1);
    BigInteger t(0);

    while (!st.empty()) {
        auto [a0, b0] = st.back();
        st.pop_back();

        auto tmp = s;
        s = t;
        t = tmp - a0 / b0 * t;
    }

    return std::make_tuple(s, t, g);
}

BigInteger BigInteger::ModularExponentiation(const BigInteger &b, const BigInteger &e, const BigInteger &m) {
    // change to base `2`
    auto e2 = QuickConvertBase(e.digits, KShift, 2);

    // precalculate b^m s.t. m = 1, 2, ... 2^(e2.size()-1)
    std::vector<BigInteger> b2({ b });
    b2.resize(e2.size());
    for (auto i = 1; i < e2.size(); i++) {
        b2[i] = (b2[i-1] * b2[i-1]) % m;
    }

    // result
    BigInteger k(1);
    for (auto i = 0; i < e2.size(); i++) {
        if (e2[i]) {
            k *= b2[i];
            k %= m;
        }
    }

    return k;
}

BigInteger::BigInteger(BigInteger &&v) noexcept {
    this->sign = v.sign;
    this->digits = v.digits;
}

BigInteger & BigInteger::operator=(BigInteger &&v) noexcept {
    this->sign = v.sign;
    this->digits = v.digits;

    return *this;
}

bool BigInteger::MillerRabinTest(const BigInteger &v) {
    if (v < 0) return MillerRabinTest(-v);

    if (v < 3) return v == 2;
    if (v.IsEven()) return false;

    auto n = v - 1;
    // n = 2^s*t
    uint32_t s = 0;
    while (n.IsEven()) {
        s++;
        n /= 2;
    }

    for (auto &a : KMillerRabinTestBase) {
        if (v == a) return true;
        // a^n = m (mod v)
        auto m = ModularExponentiation(BigInteger(a), n, v);

        if (m == 1 || m == v - 1) continue;

        bool is_break = false;
        for (auto j = 0; j < s; j++) {
            auto k = m * m % v;
            if (k == 1) return false;
            if (k == v - 1) {
                is_break = true;
                break;
            }
        }
        if (!is_break) return false;
    }

    return true;
}
