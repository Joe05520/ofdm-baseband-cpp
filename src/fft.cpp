#include "fft.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace hw1 {
namespace {
bool power_of_two(std::size_t n) { return n && !(n & (n - 1)); }

// Iterative DIT radix-2 transform
void radix2(Vector& a, bool inverse) {
    const auto n = a.size();
    if (!power_of_two(n)) throw std::invalid_argument("radix2 requires power-of-two N");
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        auto bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n;) {
        for (std::size_t j = 0; j < len / 2; ++j) {
            const auto w = std::polar(1.0, (inverse ? 2.0 : -2.0) * pi *
                                           double(j) / double(len));
            for (std::size_t start = 0; start < n; start += len) {
                const auto u = a[start + j];
                const auto v = w * a[start + j + len / 2];
                a[start + j] = u + v;
                a[start + j + len / 2] = u - v;
            }
        }
        if (len == n) break; // Avoid size_t overflow after the final stage.
        len *= 2;
    }
    if (inverse) for (auto& z : a) z /= double(n);
}

} // namespace

// Forward FFT for any positive N. Bluestein reduces other lengths to a
// zero-padded convolution computed with our Radix-2 FFT.
Vector fft(Vector a) {
    const auto n = a.size();
    if (n == 0) throw std::invalid_argument("FFT input must be nonempty");
    if (power_of_two(n)) { radix2(a, false); return a; }
    if (n > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("FFT size overflow");
    std::size_t m = 1;
    while (m < 2 * n - 1) {
        if (m > std::numeric_limits<std::size_t>::max() / 2)
            throw std::length_error("FFT convolution size overflow");
        m *= 2;
    }
    Vector u(m), v(m), chirp(n);
    for (std::size_t j = 0; j < n; ++j) {
        const long double phase = std::fmod(static_cast<long double>(j) * j,
                                            2.0L * n) / n;
        chirp[j] = std::polar(1.0, -pi * double(phase));
        u[j] = a[j] * chirp[j];
        v[j] = std::conj(chirp[j]);
        if (j) v[m - j] = v[j];
    }
    radix2(u, false); radix2(v, false);
    for (std::size_t j = 0; j < m; ++j) u[j] *= v[j];
    radix2(u, true);
    for (std::size_t j = 0; j < n; ++j) a[j] = u[j] * chirp[j];
    return a;
}
Vector ifft(Vector a) {
    const auto n = a.size();
    if (n == 0) throw std::invalid_argument("IFFT input must be nonempty");
    if (power_of_two(n)) { radix2(a, true); return a; }
    for (auto& z : a) z = std::conj(z);
    a = fft(a);
    for (auto& z : a) z = std::conj(z) / double(n);
    return a;
}
} // namespace hw1
