#pragma once
#include <complex>
#include <cstddef>
#include <vector>

namespace hw1 {
using Complex = std::complex<double>;
using Vector = std::vector<Complex>;
inline constexpr double pi = 3.141592653589793238462643383279502884;

// Negative exponent, no normalization. Any positive input length.
Vector fft(Vector input);
// Positive exponent, normalized by 1/N. Any positive input length.
Vector ifft(Vector input);
} // namespace hw1
