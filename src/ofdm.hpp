#pragma once
#include "fft.hpp"
#include <array>
#include <cstdint>
#include <random>
#include <string>

namespace hw2 {
using hw1::Complex;
using hw1::Vector;
constexpr std::size_t N = 64, CP = 16, L = N + CP;
using Bits = std::vector<std::uint8_t>;
enum class Modulation { QPSK, QAM16 };
unsigned bits_per_symbol(Modulation m);
std::string name(Modulation m);
Bits random_bits(std::size_t count, std::uint64_t seed);
Complex map_symbol(const std::uint8_t* bits, Modulation m);
std::array<std::uint8_t,4> demap_symbol(Complex z, Modulation m);
Vector modulate(const Bits& bits, Modulation m);
Vector add_delay(const Vector& tx, std::size_t delay);
class Gaussian {
    std::mt19937_64 rng_;
    double uniform_open();
public:
    explicit Gaussian(std::uint64_t seed) : rng_(seed) {}
    Complex sample(double complex_variance);
};
struct NoiseStats { double power = 0, mean_i = 0, mean_q = 0, iq = 0; };
Vector awgn(const Vector& tx, double variance, std::uint64_t seed, NoiseStats* stats = nullptr);
struct Timing { std::size_t cp_start = 0; std::array<double,L> metric{}; };
Timing synchronize(const Vector& rx, std::size_t symbols);
std::vector<Timing> synchronize_sweep(const Vector& rx, const std::vector<std::size_t>& counts);
struct Decode { std::uint64_t errors = 0, compared_bits = 0; Vector first_symbols; };
Decode demodulate(const Vector& rx, const Bits& bits, Modulation m, std::size_t estimated_cp);
double mean_power(const Vector& x);
}
