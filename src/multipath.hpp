#pragma once
#include "ofdm.hpp"
namespace hw2::multipath {
Vector transmit(const Bits& bits, Modulation mode, std::size_t cp);
// A packet is one continuous stream; channel memory crosses symbol boundaries.
Vector filter(const Vector& tx, const Vector& h);
Vector rayleigh(std::size_t delay_step, std::uint64_t seed);
Vector response(const Vector& h);
Vector window_fft(const Vector& rx, std::size_t frame, std::size_t cp);
Vector zero_force(const Vector& y, const Vector& H);
std::uint64_t errors(const Vector& y, const Bits& bits, Modulation mode, std::size_t frame);
}
