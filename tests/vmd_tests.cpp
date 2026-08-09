#include <vmd/vmd.hpp>
#include "VMD.h"

#include <unsupported/Eigen/FFT>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

Eigen::VectorXd make_signal(std::size_t n) {
    Eigen::VectorXd value(static_cast<Eigen::Index>(n));
    for (std::size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i + 1) / static_cast<double>(n);
        value(static_cast<Eigen::Index>(i)) = std::cos(2.0 * kPi * 2.0 * t) +
            0.1 * std::cos(2.0 * kPi * 24.0 * t) +
            0.005 * std::cos(2.0 * kPi * 288.0 * t);
    }
    return value;
}

vmd::Options make_options(std::size_t k = 8) {
    vmd::Options value;
    value.mode_count = k;
    value.alpha = {50.0};
    value.tau = 0.0;
    value.initialization = vmd::Initialization::Uniform;
    value.tolerance = 1e-7;
    value.max_iterations = 500;
    value.sample_rate_hz = 1200.0;
    return value;
}

void test_reference_signal() {
    const Eigen::VectorXd input = make_signal(1200);
    const vmd::Result result = vmd::decompose(input, make_options());
    check(result.modes.rows() == 8 && result.modes.cols() == 1200, "mode dimensions are K x N");
    check(result.spectra.rows() == 1200 && result.spectra.cols() == 8, "spectrum dimensions are N x K");
    check(result.diagnostics.reconstruction_error <= 5e-5, "reference reconstruction error <= 5e-5");
    check(result.modes.array().isFinite().all(), "reference modes are finite");
    check(result.spectra.real().array().isFinite().all(), "reference spectra are finite");
}

void test_spectrum_contract() {
    const Eigen::VectorXd input = make_signal(257);
    const vmd::Result result = vmd::decompose(input, make_options(5));
    Eigen::FFT<double> fft;
    double squared_error = 0.0;
    double squared_reference = 0.0;
    for (Eigen::Index k = 0; k < result.modes.rows(); ++k) {
        std::vector<double> mode(static_cast<std::size_t>(input.size()));
        for (Eigen::Index i = 0; i < input.size(); ++i) mode[static_cast<std::size_t>(i)] = result.modes(k, i);
        std::vector<std::complex<double>> raw;
        fft.fwd(raw, mode);
        const Eigen::Index shift = (input.size() + 1) / 2;
        for (Eigen::Index i = 0; i < input.size(); ++i) {
            const auto expected = raw[static_cast<std::size_t>((i + shift) % input.size())];
            squared_error += std::norm(result.spectra(i, k) - expected);
            squared_reference += std::norm(expected);
        }
    }
    const double relative = std::sqrt(squared_error / std::max(squared_reference, 1e-300));
    check(relative < 1e-12, "spectrum equals fftshift(fft(mode))");
}

void test_degenerate_and_lengths() {
    for (Eigen::Index n : {Eigen::Index{1}, Eigen::Index{2}, Eigen::Index{31}, Eigen::Index{32}}) {
        Eigen::VectorXd zero = Eigen::VectorXd::Zero(n);
        vmd::Options value = make_options(3);
        value.dc_mode = true;
        const vmd::Result zero_result = vmd::decompose(zero, value);
        check(zero_result.modes.cols() == n, "zero signal preserves length");
        check(zero_result.omega_normalized.array().isFinite().all(), "zero frequencies are finite");

        Eigen::VectorXd constant = Eigen::VectorXd::Constant(n, 3.5);
        const vmd::Result constant_result = vmd::decompose(constant, value);
        check(constant_result.modes.cols() == n, "constant signal preserves odd/even length");
        check(constant_result.modes.array().isFinite().all(), "constant modes are finite");
    }
}

void test_random_seed() {
    const Eigen::VectorXd input = make_signal(129);
    vmd::Options value = make_options(4);
    value.initialization = vmd::Initialization::Random;
    value.random_seed = 42;
    const vmd::Result a = vmd::decompose(input, value);
    const vmd::Result b = vmd::decompose(input, value);
    check((a.modes - b.modes).norm() == 0.0, "fixed random seed is reproducible");
    check((a.omega_normalized - b.omega_normalized).norm() == 0.0, "frequency history is reproducible");
}

void test_per_mode_alpha_and_warm_start() {
    const Eigen::VectorXd input = make_signal(128);
    vmd::Options value = make_options(3);
    value.alpha = {40.0, 50.0, 60.0};
    const vmd::Result first = vmd::decompose(input, value);
    value.warm_start_modes = first.modes;
    value.warm_start_omega_normalized = first.omega_normalized.bottomRows(1).transpose();
    const vmd::Result warm = vmd::decompose(input, value);
    check(warm.modes.array().isFinite().all(), "warm start produces finite modes");
}


void test_single_frequency() {
    constexpr Eigen::Index n = 256;
    constexpr double frequency_hz = 16.0;
    Eigen::VectorXd input(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        input(i) = std::cos(2.0 * kPi * frequency_hz * static_cast<double>(i) / static_cast<double>(n));
    }
    vmd::Options value = make_options(1);
    value.alpha = {2000.0};
    value.sample_rate_hz = static_cast<double>(n);
    const vmd::Result result = vmd::decompose(input, value);
    check(result.diagnostics.reconstruction_error < 5e-2, "single-mode sinusoid reconstruction stays within boundary tolerance");
    check(std::abs(result.omega_hz(result.omega_hz.rows() - 1, 0) - frequency_hz) < 0.25,
          "single-mode center frequency matches the sinusoid");
}

void test_noisy_and_chirp_signals() {
    constexpr Eigen::Index n = 257;
    std::mt19937_64 generator(1234);
    std::normal_distribution<double> noise(0.0, 0.03);
    Eigen::VectorXd noisy(n);
    Eigen::VectorXd chirp(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        noisy(i) = std::cos(2.0 * kPi * 13.0 * t) + 0.2 * std::cos(2.0 * kPi * 47.0 * t) + noise(generator);
        chirp(i) = std::cos(2.0 * kPi * (5.0 * t + 35.0 * t * t / 2.0));
    }
    vmd::Options value = make_options(4);
    value.alpha = {500.0};
    value.sample_rate_hz = static_cast<double>(n);
    const vmd::Result noisy_result = vmd::decompose(noisy, value);
    const vmd::Result chirp_result = vmd::decompose(chirp, value);
    check(noisy_result.modes.array().isFinite().all() && noisy_result.spectra.real().array().isFinite().all(),
          "noisy signal produces finite outputs");
    check(chirp_result.modes.array().isFinite().all() && chirp_result.omega_hz.array().isFinite().all(),
          "chirp signal produces finite outputs");
}

void test_dc_and_advanced_options() {
    constexpr Eigen::Index n = 192;
    Eigen::VectorXd input(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        input(i) = 2.0 + 0.25 * std::cos(2.0 * kPi * 18.0 * t);
    }
    vmd::Options value = make_options(2);
    value.alpha = {1000.0};
    value.dc_mode = true;
    value.keep_omega_history = false;
    value.sample_rate_hz = static_cast<double>(n);
    const vmd::Result dc_result = vmd::decompose(input, value);
    check(dc_result.omega_normalized.rows() == 1, "omega history can be disabled");
    check(dc_result.omega_normalized(0, 0) == 0.0, "DC mode remains fixed at zero frequency");

    value.keep_omega_history = true;
    value.adaptive_tau = true;
    value.tau = 0.1;
    value.max_iterations = 30;
    const vmd::Result adaptive = vmd::decompose(input, value);
    check(std::isfinite(adaptive.diagnostics.final_tau) && adaptive.diagnostics.final_tau >= value.tau_min &&
              adaptive.diagnostics.final_tau <= value.tau_max,
          "adaptive tau remains finite and bounded");
}

void test_validation() {
    bool threw = false;
    try {
        Eigen::VectorXd empty;
        (void)vmd::decompose(empty, make_options());
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "empty signal is rejected");

    threw = false;
    try {
        Eigen::VectorXd invalid(2);
        invalid << 0.0, std::numeric_limits<double>::quiet_NaN();
        (void)vmd::decompose(invalid, make_options());
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "non-finite signal is rejected");

    threw = false;
    try {
        Eigen::VectorXd input = Eigen::VectorXd::Zero(4);
        vmd::Options invalid = make_options();
        invalid.initialization = static_cast<vmd::Initialization>(99);
        (void)vmd::decompose(input, invalid);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "invalid initialization enum is rejected");
}

void test_legacy_wrapper() {
    const Eigen::VectorXd input = make_signal(96);
    vectord values(input.data(), input.data() + input.size());
    Eigen::MatrixXd modes;
    Eigen::MatrixXcd spectra;
    Eigen::MatrixXd omega;
    VMD(modes, spectra, omega, values, 50.0, 0.0, 4, 0, 1, 1e-7, 2.2204e-16);
    check(modes.rows() == 4 && modes.cols() == 96, "legacy wrapper preserves mode orientation");
    check(spectra.rows() == 96 && spectra.cols() == 4, "legacy wrapper returns corrected spectra");
}

}  // namespace

int main() {
    test_reference_signal();
    test_spectrum_contract();
    test_degenerate_and_lengths();
    test_random_seed();
    test_per_mode_alpha_and_warm_start();
    test_single_frequency();
    test_noisy_and_chirp_signals();
    test_dc_and_advanced_options();
    test_validation();
    test_legacy_wrapper();
    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "All VMD tests passed\n";
    return EXIT_SUCCESS;
}
