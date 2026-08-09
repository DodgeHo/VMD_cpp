#include <vmd/vmd.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

template <typename Integer>
Integer positive_integer(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0;
    if (value.empty() || value.front() == '-') {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    const unsigned long long parsed = std::stoull(value, &consumed);
    if (consumed != value.size() || parsed == 0) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    return static_cast<Integer>(parsed);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Eigen::Index n = argc > 1 ? positive_integer<Eigen::Index>(argv[1], "N") : 1200;
        const std::size_t k = argc > 2 ? positive_integer<std::size_t>(argv[2], "K") : 8;
        const std::size_t repetitions = argc > 3 ? positive_integer<std::size_t>(argv[3], "repetitions") : 1;
        constexpr double pi = 3.14159265358979323846;
        Eigen::VectorXd signal(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(n);
            signal(i) = std::cos(2.0 * pi * 3.0 * t) + 0.2 * std::cos(2.0 * pi * 71.0 * t);
        }

        vmd::Options options;
        options.mode_count = k;
        options.alpha = {2000.0};
        options.max_iterations = 500;
        options.keep_omega_history = true;
        options.sample_rate_hz = static_cast<double>(n);

        vmd::Result result;
        double mirror_fft_total = 0.0;
        double iteration_total = 0.0;
        double reconstruction_total = 0.0;
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t repeat = 0; repeat < repetitions; ++repeat) {
            result = vmd::decompose(signal, options);
            mirror_fft_total += result.timings.mirror_and_fft_ms;
            iteration_total += result.timings.iteration_ms;
            reconstruction_total += result.timings.reconstruction_ms;
        }
        const double total = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();

        const double working_set_estimate = static_cast<double>(
            2 * k * static_cast<std::size_t>(n) * sizeof(std::complex<double>) +
            7 * static_cast<std::size_t>(n) * sizeof(std::complex<double>) +
            5 * static_cast<std::size_t>(n) * sizeof(double));
        const double divisor = static_cast<double>(repetitions);

        std::cout << std::fixed << std::setprecision(3)
                  << "N=" << n << " K=" << k << " repetitions=" << repetitions
                  << " iterations=" << result.diagnostics.iterations << '\n'
                  << "mirror_fft_mean_ms=" << mirror_fft_total / divisor << '\n'
                  << "iteration_mean_ms=" << iteration_total / divisor << '\n'
                  << "reconstruction_mean_ms=" << reconstruction_total / divisor << '\n'
                  << "total_mean_ms=" << total / divisor << '\n'
                  << "estimated_core_workspace_mib=" << working_set_estimate / (1024.0 * 1024.0) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "benchmark error: " << error.what() << '\n';
        return 2;
    }
}
