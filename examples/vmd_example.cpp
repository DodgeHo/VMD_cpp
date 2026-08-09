#include <vmd/vmd.hpp>

#include <cmath>
#include <iostream>

int main() {
    constexpr double pi = 3.14159265358979323846;
    constexpr Eigen::Index n = 1200;
    Eigen::VectorXd signal(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        const double t = static_cast<double>(i + 1) / static_cast<double>(n);
        signal(i) = std::cos(2.0 * pi * 2.0 * t) +
                    std::cos(2.0 * pi * 24.0 * t) / 10.0 +
                    std::cos(2.0 * pi * 288.0 * t) / 200.0;
    }

    vmd::Options options;
    options.mode_count = 8;
    options.alpha = {50.0};
    options.sample_rate_hz = 1200.0;
    const vmd::Result result = vmd::decompose(signal, options);

    std::cout << "modes: " << result.modes.rows() << " x " << result.modes.cols() << '\n'
              << "iterations: " << result.diagnostics.iterations << '\n'
              << "converged: " << (result.diagnostics.converged ? "yes" : "no") << '\n'
              << "reconstruction error: " << result.diagnostics.reconstruction_error << '\n'
              << "center frequencies (Hz): " << result.omega_hz.bottomRows(1) << '\n';
    return result.diagnostics.converged ? 0 : 4;
}
