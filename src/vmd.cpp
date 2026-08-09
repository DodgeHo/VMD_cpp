#include <vmd/vmd.hpp>

#include "fft_backend.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace vmd {
namespace {

using Clock = std::chrono::steady_clock;
using Complex = std::complex<double>;
using RowMatrixXcd = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

constexpr double kEnergyFloor = 64.0 * std::numeric_limits<double>::epsilon();

void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}

std::vector<double> expanded_alpha(const Options& options) {
    require(options.alpha.size() == 1 || options.alpha.size() == options.mode_count,
            "alpha must contain one value or exactly mode_count values");
    std::vector<double> values(options.mode_count, options.alpha.front());
    if (options.alpha.size() == options.mode_count) values = options.alpha;
    for (double value : values) {
        require(std::isfinite(value) && value > 0.0, "alpha values must be finite and positive");
    }
    return values;
}

void validate(Eigen::Ref<const Eigen::VectorXd> signal, const Options& options) {
    require(signal.size() > 0, "signal must not be empty");
    require(static_cast<std::uint64_t>(signal.size()) <=
                static_cast<std::uint64_t>(std::numeric_limits<Eigen::Index>::max() / 2),
            "signal is too large");
    require(signal.array().isFinite().all(), "signal must contain only finite values");
    require(options.mode_count > 0, "mode_count must be positive");
    require(options.mode_count <= static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()),
            "mode_count is too large");
    (void)expanded_alpha(options);
    require(std::isfinite(options.tau) && options.tau >= 0.0, "tau must be finite and non-negative");
    require(std::isfinite(options.tolerance) && options.tolerance > 0.0,
            "tolerance must be finite and positive");
    require(options.max_iterations > 0, "max_iterations must be positive");
    require(std::isfinite(options.sample_rate_hz) && options.sample_rate_hz > 0.0,
            "sample_rate_hz must be finite and positive");
    switch (options.initialization) {
        case Initialization::Zero:
        case Initialization::Uniform:
        case Initialization::Random:
            break;
        default:
            throw std::invalid_argument("initialization is not a valid enum value");
    }
    require(std::isfinite(options.tau_min) && std::isfinite(options.tau_max) &&
                options.tau_min >= 0.0 && options.tau_max >= options.tau_min,
            "adaptive tau bounds must be finite, non-negative, and ordered");
    require(std::isfinite(options.adaptive_tau_balance) &&
                std::isfinite(options.adaptive_tau_increase) &&
                std::isfinite(options.adaptive_tau_decrease) &&
                options.adaptive_tau_balance > 1.0 && options.adaptive_tau_increase > 1.0 &&
                options.adaptive_tau_decrease > 1.0,
            "adaptive tau factors must be finite and greater than one");
    if (options.adaptive_tau) require(options.tau > 0.0, "adaptive_tau requires tau > 0");

    const Eigen::Index k = static_cast<Eigen::Index>(options.mode_count);
    if (options.warm_start_modes.size() != 0) {
        require(options.warm_start_modes.rows() == k && options.warm_start_modes.cols() == signal.size(),
                "warm_start_modes must be K x N");
        require(options.warm_start_modes.array().isFinite().all(),
                "warm_start_modes must contain only finite values");
    }
    if (options.warm_start_omega_normalized.size() != 0) {
        require(options.warm_start_omega_normalized.size() == k,
                "warm_start_omega_normalized must contain K values");
        require(options.warm_start_omega_normalized.array().isFinite().all(),
                "warm_start_omega_normalized must contain only finite values");
        require((options.warm_start_omega_normalized.array() >= 0.0).all() &&
                    (options.warm_start_omega_normalized.array() <= 0.5).all(),
                "warm-start frequencies must lie in [0, 0.5]");
    }
}

std::vector<double> mirror_signal(Eigen::Ref<const Eigen::VectorXd> signal) {
    const Eigen::Index n = signal.size();
    const Eigen::Index left = n / 2;
    const Eigen::Index right = n - left;
    std::vector<double> mirrored(static_cast<std::size_t>(2 * n));
    for (Eigen::Index i = 0; i < left; ++i) mirrored[static_cast<std::size_t>(i)] = signal(left - 1 - i);
    for (Eigen::Index i = 0; i < n; ++i) mirrored[static_cast<std::size_t>(left + i)] = signal(i);
    for (Eigen::Index i = 0; i < right; ++i) {
        mirrored[static_cast<std::size_t>(left + n + i)] = signal(n - 1 - i);
    }
    return mirrored;
}

Eigen::RowVectorXd initialize_omega(const Options& options, Eigen::Index signal_size) {
    const Eigen::Index k_count = static_cast<Eigen::Index>(options.mode_count);
    Eigen::RowVectorXd omega = Eigen::RowVectorXd::Zero(k_count);
    if (options.warm_start_omega_normalized.size() != 0) {
        omega = options.warm_start_omega_normalized.transpose();
    } else if (options.initialization == Initialization::Uniform) {
        for (Eigen::Index k = 0; k < k_count; ++k) {
            omega(k) = 0.5 * static_cast<double>(k) / static_cast<double>(k_count);
        }
    } else if (options.initialization == Initialization::Random && signal_size > 2) {
        const double low = std::min(0.5, 1.0 / static_cast<double>(signal_size));
        std::mt19937_64 generator(options.random_seed);
        std::uniform_real_distribution<double> distribution(std::log(low), std::log(0.5));
        for (Eigen::Index k = 0; k < k_count; ++k) omega(k) = std::exp(distribution(generator));
        std::sort(omega.data(), omega.data() + omega.size());
    }
    if (options.dc_mode) omega(0) = 0.0;
    return omega;
}

void initialize_modes_from_warm_start(
    const Options& options,
    detail::FFTBackend& fft,
    RowMatrixXcd& modes) {
    if (options.warm_start_modes.size() == 0) return;
    const Eigen::Index n = options.warm_start_modes.cols();
    const Eigen::Index half = n;
    for (Eigen::Index k = 0; k < options.warm_start_modes.rows(); ++k) {
        Eigen::VectorXd mode = options.warm_start_modes.row(k).transpose();
        auto mirrored = mirror_signal(mode);
        std::vector<Complex> spectrum;
        fft.forward_real(mirrored, spectrum);
        for (Eigen::Index j = 0; j < half; ++j) modes(k, j) = spectrum[static_cast<std::size_t>(j)];
    }
}

Eigen::MatrixXd build_history(
    const std::vector<Eigen::RowVectorXd>& history,
    const Eigen::RowVectorXd& final_omega,
    bool keep_history) {
    if (!keep_history) {
        Eigen::MatrixXd result(1, final_omega.size());
        result.row(0) = final_omega;
        return result;
    }
    Eigen::MatrixXd result(static_cast<Eigen::Index>(history.size()), final_omega.size());
    for (std::size_t i = 0; i < history.size(); ++i) {
        result.row(static_cast<Eigen::Index>(i)) = history[i];
    }
    return result;
}

void verify_finite(const Result& result) {
    if (!result.modes.array().isFinite().all() || !result.spectra.real().array().isFinite().all() ||
        !result.spectra.imag().array().isFinite().all() ||
        !result.omega_normalized.array().isFinite().all() || !result.omega_hz.array().isFinite().all()) {
        throw std::runtime_error("VMD produced a non-finite result");
    }
}

}  // namespace

Result decompose(Eigen::Ref<const Eigen::VectorXd> signal, const Options& options) {
    validate(signal, options);
    const auto alphas = expanded_alpha(options);
    const Eigen::Index n = signal.size();
    const Eigen::Index mirrored_size = 2 * n;
    const Eigen::Index half = n;
    const Eigen::Index k_count = static_cast<Eigen::Index>(options.mode_count);

    Result result;
    detail::FFTBackend fft;

    const auto preparation_start = Clock::now();
    Eigen::RowVectorXcd input_positive(half);
    Eigen::RowVectorXd frequencies(half);
    {
        const auto mirrored = mirror_signal(signal);
        std::vector<Complex> mirrored_spectrum;
        fft.forward_real(mirrored, mirrored_spectrum);
        for (Eigen::Index j = 0; j < half; ++j) {
            input_positive(j) = mirrored_spectrum[static_cast<std::size_t>(j)];
            frequencies(j) = static_cast<double>(j) / static_cast<double>(mirrored_size);
        }
    }

    RowMatrixXcd previous = RowMatrixXcd::Zero(k_count, half);
    RowMatrixXcd current = RowMatrixXcd::Zero(k_count, half);
    initialize_modes_from_warm_start(options, fft, previous);

    Eigen::RowVectorXd omega_previous = initialize_omega(options, n);
    Eigen::RowVectorXd omega_current = omega_previous;
    std::vector<Eigen::RowVectorXd> omega_history;
    if (options.keep_omega_history) {
        omega_history.reserve(options.max_iterations + 1);
        omega_history.push_back(omega_previous);
    }
    result.timings.mirror_and_fft_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - preparation_start).count();

    const auto iteration_start = Clock::now();
    Eigen::RowVectorXcd lambda = Eigen::RowVectorXcd::Zero(half);
    Eigen::RowVectorXcd remaining_previous(half);
    Eigen::RowVectorXcd current_sum(half);
    Eigen::RowVectorXcd numerator(half);
    Eigen::RowVectorXcd spectral_residual(half);
    Eigen::RowVectorXd denominator(half);
    Eigen::RowVectorXd frequency_delta(half);
    Eigen::RowVectorXd energy(half);

    const double input_scale = std::max(input_positive.norm(), kEnergyFloor);
    const double energy_floor = kEnergyFloor * std::max(1.0, input_positive.squaredNorm());
    double tau = options.tau;
    double final_delta = std::numeric_limits<double>::infinity();
    double primal_residual = std::numeric_limits<double>::infinity();
    double dual_residual = std::numeric_limits<double>::infinity();
    bool converged = false;
    std::size_t iterations = 0;

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        remaining_previous = previous.colwise().sum();
        current_sum.setZero();

        for (Eigen::Index k = 0; k < k_count; ++k) {
            remaining_previous -= previous.row(k);
            numerator = input_positive - current_sum - remaining_previous - 0.5 * lambda;
            frequency_delta = frequencies.array() - omega_previous(k);
            denominator = 1.0 + alphas[static_cast<std::size_t>(k)] * frequency_delta.array().square();
            current.row(k) = numerator.array() / denominator.array();

            if (options.dc_mode && k == 0) {
                omega_current(k) = 0.0;
            } else {
                energy = current.row(k).cwiseAbs2();
                const double mode_energy = energy.sum();
                omega_current(k) = mode_energy > energy_floor
                    ? frequencies.dot(energy) / mode_energy
                    : omega_previous(k);
            }
            current_sum += current.row(k);
        }

        spectral_residual = current_sum - input_positive;
        primal_residual = spectral_residual.norm() / input_scale;
        const double iterate_norm = (current - previous).norm();
        dual_residual = tau * iterate_norm / input_scale;
        lambda += tau * spectral_residual;
        final_delta = std::numeric_limits<double>::epsilon() +
                      (current - previous).squaredNorm() / static_cast<double>(mirrored_size);
        iterations = iteration + 1;

        if (options.keep_omega_history) omega_history.push_back(omega_current);

        if (options.adaptive_tau) {
            if (primal_residual > options.adaptive_tau_balance * std::max(dual_residual, kEnergyFloor)) {
                tau = std::min(options.tau_max, tau * options.adaptive_tau_increase);
            } else if (dual_residual > options.adaptive_tau_balance * std::max(primal_residual, kEnergyFloor)) {
                tau = std::max(options.tau_min, tau / options.adaptive_tau_decrease);
            }
        }

        previous.swap(current);
        omega_previous.swap(omega_current);
        if (final_delta <= options.tolerance) {
            converged = true;
            break;
        }
    }

    result.timings.iteration_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - iteration_start).count();

    const auto reconstruction_start = Clock::now();
    result.modes.resize(k_count, n);
    result.spectra.resize(n, k_count);
    std::vector<Complex> full_spectrum(static_cast<std::size_t>(mirrored_size));
    std::vector<Complex> mirrored_mode;
    std::vector<double> cropped_mode(static_cast<std::size_t>(n));
    std::vector<Complex> cropped_half_spectrum;
    const Eigen::Index crop_start = n / 2;

    for (Eigen::Index k = 0; k < k_count; ++k) {
        std::fill(full_spectrum.begin(), full_spectrum.end(), Complex{});
        for (Eigen::Index j = 0; j < half; ++j) {
            full_spectrum[static_cast<std::size_t>(j)] = previous(k, j);
        }
        for (Eigen::Index j = 1; j < half; ++j) {
            full_spectrum[static_cast<std::size_t>(mirrored_size - j)] = std::conj(previous(k, j));
        }
        fft.inverse(full_spectrum, mirrored_mode);
        for (Eigen::Index i = 0; i < n; ++i) {
            const double value = mirrored_mode[static_cast<std::size_t>(crop_start + i)].real();
            result.modes(k, i) = value;
            cropped_mode[static_cast<std::size_t>(i)] = value;
        }

        fft.forward_real(cropped_mode, cropped_half_spectrum);
        const Eigen::Index shift = (n + 1) / 2;
        for (Eigen::Index i = 0; i < n; ++i) {
            const Eigen::Index source = (i + shift) % n;
            result.spectra(i, k) = source <= n / 2
                ? cropped_half_spectrum[static_cast<std::size_t>(source)]
                : std::conj(cropped_half_spectrum[static_cast<std::size_t>(n - source)]);
        }
    }

    result.omega_normalized = build_history(omega_history, omega_previous, options.keep_omega_history);
    result.omega_hz = result.omega_normalized * options.sample_rate_hz;
    result.diagnostics.converged = converged;
    result.diagnostics.iterations = iterations;
    result.diagnostics.final_mode_delta = final_delta;
    const Eigen::VectorXd reconstruction = result.modes.colwise().sum().transpose();
    const double signal_norm = signal.norm();
    result.diagnostics.reconstruction_error = (reconstruction - signal).norm() /
        (signal_norm > kEnergyFloor ? signal_norm : 1.0);
    result.diagnostics.primal_residual = primal_residual;
    result.diagnostics.dual_residual = dual_residual;
    result.diagnostics.final_tau = tau;
    result.timings.reconstruction_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - reconstruction_start).count();

    verify_finite(result);
    return result;
}

}  // namespace vmd
