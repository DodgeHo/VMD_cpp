#pragma once

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace vmd {

/** Center-frequency initialization strategy. */
enum class Initialization {
    Zero,     ///< Initialize every mode at zero cycles per sample.
    Uniform,  ///< Distribute modes uniformly from DC toward Nyquist.
    Random    ///< Use sorted, seeded log-uniform center frequencies.
};

/** Timings for the core numerical stages, in milliseconds. */
struct Timings {
    double mirror_and_fft_ms = 0.0;  ///< Mirror extension and input FFT.
    double iteration_ms = 0.0;       ///< Main alternating-update loop.
    double reconstruction_ms = 0.0;  ///< IMF reconstruction and public spectra.
};

/** Configuration for variational mode decomposition. */
struct Options {
    std::size_t mode_count = 3;  ///< Number of modes K; must be positive.

    /** One positive bandwidth penalty or exactly K positive penalties. */
    std::vector<double> alpha{2000.0};

    double tau = 0.0;  ///< Non-negative dual-ascent step.
    bool dc_mode = false;  ///< Keep the first mode centered at zero frequency.
    Initialization initialization = Initialization::Uniform;  ///< Frequency initializer.
    double tolerance = 1e-7;  ///< Positive stopping threshold for mode delta.
    std::size_t max_iterations = 500;  ///< Positive iteration limit.
    std::uint64_t random_seed = 0;  ///< Seed for Random initialization.
    double sample_rate_hz = 1.0;  ///< Positive sample rate used for omega_hz.
    bool keep_omega_history = true;  ///< Keep all center-frequency iterates.

    /** Optional finite K x N time-domain modes used as the initial iterate. */
    Eigen::MatrixXd warm_start_modes;

    /** Optional finite K-vector of initial frequencies in cycles per sample. */
    Eigen::VectorXd warm_start_omega_normalized;

    /** Enable residual-balanced tau adaptation. Disabled by default. */
    bool adaptive_tau = false;
    double adaptive_tau_balance = 10.0;  ///< Residual imbalance threshold (> 1).
    double adaptive_tau_increase = 2.0;  ///< Multiplicative increase factor (> 1).
    double adaptive_tau_decrease = 2.0;  ///< Multiplicative decrease divisor (> 1).
    double tau_min = 1e-8;  ///< Finite non-negative adaptive lower bound.
    double tau_max = 1e3;  ///< Finite adaptive upper bound, at least tau_min.
};

/** Numerical convergence and reconstruction information. */
struct Diagnostics {
    bool converged = false;  ///< True when final_mode_delta <= tolerance.
    std::size_t iterations = 0;  ///< Number of completed mode-update passes.
    double final_mode_delta = 0.0;  ///< Compatibility squared-iterate metric.
    double reconstruction_error = 0.0;  ///< Relative time-domain residual norm.
    double primal_residual = 0.0;  ///< Relative spectral constraint residual.
    double dual_residual = 0.0;  ///< Tau-scaled relative iterate change.
    double final_tau = 0.0;  ///< Final dual-ascent step after optional adaptation.
};

/** VMD outputs. Matrix orientations and units are part of the API contract. */
struct Result {
    Eigen::MatrixXd modes;             ///< K x N time-domain IMFs.
    Eigen::MatrixXcd spectra;          ///< N x K unnormalized fftshift(fft(IMF)).
    Eigen::MatrixXd omega_normalized;  ///< I x K center frequencies, cycles/sample.
    Eigen::MatrixXd omega_hz;          ///< I x K center frequencies, hertz.
    Diagnostics diagnostics;           ///< Convergence and residual metrics.
    Timings timings;                   ///< Core-stage durations in milliseconds.
};

/**
 * Decompose a finite, non-empty real signal into band-limited modes.
 *
 * @param signal Borrowed N-sample input vector. The function does not retain it.
 * @param options Decomposition settings. A one-element alpha vector broadcasts
 *        to all modes; otherwise alpha must contain exactly mode_count values.
 * @return Owning mode, spectrum, frequency, diagnostic, and timing data.
 *
 * @throws std::invalid_argument for invalid options, non-finite samples, or
 *         incompatible warm-start dimensions.
 * @throws std::runtime_error if the numerical backend produces non-finite output.
 *
 * When frequency history is enabled, I equals iterations + 1 and row zero is
 * the initialization. Otherwise I is one and contains only the final values.
 * Calls own all temporary state and are reentrant across independent threads.
 */
Result decompose(Eigen::Ref<const Eigen::VectorXd> signal, const Options& options);

}  // namespace vmd
