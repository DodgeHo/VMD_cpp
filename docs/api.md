# C++ API

Include the public header:

```cpp
#include <vmd/vmd.hpp>
```

All new API types are in namespace `vmd`. The library requires C++17 and exposes
Eigen matrix types in its public interface.

## `Initialization`

| Value | Behavior |
|---|---|
| `Zero` | All initial center frequencies are zero. |
| `Uniform` | Frequencies are distributed from zero toward Nyquist. |
| `Random` | Sorted log-uniform frequencies generated from `random_seed`. |

## `Options`

| Field | Default | Meaning |
|---|---:|---|
| `mode_count` | `3` | Number of modes `K`; must be positive. |
| `alpha` | `{2000}` | One positive value or exactly `K` positive values. |
| `tau` | `0` | Non-negative dual ascent step. |
| `dc_mode` | `false` | Fix the first center frequency at zero. |
| `initialization` | `Uniform` | Initial frequency strategy. |
| `tolerance` | `1e-7` | Positive stopping threshold for mode delta. |
| `max_iterations` | `500` | Positive iteration limit. |
| `random_seed` | `0` | Seed used by random initialization. |
| `sample_rate_hz` | `1` | Positive sample rate used to derive `omega_hz`. |
| `keep_omega_history` | `true` | Keep every frequency row instead of only the final row. |
| `warm_start_modes` | empty | Optional `K x N` finite time-domain modes. |
| `warm_start_omega_normalized` | empty | Optional finite `K`-vector in `[0, 0.5]`. |
| `adaptive_tau` | `false` | Enable residual-balanced step adjustment. |
| `adaptive_tau_balance` | `10` | Residual imbalance threshold, greater than one. |
| `adaptive_tau_increase` | `2` | Multiplicative increase, greater than one. |
| `adaptive_tau_decrease` | `2` | Multiplicative decrease divisor, greater than one. |
| `tau_min` | `1e-8` | Finite non-negative lower bound. |
| `tau_max` | `1e3` | Finite upper bound, not below `tau_min`. |

## `Result`

| Field | Shape/unit | Meaning |
|---|---|---|
| `modes` | `K x N` | Time-domain intrinsic mode functions. |
| `spectra` | `N x K` | `fftshift(fft(mode))`, complex and unnormalized. |
| `omega_normalized` | `I x K` | Center frequencies in cycles per sample. |
| `omega_hz` | `I x K` | Center frequencies in hertz. |
| `diagnostics` | scalar fields | Convergence and residual metrics. |
| `timings` | milliseconds | Mirror/FFT, iteration, and reconstruction durations. |

`I` is `iterations + 1` when history is enabled because row zero contains the
initial frequencies. It is one when history is disabled.

## Diagnostics

- `converged`: the final mode delta met `tolerance`.
- `iterations`: completed update passes.
- `final_mode_delta`: compatibility stopping metric.
- `reconstruction_error`: relative time-domain reconstruction norm.
- `primal_residual`: relative frequency-domain constraint residual.
- `dual_residual`: scaled relative iterate change.
- `final_tau`: final dual ascent step after optional adaptation.

## Function

```cpp
vmd::Result vmd::decompose(
    Eigen::Ref<const Eigen::VectorXd> signal,
    const vmd::Options& options);
```

The signal must be non-empty and finite. Invalid arguments throw
`std::invalid_argument`. A non-finite internal result throws `std::runtime_error`.
Nonconvergence is not an exception; inspect `Result::diagnostics.converged`.

Inputs are borrowed only for the duration of the call. Results and work buffers
own their storage. The implementation has no mutable decomposition globals, so
independent calls are reentrant. Backend library initialization rules still apply
when an application selects FFTW or MKL.

## Warm Start Example

```cpp
vmd::Result first = vmd::decompose(signal, options);

vmd::Options next = options;
next.warm_start_modes = first.modes;
next.warm_start_omega_normalized =
    first.omega_normalized.bottomRows(1).transpose();

vmd::Result refined = vmd::decompose(signal, next);
```

Warm starts must use the same sample count and mode count.

## Legacy API

`VMD.h` retains the original global wrapper:

```cpp
void VMD(Eigen::MatrixXd& u,
         Eigen::MatrixXcd& u_hat,
         Eigen::MatrixXd& omega,
         const std::vector<double>& signal,
         double alpha, double tau, int K, int DC, int init,
         double tol, double eps);
```

The wrapper uses at most 500 iterations and maps `init` values 0, 1, and 2 to
zero, uniform, and random initialization. `eps` is accepted for source
compatibility. Prefer `vmd::decompose` for explicit seeds, sample rates, per-mode
alpha, warm starts, adaptive tau, and diagnostics.
