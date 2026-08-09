# Algorithm and Numerical Behavior

## Objective

VMD represents a real signal `f(t)` as a sum of band-limited modes `u_k(t)`.
For center frequencies `omega_k`, the constrained objective minimizes the summed
bandwidth of each demodulated analytic mode while enforcing reconstruction. The
implementation follows the standard augmented-Lagrangian/ADMM update in the
frequency domain:

```text
u_k <- (f_plus - sum_{i != k} u_i - lambda / 2)
       / (1 + alpha_k * (frequency - omega_k)^2)

omega_k <- dot(frequency, abs(u_k)^2) / sum(abs(u_k)^2)

lambda <- lambda + tau * (sum_k u_k - f_plus)
```

`alpha_k` may differ by mode. A single `Options::alpha` value is broadcast to all
modes.

## Frequency Layout

The mirrored signal has length `2N`. A real-to-complex FFT produces its
`N + 1` non-negative bins; the solver stores only bins `j = 0 ... N-1` because
the Nyquist bin is outside the analytic update interval. Normalized frequencies
are `j / (2N)`, in cycles per sample, and therefore lie in `[0, 0.5)`. This is
equivalent to the positive half of the shifted spectrum used by the reference
formulation, but avoids allocating or transforming a redundant negative half.

Public spectra have a different, explicit contract: after every time-domain IMF
is cropped back to `N` samples, the implementation performs a real-to-complex
FFT, reconstructs the Hermitian half, and applies `fftshift`. `Result::spectra`
is always `N x K`. It never exposes the mirrored-domain iterate.

For an output row `i`, the corresponding shifted frequency in hertz is:

```text
(i - floor(N / 2)) * sample_rate_hz / N
```

## Mirror Boundary

For input length `N`, the extension contains:

```text
reverse(first floor(N/2) samples)
original N samples
reverse(last ceil(N/2) samples)
```

The total mirrored length is exactly `2N` for both odd and even inputs. Cropping
starts at `floor(N/2)` and therefore returns exactly the original sample count.

## Initialization

- `Zero`: all initial center frequencies are zero.
- `Uniform`: mode `k` starts at `0.5 * k / K`.
- `Random`: log-uniform values are drawn in the resolvable positive frequency
  interval with `std::mt19937_64`, sorted, and controlled by `random_seed`.

If `dc_mode` is enabled, mode zero is fixed at frequency zero. Optional warm-start
modes are mirrored and transformed into the internal positive-frequency layout.
Optional warm-start frequencies are normalized cycles per sample in `[0, 0.5]`.

## Degenerate Inputs

The center-frequency denominator can be zero for a zero signal or an empty mode.
When mode energy is below a scale-aware floor, the previous frequency is retained
instead of dividing by zero. Length-one FFTs are handled as identity transforms
because some FFT backends do not support that plan size safely. Inputs and all
public results are checked for finite values.

## Stopping and Diagnostics

The compatibility stopping metric is retained:

```text
mode_delta = epsilon + squaredNorm(current - previous) / mirrored_length
```

The solver stops when `mode_delta <= tolerance` or `max_iterations` is reached.
It also reports:

- `reconstruction_error`: time-domain relative norm of `sum(modes) - signal`.
- `primal_residual`: relative spectral norm of `sum(mode_iterates) - input_plus`.
- `dual_residual`: `tau` times the relative iterate-change norm.

For a zero input, reconstruction error is normalized by 1 rather than zero.

## Adaptive Tau

Adaptive `tau` is disabled by default. When enabled, residual balancing increases
or decreases `tau` if one residual exceeds the other by
`adaptive_tau_balance`, while clamping to `[tau_min, tau_max]`. This is an
advanced convergence heuristic and may change the decomposition path; use a fixed
`tau` when reproducing historical results.

## Memory Model

The main mode iterate uses two `K x N` complex buffers. The dual variable is
updated in place. Residual, numerator, denominator, frequency, and energy vectors
are preallocated. The core workspace is therefore `O(K * N)` and does not grow
with `max_iterations`.

Center-frequency history is the only optional iteration-sized output. Set
`keep_omega_history = false` to retain only the final `1 x K` row.

## FFT Backends

- `Eigen`: default and dependency-free beyond Eigen.
- `FFTW`: selected with `-DVMD_FFT_BACKEND=FFTW`.
- `MKL`: selected with `-DVMD_FFT_BACKEND=MKL`.

Each decomposition owns one backend object and reuses its plans and buffers across
mirror, warm-start, reconstruction, and public-spectrum transforms. Real inputs
use the backend's dedicated real-to-complex half-spectrum path.
