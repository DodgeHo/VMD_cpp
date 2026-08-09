# Validation and Performance

## Correctness Tests

The `vmd_tests` executable covers:

- The historical 1,200-sample multifrequency example
- Single-frequency center estimation
- Deterministic noisy input and a chirp
- DC-constrained decomposition
- Zero and constant inputs
- Length 1, length 2, odd, and even signals
- Reproducible random initialization
- Per-mode alpha and warm starts
- Optional frequency-history suppression and adaptive tau bounds
- Invalid and non-finite arguments
- The legacy `VMD(...)` wrapper

The spectrum contract is checked independently with Eigen FFT. The relative norm
between every returned spectrum and `fftshift(fft(mode))` must be below `1e-12`.
The historical example must have relative reconstruction error no greater than
`5e-5`, and all public numeric matrices must remain finite.

Run:

```sh
ctest --test-dir build -C Release --output-on-failure
```

## CLI End-to-End Test

`vmd.cli.e2e` executes the demo and verifies that all CSV, JSON, and HTML files
exist and contain their expected structural markers. It also verifies:

- Indexed input does not lose its first row when another column contains text.
- Invalid numeric input returns code 3 with line and column information.
- A deliberately nonconverged run returns code 4 and preserves output files.

## Benchmark

Build and run:

```sh
cmake --build build --config Release --target vmd_benchmark
./build/vmd_benchmark 1200 8 20
```

The optional third argument repeats the same deterministic decomposition and the
reported stage and wall times are per-run means. Repetitions make short runs easier
to profile without changing their workspace shape. The benchmark reports
mirror/FFT, main iteration, output reconstruction, total wall time, and an
estimated core workspace. CLI `summary.json` additionally
reports input parsing, CSV serialization, and HTML report generation time.

Performance numbers should always record:

- CPU model and logical core count
- Operating system and compiler version
- CMake build type and FFT backend
- Signal length, mode count, tolerance, and completed iterations

### Local Release Snapshot

The final local regression on July 31, 2026 used Windows 11 64-bit,
an Intel Core i5-13500HX (20 logical processors), MinGW GCC 15.2.0, Release
mode, and the Eigen real-to-complex FFT backend. For `N = 1200`, `K = 8`,
500 repeated decompositions, and 452 iterations per decomposition, it measured:

| Stage | Mean time |
| --- | ---: |
| Mirror and input FFT | 0.154 ms |
| Main iteration | 32.029 ms |
| Mode and spectrum reconstruction | 0.480 ms |
| Total | 32.847 ms |

The estimated core workspace was 0.467 MiB and the sampled process working-set
peak was 6.023 MiB. These values are a reproducibility snapshot, not a portable
performance guarantee.

## Complexity

The iteration keeps two complex `K x N` mode buffers and a fixed number of
length-`N` work vectors. Core working memory is `O(K * N)` and does not scale with
`max_iterations`. Full omega history is `O(iterations * K)` and may be disabled.

For `N = 1200` and `K = 8`, the benchmark estimate is far below the 32 MiB core
workspace target. Process peak memory is higher because it also includes the C++
runtime, Eigen/FFT backend, executable image, and serialized result matrices.

## Sanitizers and Warnings

High warnings are enabled by default. On GCC or Clang:

```sh
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug \
  -DVMD_ENABLE_SANITIZERS=ON -DVMD_BUILD_BENCHMARKS=OFF
cmake --build build-asan --parallel
ctest --test-dir build-asan --output-on-failure
```

The CI workflow covers MSVC, GCC, and Clang plus an ASan/UBSan Linux job.

## Reference Comparison

The original implementation and standard VMD formulation use a shifted
positive-frequency signal and the same sequential Wiener-filter update. This
implementation stores the equivalent unshifted positive half to reduce memory.
Behavioral comparison should use time-domain modes, final normalized center
frequencies, and reconstruction residual rather than the historical `u_hat`,
which mixed domains and did not satisfy its documented shape/content contract.
