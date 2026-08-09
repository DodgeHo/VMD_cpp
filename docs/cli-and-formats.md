# CLI and File Formats

## Invocation

```text
vmd_cli --input signal.csv --column value --sample-rate 1000 \
        --modes 6 --alpha 2000 --output-dir results --report
```

Use `--input -` for standard input. Use `--demo` for the deterministic built-in
1,200-sample signal. `--demo` and `--input` are mutually exclusive.

## Options

| Option | Value | Default |
|---|---|---|
| `--input` | path or `-` | none |
| `--column` | name or zero-based index | `0` |
| `--demo` | flag | off |
| `--sample-rate` | positive hertz | `1`; demo uses `1200` |
| `--modes` | positive integer | `3` |
| `--alpha` | one value or comma-separated `K` values | `2000` |
| `--tau` | non-negative real | `0` |
| `--dc` | flag | off |
| `--init` | `zero`, `uniform`, or `random` | `uniform` |
| `--tol` | positive real | `1e-7` |
| `--max-iterations` | positive integer | `500` |
| `--seed` | unsigned integer | `0` |
| `--adaptive-tau` | flag | off |
| `--no-omega-history` | flag | off |
| `--output-dir` | path | `vmd-output` |
| `--report` | flag | off |
| `--help`, `-h` | flag | off |

## Input Rules

- Comma and tab delimiters are detected from the first non-empty row.
- UTF-8 BOM, CRLF, and blank physical lines are accepted.
- The selected column may be a header name or a zero-based index.
- Header detection for an indexed column uses that selected field, allowing text
  metadata in unselected columns.
- Quoted fields and doubled quotes are supported.
- Every non-empty row must have the same number of fields.
- The selected sample must parse completely as a finite `double`.
- Missing values, `NaN`, infinity, malformed numbers, and unterminated quotes are
  rejected. Errors include source, physical line, and one-based column where
  applicable.

Input does not carry time stamps into the solver. Output time is reconstructed as
`index / sample_rate_hz`.

## `modes.csv`

One row per input sample:

```text
index,time_s,signal,imf_1,...,imf_K,reconstruction,residual
```

`reconstruction` is the sum of all IMFs. `residual` is
`signal - reconstruction`.

## `spectra.csv`

One row per shifted FFT bin:

```text
index,frequency_hz,
imf_1_real,imf_1_imag,imf_1_magnitude,...
```

Frequency is `(index - floor(N/2)) * sample_rate_hz / N`. Complex values are the
unnormalized forward FFT of each cropped IMF after `fftshift`.

## `omega.csv`

```text
iteration,imf_1_hz,...,imf_K_hz
```

Row zero is initialization when history is enabled. With
`--no-omega-history`, the file contains a header and one final-frequency row.

## `summary.json`

The summary contains:

- `schema_version`
- `input`: source, selected column, sample count, rate, delimiter, and header flag
- `options`: effective algorithm settings
- `diagnostics`: convergence and residual values
- `timings_ms`: input, core stages, serialization, and report generation
- `outputs`: generated file names

The machine-readable schema is [summary.schema.json](summary.schema.json). Paths
in the summary are file names relative to the output directory. `outputs.report`
is `null` when `--report` is not requested.

## `report.html`

The report is one offline file. CSS, JavaScript, plot data, and SVG generation are
inline; there are no CDN or runtime dependencies. It provides:

- Input and reconstruction comparison
- Reconstruction residual
- Per-mode time plots and visibility controls
- One-sided magnitude spectra
- Center-frequency convergence
- Hover readouts and SVG download buttons
- Links to the CSV and JSON files

Plots use a min/max envelope capped at 2,000 plotted points per series. The CSV
files are never downsampled.

## Exit Codes

| Code | Meaning |
|---:|---|
| `0` | Successful convergence or help output. |
| `2` | CLI argument or algorithm option error. |
| `3` | Input, output, parsing, or runtime error. |
| `4` | Maximum iterations reached before convergence. Outputs are retained. |
