#include "cli_io.hpp"
#include "report.hpp"

#include <vmd/vmd.hpp>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using Clock = std::chrono::steady_clock;

struct Arguments {
    std::string input;
    std::string column;
    std::filesystem::path output_directory = "vmd-output";
    vmd::Options options;
    bool sample_rate_set = false;
    bool demo = false;
    bool report = false;
    bool help = false;
};

void print_help(std::ostream& output) {
    output << R"HELP(VMD_cpp command-line interface

Usage:
  vmd_cli --input signal.csv --column value --sample-rate 1000 \
          --modes 6 --alpha 2000 --output-dir results --report
  vmd_cli --demo --sample-rate 1200 --modes 8 --alpha 50 --report

Input:
  --input PATH              CSV/TSV path, or - for standard input
  --column NAME|INDEX       Column name or zero-based column index (default: 0)
  --demo                    Use the built-in deterministic example signal

Algorithm:
  --sample-rate HZ          Sampling rate in hertz (default: 1; demo: 1200)
  --modes K                 Number of modes (default: 3)
  --alpha A[,A...]          One alpha or one value per mode (default: 2000)
  --tau VALUE               Dual ascent step, non-negative (default: 0)
  --dc                      Keep the first mode at zero frequency
  --init zero|uniform|random
  --tol VALUE               Convergence tolerance (default: 1e-7)
  --max-iterations N        Maximum iterations (default: 500)
  --seed N                  Random initialization seed (default: 0)
  --adaptive-tau            Enable residual-balanced tau adaptation
  --no-omega-history        Keep only final center frequencies

Output:
  --output-dir PATH         Output directory (default: vmd-output)
  --report                  Generate self-contained report.html
  --help                    Show this help

Exit codes: 0 success, 2 argument error, 3 I/O error, 4 nonconvergence.
)HELP";
}

template <typename Integer>
Integer integer_value(const std::string& text, const std::string& option) {
    std::size_t consumed = 0;
    try {
        if (text.empty() || text.front() == '-') throw std::invalid_argument("");
        const unsigned long long value = std::stoull(text, &consumed);
        if (consumed != text.size() || value > static_cast<unsigned long long>(std::numeric_limits<Integer>::max())) {
            throw std::invalid_argument("");
        }
        return static_cast<Integer>(value);
    } catch (const std::exception&) {
        throw std::invalid_argument(option + " expects a non-negative integer");
    }
}

double real_value(const std::string& text, const std::string& option) {
    std::size_t consumed = 0;
    try {
        const double value = std::stod(text, &consumed);
        if (consumed != text.size() || !std::isfinite(value)) throw std::invalid_argument("");
        return value;
    } catch (const std::exception&) {
        throw std::invalid_argument(option + " expects a finite number");
    }
}

Arguments parse_arguments(int argc, char** argv) {
    Arguments arguments;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        const auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(option + " requires a value");
            return argv[++i];
        };

        if (option == "--help" || option == "-h") arguments.help = true;
        else if (option == "--input") arguments.input = value();
        else if (option == "--column") arguments.column = value();
        else if (option == "--sample-rate") { arguments.options.sample_rate_hz = real_value(value(), option); arguments.sample_rate_set = true; }
        else if (option == "--modes") arguments.options.mode_count = integer_value<std::size_t>(value(), option);
        else if (option == "--alpha") arguments.options.alpha = vmd::cli::parse_alpha_list(value());
        else if (option == "--tau") arguments.options.tau = real_value(value(), option);
        else if (option == "--dc") arguments.options.dc_mode = true;
        else if (option == "--init") {
            const std::string init = value();
            if (init == "zero") arguments.options.initialization = vmd::Initialization::Zero;
            else if (init == "uniform") arguments.options.initialization = vmd::Initialization::Uniform;
            else if (init == "random") arguments.options.initialization = vmd::Initialization::Random;
            else throw std::invalid_argument("--init expects zero, uniform, or random");
        } else if (option == "--tol") arguments.options.tolerance = real_value(value(), option);
        else if (option == "--max-iterations") arguments.options.max_iterations = integer_value<std::size_t>(value(), option);
        else if (option == "--seed") arguments.options.random_seed = integer_value<std::uint64_t>(value(), option);
        else if (option == "--adaptive-tau") arguments.options.adaptive_tau = true;
        else if (option == "--no-omega-history") arguments.options.keep_omega_history = false;
        else if (option == "--output-dir") arguments.output_directory = value();
        else if (option == "--report") arguments.report = true;
        else if (option == "--demo") arguments.demo = true;
        else throw std::invalid_argument("unknown argument: " + option);
    }
    if (arguments.demo && !arguments.input.empty()) throw std::invalid_argument("--demo and --input cannot be used together");
    if (arguments.demo && !arguments.sample_rate_set) arguments.options.sample_rate_hz = 1200.0;
    return arguments;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Arguments arguments = parse_arguments(argc, argv);
        if (arguments.help) {
            print_help(std::cout);
            return 0;
        }
        if (!arguments.demo && arguments.input.empty()) {
            print_help(std::cerr);
            return 2;
        }

        vmd::cli::CliTimings timings;
        const auto input_start = Clock::now();
        vmd::cli::InputSignal input = arguments.demo
            ? vmd::cli::demo_signal(arguments.options.sample_rate_hz)
            : vmd::cli::read_signal(arguments.input, arguments.column);
        timings.input_ms = std::chrono::duration<double, std::milli>(Clock::now() - input_start).count();

        const vmd::Result result = vmd::decompose(input.values, arguments.options);

        const auto serialization_start = Clock::now();
        vmd::cli::OutputPaths outputs = vmd::cli::write_csv_outputs(
            arguments.output_directory, input, result, arguments.options.sample_rate_hz);
        timings.serialization_ms = std::chrono::duration<double, std::milli>(Clock::now() - serialization_start).count();

        if (arguments.report) {
            const auto report_start = Clock::now();
            vmd::cli::write_html_report(outputs.report, input, arguments.options, result);
            timings.report_ms = std::chrono::duration<double, std::milli>(Clock::now() - report_start).count();
        }
        vmd::cli::write_summary_json(
            outputs.summary, input, arguments.options, result, timings, outputs, arguments.report);

        std::cout << std::setprecision(8)
                  << "VMD " << (result.diagnostics.converged ? "converged" : "did not converge")
                  << " after " << result.diagnostics.iterations << " iterations\n"
                  << "reconstruction error: " << result.diagnostics.reconstruction_error << '\n'
                  << "outputs: " << std::filesystem::absolute(arguments.output_directory).string() << '\n';
        return result.diagnostics.converged ? 0 : 4;
    } catch (const std::invalid_argument& error) {
        std::cerr << "argument error: " << error.what() << '\n';
        return 2;
    } catch (const vmd::cli::IoError& error) {
        std::cerr << "I/O error: " << error.what() << '\n';
        return 3;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 3;
    }
}
