#include "cli_io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string_view>

namespace vmd::cli {
namespace {

std::string trim(std::string value) {
    const auto not_space = [](unsigned char character) { return std::isspace(character) == 0; };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
}

std::vector<std::string> parse_row(const std::string& line, char delimiter, const std::string& source, std::size_t line_number) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char character = line[i];
        if (character == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (character == delimiter && !quoted) {
            fields.push_back(trim(field));
            field.clear();
        } else {
            field.push_back(character);
        }
    }
    if (quoted) {
        throw IoError(source + ":" + std::to_string(line_number) + ": unterminated quoted field");
    }
    fields.push_back(trim(field));
    return fields;
}

std::size_t delimiter_count(const std::string& line, char delimiter) {
    bool quoted = false;
    std::size_t count = 0;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') ++i;
            else quoted = !quoted;
        } else if (!quoted && line[i] == delimiter) {
            ++count;
        }
    }
    return count;
}

std::optional<double> number(const std::string& field) {
    if (field.empty()) return std::nullopt;
    std::size_t consumed = 0;
    try {
        const double value = std::stod(field, &consumed);
        if (consumed != field.size() || !std::isfinite(value)) return std::nullopt;
        return value;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

bool is_index(const std::string& value) {
    return !value.empty() && std::all_of(value.begin(), value.end(), [](unsigned char character) {
        return std::isdigit(character) != 0;
    });
}

std::string json_escape(const std::string& value) {
    std::ostringstream output;
    for (unsigned char character : value) {
        switch (character) {
        case '"': output << "\\\""; break;
        case '\\': output << "\\\\"; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (character < 0x20) {
                output << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(character);
            } else {
                output << static_cast<char>(character);
            }
        }
    }
    return output.str();
}

std::ofstream open_output(const std::filesystem::path& path) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw IoError("cannot open output file: " + path.string());
    output << std::setprecision(17);
    return output;
}

std::string output_name(const std::filesystem::path& path) {
    return path.filename().string();
}

}  // namespace

InputSignal read_signal(const std::string& path, const std::string& column) {
    std::ifstream file;
    std::istream* stream = nullptr;
    const std::string source = path == "-" ? "<stdin>" : path;
    if (path == "-") {
        stream = &std::cin;
    } else {
        file.open(path, std::ios::binary);
        if (!file) throw IoError("cannot open input file: " + path);
        stream = &file;
    }

    std::vector<std::pair<std::size_t, std::string>> lines;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(*stream, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (trim(line).empty()) continue;
        lines.emplace_back(line_number, line);
    }
    if (lines.empty()) throw IoError(source + ": no data rows found");
    if (lines.front().second.size() >= 3 &&
        static_cast<unsigned char>(lines.front().second[0]) == 0xEF &&
        static_cast<unsigned char>(lines.front().second[1]) == 0xBB &&
        static_cast<unsigned char>(lines.front().second[2]) == 0xBF) {
        lines.front().second.erase(0, 3);
    }

    const std::size_t commas = delimiter_count(lines.front().second, ',');
    const std::size_t tabs = delimiter_count(lines.front().second, '\t');
    const char delimiter = tabs > commas ? '\t' : ',';
    auto first = parse_row(lines.front().second, delimiter, source, lines.front().first);
    if (first.empty()) throw IoError(source + ": first row is empty");

    std::size_t selected = 0;
    std::string selected_name;
    bool has_header = false;
    if (!column.empty() && is_index(column)) {
        try {
            selected = static_cast<std::size_t>(std::stoull(column));
        } catch (const std::exception&) {
            throw IoError("column index is out of range: " + column);
        }
        if (selected >= first.size()) {
            throw IoError(source + ": selected column " + std::to_string(selected) +
                          " exceeds row width " + std::to_string(first.size()));
        }
        has_header = !number(first[selected]).has_value();
    } else if (!column.empty()) {
        has_header = true;
        const auto found = std::find(first.begin(), first.end(), column);
        if (found == first.end()) throw IoError(source + ": header has no column named '" + column + "'");
        selected = static_cast<std::size_t>(std::distance(first.begin(), found));
        if (std::find(found + 1, first.end(), column) != first.end()) {
            throw IoError(source + ": duplicate header column '" + column + "'");
        }
    } else {
        has_header = !number(first.front()).has_value();
    }

    selected_name = has_header ? first[selected] : std::to_string(selected);

    const std::size_t data_begin = has_header ? 1 : 0;
    std::vector<double> values;
    values.reserve(lines.size() - data_begin);
    const std::size_t expected_width = first.size();
    for (std::size_t row_index = data_begin; row_index < lines.size(); ++row_index) {
        const auto& [physical_line, text] = lines[row_index];
        const auto fields = parse_row(text, delimiter, source, physical_line);
        if (fields.size() != expected_width) {
            throw IoError(source + ":" + std::to_string(physical_line) + ": mixed column width; expected " +
                          std::to_string(expected_width) + ", found " + std::to_string(fields.size()));
        }
        const auto parsed = number(fields[selected]);
        if (!parsed.has_value()) {
            throw IoError(source + ":" + std::to_string(physical_line) + ":" +
                          std::to_string(selected + 1) + ": missing, non-finite, or invalid number '" + fields[selected] + "'");
        }
        values.push_back(*parsed);
    }
    if (values.empty()) throw IoError(source + ": no numeric samples found");

    InputSignal result;
    result.values = Eigen::Map<Eigen::VectorXd>(values.data(), static_cast<Eigen::Index>(values.size()));
    result.source = source;
    result.column_name = selected_name;
    result.delimiter = delimiter;
    result.has_header = has_header;
    return result;
}

InputSignal demo_signal(double sample_rate_hz, std::size_t sample_count) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    InputSignal input;
    input.values.resize(static_cast<Eigen::Index>(sample_count));
    for (std::size_t i = 0; i < sample_count; ++i) {
        const double time = static_cast<double>(i) / sample_rate_hz;
        input.values(static_cast<Eigen::Index>(i)) =
            std::cos(2.0 * pi * 2.0 * time) +
            0.1 * std::cos(2.0 * pi * 24.0 * time) +
            0.005 * std::cos(2.0 * pi * 288.0 * time);
    }
    input.source = "demo";
    input.column_name = "value";
    input.has_header = true;
    return input;
}

std::vector<double> parse_alpha_list(const std::string& text) {
    std::vector<double> values;
    std::stringstream stream(text);
    std::string field;
    while (std::getline(stream, field, ',')) {
        const auto parsed = number(trim(field));
        if (!parsed.has_value() || *parsed <= 0.0) throw std::invalid_argument("--alpha expects positive comma-separated numbers");
        values.push_back(*parsed);
    }
    if (values.empty()) throw std::invalid_argument("--alpha must not be empty");
    return values;
}

std::string initialization_name(Initialization initialization) {
    switch (initialization) {
    case Initialization::Zero: return "zero";
    case Initialization::Uniform: return "uniform";
    case Initialization::Random: return "random";
    }
    return "unknown";
}

OutputPaths write_csv_outputs(
    const std::filesystem::path& output_directory,
    const InputSignal& input,
    const Result& result,
    double sample_rate_hz) {
    std::error_code error;
    std::filesystem::create_directories(output_directory, error);
    if (error) throw IoError("cannot create output directory '" + output_directory.string() + "': " + error.message());

    OutputPaths paths;
    paths.modes = output_directory / "modes.csv";
    paths.spectra = output_directory / "spectra.csv";
    paths.omega = output_directory / "omega.csv";
    paths.summary = output_directory / "summary.json";
    paths.report = output_directory / "report.html";

    {
        auto output = open_output(paths.modes);
        output << "index,time_s,signal";
        for (Eigen::Index k = 0; k < result.modes.rows(); ++k) output << ",imf_" << (k + 1);
        output << ",reconstruction,residual\n";
        for (Eigen::Index i = 0; i < input.values.size(); ++i) {
            const double reconstruction = result.modes.col(i).sum();
            output << i << ',' << static_cast<double>(i) / sample_rate_hz << ',' << input.values(i);
            for (Eigen::Index k = 0; k < result.modes.rows(); ++k) output << ',' << result.modes(k, i);
            output << ',' << reconstruction << ',' << input.values(i) - reconstruction << '\n';
        }
    }

    {
        auto output = open_output(paths.spectra);
        output << "index,frequency_hz";
        for (Eigen::Index k = 0; k < result.spectra.cols(); ++k) {
            output << ",imf_" << (k + 1) << "_real,imf_" << (k + 1) << "_imag,imf_" << (k + 1) << "_magnitude";
        }
        output << '\n';
        const Eigen::Index n = result.spectra.rows();
        const Eigen::Index center = n / 2;
        for (Eigen::Index i = 0; i < n; ++i) {
            const double frequency = static_cast<double>(i - center) * sample_rate_hz / static_cast<double>(n);
            output << i << ',' << frequency;
            for (Eigen::Index k = 0; k < result.spectra.cols(); ++k) {
                const auto value = result.spectra(i, k);
                output << ',' << value.real() << ',' << value.imag() << ',' << std::abs(value);
            }
            output << '\n';
        }
    }

    {
        auto output = open_output(paths.omega);
        output << "iteration";
        for (Eigen::Index k = 0; k < result.omega_hz.cols(); ++k) output << ",imf_" << (k + 1) << "_hz";
        output << '\n';
        for (Eigen::Index i = 0; i < result.omega_hz.rows(); ++i) {
            output << i;
            for (Eigen::Index k = 0; k < result.omega_hz.cols(); ++k) output << ',' << result.omega_hz(i, k);
            output << '\n';
        }
    }
    return paths;
}

void write_summary_json(
    const std::filesystem::path& path,
    const InputSignal& input,
    const Options& options,
    const Result& result,
    const CliTimings& timings,
    const OutputPaths& outputs,
    bool report_generated) {
    auto output = open_output(path);
    output << "{\n"
           << "  \"schema_version\": 1,\n"
           << "  \"input\": {\n"
           << "    \"source\": \"" << json_escape(input.source) << "\",\n"
           << "    \"column\": \"" << json_escape(input.column_name) << "\",\n"
           << "    \"samples\": " << input.values.size() << ",\n"
           << "    \"sample_rate_hz\": " << options.sample_rate_hz << ",\n"
           << "    \"delimiter\": \"" << (input.delimiter == '\t' ? "tab" : "comma") << "\",\n"
           << "    \"has_header\": " << (input.has_header ? "true" : "false") << "\n"
           << "  },\n"
           << "  \"options\": {\n"
           << "    \"mode_count\": " << options.mode_count << ",\n"
           << "    \"alpha\": [";
    for (std::size_t i = 0; i < options.alpha.size(); ++i) {
        if (i != 0) output << ", ";
        output << options.alpha[i];
    }
    output << "],\n"
           << "    \"tau\": " << options.tau << ",\n"
           << "    \"dc_mode\": " << (options.dc_mode ? "true" : "false") << ",\n"
           << "    \"initialization\": \"" << initialization_name(options.initialization) << "\",\n"
           << "    \"tolerance\": " << options.tolerance << ",\n"
           << "    \"max_iterations\": " << options.max_iterations << ",\n"
           << "    \"random_seed\": " << options.random_seed << ",\n"
           << "    \"keep_omega_history\": " << (options.keep_omega_history ? "true" : "false") << ",\n"
           << "    \"adaptive_tau\": " << (options.adaptive_tau ? "true" : "false") << "\n"
           << "  },\n"
           << "  \"diagnostics\": {\n"
           << "    \"converged\": " << (result.diagnostics.converged ? "true" : "false") << ",\n"
           << "    \"iterations\": " << result.diagnostics.iterations << ",\n"
           << "    \"final_mode_delta\": " << result.diagnostics.final_mode_delta << ",\n"
           << "    \"reconstruction_error\": " << result.diagnostics.reconstruction_error << ",\n"
           << "    \"primal_residual\": " << result.diagnostics.primal_residual << ",\n"
           << "    \"dual_residual\": " << result.diagnostics.dual_residual << ",\n"
           << "    \"final_tau\": " << result.diagnostics.final_tau << "\n"
           << "  },\n"
           << "  \"timings_ms\": {\n"
           << "    \"input\": " << timings.input_ms << ",\n"
           << "    \"mirror_fft\": " << result.timings.mirror_and_fft_ms << ",\n"
           << "    \"iteration\": " << result.timings.iteration_ms << ",\n"
           << "    \"spectrum_reconstruction\": " << result.timings.reconstruction_ms << ",\n"
           << "    \"serialization\": " << timings.serialization_ms << ",\n"
           << "    \"report\": " << timings.report_ms << "\n"
           << "  },\n"
           << "  \"outputs\": {\n"
           << "    \"modes\": \"" << json_escape(output_name(outputs.modes)) << "\",\n"
           << "    \"spectra\": \"" << json_escape(output_name(outputs.spectra)) << "\",\n"
           << "    \"omega\": \"" << json_escape(output_name(outputs.omega)) << "\",\n"
           << "    \"summary\": \"" << json_escape(output_name(outputs.summary)) << "\",\n"
           << "    \"report\": " << (report_generated ? "\"report.html\"" : "null") << "\n"
           << "  }\n"
           << "}\n";
}

}  // namespace vmd::cli
