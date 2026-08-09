#pragma once

#include <vmd/vmd.hpp>

#include <Eigen/Core>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace vmd::cli {

class IoError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct InputSignal {
    Eigen::VectorXd values;
    std::string source;
    std::string column_name;
    char delimiter = ',';
    bool has_header = false;
};

struct OutputPaths {
    std::filesystem::path modes;
    std::filesystem::path spectra;
    std::filesystem::path omega;
    std::filesystem::path summary;
    std::filesystem::path report;
};

struct CliTimings {
    double input_ms = 0.0;
    double serialization_ms = 0.0;
    double report_ms = 0.0;
};

InputSignal read_signal(const std::string& path, const std::string& column);
InputSignal demo_signal(double sample_rate_hz, std::size_t sample_count = 1200);
std::vector<double> parse_alpha_list(const std::string& text);
std::string initialization_name(Initialization initialization);

OutputPaths write_csv_outputs(
    const std::filesystem::path& output_directory,
    const InputSignal& input,
    const Result& result,
    double sample_rate_hz);

void write_summary_json(
    const std::filesystem::path& path,
    const InputSignal& input,
    const Options& options,
    const Result& result,
    const CliTimings& timings,
    const OutputPaths& outputs,
    bool report_generated);

}  // namespace vmd::cli
