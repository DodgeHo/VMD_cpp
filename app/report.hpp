#pragma once

#include "cli_io.hpp"

#include <filesystem>

namespace vmd::cli {

void write_html_report(
    const std::filesystem::path& path,
    const InputSignal& input,
    const Options& options,
    const Result& result);

}  // namespace vmd::cli
