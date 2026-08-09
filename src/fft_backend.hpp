#pragma once

#include <complex>
#include <memory>
#include <vector>

namespace vmd::detail {

class FFTBackend {
public:
    FFTBackend();
    ~FFTBackend();
    FFTBackend(FFTBackend&&) noexcept;
    FFTBackend& operator=(FFTBackend&&) noexcept;
    FFTBackend(const FFTBackend&) = delete;
    FFTBackend& operator=(const FFTBackend&) = delete;

    void forward_real(
        const std::vector<double>& input,
        std::vector<std::complex<double>>& output);
    void forward(
        const std::vector<std::complex<double>>& input,
        std::vector<std::complex<double>>& output);
    void inverse(
        const std::vector<std::complex<double>>& input,
        std::vector<std::complex<double>>& output);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace vmd::detail
