#include "fft_backend.hpp"

#include <new>
#include <stdexcept>

#if defined(VMD_FFT_BACKEND_FFTW)
#include <fftw3.h>
#elif defined(VMD_FFT_BACKEND_MKL)
#include <mkl_dfti.h>
#else
#include <unsupported/Eigen/FFT>
#endif

namespace vmd::detail {

struct FFTBackend::Impl {
#if defined(VMD_FFT_BACKEND_FFTW)
    std::size_t complex_size = 0;
    fftw_complex* complex_input = nullptr;
    fftw_complex* complex_output = nullptr;
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;

    std::size_t real_size = 0;
    double* real_input = nullptr;
    fftw_complex* real_output = nullptr;
    fftw_plan real_forward_plan = nullptr;

    ~Impl() {
        reset_complex();
        reset_real();
    }

    void reset_complex() {
        if (forward_plan != nullptr) fftw_destroy_plan(forward_plan);
        if (inverse_plan != nullptr) fftw_destroy_plan(inverse_plan);
        if (complex_input != nullptr) fftw_free(complex_input);
        if (complex_output != nullptr) fftw_free(complex_output);
        complex_size = 0;
        complex_input = nullptr;
        complex_output = nullptr;
        forward_plan = nullptr;
        inverse_plan = nullptr;
    }

    void reset_real() {
        if (real_forward_plan != nullptr) fftw_destroy_plan(real_forward_plan);
        if (real_input != nullptr) fftw_free(real_input);
        if (real_output != nullptr) fftw_free(real_output);
        real_size = 0;
        real_input = nullptr;
        real_output = nullptr;
        real_forward_plan = nullptr;
    }

    void ensure_complex(std::size_t requested) {
        if (requested == complex_size) return;
        reset_complex();
        complex_size = requested;
        complex_input = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * complex_size));
        complex_output = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * complex_size));
        if (complex_input == nullptr || complex_output == nullptr) throw std::bad_alloc();
        forward_plan = fftw_plan_dft_1d(
            static_cast<int>(complex_size), complex_input, complex_output, FFTW_FORWARD, FFTW_ESTIMATE);
        inverse_plan = fftw_plan_dft_1d(
            static_cast<int>(complex_size), complex_input, complex_output, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (forward_plan == nullptr || inverse_plan == nullptr) {
            throw std::runtime_error("failed to create FFTW complex plan");
        }
    }

    void ensure_real(std::size_t requested) {
        if (requested == real_size) return;
        reset_real();
        real_size = requested;
        real_input = static_cast<double*>(fftw_malloc(sizeof(double) * real_size));
        real_output = static_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * (real_size / 2 + 1)));
        if (real_input == nullptr || real_output == nullptr) throw std::bad_alloc();
        real_forward_plan = fftw_plan_dft_r2c_1d(
            static_cast<int>(real_size), real_input, real_output, FFTW_ESTIMATE);
        if (real_forward_plan == nullptr) {
            throw std::runtime_error("failed to create FFTW real plan");
        }
    }
#elif defined(VMD_FFT_BACKEND_MKL)
    std::size_t complex_size = 0;
    DFTI_DESCRIPTOR_HANDLE complex_descriptor = nullptr;
    std::size_t real_size = 0;
    DFTI_DESCRIPTOR_HANDLE real_descriptor = nullptr;

    ~Impl() {
        if (complex_descriptor != nullptr) DftiFreeDescriptor(&complex_descriptor);
        if (real_descriptor != nullptr) DftiFreeDescriptor(&real_descriptor);
    }

    void ensure_complex(std::size_t requested) {
        if (requested == complex_size) return;
        if (complex_descriptor != nullptr) DftiFreeDescriptor(&complex_descriptor);
        complex_size = requested;
        if (DftiCreateDescriptor(
                &complex_descriptor, DFTI_DOUBLE, DFTI_COMPLEX, 1,
                static_cast<MKL_LONG>(complex_size)) != 0 ||
            DftiSetValue(complex_descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE) != 0 ||
            DftiSetValue(
                complex_descriptor, DFTI_BACKWARD_SCALE,
                1.0 / static_cast<double>(complex_size)) != 0 ||
            DftiCommitDescriptor(complex_descriptor) != 0) {
            throw std::runtime_error("failed to create MKL complex DFT descriptor");
        }
    }

    void ensure_real(std::size_t requested) {
        if (requested == real_size) return;
        if (real_descriptor != nullptr) DftiFreeDescriptor(&real_descriptor);
        real_size = requested;
        if (DftiCreateDescriptor(
                &real_descriptor, DFTI_DOUBLE, DFTI_REAL, 1,
                static_cast<MKL_LONG>(real_size)) != 0 ||
            DftiSetValue(real_descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE) != 0 ||
            DftiSetValue(real_descriptor, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX) != 0 ||
            DftiCommitDescriptor(real_descriptor) != 0) {
            throw std::runtime_error("failed to create MKL real DFT descriptor");
        }
    }
#else
    Eigen::FFT<double> fft;

    Impl() { fft.SetFlag(Eigen::FFT<double>::HalfSpectrum); }
#endif
};

FFTBackend::FFTBackend() : impl_(std::make_unique<Impl>()) {}
FFTBackend::~FFTBackend() = default;
FFTBackend::FFTBackend(FFTBackend&&) noexcept = default;
FFTBackend& FFTBackend::operator=(FFTBackend&&) noexcept = default;

void FFTBackend::forward_real(
    const std::vector<double>& input,
    std::vector<std::complex<double>>& output) {
    if (input.size() <= 1) {
        output.assign(input.size(), input.empty() ? std::complex<double>{} : std::complex<double>{input.front(), 0.0});
        return;
    }
#if defined(VMD_FFT_BACKEND_FFTW)
    impl_->ensure_real(input.size());
    output.resize(input.size() / 2 + 1);
    for (std::size_t i = 0; i < input.size(); ++i) impl_->real_input[i] = input[i];
    fftw_execute(impl_->real_forward_plan);
    for (std::size_t i = 0; i < output.size(); ++i) {
        output[i] = {impl_->real_output[i][0], impl_->real_output[i][1]};
    }
#elif defined(VMD_FFT_BACKEND_MKL)
    impl_->ensure_real(input.size());
    output.resize(input.size() / 2 + 1);
    if (DftiComputeForward(
            impl_->real_descriptor, const_cast<double*>(input.data()), output.data()) != 0) {
        throw std::runtime_error("MKL real forward FFT failed");
    }
#else
    impl_->fft.fwd(output, input);
#endif
}

void FFTBackend::forward(
    const std::vector<std::complex<double>>& input,
    std::vector<std::complex<double>>& output) {
    if (input.size() <= 1) {
        output = input;
        return;
    }
#if defined(VMD_FFT_BACKEND_FFTW)
    impl_->ensure_complex(input.size());
    output.resize(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        impl_->complex_input[i][0] = input[i].real();
        impl_->complex_input[i][1] = input[i].imag();
    }
    fftw_execute(impl_->forward_plan);
    for (std::size_t i = 0; i < output.size(); ++i) {
        output[i] = {impl_->complex_output[i][0], impl_->complex_output[i][1]};
    }
#elif defined(VMD_FFT_BACKEND_MKL)
    impl_->ensure_complex(input.size());
    output.resize(input.size());
    if (DftiComputeForward(
            impl_->complex_descriptor, const_cast<std::complex<double>*>(input.data()), output.data()) != 0) {
        throw std::runtime_error("MKL forward FFT failed");
    }
#else
    impl_->fft.fwd(output, input);
#endif
}

void FFTBackend::inverse(
    const std::vector<std::complex<double>>& input,
    std::vector<std::complex<double>>& output) {
    if (input.size() <= 1) {
        output = input;
        return;
    }
#if defined(VMD_FFT_BACKEND_FFTW)
    impl_->ensure_complex(input.size());
    output.resize(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        impl_->complex_input[i][0] = input[i].real();
        impl_->complex_input[i][1] = input[i].imag();
    }
    fftw_execute(impl_->inverse_plan);
    const double scale = 1.0 / static_cast<double>(input.size());
    for (std::size_t i = 0; i < output.size(); ++i) {
        output[i] = {
            impl_->complex_output[i][0] * scale,
            impl_->complex_output[i][1] * scale};
    }
#elif defined(VMD_FFT_BACKEND_MKL)
    impl_->ensure_complex(input.size());
    output.resize(input.size());
    if (DftiComputeBackward(
            impl_->complex_descriptor, const_cast<std::complex<double>*>(input.data()), output.data()) != 0) {
        throw std::runtime_error("MKL inverse FFT failed");
    }
#else
    impl_->fft.inv(output, input);
#endif
}

}  // namespace vmd::detail
