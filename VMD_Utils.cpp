#include "VMD.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

void VMD(
    Eigen::MatrixXd& u,
    Eigen::MatrixXcd& u_hat,
    Eigen::MatrixXd& omega,
    const vectord& signal,
    double alpha,
    double tau,
    int K,
    int DC,
    int init,
    double tol,
    double eps) {
    if (K <= 0) throw std::invalid_argument("K must be positive");
    if (init < 0 || init > 2) throw std::invalid_argument("init must be 0, 1, or 2");
    vmd::Options options;
    options.mode_count = static_cast<std::size_t>(K);
    options.alpha = {alpha};
    options.tau = tau;
    options.dc_mode = DC != 0;
    options.tolerance = tol;
    options.max_iterations = 500;
    options.keep_omega_history = true;
    options.initialization = init == 1 ? vmd::Initialization::Uniform
                           : init == 2 ? vmd::Initialization::Random
                                       : vmd::Initialization::Zero;
    (void)eps;

    const Eigen::Map<const Eigen::VectorXd> mapped(signal.data(), static_cast<Eigen::Index>(signal.size()));
    vmd::Result result = vmd::decompose(mapped, options);
    u = std::move(result.modes);
    u_hat = std::move(result.spectra);
    omega = std::move(result.omega_normalized);
}

vectorcd circshift(const vectorcd& data, int offset) {
    if (data.empty()) return {};
    const int size = static_cast<int>(data.size());
    int shift = offset % size;
    if (shift < 0) shift += size;
    vectorcd result(data.size());
    for (int i = 0; i < size; ++i) result[static_cast<std::size_t>((i + shift) % size)] = data[static_cast<std::size_t>(i)];
    return result;
}
