#pragma once

#include <vmd/vmd.hpp>

#include <complex>
#include <vector>

using vectord = std::vector<double>;
using vectorcd = std::vector<std::complex<double>>;

/**
 * Compatibility wrapper for the original VMD_cpp API.
 *
 * `u` is K x N, `u_hat` is N x K, and `omega` contains normalized
 * center frequencies in cycles per sample. New code should prefer
 * `vmd::decompose()`.
 */
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
    double eps);

vectorcd circshift(const vectorcd& data, int offset);
