#include <vmd/vmd.hpp>

#include <Eigen/Core>

#include <iostream>

int main() {
    Eigen::VectorXd signal(8);
    signal << 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0;

    vmd::Options options;
    options.mode_count = 1;
    options.alpha = {100.0};
    options.max_iterations = 50;

    const vmd::Result result = vmd::decompose(signal, options);
    if (result.modes.rows() != 1 || result.modes.cols() != signal.size() ||
        result.spectra.rows() != signal.size() || result.spectra.cols() != 1) {
        std::cerr << "unexpected installed-package result shape\n";
        return 1;
    }
    std::cout << "installed package consumer passed\n";
    return 0;
}
