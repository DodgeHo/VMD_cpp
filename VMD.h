#pragma once
#include <vector>
#include <cmath>
#include <ctime>
#include "eigen/Eigen/Eigen"
#include "eigen/unsupported/Eigen/FFT"
//#include <eigen/Eigen/Core>
//#include <eigen/unsupported/Eigen/FFT>

#define pI acos(-1)
using namespace Eigen;
typedef std::vector<double> vectord;
typedef std::vector<std::complex<double> > vectorcd;
typedef std::vector<MatrixXcd> Matrix3DXd;

void VMD(MatrixXd& u, MatrixXcd& u_hat, MatrixXd& omega,
	vectord& signal, const double alpha, const double tau,
	const int K, const int DC, const int init, const double tol, const double eps);

vectorcd circshift(vectorcd& data, int offset);
vectord omega_init_method2(int K, const double fs);
vectorcd ExtractColFromMatrixXcd(MatrixXcd& Input, const int k, const int T);
vectorcd ExtractRowFromMatrixXd(MatrixXd& Input, const int k, const int T);