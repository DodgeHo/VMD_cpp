#pragma once
#include <vector>
#include <cmath>
#include "eigen3/Eigen/Eigen"
#include "eigen3/unsupported/Eigen/FFT"

#define pI acos(-1)
using namespace Eigen;
typedef std::vector<double> vectord;
typedef std::vector<std::complex<double> > vectorcd;
typedef std::vector<MatrixXcd> Matrix3DXd;

void VMD(MatrixXd& u, MatrixXcd& u_hat, MatrixXd& omega,
	vectord& signal, const double alpha, const double tau,
	const int K, const int DC, const int init, const double tol, const double eps);

vectorcd circshift(vectorcd& data, int offset);
template<typename T> void reverse(T& v, const int s, const int l);
vectord omega_init_method2(int K, const double fs);
MatrixXcd vector_to_MatrixXcd_in_row(vectorcd& Input);
MatrixXcd vector_to_MatrixXcd_in_col(vectorcd& Input);
vectorcd ExtractColFromMatrixXcd(MatrixXcd& Input, const int k, const int T);
vectorcd ExtractRowFromMatrixXd(MatrixXd& Input, const int k, const int T);
MatrixXcd sum(Matrix3DXd& u_hat_plus, const int n);


