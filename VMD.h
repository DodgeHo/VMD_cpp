#pragma once
//#pragma warning(disable: 4267)
//#pragma warning(disable: 4244)
//#pragma warning(disable: 26451)
#include <iostream>
#include <vector>
#include <ctime>
#include <exception>
#include <math.h>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"
#include "eigen3/unsupported/Eigen/FFT"


#define eps_for_VMD 2.2204e-16
using namespace Eigen;
typedef std::vector<double> vectord;
typedef std::vector<std::complex<double> > vectorcd;
typedef std::vector<MatrixXcd> Matrix3DXd;

void VMD(MatrixXd& u, MatrixXcd& u_hat, MatrixXd& omega,
	vectord& signal,
	double alpha, double tau, const int K, const int DC, const int init, double tol);

vectorcd circshift(vectorcd& data, int offset);

template<typename T>
void reverse(T& v, int s, int l);

vectord omega_init_method2(int K, double fs);

MatrixXcd vector_to_MatrixXcd_in_row(vectorcd& Input);
MatrixXcd vector_to_MatrixXcd_in_col(vectorcd& Input);

vectorcd ExtractColFromMatrixXcd(MatrixXcd& Input, int k, int T);
vectorcd ExtractRowFromMatrixXd(MatrixXd& Input, int k, int T);

MatrixXcd sum(Matrix3DXd& u_hat_plus, int n);


