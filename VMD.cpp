#include "VMD.h"
#include <fstream>
#include <iostream>
using namespace Eigen;
using namespace std;
/*
<VMD_CPP: C++ implementation of Variational Mode Decomposition using Eigen.>
Copyright (C) <2019>  <Lang He: asdsay@gmail.com>
Mozilla Public License v. 2.0.
*/

int main() {

	// create a signal to simulation the procedure.
	double f_1 = 2.0, f_2 = 24.0, f_3 = 288.0;
	int T = 1000;
	vectord t(T), v_1(T), v_2(T), v_3(T),signal(T);
	for (int i = 0; i < T; i++) {
		t[i] = double(i + 1) / T;
		v_1[i] = cos(2 * pI * f_1 * t[i]);
		v_2[i] = cos(2 * pI * f_2 * t[i]) / 4.0;
		v_3[i] = cos(2 * pI * f_3 * t[i]) / 16.0;
		signal[i] = v_1[i] + v_2[i] + v_3[i];
	}

	// initial some input parameters
	const double alpha = 2000.0, tau = 0, tol = 1e-7, eps = 2.2204e-16;
	const int K = 3, DC = 0, init = 1;
	const static double CSVFormat(4);

	// Example 1: If you want to get the full results as a 2D matrix of VMD. 	
	MatrixXd u, omega;
	MatrixXcd u_hat;
	VMD(u, u_hat, omega, signal, alpha, tau, K, DC, init, tol, eps);

	//Example 2: If you only wants to get sum result of the first n mode of signals.
	const double hp_cut_off = 200, lp_cut_off = 50; // Hz
	const double Fs = 50;
	/* Same as Exmaple 1
	MatrixXd u, omega;
	MatrixXcd u_hat;
	VMD(u, u_hat, omega, signal, alpha, tau, K, DC, init, tol, eps);*/
	double center_freq;
	MatrixXd p_data = u.row(0); p_data.fill(0);
	for (int k = 0; k < K; k++) {
		center_freq = Fs * omega(omega.rows() - 1, k);
		if (hp_cut_off < center_freq && center_freq < lp_cut_off)
			p_data = p_data + u.row(k);
	}


	for (int j = 0; j < T; j++) {
		if (p_data(0, j) < 0) {
			p_data = p_data.array() - p_data.minCoeff();
			break;
		}
	}

	/*
	// output example 
	cout << "Decomposition results" << endl;
	for (int i = 0; i < u.rows(); i++) {
		for (int j = 0; j < u.cols(); j++)
			cout << u(i, j) << ' ';
		cout << endl << endl;
	}
	string name = "test_result_1.csv";
	ofstream file1(name.c_str());
	file1 << u.format(CSVFormat);
	file1.close();
	*/
	return 0;
}; 