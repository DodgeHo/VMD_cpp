﻿#include "VMD.h"
#include <fstream>
#include <iostream>
#include <thread>
using namespace Eigen;
using namespace std;
/*
<VMD_CPP: C++ implementation of Variational Mode Decomposition using Eigen.>
Copyright (C) <2019>  <Lang HE: asdsay@gmail.com>
Mozilla Public License v. 2.0.
*/

void printMatrix(const MatrixXd& u) {
	std::ostringstream out; // use ostringstream to accumulate output
	for (int i = 0; i < u.rows(); i++) {
		for (int j = 0; j < u.cols(); j++)
			out << u(i, j) << ' ';
		out << "\n\n";
	}
	std::cout << out.str(); // output once
}

int main() {

	// create a signal to simulation the procedure.
	double f_1 = 2.0, f_2 = 24.0, f_3 = 288.0;
	int T = 1200;
	vectord t(T), v_1(T), v_2(T), v_3(T),signal(T);
	for (int i = 0; i < T; i++) {
		t[i] = double(i + 1) / T;
		v_1[i] = cos(2 * pI * f_1 * t[i]);
		v_2[i] = cos(2 * pI * f_2 * t[i]) / 10.0;
		v_3[i] = cos(2 * pI * f_3 * t[i]) / 200.0;
		signal[i] = v_1[i] + v_2[i] + v_3[i];
	}

	// initial some input parameters
	const double alpha = 50.0, tau = 0, tol = 1e-7, eps = 2.2204e-16;
	const int K = 8, DC = 0, init = 1;
	const static double CSVFormat(4);
	Eigen::setNbThreads(std::thread::hardware_concurrency()); // Set the numbers of threads that Eigen uses

	// Example 1: If you want to get the full results as a 2D matrix of VMD. 	
	MatrixXd u, omega;
	MatrixXcd u_hat;
	VMD(u, u_hat, omega, signal, alpha, tau, K, DC, init, tol, eps);

	//Example 2: If you only wants to get sum result of the first n mode of signals.
	const double hp_cut_off = 1, lp_cut_off = 15; // Hz
	const double Fs = 50;
	/* Same as Example 1
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

	
	// Output results
	cout << "VMD Decomposition Results:" << endl;
	cout << "Number of modes: " << K << endl;
	cout << "Matrix dimensions: " << u.rows() << " x " << u.cols() << endl;
	
	// Save decomposition results to CSV
	const string output_filename = "vmd_decomposition_results.csv";
	ofstream output_file(output_filename);
	if (!output_file) {
		cerr << "Error: Could not open file " << output_filename << endl;
		return 1;
	}
	output_file << u.format(CSVFormat);
	output_file.close();
	cout << "Results saved to: " << output_filename << endl;

	return 0;
}; 