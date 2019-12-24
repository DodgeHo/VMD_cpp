#include "VMD1DProcess.h"
MatrixXd VMD1DProcess(vectord data, double Fs, double hp_cut_off, double lp_cut_off, int mode_num, int init) {
	int L = data.size();
	double alpha = 2000;// moderate bandwidth constraint 5x5: 2000default: 2000
	double 	tau = 0;// noise - tolerance
	int DC = 0;// no DC part imposed
	double tol = 1e-7;
	int K = mode_num;//  mode_num;

	MatrixXd u, omega;
	MatrixXcd u_hat;
	VMD(u, u_hat, omega, data, alpha, tau, K, DC, init, tol);
	double center_freq;
	MatrixXd p_data = u.row(0); p_data.fill(0);
	for (int k = 0; k < K; k++) {
		center_freq = Fs * omega(omega.rows() - 1, k);
		if (hp_cut_off < center_freq && center_freq < lp_cut_off)
			p_data = p_data + u.row(k);
	}


	for (int j = 0; j < L; j++) {
		if (p_data(0, j) < 0) {
			p_data = p_data.array() - p_data.minCoeff();
			break;
		}
	}
	return p_data;
}