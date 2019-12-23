#include "VMD.h"
#include <fstream>
using namespace Eigen;
using namespace std;
const static double CSVFormat(4);
#define pi 3.1415926
int main() {
	double f_1 = 2.0, f_2 = 24.0, f_3 = 288.0;
	int T = 1000;
	vectord t(T), v_1(T), v_2(T), v_3(T),signal(T);
	for (int i = 0; i < T; i++) {
		t[i] = double(i + 1) / T;
		v_1[i] = cos(2 * pi * f_1 * t[i]);
		v_2[i] = cos(2 * pi * f_2 * t[i]) / 4.0;
		v_3[i] = cos(2 * pi * f_3 * t[i]) / 16.0;
		signal[i] = v_1[i] + v_2[i] + v_3[i];
	}


	double alpha = 2000.0, tau = 0, tol = 1e-7;
	int K = 3, DC = 0, init = 1;
		
	MatrixXd u, omega;
	MatrixXcd u_hat;
	VMD(u, u_hat, omega, signal, alpha, tau, K, DC, init, tol);
	/*
	for (int i = 0; i < u.rows(); i++) {
		for (int j = 0; j < u.cols(); j++)
			cout << u(i, j) << ' ';
		cout << endl;
	}*/
		
	std::string name = "test_result.csv";
	ofstream file(name.c_str());
	file << u.format(CSVFormat);
	file.close();
	system("PAUSE");
	return 0;

}; 