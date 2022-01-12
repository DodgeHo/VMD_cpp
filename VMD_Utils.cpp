#include "VMD.h"
using namespace Eigen;
using namespace std;

void VMD
	(MatrixXd& u, MatrixXcd& u_hat, MatrixXd& omega,
	vectord& signal, const double alpha, const double tau,
	const int K, const int DC, const int init, const double tol, const double eps) {
	/* ---------------------
	    
	Output:
	-------
	u - the collection of decomposed modes (2D double Matrix in Eigen -MatrixXd)
	u_hat - spectra of the modes (2D complex<double> Matrix in Eigen -MatrixXd)
	omega - estimated mode center - frequencies (2D double Matrix in Eigen -MatrixXd)
	-------
	Input:
	-------
	signal - the time domain signal(1D vector) to be decomposed
	alpha - the balancing parameter of the data - fidelity constraint
	tau - time - step of the dual ascent(pick 0 for noise - slack)
	K - the number of modes to be recovered
	DC - true if the first mode is putand kept at DC(0 - freq)
	init - 0 = all omegas start at 0
		                1 = all omegas start uniformly distributed
		                2 = all omegas initialized randomly
	tol - tolerance of convergence criterion; typically around 1e-6
		
	*/

	// ----------Preparations
	// Periodand sampling frequency of input signal
	int T = int(signal.size());
	int saveT = T;
	double fs = 1.0 / T;
	
	//extend the signal by mirroring
	vectord  f(2 * T, 0.0);
	copy(signal.begin(), signal.end(), f.begin() + T / 2);
	for (int i = 0; i < T / 2; i++)
		f[i] = signal[T / 2 - 1 - i];
	for (int i = 3 * T / 2; i < 2 * T; i++)
		f[i] = signal[T + 3 * T / 2 - 1 - i];

	// Time Domain 0 to T (of mirrored signal)
	// Spectral Domain discretization
	T = int(f.size());
	vectorcd freqs(T, 0.0);
	vectord timevec(T, 0.0);
	for (int i = 0; i < T; i ++) {
		timevec[i] = double(i + 1.0) / T;
		freqs[i] = (timevec[i] - 0.5) - double(1 / T);
	}

	// Maximum number of iterations(if not converged yet, then it won't anyway)
	int N = 500;

	// Construct and center f_hat
	vectorcd freqvec(T, 0.0);
	FFT<double> fft; fft.fwd(freqvec,f);
	vectorcd f_hat = circshift(freqvec, T / 2);
	vectorcd f_hat_plus(f_hat.size(), 0.0);
	copy(f_hat.begin() + T / 2, f_hat.end(), f_hat_plus.begin() + T / 2);

	// matrix keeping track of every iterant // could be discarded for mem
	Matrix3DXd u_hat_plus(N, MatrixXcd::Zero(K, T));

	// Initialization of omega_k
	MatrixXcd omega_plus = MatrixXcd::Zero(N, K);
	vectord tmp;
	switch (init) {
	case 1:
		for (int i = 0; i < K; i++){
			omega_plus(0, i) = double(0.5 / K) * (i);
			for (int j = 1; j < N; j++)
				omega_plus(j, i) = 0.0;
		}
		break;
	case 2:
		tmp = omega_init_method2(K, fs);
		for (int i = 0; i < K; i++) {
			omega_plus(0, i) = tmp[i];
			for (int j = 1; j < N; j++)
				omega_plus(j, i) = 0.0;
		}
		break;
	default:
		break;
	}

	//% if DC mode imposed, set its omega to 0
	if (DC)
		omega_plus(0, 0) = 0;

	// start with empty dual variables
	MatrixXcd lambda_hat = MatrixXcd::Zero(N, T);

	// other inits
	double uDiff = tol + eps;//% update step
	int n = 1;// loop counter
	MatrixXcd sum_uk = MatrixXcd::Zero(1, T);
	// accumulator
	int k ;
	//vectord sum_uk(freqs.size());


	// ----------- Main loop for iterative updates
	while (uDiff > tol && n < N) {

		//update first mode accumulator
		k = 0;
		sum_uk = u_hat_plus[n - 1].row(K - 1) + sum_uk - u_hat_plus[n-1].row(0);
		
		//update spectrum of first mode through Wiener filter of residuals
		MatrixXcd Dividend_vec = vector_to_MatrixXcd_in_col(f_hat_plus) - sum_uk - (lambda_hat.row(n - 1) / 2.0);
		MatrixXcd Divisor_vec = (1 + alpha *
			((vector_to_MatrixXcd_in_col(freqs).array() - omega_plus(n - 1, k))).array().square());
		u_hat_plus[n].row(k) = Dividend_vec.cwiseQuotient(Divisor_vec);

		//update first omega if not held at 0
		if (!DC) {
			std::complex<double> Dividend{ 0,0 }, Divisor{ 0, 0 }, Addend{ 0, 0 };
			for (int i = 0; i < T - T / 2; i++) {
				Addend = abs(u_hat_plus[n](k, T / 2 + i))* abs(u_hat_plus[n](k, T / 2 + i));
				Divisor += Addend;
				Dividend += freqs[T / 2 + i] * Addend;
			}
			omega_plus(n, k) = Dividend/ Divisor;
			
		}
		// Dual ascent

		for (k = 1; k < K ; k++) {
			//accumulator
			sum_uk = u_hat_plus[n].row(k - 1) + sum_uk - u_hat_plus[n - 1].row(k);

			//mode spectrum
			MatrixXcd Dividend_vec = vector_to_MatrixXcd_in_col(f_hat_plus) - sum_uk - (lambda_hat.row(n - 1) / 2.0);
			MatrixXcd Divisor_vec = (1 + alpha *
				((vector_to_MatrixXcd_in_col(freqs).array() - omega_plus(n - 1, k))).array().square());
			u_hat_plus[n].row(k) = Dividend_vec.cwiseQuotient(Divisor_vec);

			//center frequencies
			std::complex<double> Dividend{ 0,0 }, Divisor{ 0, 0 }, Addend{ 0, 0 };
			for (int i = 0; i < T - T / 2; i++) {
				Addend = abs(u_hat_plus[n](k, T / 2 + i))* abs(u_hat_plus[n](k, T / 2 + i));
				Divisor += Addend;
				Dividend += freqs[T / 2 + i] * Addend;
			}
			omega_plus(n, k) = Dividend/ Divisor;
		}
	
		lambda_hat.row(n) = lambda_hat.row(n - 1) +	tau * 
			(sum(u_hat_plus, n) - vector_to_MatrixXcd_in_col(f_hat_plus));
		n++;
		//uDiff = eps;

		std::complex<double> acc{ eps, 0 };
		for (int i = 0; i < K; i++) {
			MatrixXcd tmp = u_hat_plus[n-1].row(i) - u_hat_plus[n-2].row(i);
			tmp =  (tmp * (tmp.adjoint()));
			acc = acc + tmp(0,0) / double(T);

		}
		uDiff = abs(acc);

	}

	// ------ Postprocessing and cleanup

	// discard empty space if converged early
	N = std::min(N, n);
	omega = omega_plus.topRows(N).real();

	//Signal reconstruction
	u_hat = MatrixXcd::Zero(T, K);
	for (int i = T / 2; i < T; i++)
		for (int k = 0; k < K; k++)
			u_hat(i, k) = u_hat_plus[N-1](k, i);
	
	for (int i = T / 2; i >= 0; i--)
		for (int k = 0; k < K; k++) 
			u_hat(i, k) = conj(u_hat_plus[N - 1](k, T - i - 1));

			
	u_hat.row(0) = u_hat.row(T - 1).transpose().adjoint();
	u.resize(K, saveT);
	vectord result_col;
	for (int k = 0; k < K; k++) {
		vectorcd u_hat_col = ExtractColFromMatrixXcd(u_hat, k, T);
		u_hat_col = circshift(u_hat_col, int(floor(T / 2)));
		fft.inv(result_col, u_hat_col);
		for (int t = 0; t < saveT; t++)
			u(k, t) = result_col[t + T / 4];
	}


	u_hat.fill(0);
	vectord result_timevec(saveT, 0);
	for (int i = 0; i < saveT; i += 1) {
		result_timevec[i] = double(i + 1) / saveT;
	}

	for (int k = 0; k < K; k++) {
		vectorcd u_row = ExtractRowFromMatrixXd(u, k, saveT);
		fft.inv(result_timevec, u_row);
		u_row = circshift(u_row, saveT / 2);
		for (int t = 0; t < saveT; t++)
			u_hat(t,k) = u_row[t].real();
	}


	return;
}

#pragma region Ancillary Functions 

vectorcd circshift(vectorcd& data, int offset){
	int n = int(data.size());
	if (offset == 0) {
		vectorcd out_data(data);
		return out_data;
	}
	else{
		if (offset > 0) offset = n - offset;		// move to right by offset positions
		else              offset = -offset;			// move to left by offset positions
		vectorcd out_data(data.begin() + offset, data.end());
		out_data.insert(out_data.end(), data.begin(), data.begin() + offset);
		return out_data;
	}
}

vectord omega_init_method2(int K, const double fs) {
	vectord res(K, 0);
	int N = INT_MAX/2;
	srand(int(time(NULL)));
	for (int i = 0; i < K; i++) {
		res[i] = exp(log(fs) + (log(0.5) - log(fs)) *
			(rand() % (N + 1) / (float)(N + 1))
		);
	}
	sort(res.begin(), res.end());
	return res;
}

MatrixXcd vector_to_MatrixXcd_in_col(vectorcd& Input) {
	std::complex<double>* dataPtr = &Input[0];
	Eigen::MatrixXcd copiedMatrix = Eigen::Map<Eigen::MatrixXcd>(dataPtr, 1, int(Input.size()));
	return copiedMatrix;
}

vectorcd ExtractColFromMatrixXcd(MatrixXcd& Input, const int ColIdx, const int RowNum) {
	vectorcd Output(RowNum, 0);
	for (int i = 0; i < RowNum; ++i)
		Output[i] = Input(i, ColIdx);
	return Output;
}

vectorcd ExtractRowFromMatrixXd(MatrixXd& Input, const int RowIdx, const int ColNum) {
	vectorcd Output(ColNum, 0);
	for (int i = 0; i < ColNum; ++i)
		Output[i] = Input(RowIdx, i);
	return Output;
}

MatrixXcd sum(Matrix3DXd& u_hat_plus, const int n) {
	MatrixXcd cov = u_hat_plus[n];
	return cov.colwise().sum();
}

#pragma endregion

