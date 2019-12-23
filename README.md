# VMD_Cpp

Variational Mode Decomposition for Cpp using Eigen.

This is cpp realization for Variatioanl Mode Decomposition

Coding by: Lang He (asdsay@gmail.com)

Refering to the paper:
K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans. on Signal Processing (in press) 



##  Output:

u - the collection of decomposed modes (2D double Matrix in Eigen -MatrixXd)

u_hat - spectra of the modes (2D complex<double> Matrix in Eigen -MatrixXd)

omega - estimated mode center - frequencies (2D double Matrix in Eigen -MatrixXd)

##  Input:

signal - the time domain signal(1D vector<double>) to be decomposed

alpha - the balancing parameter of the data - fidelity constraint

tau - time - step of the dual ascent(pick 0 for noise - slack)

K - the number of modes to be recovered

DC - true if the first mode is putand kept at DC(0 - freq)

init - 0 = all omegas start at 0
         1 = all omegas start uniformly distributed
         2 = all omegas initialized randomly

tol - tolerance of convergence criterion; typically around 1e-6

