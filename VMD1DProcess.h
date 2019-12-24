#pragma once
#include "VMD.h"
MatrixXd VMD1DProcess(vectord data, double Fs, double hp_cut_off, double lp_cut_off, int mode_num, int init);