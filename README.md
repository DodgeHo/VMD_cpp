# VMD_cpp
# C++ implementation of Variational Mode Decomposition using Eigen3
Written by: Dodge(Lang HE) asdsay@gmail.com
Updated date: 2022-06-23

VMD, aka Variational Mode Decomposition, is a signal processing tool that decompse the input signal into different band-limited IMFs. 
Project **VMD_cpp**  is an imitation of [that in MATLAB](https://ww2.mathworks.cn/help/wavelet/ref/vmd.html). In this project, I used eigen3 to refactor VMD in C++, so that we can use it without MATLAB. 
Detail input and output please check out function **VMD** in file [VMD_Utils.cpp](https://github.com/DodgeHo/VMD_cpp/blob/master/VMD_Utils.cpp).

This sample code was written in MSBuild. Please run in Visual Studio 2015 or newer version. If you want to run in MSVC or CMAKE/GCC, please change the "include" method from CPP file to Eigen according to the version in [20211128](https://github.com/DodgeHo/VMD_cpp/blob/7bdc0ea6702175bc81c7a48027a706dd6d369b6a/VMD.cpp).

If you are looking for document to describe Variational mode decomposition, please turn to the original paper [Variational Mode Decomposition](https://ieeexplore.ieee.org/document/6655981). You can also find the MATLAB codes here.


Updated 2022-06-23: This VMD runs too slow. I tried to use OpenMP to make Eigen in parallel computing, but it couldn't work. Need to find another way in future.
