# VMD_cpp
# C++ implementation of Variational Mode Decomposition using Eigen3
Written by: Dodge(Lang HE) asdsay@gmail.com
Updated date: 2023-03-24

VMD, aka Variational Mode Decomposition, is a signal processing tool that decompse the input signal into different band-limited IMFs. 
Project **VMD_cpp**  is an imitation of [that in MATLAB](https://ww2.mathworks.cn/help/wavelet/ref/vmd.html). In this project, I used eigen3 to refactor VMD in C++, so that we can use it without MATLAB. 
Detail input and output please check out function **VMD** in file [VMD_Utils.cpp](https://github.com/DodgeHo/VMD_cpp/blob/master/VMD_Utils.cpp).

This sample code was written in MSBuild. You can both run in Visual Studio 2022 or MSVC or CMAKE/GCC, you can use either the sln project file, or the CMakeList.txt, they both work.

If you are looking for document to describe Variational mode decomposition, please turn to the original paper [Variational Mode Decomposition](https://ieeexplore.ieee.org/document/6655981). You can also find the MATLAB codes here.


Updated 2022-06-23: This VMD runs too slow. I tried to use OpenMP to make Eigen in parallel computing, but it couldn't work. Need to find another way in future.
