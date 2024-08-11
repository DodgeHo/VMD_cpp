# VMD_cpp

# C++ implementation of Variational Mode Decomposition using Eigen3
Written by: Dodge(Lang HE) asdsay@gmail.com
Updated date: 2024-08-11

VMD, aka Variational Mode Decomposition, is a signal processing tool that decompse the input signal into different band-limited IMFs. 
Project **VMD_cpp**  is an imitation of [that in MATLAB](https://ww2.mathworks.cn/help/wavelet/ref/vmd.html). In this project, I used eigen3 to refactor VMD in C++, so that we can use it without MATLAB. 
Detail input and output please check out function **VMD** in file [VMD_Utils.cpp](https://github.com/DodgeHo/VMD_cpp/blob/master/VMD_Utils.cpp).

This sample code was written in MSBuild. You can both run in Visual Studio 2022 or MSVC or CMAKE/GCC, you can use either the sln project file, or the CMakeList.txt, they both work.

If you are looking for document to describe Variational mode decomposition, please turn to the original paper [Variational Mode Decomposition](https://ieeexplore.ieee.org/document/6655981). You can also find the MATLAB codes here.

This VMD runs too slow. I tried to use Multiple Threads in Eigen, and I did everything I can. Hopefully this is good enough for someone using it.


# VMD（变分模态分解）的C++实现，使用了Eigen3

作者：Dodge asdsay@gmail.com 
更新日期：2024-08-11

VMD（变分模态分解）是一种信号处理算法，可以将输入信号分解为不同带限的内禀模态函数（IMFs）。
本项目**VMD_cpp** 是参考于[其在MATLAB中的实现](https://ww2.mathworks.cn/help/wavelet/ref/vmd.html)。在项目中，借助eigen3来实现C++中的VMD，从而无须MATLAB的计算。详细的输入输出，可以查看[VMD_Utils.cpp](https://github.com/DodgeHo/VMD_cpp/blob/master/VMD_Utils.cpp)文件中的 **VMD**函数。

本项目用MSBuild编写，也可在Visual Studio 2022、MSVC或CMAKE/GCC中运行，sln项目文件或CMakeList.txt都能用。

如果需要描述变分模态分解的文档，可参阅原始论文[Variational Mode Decomposition](https://ieeexplore.ieee.org/document/6655981)。

这个项目稍慢。我尝试加入了Eigen自带的多线程计算，并尽可能在原代码上修改，尽可能让Eigen加快，希望效果足够好吧。
