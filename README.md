# Implicit time-space-domain staggered-grid finite-difference scheme for 2D acoustic wave equation  
## Introduction
  &nbsp;&nbsp;&nbsp;&nbsp;This GitHub project contains the core computational code for the optimized implicit time-space-domain staggered-grid finite-difference scheme. I wrote it in C/C++ on the Visual Studio 2022 platform. Each folder includes a `time_space_domain.cpp` and `time_space_domain.h` to store the corresponding functions. Additionally, each folder contains a `main.cpp` file for running the program.
  The dispersion folder contains C/C++ code for calculating phase velocity dispersion and relative error. The default parameters are set as M=6, N=3, and Courant number r=0.15. The generated data is stored in txt files.
  In the stability folder, the stability factor is directly calculated in main.cpp and saved into a txt file.
  In ITSDSGFD_Variable_Operator_length and ITSDSGFD_Fixed_Order, the program first calculates the difference coefficient table and operator length table, and then performs forward modeling. The simulation results are stored in binary .dat files with the data type of float. Each folder contains all the files for implementing the corresponding functions, and the program can be run by calling the files on any platform that supports C++ operation.
## Demo
 &nbsp;&nbsp;&nbsp;&nbsp;Here, taking the Demo as an example, I will introduce in as much detail as possible the specific steps and final results of the implicit space-time domain staggered grid finite difference method and its optimized solution for the two-dimensional acoustic wave equation.
### File import
My code compilation platform is Visual Studio 2022. First, I need to create an empty project. Then, I should import the `main.cpp` and `time_space_domain.cpp` files into the source files, and import the `time_space_domain.h` into the header files, so as to ensure that all necessary functions have been imported.

```cpp{: .line-numbers}
#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>
#include<time.h>
#include"time_space_domain.h"
using namespace std;
```
### Basic parameters
```cpp{: .line-numbers}
	//basic parameters for modeling
	int M = 6;							//Operator length M for spatial derivative
	int N = 3;							//Operator length N for temporal derivative
	int NX = 201;						//Grid number in x direction
	int NZ = 201;						//Grid number in z direction
	int NT = 1001;						//Moment number for time,starting from time zero.
	double fm = 20.0;					//Dominant frequency of Ricker wavelet(Hz)		
	double esp0 = 1e-6;					//Relative error limitation of Remez exchange algorithm optimization
	double dt = 0.001;					//Time step(s)
	double dh = 10.0;					//Grid spacing(m)
```
&nbsp;&nbsp;&nbsp;&nbsp;The program contains the following basic parameters. The accuracy can be controlled by modifying the operator length M of the spatial partial derivative, the operator length N of the temporal partial derivative, and the error limit esp0 of the Remez exchange algorithm optimization. It should be noted that N needs to be less than or equal to M, and N is preferably no more than 4. As derived in the paper, when N ≥ M, the difference coefficients of the temporal partial derivative cannot be solved. Moreover, an excessively large N will severely reduce the calculation efficiency.Other parameters such as the time step and grid spacing can be modified experimentally, but it is necessary to ensure that they comply with the stability conditions.

### PML
```cpp{: .line-numbers}
	//Parameters required for Perfectly Matched Layer (PML)
	int Nlayer = 30;					//The number of PML layers
	NX = NX + 2 * Nlayer;				//Extend the boundary in x direction
	NZ = NZ + 2 * Nlayer;				//Extend the boundary in z direction
	double R = 0.001;					//Theoretical reflection coefficient
	int n = 2;							//The order of PML
	double vmax = 5000.0;				//The maximum P-wave velocity
	double** d1 = area_2d(NX, NZ);
	double** d2 = area_2d(NX, NZ);
	double** d3 = area_2d(NX, NZ);
	double** d4 = area_2d(NX, NZ);
	for (int i = 0; i < Nlayer; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			d1[i][j] = pow((Nlayer * 1.0 - i*1.0) / Nlayer, n) * log(1 / R) * 1.5 * vmax / (Nlayer * dh);
			d1[NX - 1 - i][j] = pow((Nlayer*1.0 - i*1.0) / Nlayer, n) * log(1 / R) * 1.5 * vmax / (Nlayer * dh);
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < Nlayer; j++)
		{
			d2[i][j] = pow((Nlayer * 1.0 - 1.0 * j ) / Nlayer, n) * log(1.0 / R) * 1.5 * vmax / (Nlayer * dh);
			d2[i][NZ - 1 - j] = pow((Nlayer*1.0 - j*1.0) / Nlayer, n) * log(1 / R) * 1.5 * vmax / (Nlayer * dh);
		}
	}
	for (int i = 0; i < Nlayer; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			d3[i][j] = pow((Nlayer * 1.0 - i * 1.0) / Nlayer, n) * log(1.0 / R) * 1.5 * vmax / (Nlayer * dh);
			d3[NX - i - 1][j] = pow((Nlayer * 1.0 - i * 1.0) / Nlayer, n) * log(1.0 / R) * 1.5 * vmax / (Nlayer * dh);
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < Nlayer; j++)
		{
			d4[i][j] = pow((Nlayer * 1.0 - j * 1.0) / Nlayer, n) * log(1.0 / R) * 1.5 * vmax / (Nlayer * dh);
			d4[i][NZ - j - 1] = pow((Nlayer * 1.0 - j * 1.0) / Nlayer, n) * log(1.0 / R) * 1.5 * vmax / (Nlayer * dh);
		}
	}
```
You can modify the number of PML absorption layers Nlayer, the theoretical reflection coefficient R, and the value of the maximum P-wave velocity vmax. If Nlayer=0, it is equivalent to not adding a PML absorption boundary. The d1, d2, d3, and d4 arrays are specifically used to store the attenuation coefficients in different directions.

