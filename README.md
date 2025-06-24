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
&nbsp;&nbsp;&nbsp;&nbsp;You can modify the number of PML absorption layers Nlayer, the theoretical reflection coefficient R, and the value of the maximum P-wave velocity vmax. If Nlayer=0, it is equivalent to not adding a PML absorption boundary. The d1, d2, d3, and d4 arrays are specifically used to store the attenuation coefficients in different directions.


### FD coefficient table
```cpp{: .line-numbers}
	//Parameters required for solving the difference coefficient table.
	int Length = M + 1 + N * (N - 1) / 2;
	int v_area = 1;
	double** Parameters = area_2d(v_area, Length);
	double* x1 = linspace(N * (N - 1) / 2);
	double* x2 = linspace(M + 1);
	double* v_parameters = linspace(Length);
	double* du = linspace(M);
	double b = 0.0;
	double** duw = area_2d(N, N);
	double* bak = linspace(2);
for (int i = 0; i < v_area; i++)
{
	double r = (1500.0 + i * 1.0) * dt / dh;
	////////////////////////
	int L1 = N;
	int L2 = N * (N - 1) / 2;
	if (N % 2 == 0)
	{
		L1 = N * N / 4;
	}
	if (N % 2 != 0)
	{
		L1 = (N * N - 1) / 4;
	}
	double** supple = area_2d(N, N);
	double* fn = linspace(M + 1);
	double** qmn_equation = area_2d(L2, L2 + 1);
	double** A1 = area_2d(L2, L2);
	double* b1 = linspace(L2);
	double* gn = linspace(M + 1);
	double** Mr = area_2d(M + 1, M + 1);
	//////////////////
	solve_fn(M, r, fn);
	fp = fopen("qmn_equation.txt", "w");
	if (fp != NULL)
	{
		for (int m = 1; m <= int(N / 2); m++)
		{
			for (int n = m; n <= N - m; n++)
			{
				for (int u = 1; u <= N - 1; u++)
				{
					for (int w = 1; w <= N - u; w++)
					{
						fprintf(fp, "%lf\n", q_coe(m, n, u, w));
					}
				}
				fprintf(fp, "%lf\n", mnfmn(m, n, fn));
			}
		}
		for (int m = 1; m <= N - 1; m++)
		{
			for (int n = 1; n <= N - m; n++)
			{
				if (m >= 1 && m <= int(N / 2) && n >= m && n <= N - m)
				{
				}
				else
				{
					supple[m - 1][n - 1] = -1.0;
					supple[n - 1][m - 1] = 1.0;
					for (int u = 1; u <= N - 1; u++)
					{
						for (int w = 1; w <= N - u; w++)
						{
							fprintf(fp, "%lf\n", supple[u - 1][w - 1]);
						}
					}
					fprintf(fp, "%lf\n", 0.0);
				}
			}
		}
		fclose(fp);
	}
	fp = fopen("qmn_equation.txt", "r");
	if (fp != NULL)
	{
		for (int i = 0; i < L2; i++)
		{
			for (int j = 0; j < L2 + 1; j++)
			{
				fscanf_s(fp, "%lf", &qmn_equation[i][j]);
			}
		}
		fclose(fp);
	}
	for (int i = 0; i < L2; i++)
	{
		for (int j = 0; j < L2; j++)
		{
			A1[i][j] = qmn_equation[i][j];
		}
		b1[i] = qmn_equation[i][L2];
	}
	LU(A1, b1, x1, L2);
	int flag2 = 0;
	for (int m = 1; m <= N - 1; m++)
	{
		for (int n = 1; n <= N - m; n++)
		{
			duw[m - 1][n - 1] = x1[flag2];
			flag2 += 1;
		}
	}
	fn_to_gn(fn, gn, duw, M, N);
	for (int n = 1; n <= M + 1; n++)
	{
		for (int u = 1; u <= M; u++)
		{
			Mr[n - 1][u - 1] = alpha(n, u);
		}
	}
	for (int n = 2; n <= M + 1; n++)
	{
		Mr[n - 1][M] = beta(M, gn, n);
	}

	LU(Mr, gn, x2, M + 1);
	/////////////////
	for (int i = 0; i < M; i++)
	{
		du[i] = x2[i];
	}
	b = x2[M];
	int flag = 0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			duw[u - 1][w - 1] = x1[flag];
			flag += 1;
		}
	}
	//////////////
	Remez_Optimize_eta(esp0, M, N, r, du, b, duw, bak);		//Remez exchange algorithm (RA) optimization
	b = bak[1];
	printf("r=%lf betamax=%lf\n", r, bak[0]);
	
	//LS_OptimizeR_eta(esp0, M, N, r, du, b, duw, bak);			//Least square (LS) optimization
	//b = bak[1];
	/////////////////
	for (int i = 0; i < M; i++)
	{
		x2[i] = du[i];
	}
	x2[M] = b;
	flag = 0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			x1[flag] = duw[u - 1][w - 1];
			flag += 1;
		}
	}
	/////////////
	for (int j = 0; j < N * (N - 1) / 2; j++)
	{
		Parameters[i][j] = x1[j];
	}
	for (int j = 0; j < M + 1; j++)
	{
		Parameters[i][j + N * (N - 1) / 2] = x2[j];
	}
	delete_2d(supple, N);
	delete_2d(Mr, M + 1);
	delete[]fn;
	delete[]gn;
	delete_2d(qmn_equation, L2);
	delete_2d(A1, L2);
	delete[]b1;
}
```
&nbsp;&nbsp;&nbsp;&nbsp;Since the time step and grid spacing commonly used in forward modeling are generally fixed, the Courant number is determined solely by the P-wave velocity. By solving the Courant number for different velocities and then specifying the operator length and error limit, the ITSDSGFD difference coefficients optimized by the Remez algorithm can be obtained. The velocity range calculated in this paper is from 1500 m/s to 5000 m/s, so the difference coefficient table is stored in a two-dimensional array.

### Model parameters
```cpp{: .line-numbers}
	double** rho = area_2d(NX, NZ);
	double** vp = area_2d(NX, NZ);
	/*****************Homogeneous model**********************/
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			vp[i][j] = 1500.0;//(m/s)
			rho[i][j] = 1000.0;//(kg/m^3)

		}
	}
```
&nbsp;&nbsp;&nbsp;&nbsp;This demo uses a homogeneous model for testing, with a P-wave velocity of 1500 m/s and a density of 1000 kg/m³.


### Observation system
```cpp{: .line-numbers}
	double** DT = area_2d(NX, NZ);				    //Unit impulse function(defines the spatial position of the source)
	double** wavefront_p1 = area_2d(NX, NZ);		//The wavefield snapshot of stress p captured at time moment k0
	double** wavefront_p2 = area_2d(NX, NZ);		//The wavefield snapshot of stress p captured at time moment k1.
	double** record_p = area_2d(NX, NT);			//Seismic record or seismogram
	double* trace = linspace(NT);					//The waveform record of a certain point.
	int k0 = 400;									//The moment to capture the first wavefield snapshot.
	int k1 = 800;									//The moment to capture the second wavefield snapshot.
	double shotoffset = 20.0;						//The x-coordinate of the source point.
	double shotdepth = 20.0;						//The z-coordinate of the source point.
	DT[Nlayer + int(shotoffset / dh)][Nlayer + int(shotdepth / dh)] = 1.0;         //Source
	double offset0 = 100.0;
	double depth0 = 400.0;
	int NX0 = int(offset0 / dh);
	int NZ0 = int(depth0 / dh);
```
&nbsp;&nbsp;&nbsp;&nbsp;You can modify the parameters of different observation systems to obtain seismic wave data under various observation systems.

### Forward modeling
```cpp{: .line-numbers}
for (int k = 0; k < NT; k++)
{
	if (k % 50 == 0)
	{
		time2 = clock();
		duration = (time2 - time1) / CLOCKS_PER_SEC;
		printf("k=%d,p[NX/2][NZ/2]=%lf t=%lfs\n", k, p[NX / 2][NZ / 2],duration);
	}
	for (int j = 0; j < NZ; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			for (int u = 1; u <= M; u++)
			{
				du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
			}
			b = v_parameters[Length - 1];
			double sumvx = 0.0;
			for (int u = 1; u <= M; u++)
			{
				if ((i + u - 1) < (NX))
				{
					sumvx += du[u - 1] * vx[i + u - 1][j];
				}
				if ((i - u) >= 0)
				{
					sumvx -= du[u - 1] * vx[i - u][j];
				}
			}
			u_before1[i] = sumvx / dh;
			b_alpha1[i][i] = 1.0 - 2.0 * b;
			if (i + 1 < NX)
			{
				b_alpha1[i][i + 1] = b;
			}
			if (i - 1 >= 0)
			{
				b_alpha1[i][i - 1] = b;
			}
		}
		CatchUp(b_alpha1, u_before1, u_after1, NX);
		for (int i = 0; i < NX; i++)
		{
			dvxdx[i][j] = u_after1[i];
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			int flag1 = 0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					duw[u - 1][w - 1] = v_parameters[flag1];
					flag1 += 1;
				}
			}
			double sum = 0.0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					if ((i + u - 1) < NX && j + w < NZ)
					{
						sum += duw[u - 1][w - 1] * vx[i + u - 1][j + w];
					}
					if ((i - u) >= 0 && j + w < NZ)
					{
						sum -= duw[u - 1][w - 1] * vx[i - u][j + w];
					}
					if ((i + u - 1) < NX && j - w >= 0)
					{
						sum += duw[u - 1][w - 1] * vx[i + u - 1][j - w];
					}
					if ((i - u) >= 0 && j - w >= 0)
					{
						sum -= duw[u - 1][w - 1] * vx[i - u][j - w];
					}
				}
			}
			dvxdxt[i][j] = sum / dh;
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			for (int u = 1; u <= M; u++)
			{
				du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
			}
			b = v_parameters[Length - 1];
			double sumvz = 0.0;
			for (int u = 1; u <= M; u++)
			{
				if ((j + u - 1) < (NZ))
				{
					sumvz += du[u - 1] * vz[i][j + u - 1];
				}
				if ((j - u) >= 0)
				{
					sumvz -= du[u - 1] * vz[i][j - u];
				}
			}
			u_before2[j] = sumvz / dh;
			b_alpha2[j][j] = 1.0 - 2.0 * b;
			if (j + 1 < NZ)
			{
				b_alpha2[j][j + 1] = b;
			}
			if (j - 1 >= 0)
			{
				b_alpha2[j][j - 1] = b;
			}
		}
		CatchUp(b_alpha2, u_before2, u_after2, NZ);
		for (int j = 0; j < NZ; j++)
		{
			dvzdz[i][j] = u_after2[j];
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			int flag1 = 0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					duw[u - 1][w - 1] = v_parameters[flag1];
					flag1 += 1;
				}
			}
			double sum = 0.0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					if ((j + u - 1) < NZ && i + w < NX)
					{
						sum += duw[u - 1][w - 1] * vz[i + w][j + u - 1];
					}
					if ((j - u) >= 0 && i + w < NX)
					{
						sum -= duw[u - 1][w - 1] * vz[i + w][j - u];
					}
					if ((j + u - 1) < NZ && i - w >= 0)
					{
						sum += duw[u - 1][w - 1] * vz[i - w][j + u - 1];
					}
					if ((j - u) >= 0 && i - w >= 0)
					{
						sum -= duw[u - 1][w - 1] * vz[i - w][j - u];
					}
				}
			}
			dvzdzt[i][j] = sum / dh;
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			px[i][j] = (2.0 - d3[i][j] * dt) * px[i][j] / (2.0 + d3[i][j] * dt) - (2.0 * rho[i][j] * vp[i][j] * vp[i][j] * dt) * (dvxdx[i][j] + dvxdxt[i][j]) / (2.0 + d3[i][j] * dt);
			pz[i][j] = (2.0 - d4[i][j] * dt) * pz[i][j] / (2.0 + d4[i][j] * dt) - (2.0 * rho[i][j] * vp[i][j] * vp[i][j] * dt) * (dvzdz[i][j] + dvzdzt[i][j]) / (2.0 + d4[i][j] * dt);
			p[i][j] = px[i][j] + pz[i][j] + DT[i][j] * ricker(k * dt, fm) ;

			

		}
	}
	for (int j = 0; j < NZ; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			for (int u = 1; u <= M; u++)
			{
				du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
			}
			b = v_parameters[Length - 1];
			double sumpx = 0.0;
			for (int u = 1; u <= M; u++)
			{
				if ((i + u) < NX)
				{
					sumpx += du[u - 1] * p[i + u][j];
				}
				if ((i - u + 1) >= 0)
				{
					sumpx -= du[u - 1] * p[i - u + 1][j];
				}
			}
			u_before3[i] = sumpx / dh;
			b_alpha3[i][i] = 1.0 - 2.0 * b;
			if (i + 1 < NX)
			{
				b_alpha3[i][i + 1] = b;
			}
			if (i - 1 >= 0)
			{
				b_alpha3[i][i - 1] = b;
			}
		}
		CatchUp(b_alpha3, u_before3, u_after3, NX);
		for (int i = 0; i < NX; i++)
		{
			dpdx[i][j] = u_after3[i];
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			int flag1 = 0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					duw[u - 1][w - 1] = v_parameters[flag1];
					flag1 += 1;
				}
			}
			double sum = 0.0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					if ((i + u) < NX && j + w < NZ)
					{
						sum += duw[u - 1][w - 1] * p[i + u][j + w];
					}
					if ((i - u + 1) >= 0 && j + w < NZ)
					{
						sum -= duw[u - 1][w - 1] * p[i - u + 1][j + w];
					}
					if ((i + u) < NX && j - w >= 0)
					{
						sum += duw[u - 1][w - 1] * p[i + u][j - w];
					}
					if ((i - u + 1) >= 0 && j - w >= 0)
					{
						sum -= duw[u - 1][w - 1] * p[i - u + 1][j - w];
					}
				}
			}
			dpdxt[i][j] = sum / dh;
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			for (int u = 1; u <= M; u++)
			{
				du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
			}
			b = v_parameters[Length - 1];
			double sumpz = 0.0;
			for (int u = 1; u <= M; u++)
			{
				if ((j + u) < NZ)
				{
					sumpz += du[u - 1] * p[i][j + u];
				}
				if ((j - u + 1) >= 0)
				{
					sumpz -= du[u - 1] * p[i][j - u + 1];
				}
			}
			u_before4[j] = sumpz / dh;
			b_alpha4[j][j] = 1.0 - 2.0 * b;
			if (j + 1 < NZ)
			{
				b_alpha4[j][j + 1] = b;
			}
			if (j - 1 >= 0)
			{
				b_alpha4[j][j - 1] = b;
			}
		}
		CatchUp(b_alpha4, u_before4, u_after4, NZ);
		for (int j = 0; j < NZ; j++)
		{
			dpdz[i][j] = u_after4[j];
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			int flag = int(vp[i][j]) - 1500;
			for (int v = 0; v < Length; v++)
			{
				v_parameters[v] = Parameters[flag][v];
			}
			int flag1 = 0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					duw[u - 1][w - 1] = v_parameters[flag1];
					flag1 += 1;
				}
			}
			double sum = 0.0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					if ((j + u) < NZ && i + w < NX)
					{
						sum += duw[u - 1][w - 1] * p[i + w][j + u];
					}
					if ((j - u + 1) >= 0 && i + w < NX)
					{
						sum -= duw[u - 1][w - 1] * p[i + w][j - u + 1];
					}
					if ((j + u) < NZ && i - w >= 0)
					{
						sum += duw[u - 1][w - 1] * p[i - w][j + u];
					}
					if ((j - u + 1) >= 0 && i - w >= 0)
					{
						sum -= duw[u - 1][w - 1] * p[i - w][j - u + 1];
					}
				}
			}
			dpdzt[i][j] = sum / dh;
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			vx[i][j] = (2.0 - d1[i][j] * dt) * vx[i][j] / (2.0 + d1[i][j] * dt) - (2.0 * dt) * (dpdx[i][j] + dpdxt[i][j]) / (rho[i][j] * (2.0 + d1[i][j] * dt));
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			vz[i][j] = (2.0 - d2[i][j] * dt) * vz[i][j] / (2.0 + d2[i][j] * dt) - (2.0 * dt) * (dpdz[i][j] + dpdzt[i][j]) / (rho[i][j] * (2.0 + d2[i][j] * dt));
		}
	}
	if (k == k0)
	{
		for (int i = 0; i < NX; i++)
		{
			for (int j = 0; j < NZ; j++)
			{
				wavefront_p1[i][j] = p[i][j];
			}
		}
	}
	if (k == k1)
	{
		for (int i = 0; i < NX; i++)
		{
			for (int j = 0; j < NZ; j++)
			{
				wavefront_p2[i][j] = p[i][j];
			}
		}
	}
	//snapshot
	for (int i = 0; i < NX; i++)
	{
		record_p[i][k] = p[i][Nlayer];
	}
	//seismic record
	trace[k] = p[NX0 + Nlayer][NZ0 + Nlayer];
	//The waveform record of a certain point.
}
```
&nbsp;&nbsp;&nbsp;&nbsp;The forward modeling process is a time loop. The source is loaded onto the stress. Since the calculation process of the implicit space-time domain staggered grid finite difference method is more complex than that of the conventional explicit staggered grid finite difference method, the spatial partial derivative values are first calculated, and then the time recursion is performed.


### Output
```cpp{: .line-numbers}
	fp = fopen("wavefront1_tsd.dat", "wb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int j = Nlayer; j < NZ - Nlayer; j++)
			{
				float wpt = 0.0;
				wpt = float(wavefront_p1[i][j]);
				fwrite(&wpt, sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}
	fp = fopen("wavefront2_tsd.dat", "wb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int j = Nlayer; j < NZ - Nlayer; j++)
			{
				float wpt = 0.0;
				wpt = float(wavefront_p2[i][j]);
				fwrite(&wpt, sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}
	fp = fopen("recordra.dat", "wb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int k = 0; k < NT; k++)
			{
				float rd = 0.0;
				rd = float(record_p[i][k]);
				fwrite(&rd, sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}
	fp = fopen("recordra.txt", "wb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int k = 0; k < NT; k++)
			{
				fprintf(fp, "%lf\n", record_p[i][k]);
			}
		}
		fclose(fp);
	}
	fp = fopen("tracera.txt", "w");
	if (fp != NULL)
	{
		for (int k = 0; k < NT; k++)
		{
			fprintf(fp, "%1.12f\n", trace[k]);
		}
		fclose(fp);
	}
```
&nbsp;&nbsp;&nbsp;&nbsp;Before simulation, this demo extends an additional PML absorption boundary beyond the model boundary, so the output results need to be cropped to retain only the internal normal wavefield. The file type output by this demo is a binary .dat file, and the data type is float. You need to use other programs that can read and write binary .dat files, or use Python, MATLAB, etc. to read and plot the data.
