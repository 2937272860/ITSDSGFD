#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<time.h>
#include"time_space_domain.h"
using namespace std;

/*********************************************/
/*
This program demonstrates the solution process of FD coefficients for the implicit time-space-domain staggered-grid finite difference (ITSDSGFD) scheme.
The inputs of the program are the operator length N for the temporal partial derivative, the operator length M for the spatial partial derivative, and the Courant number r. 
During the calculation, the program will display the FD coefficients obtained by the Taylor series expansion(TE) method on the screen. 
If the Least squares(LS) method or the Remez exchange algorithm (RA) is called for optimization, the FD coefficients before and after optimization will be displayed in sequence. 
After obtaining the FD coefficients, the program will calculate the phase velocity dispersion and relative error corresponding to different FD coefficients and output them as a txt file.
Finally, the program will display the stability factor corresponding to the FD coefficients obtained in this solution, allowing us to determine whether the given operator lengths and Courant number meet the stability conditions.
*/
/***********************************************/

int main()
{
	clock_t time1, time2;
	time1 = clock();
	FILE* fp;
	/****************************************************************************************/
	//necessary parameters
	int M = 6;						//the operator length M for the spatial partial derivative
	int N = 3;						//the operator length N for the temporal partial derivative
	double v = 1500.0;				//P-wave velocity(m/s)
	double dt = 0.001;				//time step(s)
	double dh = 10.0;				//grid spacing(m)
	/***********************************************************************************/
	double r = v * dt / dh;			//Courant number
	/********************************************************/
	double* du = linspace(M);		//spatial derivative FD coefficients
	double b = 0.0;					//spatial implicit FD coefficients
	double** duw = area_2d(N, N);	//temporal derivative FD coefficients

	/*******************************************************************/
	//The process of solving the temporal partial derivative FD coefficients for  using the Taylor series expansion (TE) method
	int L1 = N;
	int L2 = N * (N - 1) / 2;
	if (N % 2 == 0)
	{
		L1 = N * N / 4;				//the number of independent temporal derivative FD coefficients
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
	double* x1 = linspace(N * (N - 1) / 2);
	double* x2 = linspace(M + 1);
	
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
	//The temporal partial derivatives FD coefficients are stored in the two-dimensional array double** duw.
	/**********************************************************************************************/



	
	/**********************************************************************************************/
	//The process of solving the spatial partial derivative implicit FD coefficients for  using the Taylor series expansion (TE) method
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
	//The spatial partial derivatives implicit FD coefficients are stored in double b and the one-dimensional array double* du
	/**********************************************************************************************/
	printf("TE-based FD coefficients\n");
	for (int i = 0; i < M; i++)
	{
		printf("du[%d]=%lf\n", i, du[i]);
	}
	printf("b=%lf\n", b);
	int flag = 0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			duw[u - 1][w - 1] = x1[flag];
			flag += 1;
			printf("duw[%d][%d]=%lf\n", u - 1, w - 1, duw[u - 1][w - 1]);
		}
	}
	//Display all the FD coefficients obtained by the TE method on the screen.
	/*********************************************************************************************************************************/



	/*********************************************************************************************************************************/
	//necessary parameters for optimization
	double betamax = 2.4;				//the maximum wavenumber betamax
	double E = 1e-6;					//error values E
	double eta = 1e-6;					//relative error limitation
	double* bak = linspace(2);			//store two parameters in the calculation process
	/**********************/

	
	/*********************************************************************************************************************************/
	//There are four choices for ptimization.

	/**********************/
	// Optimize the spatial implicit FD coefficients using the Remez exchange algorithm (RA) under the premise of a given error limitation eta. 
	// While outputting the difference coefficients, output the maximum wavenumber range that meets the error limitation eta.
	/**********************/
	Remez_Optimize_eta(eta, M, N, r, du, b, duw, bak);
	b = bak[1];								//bak[1] stores implicit FD coefficients b
	betamax = bak[0];						//bak[0] stores the maximum wavenumber
	printf("r=%lf eta=%lf betamax=%lf\n", r, eta, bak[0]);
	/**********************/

	/**********************/
	// Optimize the spatial implicit FD coefficients with the Remez exchange algorithm (RA) under the premise of a given maximum wavenumber range.
	// When outputting the difference coefficients, also output the error values E corresponding to the wavenumber range, FD operator length M,N, Courant number r and temporal FD coefficients duw.
	/**********************/
	//betamax = 2.46;
	//Remez_Optimize_beta(betamax, M, N, r, du, b, duw, bak);
	//b = bak[1];							//bak[1] stores implicit FD coefficients b
	//E = bak[0];							//bak[0] stores error values E
	//printf("r=%lf E=%1.12f betamax=%lf\n", r, E, betamax);
	/**********************/

	/**********************/
	// Optimize the spatial implicit FD coefficients using the Least square method(LS) under the premise of a given error limitation eta. 
	// While outputting the difference coefficients, output the maximum wavenumber range that meets the error limitation eta.
	/**********************/
	//LS_OptimizeR_eta(eta, M, N, r, du, b, duw, bak);
	//b = bak[1];							//bak[1] stores implicit FD coefficients b
	//betamax = bak[0];						//bak[0] stores the maximum wavenumber
	//printf("betamax=%1.12f\n", bak[0]);
	/**********************/



	/**********************/
	// Optimize the spatial implicit FD coefficients with the Least square method(LS) under the premise of a given maximum wavenumber range.
	// When outputting the difference coefficients, also output the error values corresponding to the wavenumber range, FD operator length M,N, Courant number r and temporal FD coefficients duw.
	/**********************/
	//betamax = 2.7380;
	//LS_OptimizeR_beta(betamax, M, N, r, du, b, duw, bak);
	//b = bak[1];							//bak[1] stores implicit FD coefficients b
	//E = bak[0];							//bak[0] stores error values E
	//printf("betamax=%lf E=%1.12f\n", betamax, E);
	/**********************/

	/**********************/
	//optimized FD coefficients
	printf("optimized FD coefficients\n");
	for (int i = 0; i < M; i++)
	{
		printf("du[%d]=%1.12f\n", i, du[i]);
	}
	printf("b=%1.12f\n", b);
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			printf("duw[%d][%d]=%1.12f\n", u - 1, w - 1, duw[u - 1][w - 1]);
		}
	}
	/**********************/

	/*********************************************************************************************************************************/
	//we choose the given error limitation,and we use RA to optimize the FD coefficients 
	//parameters can be checked in parameters.txt
	fp = fopen("parameters.txt", "wb");
	if (fp != NULL)
	{
		fprintf(fp, "Paramters\n");
		fprintf(fp, "M=%d\n", M);
		fprintf(fp, "N=%d\n", N);
		fprintf(fp, "r=%lf\n", r);
		fprintf(fp, "eta=%1.9f\n", eta);
		for (int i = 0; i < M; i++)
		{
			x2[i] = du[i];
			fprintf(fp,"du[%d]=%1.9f\n", i, du[i]);
		}
		x2[M] = b;
		fprintf(fp,"b=%1.9f\n", b);
		flag = 0;
		for (int u = 1; u <= N - 1; u++)
		{
			for (int w = 1; w <= N - u; w++)
			{
				x1[flag] = duw[u - 1][w - 1];
				flag += 1;
				fprintf(fp,"duw[%d][%d]=%1.9f\n", u - 1, w - 1, duw[u - 1][w - 1]);
			}
		}
		fclose(fp);
	}
	//parameters can be checked in parameters.txt

	/**********************/
	//phase velocity dispersion, defined as function named dispersin4 in time_space_domain.cpp
	fp = fopen("disp.txt", "w");
	if (fp != NULL)
	{
		//theta
		for (int j = 0; j < 1000; j++)
		{
			double kh = j * PI / 1000 + 1e-4;
			fprintf(fp, "%1.12f ", kh);
			for (int i = 0; i < 5; i++)//kh=beta
			{
				double theta = i * PI / 16;
				double delta = dispersion4(kh, theta, r, M, N, x1, x2);
				fprintf(fp, "%1.12f ", delta);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	////////
	fp = fopen("disper.txt", "w");
	if (fp != NULL)
	{
		//theta
		for (int j = 0; j < 1000; j++)
		{
			double kh = j * PI / 1000 + 1e-4;
			fprintf(fp, "%1.12f ", kh);
			double delta = dispersion4(kh, 1e-12, r, M, N, x1, x2);
			fprintf(fp, "%1.12f ", delta);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	///////////////
	fp = fopen("dispersion.txt", "w");
	if (fp != NULL)
	{
		for (int i = 0; i < 50; i++)
		{
			double beta = (i + 1) * PI / 50+1e-6;
			for (int j = 0; j < 50; j++)
			{
				double theta = j * PI / 100 + 1e-6;
				double disp = dispersion4(beta, theta, r, M, N, x1, x2);
				fprintf(fp, "%lf %lf %lf\n", beta, theta, disp);
			}
		}
		fclose(fp);
	}
	/**********************/

	/**********************/
	//relative error function, defined as function named ErrorR in time_space_domain.cpp
	fp = fopen("ErrorR.txt", "w");
	if (fp != NULL)
	{
		for (int i = 0; i < 1000; i++)
		{
			double beta = i * PI / 1000 + 1e-4;
			double delta = ErrorR(beta, M, N, r, du, b, duw);
			//fprintf(fp, "%lf %lf\n",beta, delta);
			fprintf(fp, "%1.12f %1.12f\n", beta, delta);
		}
		fclose(fp);
	}
	fp = fopen("ErrorRemez.txt", "w");
	if (fp != NULL)
	{
		for (int i = 0; i < 1000; i++)
		{
			double beta = i * PI / 1000 + 1e-4;
			double delta = ErrorR(beta, M, N, r, du, b, duw);
			//fprintf(fp, "%lf %lf\n",beta, delta);
			fprintf(fp, "%1.12f %1.12f %1.12f %1.12f\n", beta, delta, E, -E);
		}
		fclose(fp);
	}
	/**********************/

	/**********************/
	//stability factor, calculated directly in main function
	double s = 0.0;						//stability factor
	double s1 = 0.0;
	for (int u = 1; u <= M; u++)
	{
		s1 += du[u - 1] * pow(-1.0, u + 1) / (1.0 - 4.0 * b);
	}
	double s2 = 0.0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			s2 += 2.0 * duw[u - 1][w - 1] * pow(-1.0, u + w + 1);
		}
	}
	s = 1.0 / (sqrt(2.0) * (s1 + s2));//stability factor
	printf("M=%d,N=%d,r=%lf,stability factor s=%lf\n", M, N, r, s);
	/**********************/

	/*********************************************************************************************************************************/
	//Free up the memory space when the program finishes running.
	delete_2d(supple, N);
	delete_2d(Mr, M + 1);
	delete[]fn;
	delete[]gn;
	delete_2d(qmn_equation, L2);
	delete_2d(A1, L2);
	delete[]b1;
	delete[]x2;
	delete[]x1;
	delete[]du;
	delete_2d(duw, N);
	delete[]bak;



	time2 = clock();
	double duration = (time2 - time1) / CLOCKS_PER_SEC;
	printf(" Program runs in %lf seconds\n", duration);
	return 0;
}