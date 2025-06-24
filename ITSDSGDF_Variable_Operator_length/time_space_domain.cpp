#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<iostream>
using namespace std;
#define PI 3.1415926

/******************************************/
//basic function
double** area_2d(int NX, int NZ)//Generate a two-dimensional array with NX columns and NZ rows.
{
	double** b = new double* [NX];
	for (int i = 0; i < NX; i++)
		b[i] = new double[NZ];
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NZ; j++)
			b[i][j] = 0;//Initialize a two-dimensional double array with zeros.
	return b;
}

double* linspace(int L)//Define a one-dimensional double array with length L.
{
	double* l = new double[L];
	for (int i = 0; i < L; i++)
		l[i] = 0;//Initialize with zero.
	return l;
}

int* linspace_int(int L)//Define a one-dimensional int array with length L.
{
	int* l = new int[L];
	for (int i = 0; i < L; i++)
		l[i] = 0;//Initialize with zero.
	return l;
}

void delete_2d(double** b, int NX)//Release a two-dimensional array with NX columns.
{
	for (int i = 0; i < NX; i++)
		delete[]b[i];
	delete[]b;
	return;
}

double factorial(int N)//Factorial function
{
	double  p = 1.0;
	if (N == 0)
	{
		p = 1.0;
	}
	else
	{
		for (int i = 1; i <= N; i++)
		{
			p *= 1.0 * i;
		}
	}
	return p;
}

void LU(double** a, double* b, double* x, int  N)
//Solve the N-dimensional linear system of equations with coefficient matrix A and nonhomogeneous term b using LU decomposition,
// and return the solution x in the form of a one-dimensional array.
{
	double** l = area_2d(N, N);
	double** u = area_2d(N, N);
	int i, r, k;
	for (i = 0; i < N; i++)
	{
		u[0][i] = a[0][i];
	}
	for (i = 1; i < N; i++)
	{
		l[i][0] = a[i][0] / u[0][0];
	}
	for (r = 1; r < N; r++)
	{
		for (i = r; i < N; i++)
		{
			double sum1 = 0;
			for (k = 0; k < r; k++)
			{
				sum1 += l[r][k] * u[k][i];
			}
			u[r][i] = a[r][i] - sum1;
		}

		if (r != N)
			for (i = r + 1; i < N; i++)
			{
				double sum2 = 0;
				for (k = 0; k < r; k++)
				{
					sum2 += l[i][k] * u[k][r];
				}
				l[i][r] = (a[i][r] - sum2) / u[r][r];
			}

	}
	double* y = linspace(N);
	y[0] = b[0];
	for (i = 1; i < N; i++)
	{
		double sum3 = 0;
		for (k = 0; k < i; k++)
			sum3 += l[i][k] * y[k];
		y[i] = b[i] - sum3;
	}

	x[N - 1] = y[N - 1] / u[N - 1][N - 1];
	for (i = N - 2; i >= 0; i--)
	{
		double sum4 = 0;
		for (k = i + 1; k < N; k++)
			sum4 += u[i][k] * x[k];
		x[i] = (y[i] - sum4) / u[i][i];
	}
	delete_2d(l, N);
	delete_2d(u, N);
	delete[]y;
	return;
}
/******************************************/
//ITSDSGFD coefficients calculated by TE method
/**********/
void solve_fn(int M, double r, double* fn)
//Calculate the values of fn for ( n = 1, 2, ..., M, M+1 ) and store them in the array fn.
{
	fn[0] = 0.5;
	if (M == 1)
	{
		fn[1] = -r * r / (2 * factorial(4));
	}
	if (M >= 2)
	{
		fn[1] = -r * r / (2 * factorial(4));
		for (int n = 3; n <= M + 1; n++)
		{
			double caucy = 0.0;
			for (int l = 2; l <= n - 1; l++)
			{
				caucy += fn[l - 1] * fn[n - l];
			}
			double coe = pow(-1.0, n - 1) * pow(r, 2 * n - 2) / (2 * factorial(2 * n));
			fn[n - 1] = coe - caucy;
		}
	}
	return;
}
double q_coe(int m, int n, int u, int w)
{
	double c1 = pow(-1.0, m + n - 1) * pow(1.0 * u - 0.5, 2 * m - 1) * pow(1.0 * w, 2 * n);
	double c2 = 1.0 * factorial(2 * m - 1) * factorial(2 * n);
	double c3 = pow(-1.0, m + n - 1) * pow(1.0 * u - 0.5, 2 * n - 1) * pow(1.0 * w, 2 * m);
	double c4 = 1.0 * factorial(2 * n - 1) * factorial(2 * m);
	double c5 = 2.0 * c1 / c2 + 2.0 * c3 / c4;
	return c5;
}
double mnfmn(int m, int n, double* fn)
{
	double c1 = 1.0 * factorial(m + n) / (1.0 * factorial(m) * factorial(n));
	double c2 = fn[m + n - 1];
	double c3 = c1 * c2;
	return c3;
}
void fn_to_gn(double* fn, double* gn, double** duw, int M, int N)
{

	for (int n = 1; n <= M + 1; n++)
	{
		double sum = 0.0;
		for (int u = 1; u <= N - 1; u++)
		{
			for (int w = 1; w <= N - u; w++)
			{
				sum += 2.0 * duw[u - 1][w - 1] * pow(-1.0, n - 1) * pow(1.0 * u - 0.5, 2 * n - 1) / factorial(2 * n - 1);
			}
		}
		gn[n - 1] = fn[n - 1] - sum;
	}
	return;
}
double alpha(int n, int u)
{
	double c1 = 1.0;
	double c2 = pow(-1.0, n - 1);
	double c3 = pow(1.0 * u - 0.5, 2 * n - 1);
	double c4 = 1.0 / (factorial(2 * n - 1) * 1.0);
	c1 = c2 * c3 * c4;
	return c1;
}
double beta(int M, double* gn, int n)
{
	double c1 = 0.0;
	for (int p = 1; p <= n - 1; p++)
	{
		c1 += pow(-1.0, n - p) * gn[p - 1] / (factorial(2 * n - 2 * p) * 1.0);
	}
	c1 *= -2.0;
	return c1;
}
/**********/

/*************************************************/
//optimization function
double H(double beta, double r, int N, double** duw)
{
	double s1 = sin(0.5 * r * beta) / r;
	double s2 = 0.0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			s2 += 2.0 * duw[u - 1][w - 1] * sin((1.0 * u - 0.5) * beta);
		}
	}
	double s3 = s1 - s2;
	return s3;
}

double ErrorR(double beta, int M, int N, double r, double* du, double b, double** duw)//Relative error
{
	double s1 = 0.0;
	for (int i = 1; i <= M; i++)
	{
		s1 += du[i - 1] * sin((1.0 * i - 0.5) * beta);
	}
	double s2 = b * (2.0 - 2.0 * cos(beta));
	double s3 = H(beta, r, N, duw);
	double s4 = s1/s3 + s2 -1.0;
	return s4;
}

double L2R(double** RombergT, int RomberOrder, double betamax, double r, int M, int N, double* du, double b, double** duw)//The L2 norm of relative error.
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = ErrorR(xa, M, N, r, du, b, duw);
	double fb = ErrorR(xb, M, N, r, du, b, duw);
	RombergT[0][0] = h * (fa*fa + fb*fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += ErrorR(xa + (2 * i - 1) * h, M, N, r, du, b, duw)* ErrorR(xa + (2 * i - 1) * h, M, N, r, du, b, duw);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;

		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}

/*************************************************/
/*
When optimizing the difference coefficients using the least squares method,it is necessary to solve a system of linear equations. 
At this time, both the coefficient matrix and the non - homogeneous term of the system of equations are integral expressions,
and numerical integration methods need to be used to solve them.
*/
//All numerical integrations in this program group are solved using the Romberg integration formula.

double f1R(double beta, int u, int j, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = sin((1.0 * u - 0.5) * beta) * sin((1.0 * j - 0.5) * beta);
	double s3 = H(beta, r, N, duw);
	return s1 / (s3 * s3);
}

double a1R(double** RombergT, int RomberOrder, double betamax, int u, int j, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f1R(xa, u, j, r, N, duw);
	double fb = f1R(xb, u, j, r, N, duw);;
	RombergT[0][0] = h * (fa + fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f1R(xa + (2 * i - 1) * h, u, j, r, N, duw);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;

		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}

double f2R(double beta, int j, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s2 = 0.0;
	s2 = H(beta, r, N, duw);
	double s3 = 0.0;
	s3 = sin((1.0 * j - 0.5) * beta);
	double s4 = 0.0;
	s4 = s1 * s3 / s2;
	return s4;
}

double a2R(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f2R(xa, j, r, N, duw);
	double fb = f2R(xb, j, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f2R(xa + (2 * i - 1) * h, j, r, N, duw);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;

		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}

double f3R(double beta, int u, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s2 = 0.0;
	s2 = H(beta, r, N, duw);
	double s3 = 0.0;
	s3 = sin((1.0 * u - 0.5) * beta);
	double s4 = 0.0;
	s4 = s1 * s3 / s2;
	return s4;
}

double a3R(double** RombergT, int RomberOrder, double betamax, int u, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f3R(xa, u, r, N, duw);
	double fb = f3R(xb, u, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f3R(xa + (2 * i - 1) * h, u, r, N, duw);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;
		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}

double f4R(double beta, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s4 = 0.0;
	s4 = s1  * s1 ;
	return s4;
}

double a4R(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f4R(xa, r, N, duw);
	double fb = f4R(xb, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f4R(xa + (2 * i - 1) * h, r, N, duw);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;
		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}

double g1R(double beta, int j, double r, int N, double** duw)
{
	double s2 = 0.0;
	s2 = H(beta, r, N, duw);
	double s3 = 0.0;
	s3 = sin((1.0 * j - 0.5) * beta);
	double s4 = 0.0;
	s4 = s3 / s2;
	return s4;
}

double b1R(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = g1R(xa, j, r, N, duw);
	double fb = g1R(xb, j, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += g1R(xa + (2 * i - 1) * h, j, r, N, duw);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;
		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}

double g2R(double beta)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	return s1;
}

double b2R(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = g2R(xa);
	double fb = g2R(xb);
	RombergT[0][0] = h * (fa + fb);/////Replace the integrand
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += g2R(xa + (2 * i - 1) * h);///////Replace the integrand
		}
		////
		RombergT[k][0] = RombergT[k - 1][0] / 2 + h * F;
		///
		for (int m = 1; m <= k; m++)
		{
			RombergT[k - m][m] = (pow(4, m) * RombergT[k - m + 1][m - 1] - RombergT[k - m][m - 1]) / (pow(4, m) - 1);
		}
		///
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) < esp)
		{
			Rout = 1;
			I = RombergT[0][k];
		}
		if (fabs(RombergT[0][k] - RombergT[0][k - 1]) >= esp)
		{
			h *= 0.5;
			n *= 2;
			k += 1;
			Rout = 0;
		}
		//if (k == RomberOrder)
		//{
		//	Rout = 1;
		//}
	} while (Rout == 0);
	return I;
}
/*************************************************/


// Optimize the spatial implicit FD coefficients with the Least square method(LS) under the premise of a given maximum wavenumber range.
// When outputting the difference coefficients, also output the error values corresponding to the wavenumber range, FD operator length M,N, Courant number r and temporal FD coefficients duw.
void LS_OptimizeR_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB)
{
	double** Mr_new = area_2d(M + 1, M + 1);
	double* gn_new = linspace(M + 1);
	double* x2_new = linspace(M + 1);
	int RomberOrder = 15;
	double** RombergT = area_2d(RomberOrder, RomberOrder);

	for (int j = 1; j <= M; j++)
	{
		for (int u = 1; u <= M; u++)
		{
			for (int p = 0; p < RomberOrder; p++)
			{
				for (int q = 0; q < RomberOrder; q++)
				{
					RombergT[p][q] = 0.0;
				}
			}
			Mr_new[j - 1][u - 1] = a1R(RombergT, RomberOrder, betamax, u, j, r, N, duw);
		}
	}
	for (int j = 1; j <= M; j++)
	{
		for (int p = 0; p < RomberOrder; p++)
		{
			for (int q = 0; q < RomberOrder; q++)
			{
				RombergT[p][q] = 0.0;
			}
		}
		Mr_new[j - 1][M] = a2R(RombergT, RomberOrder, betamax, j, r, N, duw);
	}
	for (int u = 1; u <= M; u++)
	{
		for (int p = 0; p < RomberOrder; p++)
		{
			for (int q = 0; q < RomberOrder; q++)
			{
				RombergT[p][q] = 0.0;
			}
		}
		Mr_new[M][u - 1] = a3R(RombergT, RomberOrder, betamax, u, r, N, duw);
	}
	Mr_new[M][M] = a4R(RombergT, RomberOrder, betamax, r, N, duw);
	for (int j = 1; j <= M; j++)
	{
		for (int p = 0; p < RomberOrder; p++)
		{
			for (int q = 0; q < RomberOrder; q++)
			{
				RombergT[p][q] = 0.0;
			}
		}
		gn_new[j - 1] = b1R(RombergT, RomberOrder, betamax, j, r, N, duw);
	}
	for (int p = 0; p < RomberOrder; p++)
	{
		for (int q = 0; q < RomberOrder; q++)
		{
			RombergT[p][q] = 0.0;
		}
	}
	gn_new[M] = b2R(RombergT, RomberOrder, betamax, r, N, duw);
	LU(Mr_new, gn_new, x2_new, M + 1);
	for (int i = 0; i < M; i++)
	{
		du[i] = x2_new[i];
	}
	b = x2_new[M];
	//Solve the system of linear equations to obtain the optimized implicit difference coefficients.
	double Errorlimit = 0.0;
	for (double beta = betamax; beta >=1e-6; beta -= 1e-6)
	{
		double errori = fabs(ErrorR(beta, M, N, r, du, b, duw));
		if (errori >= Errorlimit)
		{
			Errorlimit = errori;
		}
	}
	//Find the maximum error Errorlimit.
	//printf("betamax=%lf Errorlimit=%1.12f\n", betamax, Errorlimit);
	delete_2d(Mr_new, M + 1);
	delete[]gn_new;
	delete[]x2_new;
	delete_2d(RombergT, RomberOrder);
	//printf("betamax=%lf\n", betamax);
	BetaMaxNewAndB[0] = Errorlimit;
	BetaMaxNewAndB[1] = b;
	return;
}


// Optimize the spatial implicit FD coefficients using the Least square method(LS) under the premise of a given error limitation eta. 
// While outputting the difference coefficients, output the maximum wavenumber range that meets the error limitation eta.
//Use the bisection method to find the maximum wavenumber that satisfies the error limitation.
void LS_OptimizeR_eta(double eta, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB)
{
	double betaa = 1e-12;
	double betab = PI - 1e-12;
	double Errorlimit = 1e-12;
	int lsoutmark = 0;
	double betamax = (betaa + betab) / 2;
	do {
		betamax = (betaa + betab) / 2;
		LS_OptimizeR_beta(betamax, M, N, r, du, b, duw, BetaMaxNewAndB);
		Errorlimit = BetaMaxNewAndB[0];
		b = BetaMaxNewAndB[1];
		if (Errorlimit <= eta)
		{
			betaa = betamax;
		}
		if (Errorlimit > eta)
		{
			betab = betamax;
		}
		if ((betab - betaa) <1e-9)
		{
			lsoutmark = 1;
		}
		BetaMaxNewAndB[0]=betamax;
		printf("betamx=%lf Errorlimit=%lf\n", betamax, Errorlimit);
	} while (lsoutmark == 0);
	return;
}

double dispersion4(double kh, double theta, double r, int M, int N, double* x1, double* x2)//phase velocity dispersion
//x1 stores temporal FD coefficients, and x2 stores spatial FD coefficients
{
	double d = 0.0;
	double** duw = area_2d(N, N);

	int flag1 = 0;
	for (int m = 1; m <= N - 1; m++)
	{
		for (int n = 1; n <= N - m; n++)
		{
			duw[m - 1][n - 1] = x1[flag1];
			flag1 += 1;
			//printf("duw[%d][%d]=%lf ", m - 1, n - 1, duw[m - 1][n - 1]);
		}
		//printf("\n");
	}
	double ckx = 0.0;
	double ckz = 0.0;
	for (int u = 1; u <= M; u++)
	{
		ckx += x2[u - 1] * sin((1.0 * u - 0.5) * kh * cos(theta)) / (1.0 - 2.0 * x2[M] + 2.0 * x2[M] * cos(kh * cos(theta)));
		ckz += x2[u - 1] * sin((1.0 * u - 0.5) * kh * sin(theta)) / (1.0 - 2.0 * x2[M] + 2.0 * x2[M] * cos(kh * sin(theta)));
	}
	double dkx = 0.0;
	double dkz = 0.0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			dkx += 2.0 * duw[u - 1][w - 1] * sin((1.0 * u - 0.5) * kh * cos(theta)) * cos(w * kh * sin(theta));
			dkz += 2.0 * duw[u - 1][w - 1] * sin((1.0 * u - 0.5) * kh * sin(theta)) * cos(w * kh * cos(theta));
		}
	}
	double c1 = 1.0 - 2.0 * r * r * (ckx + dkx) * (ckx + dkx) - 2.0 * r * r * (ckz + dkz) * (ckz + dkz);
	d = acos(c1) / (r * kh);
	delete_2d(duw, N);
	return d;
}


// Optimize the spatial implicit FD coefficients with the Remez exchange algorithm (RA) under the premise of a given maximum wavenumber range.
// When outputting the difference coefficients, also output the error values E corresponding to the wavenumber range, FD operator length M,N, Courant number r and temporal FD coefficients duw.
void Remez_Optimize_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* bak)
{
	double esp = 1e-9;
	double E0 = 1e-9;
	double E = 1e-9;
	int Remezoutmark = 0;
	int RemezT = 0;
	double* betai = linspace(M + 2);
	for (int i = 1; i <= M + 2; i++)
	{
		betai[i - 1] = i * betamax / (M + 2);
		//printf("betai[%d]=%lf ", i - 1, betai[i - 1]);
	}
	//printf("\n");
	double** Mr_new = area_2d(M + 2, M + 2);
	double* gn_new = linspace(M + 2);
	double* Remez_du = linspace(M + 2);
	double* Remez_root = linspace(M + 2);
	///////////
	do {
		for (int i = 0; i < M + 2; i++)
		{
			for (int j = 0; j < M; j++)
			{
				Mr_new[i][j] = sin((1.0 * (j + 1) - 0.5) * betai[i]) / H(betai[i], r, N, duw);
			}
			Mr_new[i][M] = (2.0 - 2.0 * cos(betai[i]));
			Mr_new[i][M + 1] = pow(-1.0, i + 1);
			gn_new[i] = 1.0;
		}
		LU(Mr_new, gn_new, Remez_du, M + 2);
		E = Remez_du[M + 1];
		b = Remez_du[M];
		for (int i = 0; i < M; i++)
		{
			du[i] = Remez_du[i];
		}
		for (int i = 1; i <= M + 1; i++)
		{
			double x0 = betai[i - 1];//Lower limit of interval
			double x1 = betai[i];//Upper limit of interval
			double x2 = (x1+x0)/2;//Zero point
			int roott = 0;
			do {
				roott += 1;
				x2 = (x1 + x0) / 2;
				double f0 = ErrorR(x0, M, N, r, du, b, duw);
				double f1 = ErrorR(x1, M, N, r, du, b, duw);
				double f2 = ErrorR(x2, M, N, r, du, b, duw);
				if (f0 * f2 <= 0)
				{
					x1 = x2;
				}
				if (f1 * f2 <= 0)
				{
					x0 = x2;
				}
				//printf("i=%d roott=%d x2=%lf\n", i, roott, x2);
			} while (fabs(x1 - x0) >=1e-9&&roott<40);
			Remez_root[i] = x2;
		}
		//for (int i = 0; i < M + 2; i++)
		//{
		//	printf("Remez_root[%d]=%lf ", i, Remez_root[i]);
		//}
		//printf("\n");
		for (int i = 0; i <= M; i++)
		{
			double left = Remez_root[i];//Lower limit of interval
			double right = Remez_root[i + 1];//pper limit of interval
			double midl, midr;
			//Use the trisection method to find the maximum value of a function and the point where the maximum occurs.
			while (right - left >= 1e-9)
			{
				midl = (2 * left + right) / 3;
				midr = (left + 2 * right) / 3;
				if (fabs(ErrorR(midl, M, N, r, du, b, duw)) < fabs(ErrorR(midr, M, N, r, du, b, duw)))////Find the maximum value of a function.
				{
					left = midl;
				}
				else
				{
					right = midr;
				}
			}
			betai[i] = (left + right) / 2;
		}
		E = fabs(E);
		//printf("RemezT=%d E=%1.12f\n",RemezT, E);
		if (fabs(E - E0) / fabs(E0) < 1E-10)//Error convergence
		{
			Remezoutmark = 1;
		}
		E0 = E;
		RemezT += 1;
		if (RemezT > 20)//Automatically exit the loop after too many iterations.
		{
			Remezoutmark = 1;
		}
		//for (int i = 0; i < M + 2; i++)
		//{
		//	printf("%betai[%d]=%lf ", i, betai[i]);
		//}
		//printf("\n");
	} while (Remezoutmark == 0);
	//printf("betamax=%lf E=%1.12f\n", betamax, E);
	///////////
	delete[]Remez_du;
	delete[]gn_new;
	delete_2d(Mr_new, M + 2);
	delete[]betai;
	delete[]Remez_root;
	bak[0] = E;
	bak[1] = b;
	return;
}

// Optimize the spatial implicit FD coefficients using the Least square method(LS) under the premise of a given error limitation eta. 
// While outputting the difference coefficients, output the maximum wavenumber range that meets the error limitation eta.
//Use the bisection method to find the maximum wavenumber that satisfies the error limitation.
void Remez_Optimize_eta(double eta, int M, int N, double r, double* du, double b, double** duw,double*bak)
{
	double betaa = 1e-9;
	double betab = PI - 1e-9;
	double E = 1e-9;
	int etat = 0;
	int outmark = 0;
	do {
		etat += 1;
		double betamax = (betaa + betab) / 2;
		Remez_Optimize_beta(betamax, M, N, r, du, b, duw, bak);
		b = bak[1];
		E = fabs(bak[0]);
		if (E > eta)
		{
			betab = betamax;
		}
		if (E <= eta)
		{
			betaa = betamax;
		}
		bak[0] = betamax;
		//printf("etat=%d betamax=%lf E=%1.12f b-a=%lf\n",etat, betamax, E, betab - betaa);
	} while ((betab-betaa)>=1e-9&&etat<100);
	//printf("r=%lf betamax=%lf E=%1.12f\n", r, betaa, E);

	return;
}

/************************************************/
//necessary function for modeling
double ricker(double idt, double fm)
//The source function is a Ricker wavelet, 
// which takes the time instant idt and the dominant frequency fm as inputs and outputs the wavefield value at the corresponding instant.
{
	double r = 0.0;
	r = (1 - 2 * pow(PI * fm * (idt - 1.0 / fm), 2)) * exp(-pow(PI * fm * (idt - 1.0 / fm), 2));
	return r;
}


void CatchUp(double** A, double* u_before, double* u_after, int N)
//The CatchUp function is used to solve a system of linear equations with a tridiagonal matrix as the coefficient matrix. 
// Here, double** A is a two-dimensional array storing the coefficient matrix,
// double* u_before is the non-homogeneous term, 
// double* u_after is the solution to the equations,
// and int N is the order of the equations.
{
	double* alpha = linspace(N);
	double* beta = linspace(N);
	double* gama = linspace(N);
	double* y = linspace(N);
	alpha[0] = A[0][0];
	beta[0] = A[0][1] / A[0][0];
	y[0] = u_before[0] / A[0][0];
	for (int i = 1; i < N; i++)
	{
		gama[i] = A[i][i - 1];
		alpha[i] = A[i][i] - gama[i] * beta[i - 1];
		if (i < N - 1)
		{
			beta[i] = A[i][i + 1] / alpha[i];
		}
		y[i] = (u_before[i] - A[i][i - 1] * y[i - 1]) / alpha[i];
	}
	u_after[N - 1] = y[N - 1];
	for (int i = N - 2; i >= 0; i--)
	{
		u_after[i] = y[i] - beta[i] * u_after[i + 1];
	}
	delete[]y;
	delete[]alpha;
	delete[]beta;
	delete[]gama;
	return;
}
/*******************************************/


/*****************************************************************************************************************/
//Variable-Operator-Length Strategy for RA-based ITSDSGFD
void vtoMvRemez(int N, double dt, double dh, double fmax, double esp0, double eta0, int* Mv)
{
	FILE* fp;
	int Mout = 30;
	int lloutmark = 0;
	int M = 30;
	int v_area = 3501;
	for (int i = 0; i < v_area; i++)
	{
		double v = i + 1500.0;
		double r = v * dt / dh;
		double betamax = 2.0 * PI * fmax * dh / v;
		do {
			////////////////////
			double* du = linspace(M);
			double b = 0.0;
			double** duw = area_2d(N, N);
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
			double* x1 = linspace(N * (N - 1) / 2);
			double* x2 = linspace(M + 1);
			double* bak = linspace(2);
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
			//////////////////////////////////////
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
				//printf("du[%d]=%lf\n", i, du[i]);
			}
			b = x2[M];
			//printf("b=%lf\n", b);
			int flag = 0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					duw[u - 1][w - 1] = x1[flag];
					flag += 1;
					//printf("duw[%d][%d]=%lf\n", u - 1, w - 1, duw[u - 1][w - 1]);
				}
			}
			/////////////////////////////////Remez optimization
			Remez_Optimize_eta(esp0, M, N, r, du, b, duw, bak);
			b = bak[1];
			//////////////
			for (int i = 0; i < M; i++)
			{
				x2[i] = du[i];
			}
			x2[M] = b;
			//printf("b=%lf\n", b);
			flag = 0;
			for (int u = 1; u <= N - 1; u++)
			{
				for (int w = 1; w <= N - u; w++)
				{
					x1[flag] = duw[u - 1][w - 1];
					flag += 1;
				}
			}
			double etamax = 0.0;
			double etamin = 0.0;
			for (double kh = betamax; kh >= 1e-5; kh -= (betamax / 50))
			{
				for (double theta = PI / 2; theta >= 1e-5; theta -= (PI / 100))
				{
					double eta = (dh / v) * (1.0 / dispersion4(kh, theta, r, M, N, x1, x2) - 1);
					if (eta >= etamax)
					{
						etamax = eta;
					}
					if (eta <= etamin)
					{
						etamin = eta;
					}
				}
			}
			if (etamax <= eta0)
			{
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
				Mout = M;
				M = M - 1;
				lloutmark = 0;
			}
			if (etamax > eta0)
			{
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
				lloutmark = 1;
			}
			if (M < N)
			{
				lloutmark = 1;
			}
		} while (lloutmark == 0);
		printf("v=%lf r=%lf M=%d\n", v, r, Mout);
		M = Mout;
		//printf("r=%lf M=%d\n", r, M);
		Mv[i] = M;
	}


	return;
}

///*********************************************************************/