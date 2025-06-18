#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<iostream>
using namespace std;
#define PI 3.1415926
double** area_2d(int NX, int NZ)//生成NX列NZ行二维数组
{
	double** b = new double* [NX];
	for (int i = 0; i < NX; i++)
		b[i] = new double[NZ];
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NZ; j++)
			b[i][j] = 0;//二维double数组赋初值为零
	return b;
}

double* linspace(int L)//定义L长度一维数组
{
	double* l = new double[L];
	for (int i = 0; i < L; i++)
		l[i] = 0;//一维长度赋初值为零
	return l;
}

int* linspace_int(int L)//定义L长度一维数组
{
	int* l = new int[L];
	for (int i = 0; i < L; i++)
		l[i] = 0;//一维长度赋初值为零
	return l;
}

void delete_2d(double** b, int NX)//deleteNX列的二维数组
{
	for (int i = 0; i < NX; i++)
		delete[]b[i];
	delete[]b;
	return;
}

double factorial(int N)
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

void LU(double** a, double* b, double* x, int  N)//求解系数矩阵A，非齐次项b的N维线性方程组，返回解x，解的形式为一维数组
{
	double** l = area_2d(N, N);
	double** u = area_2d(N, N);
	int i, r, k;
	//进行U的第一行的赋值
	for (i = 0; i < N; i++)
	{
		u[0][i] = a[0][i];
	}

	//进行L的第一列的赋值
	for (i = 1; i < N; i++)
	{
		l[i][0] = a[i][0] / u[0][0];
	}

	//计算U的剩下的行数和L的剩下的列数
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



void solve_fn(int M, double r, double* fn)
//求解n=1,2,……，M，M+1时fn的值并存入fn数组
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

void time_space_domain_coefficience(int M, int N, double r, double* x1, double* x2)
{
	double* fn = linspace(M + 1);
	solve_fn(M, r, fn);
	//第二步，计算qmn线性方程组，求duw
	//有效方程组个数
	int L1 = N;
	if (N % 2 == 0)
	{
		L1 = N * N / 4;
	}
	if (N % 2 != 0)
	{
		L1 = (N * N - 1) / 4;
	}
	int L2 = N * (N - 1) / 2;
	//printf("总方程数目L2=%d\n", L2);
	//printf("有效方程组个数L=%d\n", L1);
	//补充方程组个数
	double** supple = area_2d(N, N);
	FILE* fp;
	fp = fopen("qmn_equation.txt", "wb");
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
						//
						//fprintf(fp, "%lf\n", q_coe(m, n, u, w));
						double qoemnuw = q_coe(m, n, u, w);
						fwrite(&qoemnuw, sizeof(double), 1, fp);
					}
				}
				double mnfmnvalue = mnfmn(m, n, fn);
				fwrite(&mnfmnvalue, sizeof(double), 1, fp);
				//fprintf(fp, "%lf\n", mnfmn(m, n, fn));
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
							//fprintf(fp, "%lf\n", supple[u - 1][w - 1]);
							fwrite(&supple[u - 1][w - 1], sizeof(double), 1, fp);
						}
					}
					double supple0 = 0.0;
					//fprintf(fp, "%lf\n", 0.0);
					fwrite(&supple0, sizeof(double), 1, fp);
				}
			}
		}
		fclose(fp);
	}
	delete_2d(supple, N);
	double** qmn_equation = area_2d(L2, L2 + 1);
	double** A1 = area_2d(L2, L2);
	double* b1 = linspace(L2);
	double* x10 = linspace(L2);
	fp = fopen("qmn_equation.txt", "rb");
	if (fp != NULL)
	{
		for (int i = 0; i < L2; i++)
		{
			for (int j = 0; j < L2 + 1; j++)
			{
				//fscanf_s(fp, "%lf", &qmn_equation[i][j]);
				fread(&qmn_equation[i][j], sizeof(double), 1, fp);
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
	LU(A1, b1, x10, L2);
	delete_2d(qmn_equation, L2);
	delete_2d(A1, L2);
	delete[]b1;
	int flag2 = 0;
	double** duw = area_2d(N, N);
	for (int m = 1; m <= N - 1; m++)
	{
		for (int n = 1; n <= N - m; n++)
		{
			duw[m - 1][n - 1] = x1[flag2];
			flag2 += 1;
		}
	}
	double* gn = linspace(M + 1);
	fn_to_gn(fn, gn, duw, M, N);
	double** Mr = area_2d(M + 1, M + 1);
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
	double* x20 = linspace(M + 1);
	LU(Mr, gn, x20, M + 1);
	for (int i = 0; i < L2; i++)
	{
		x1[i] = x10[i];
	}
	for (int i = 0; i <= M; i++)
	{
		x2[i] = x20[i];
	}
	//printf("\n");
	//fprintf(fp1, "\n");
	// 该线性方程组为N*(N-1)/2维线性方程组
	//第三步，计算gn
	//第四步，计算M+1线性方程组，求du，b
	delete_2d(Mr, M + 1);
	delete_2d(duw, N);
	delete[]x10;
	delete[]x20;
	delete[]fn;
	delete[]gn;
}

double ricker(double idt, double fm)//震源函数雷克子波，输入时间时刻idt和主频fm，
{
	double r = 0.0;
	r = (1 - 2 * pow(PI * fm * (idt - 1.0 / fm), 2)) * exp(-pow(PI * fm * (idt - 1.0 / fm), 2));
	return r;
}

void CatchUpForImplicit(int NX, double b, double* u_before, double* u_after)
{
	double* beta = linspace(NX - 1);
	beta[0] = b / (1.0 - 2.0 * b);
	for (int i = 2; i <= NX - 1; i++)
	{
		beta[i - 1] = b / (1.0 - 2.0 * b - b * beta[i - 2]);
	}
	double* y = linspace(NX);
	y[0] = u_before[0] / (1.0 - 2.0 * b);
	for (int i = 2; i <= NX; i++)
	{
		y[i - 1] = (u_before[i - 1] - b * y[i - 2]) / (1.0 - 2.0 * b - b * beta[i - 2]);
	}
	u_after[NX - 1] = y[NX - 1];
	for (int i = NX - 1; i >= 1; i--)
	{
		u_after[i - 1] = y[i - 1] - beta[i - 1] * u_after[i];
	}
	delete[]beta;
	delete[]y;
	return;
}

void CatchUp(double** A, double* u_before, double* u_after, int N)//求解N维三对角矩阵
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

double error(double beta, int M, int N, double r, double* du, double b, double** duw)//相对误差
{
	double s1 = 0.0;
	for (int i = 1; i <= M; i++)
	{
		s1 += du[i - 1] * sin((1.0 * i - 0.5) * beta);
	}
	double s2 = b * (2.0 - 2.0 * cos(beta));
	double s3 = H(beta, r, N, duw);
	double s4 = s2 + s1 / s3 - 1.0;
	return s4;
}

double ErrorA(double beta, int M, int N, double r, double* du, double b, double** duw)//绝对误差
{
	double s1 = 0.0;
	for (int i = 1; i <= M; i++)
	{
		s1 += du[i - 1] * sin((1.0 * i - 0.5) * beta);
	}
	double s2 = b * (2.0 - 2.0 * cos(beta));
	double s3 = H(beta, r, N, duw);
	double s4 = s1 + s2 * s3 - s3;
	return s4;
}

double ErrorR(double beta, int M, int N, double r, double* du, double b, double** duw)//相对误差
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


double L2A(double** RombergT, int RomberOrder, double betamax, double r, int M, int N, double* du, double b, double** duw)//绝对误差L2范数
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = ErrorA(xa, M, N, r, du, b, duw);
	double fb = ErrorA(xb, M, N, r, du, b, duw);
	RombergT[0][0] = h * (fa*fa + fb*fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += ErrorA(xa + (2 * i - 1) * h, M, N, r, du, b, duw)* ErrorA(xa + (2 * i - 1) * h, M, N, r, du, b, duw);///////替换被积函数
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


double L2R(double** RombergT, int RomberOrder, double betamax, double r, int M, int N, double* du, double b, double** duw)//相对误差L2范数
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
	RombergT[0][0] = h * (fa*fa + fb*fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += ErrorR(xa + (2 * i - 1) * h, M, N, r, du, b, duw)* ErrorR(xa + (2 * i - 1) * h, M, N, r, du, b, duw);///////替换被积函数
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


double f1A(double beta, int u, int j, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = sin((1.0 * u - 0.5) * beta) * sin((1.0 * j - 0.5) * beta);
	return s1 ;
}

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
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f1R(xa + (2 * i - 1) * h, u, j, r, N, duw);///////替换被积函数
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

double a1A(double** RombergT, int RomberOrder, double betamax, int u, int j, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f1A(xa, u, j, r, N, duw);
	double fb = f1A(xb, u, j, r, N, duw);;
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f1A(xa + (2 * i - 1) * h, u, j, r, N, duw);///////替换被积函数
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



double f2A(double beta, int j, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s2 = 0.0;
	s2 = H(beta, r, N, duw);
	double s3 = 0.0;
	s3 = sin((1.0 * j - 0.5) * beta);
	double s4 = 0.0;
	s4 = s1 * s3 * s2;
	return s4;
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

double a2A(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f2A(xa, j, r, N, duw);
	double fb = f2A(xb, j, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f2A(xa + (2 * i - 1) * h,j,r,N,duw);///////替换被积函数
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
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f2R(xa + (2 * i - 1) * h, j, r, N, duw);///////替换被积函数
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

double f3A(double beta, int u, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s2 = 0.0;
	s2 = H(beta, r, N, duw);
	double s3 = 0.0;
	s3 = sin((1.0 * u - 0.5) * beta);
	double s4 = 0.0;
	s4 = s1 * s3 * s2;
	return s4;
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


double a3A(double** RombergT, int RomberOrder, double betamax, int u, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f3A(xa, u, r, N, duw);
	double fb = f3A(xb, u, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f3A(xa + (2 * i - 1) * h, u, r, N, duw);///////替换被积函数
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
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f3R(xa + (2 * i - 1) * h, u, r, N, duw);///////替换被积函数
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

double f4A(double beta, double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s2 = H(beta, r, N, duw);
	double s4 = 0.0;
	s4 = s1 * s1 * s2 * s2;
	return s4;
}


double a4A(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = f4A(xa,r,N,duw);
	double fb = f4A(xb,r,N,duw);
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f4A(xa + (2 * i - 1) * h, r, N, duw);///////替换被积函数
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
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += f4R(xa + (2 * i - 1) * h, r, N, duw);///////替换被积函数
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

double g1A(double beta, int j, double r, int N, double** duw)
{
	double s2 = 0.0;
	s2 = H(beta, r, N, duw);
	double s3 = 0.0;
	s3 = sin((1.0 * j - 0.5) * beta);
	return s3 * s2;
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
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += g1R(xa + (2 * i - 1) * h, j, r, N, duw);///////替换被积函数
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

double b1A(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = g1A(xa, j, r, N, duw);
	double fb = g1A(xb, j, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += g1A(xa + (2 * i - 1) * h, j, r, N, duw);///////替换被积函数
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

double g2A(double beta,  double r, int N, double** duw)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	double s2 = H(beta, r, N, duw);
	return s1 * s2 * s2;
}

double g2R(double beta)
{
	double s1 = 0.0;
	s1 = 2.0 - 2.0 * cos(beta);
	return s1;
}

double b2A(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw)
{
	double xa=1e-12;
	double xb = betamax;
	double esp=1e-12;
	int Rout = 0;
	double I = 0.0;
	///////
	double h = (xb - xa) / 2.0;
	double fa = g2A(xa, r, N, duw);
	double fb = g2A(xb, r, N, duw);
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += g2A(xa + (2 * i - 1) * h, r, N, duw);///////替换被积函数
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
	RombergT[0][0] = h * (fa + fb);/////替换被积函数
	int k = 1;
	int n = 1;
	/////
	do {
		double F = 0.0;
		for (int i = 1; i <= n; i++)
		{
			F += g2R(xa + (2 * i - 1) * h);///////替换被积函数
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



void LS_OptimizeR_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB)//给定最大波数和阶数，获得相应波数下的L1误差和差分系数
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
	//////////////////////求解线性方程组，获得优化隐式差分系数
	double Errorlimit = 0.0;
	for (double beta = betamax; beta >=1e-6; beta -= 1e-6)
	{
		double errori = fabs(ErrorR(beta, M, N, r, du, b, duw));
		if (errori >= Errorlimit)
		{
			Errorlimit = errori;
		}
	}
	///////////////////////////遍历
	//printf("betamax=%lf Errorlimit=%1.12f\n", betamax, Errorlimit);
	//优化完毕，释放优化内存
	delete_2d(Mr_new, M + 1);
	delete[]gn_new;
	delete[]x2_new;
	delete_2d(RombergT, RomberOrder);
	//printf("betamax=%lf\n", betamax);
	BetaMaxNewAndB[0] = Errorlimit;
	BetaMaxNewAndB[1] = b;
	return;
}

void LS_OptimizeR_eta(double eta, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB)//给定误差限和阶数，获得最大波数和差分系数
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


double dispersion4(double kh, double theta, double r, int M, int N, double* x1, double* x2)//时空域频散特性
//x1为N*（N-1）/2个duw系数，x2为M+1个du和b系数
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
			double x0 = betai[i - 1];//区间下限
			double x1 = betai[i];//区间上限
			double x2 = (x1+x0)/2;//零点
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
			double left = Remez_root[i];//区间下限
			double right = Remez_root[i + 1];//区间上限
			double midl, midr;
			////////////////////////三分法求函数极大值和极大值点
			while (right - left >= 1e-9)
			{
				midl = (2 * left + right) / 3;
				midr = (left + 2 * right) / 3;
				if (fabs(ErrorR(midl, M, N, r, du, b, duw)) < fabs(ErrorR(midr, M, N, r, du, b, duw)))////求函数极大值
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
		if (fabs(E - E0) / fabs(E0) < 1E-10)//误差收敛
		{
			Remezoutmark = 1;
		}
		E0 = E;
		RemezT += 1;
		if (RemezT > 20)//循环次数太多
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
		Remez_Optimize_beta(betamax, M, N, r, du, b, duw, bak);//二分法获得最大波数
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
		printf("etat=%d betamax=%lf E=%1.12f b-a=%lf\n",etat, betamax, E, betab - betaa);
	} while ((betab-betaa)>=1e-9&&etat<100);
	//printf("r=%lf betamax=%lf E=%1.12f\n", r, betaa, E);

	return;
}

