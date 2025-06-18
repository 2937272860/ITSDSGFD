#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<time.h>
#include"time_space_domain.h"
using namespace std;



int main()
{
	clock_t time1, time2;
	time1 = clock();
	FILE* fp;
	int M = 6;
	int N = 3;
	double v = 1500.0;
	double dt = 0.001;
	double dh = 10.0;
	double r = v * dt / dh;
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
	printf("空间差分系数个数M+1=%d\n", M+1);
	printf("有效方程组个数L=%d\n", L1);
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
	//////////////////////////////////////TE常规方法求解隐式差分系数du，b
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
		printf("du[%d]=%lf\n", i, du[i]);
	}
	b = x2[M];
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
	//////////////////////////////////////优化方法
	double betamax = 2.4;
	double E = 1e-6;
	double eta = 1e-6;
	double* bak = linspace(2);
	/////////////////////////////////Remez给定误差限求最大波数和差分系数
	Remez_Optimize_eta(eta, M, N, r, du, b, duw, bak);
	b = bak[1];
	betamax = bak[0];
	printf("r=%lf eta=%lf betamax=%lf\n", r, eta, bak[0]);



	///////////////////////Remez给定最大波数求误差限和差分系数
	//betamax = 2.46;
	//Remez_Optimize_beta(betamax, M, N, r, du, b, duw, bak);
	//b = bak[1];
	//E = bak[0];
	//printf("r=%lf E=%1.12f betamax=%lf\n", r, E, betamax);
	// 
	//////////////////////////////////////////////LS给定误差限,求最大波数和差分系数
	//LS_OptimizeR_eta(eta, M, N, r, du, b, duw, bak);
	//b = bak[1];
	//betamax = bak[0];
	//printf("betamax=%1.12f\n", bak[0]);




	/////////////////////////////////////////LS给定最大波数，求误差和差分系数
	//betamax = 2.7380;
	//LS_OptimizeR_beta(betamax, M, N, r, du, b, duw, bak);
	//b = bak[1];
	//E = bak[0];
	//printf("betamax=%lf E=%1.12f\n", betamax, E);







	//////////////
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

	for (int i = 0; i < M; i++)
	{
		x2[i] = du[i];
		printf("du[%d]=%1.12f\n", i, du[i]);
	}
	x2[M] = b;
	printf("b=%1.12f\n", b);
	flag = 0;
	for (int u = 1; u <= N - 1; u++)
	{
		for (int w = 1; w <= N - u; w++)
		{
			x1[flag] = duw[u - 1][w - 1];
			flag += 1;
			printf("duw[%d][%d]=%1.12f\n", u - 1, w - 1, duw[u - 1][w - 1]);
		}
	}
	fp = fopen("coef.txt", "w");
	if (fp != NULL)
	{
		for (int i = 0; i < M; i++)
		{
			fprintf(fp,"%1.12f\n",  du[i]);
		}
		x2[M] = b;
		fprintf(fp,"%1.12f\n", b);
		flag = 0;
		for (int u = 1; u <= N - 1; u++)
		{
			for (int w = 1; w <= N - u; w++)
			{
				x1[flag] = duw[u - 1][w - 1];
				flag += 1;
				fprintf(fp,"%1.12f\n",duw[u - 1][w - 1]);
			}
		}
		fclose(fp);
	}
	/////////////////

	fp = fopen("disp.txt", "w");
	if (fp != NULL)
	{
		//theta入射角
		for (int j = 0; j < 1000; j++)
		{
			double kh = j * PI / 1000 + 1e-4;
			fprintf(fp, "%1.12f ", kh);
			for (int i = 0; i < 5; i++)//kh波数
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
		//theta入射角
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
	/////////////
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
	////////

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
	////////////
	double s = 0.0;
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
	s = 1.0 / (sqrt(2.0) * (s1 + s2));
	printf("M=%d,N=%d,r=%lf,s=%lf\n", M, N, r, s);
	////////////////
	/////////
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
	printf("参数计算时间%lf秒\n", duration);
	return 0;
}