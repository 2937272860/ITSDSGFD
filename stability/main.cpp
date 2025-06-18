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
	int M = 10;
	FILE* sfp;
	sfp = fopen("sf.txt", "w");
	if (sfp != NULL)
	{
		for (M = 4; M <= 10; M++)
		{
			int N = 4;
			double v = 4500.0;
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
			//////////////////////////////////////优化方法
			double betamax = 2.4;
			double E = 1e-6;
			double eta = 1e-6;
			/////////////////////////////////RA给定误差限,求最大波数和差分系数
			E = eta;
			Remez_Optimize_eta(eta, M, N, r, du, b, duw, bak);
			b = bak[1];
			betamax = bak[0];
			printf("r=%lf eta=%lf betamax=%lf\n", r, eta, bak[0]);

			///////////////////////RA给定最大波数,求误差限和差分系数
			//Remez_Optimize_beta(betamax, M, N, r, du, b, duw, bak);
			//b = bak[1];
			//E = bak[0];
			//printf("r=%lf E=%1.12f betamax=%lf\n", r, E, betamax);

			////////////////////////////////////////////////LS给定误差限,求最大波数和差分系数
			//LS_OptimizeR_eta(eta, M, N, r, du, b, duw, bak);
			//b = bak[1];
			//betamax = bak[0];
			//printf("betamax=%1.12f\n", bak[0]);

			/////////////////////////////////////////给定最大波数，求误差和差分系数
			//double betamax = 2.4;
			//LS_OptimizeR_beta(betamax, M, N, r, du, b, duw, bak);
			//b = bak[1];
			//double E = bak[0];
			//printf("betamax=%lf E=%1.12f\n", betamax, E);

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
			fprintf(sfp, "%d %lf\n", M, s);
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
		}
		fclose(sfp);
	}

	time2 = clock();
	double duration = (time2 - time1) / CLOCKS_PER_SEC;
	printf("参数计算时间%lf秒\n", duration);
	return 0;
}