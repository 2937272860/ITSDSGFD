#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<time.h>
#include"time_space_domain.h"
using namespace std;




int main()
{
	clock_t time1, time2, time3;
	time3 = clock();
	FILE* fp;
	/********************************************************/
	//basic parameter
	int N = 3;
	double dt = 0.001;
	double dh = 10.0;
	double esp0 = 1e-6;
	double fmax = 65.0;
	double eta0 = 1e-6;
	////////
	int NX = 501;
	int NZ = 353;
	int NT = 4001;
	double fm = 23.0;
	/*****************************************************/
	int k0 = 1400;
	int Nlayer = 30;
	NX = NX + 2 * Nlayer;
	NZ = NZ + 2 * Nlayer;
	double R = 0.001;
	int n = 2;
	double vmax = 5000.0;
	//////
	/************************************************/
	int v_area = 3501;
	int* Mv = linspace_int(v_area);
	vtoMvRemez(N, dt, dh, fmax, esp0, eta0, Mv);
	fp = fopen("Mv.txt", "wb");
	if (fp != NULL)
	{
		for (int i = 0; i < v_area; i++)
		{
			fprintf(fp, "%d\n", Mv[i]);
		}
		fclose(fp);
	}
	int Mmax = Mv[0];
	//////
	printf("calculate FD coefficients table now\n");
	int Length = Mmax + 1 + N * (N - 1) / 2;
	double** Parameters = area_2d(v_area, Length);
	double* v_parameters = linspace(Length);

	for (int i = 0; i < v_area; i++)
	{
		double v = i * 1.0 + 1500.0;
		double r = v * dt / dh;
		double betamax = 2.0 * PI * fmax * dh / v;
		int M = Mv[i];
		/////////////////
		double* x1 = linspace(N * (N - 1) / 2);
		double* x2 = linspace(M + 1);
		double* du = linspace(M);
		double b = 0.0;
		double** duw = area_2d(N, N);
		//////////////////////////
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
		Remez_Optimize_eta(esp0, M, N, r, du, b, duw, bak);
		b = bak[1];

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
		////////////////////
		for (int j = 0; j < N * (N - 1) / 2; j++)
		{
			Parameters[i][j] = x1[j];
		}
		for (int j = N * (N - 1) / 2; j < M + 1 + N * (N - 1) / 2; j++)
		{
			Parameters[i][j] = x2[j - N * (N - 1) / 2];
		}
		////////////////////
		delete[]bak;
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
	}
	double** duw = area_2d(N, N);
	double* du = linspace(Mmax);
	double b = 0.0;
	fp = fopen("parameter.dat", "wb");
	if (fp != NULL)
	{
		for (int i = 0; i < v_area; i++)
		{
			for (int j = 0; j < Length; j++)
			{
				float rd = 0.0;
				rd = float(Parameters[i][j]);
				fwrite(&rd, sizeof(float), 1, fp);
				//fprintf(fp, "%lf ", Parameters[i][j]);
			}
			//fprintf(fp, "\n");
		}
		fclose(fp);
	}
	printf("calculate FD coefficients table finished\n");
	/********************Marmousi P wave model***********************/
	double** vp = area_2d(NX, NZ);
	double** rho = area_2d(NX, NZ);
	fp = fopen("vp501_353.dat", "rb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int j = Nlayer; j < NZ - Nlayer; j++)
			{
				float vpt = 0.0;
				fread(&vpt, sizeof(float), 1, fp);
				vp[i][j] = double(vpt);
				rho[i][j] = 0.23 * pow(vp[i][j], 0.25) * 1000.0;
			}
		}
		fclose(fp);
	}

	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			if (i < Nlayer && j < Nlayer)
			{
				vp[i][j] = vp[Nlayer][Nlayer];
			}
			if (i >= NX - Nlayer && j < Nlayer)
			{
				vp[i][j] = vp[NX - Nlayer - 1][Nlayer];
			}
			if (i < Nlayer && j >= NZ - Nlayer)
			{
				vp[i][j] = vp[Nlayer][NZ - Nlayer - 1];
			}
			if (i >= NX - Nlayer && j >= NZ - Nlayer)
			{
				vp[i][j] = vp[NX - Nlayer - 1][NZ - Nlayer - 1];
			}
			if (i < NX - Nlayer && i >= Nlayer && j < Nlayer)
			{
				vp[i][j] = vp[i][Nlayer];
			}
			if (i < Nlayer && j >= Nlayer && j < NZ - Nlayer)
			{
				vp[i][j] = vp[Nlayer][j];
			}
			if (i < NX - Nlayer && i >= Nlayer && j >= NZ - Nlayer)
			{
				vp[i][j] = vp[i][NZ - Nlayer - 1];
			}
			if (i >= NX - Nlayer && j >= Nlayer && j < NZ - Nlayer)
			{
				vp[i][j] = vp[NX - Nlayer - 1][j];
			}
		}
	}
	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NZ; j++)
		{
			if (i < Nlayer && j < Nlayer)
			{
				rho[i][j] = rho[Nlayer][Nlayer];
			}
			if (i >= NX - Nlayer && j < Nlayer)
			{
				rho[i][j] = rho[NX - Nlayer - 1][Nlayer];
			}
			if (i < Nlayer && j >= NZ - Nlayer)
			{
				rho[i][j] = rho[Nlayer][NZ - Nlayer - 1];
			}
			if (i >= NX - Nlayer && j >= NZ - Nlayer)
			{
				rho[i][j] = rho[NX - Nlayer - 1][NZ - Nlayer - 1];
			}
			if (i < NX - Nlayer && i >= Nlayer && j < Nlayer)
			{
				rho[i][j] = rho[i][Nlayer];
			}
			if (i < Nlayer && j >= Nlayer && j < NZ - Nlayer)
			{
				rho[i][j] = rho[Nlayer][j];
			}
			if (i < NX - Nlayer && i >= Nlayer && j >= NZ - Nlayer)
			{
				rho[i][j] = rho[i][NZ - Nlayer - 1];
			}
			if (i >= NX - Nlayer && j >= Nlayer && j < NZ - Nlayer)
			{
				rho[i][j] = rho[NX - Nlayer - 1][j];
			}
		}
	}


	/************Define the Perfectly Matched Layer (PML) absorbing boundary.*******************/
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
			d2[i][j] = pow((Nlayer * 1.0 - 1.0 * j) / Nlayer, n) * log(1.0 / R) * 1.5 * vmax / (Nlayer * dh);
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
	/********************Define wavefield variables.*****************/
	double** vx = area_2d(NX, NZ);
	double** p = area_2d(NX, NZ);
	double** px = area_2d(NX, NZ);
	double** pz = area_2d(NX, NZ);
	double** vz = area_2d(NX, NZ);

	double** dvxdx = area_2d(NX, NZ);
	double** dvzdz = area_2d(NX, NZ);
	double** dpdx = area_2d(NX, NZ);
	double** dpdz = area_2d(NX, NZ);

	double** dvxdxt = area_2d(NX, NZ);
	double** dvzdzt = area_2d(NX, NZ);
	double** dpdxt = area_2d(NX, NZ);
	double** dpdzt = area_2d(NX, NZ);

	double* u_before1 = linspace(NX);
	double* u_before2 = linspace(NZ);
	double* u_before3 = linspace(NX);
	double* u_before4 = linspace(NZ);

	double** b_alpha1 = area_2d(NX, NX);
	double** b_alpha2 = area_2d(NZ, NZ);
	double** b_alpha3 = area_2d(NX, NX);
	double** b_alpha4 = area_2d(NZ, NZ);

	double* u_after1 = linspace(NX);
	double* u_after2 = linspace(NZ);
	double* u_after3 = linspace(NX);
	double* u_after4 = linspace(NZ);

	/**********************************************Define observation system********************************/
	double** DT = area_2d(NX, NZ);
	double** wavefront_p = area_2d(NX, NZ);
	double** record_p = area_2d(NX, NT);
	double* trace = linspace(NT);
	double offset0 = 2500.0;
	double depth0 = 400.0;
	int NX0 = int(offset0 / dh);
	int NZ0 = int(depth0 / dh);
	//DT[NX / 2][Nlayer] = 1.0;//
	double shotoffset = 20.0;
	double shotdepth = 20.0;
	DT[Nlayer + int(shotoffset / dh)][Nlayer + int(shotdepth / dh)] = 1.0;// source
	/*******************************Forward Modeling*******************************/
	time1 = clock();
	double duration = 0.0;
	duration = (time1 - time3) / CLOCKS_PER_SEC;
	printf("The preparatory work for forward modeling took %lf seconds.\n", duration);
	for (int k = 0; k < NT; k++)
	{
		if (k % 50 == 0)
		{
			printf("k=%d,p[NX/2][NZ/2]=%lf\n", k, p[NX / 2][NZ / 2]);
		}
		for (int j = 0; j < NZ; j++)
		{
			for (int i = 0; i < NX; i++)
			{
				int flag = int(vp[i][j]) - 1500;
				int M = Mv[flag];
				for (int v = 0; v < Length; v++)
				{
					v_parameters[v] = Parameters[flag][v];
				}
				for (int u = 1; u <= M; u++)
				{
					du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
				}
				b = v_parameters[M + N * (N - 1) / 2];
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
				int M = Mv[flag];
				for (int v = 0; v < Length; v++)
				{
					v_parameters[v] = Parameters[flag][v];
				}
				for (int u = 1; u <= M; u++)
				{
					du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
				}
				b = v_parameters[M + N * (N - 1) / 2];
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
			CatchUp(b_alpha2, u_before2, u_after2, NZ);;
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
				p[i][j] = px[i][j] + pz[i][j] + DT[i][j] * ricker(k * dt, fm);
			}
		}
		for (int j = 0; j < NZ; j++)
		{
			for (int i = 0; i < NX; i++)
			{
				int flag = int(vp[i][j]) - 1500;
				int M = Mv[flag];
				for (int v = 0; v < Length; v++)
				{
					v_parameters[v] = Parameters[flag][v];
				}
				for (int u = 1; u <= M; u++)
				{
					du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
				}
				b = v_parameters[M + N * (N - 1) / 2];
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
				int M = Mv[flag];
				for (int v = 0; v < Length; v++)
				{
					v_parameters[v] = Parameters[flag][v];
				}
				for (int u = 1; u <= M; u++)
				{
					du[u - 1] = v_parameters[u - 1 + N * (N - 1) / 2];
				}
				b = v_parameters[M + N * (N - 1) / 2];
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
					wavefront_p[i][j] = p[i][j];
				}
			}
		}
		for (int i = 0; i < NX; i++)
		{
			record_p[i][k] = p[i][Nlayer];
		}
		trace[k] = p[NX0 + Nlayer][NZ0 + Nlayer];
	}
	time2 = clock();
	duration = (time2 - time1) / CLOCKS_PER_SEC;
	printf("The main loop took %lf seconds.\n", duration);

	duration = (time2 - time3) / CLOCKS_PER_SEC;
	printf("The total forward modeling time is %lf seconds.\n", duration);

	/****************Output************/
	fp = fopen("wavefront.dat", "wb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int j = Nlayer; j < NZ - Nlayer; j++)
			{
				float wp = float(wavefront_p[i][j]);
				fwrite(&wp, sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}
	fp = fopen("recordvol.dat", "wb");
	if (fp != NULL)
	{
		for (int i = Nlayer; i < NX - Nlayer; i++)
		{
			for (int k = 0; k < NT; k++)
			{
				float rd = float(record_p[i][k]);
				fwrite(&rd, sizeof(float), 1, fp);
			}
		}
		fclose(fp);
	}
	fp = fopen("recordvol.txt", "wb");
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
	fp = fopen("tracevol.txt", "w");
	if (fp != NULL)
	{
		for (int k = 0; k < NT; k++)
		{
			fprintf(fp, "%lf\n", trace[k]);
		}
		fclose(fp);
	}

	/**********free space*************/
	delete_2d(d1, NX);
	delete_2d(d2, NX);
	delete_2d(d3, NX);
	delete_2d(d4, NX);

	delete_2d(vp, NX);
	delete_2d(rho, NX);

	delete[]u_before4;
	delete[]u_before3;
	delete[]u_before2;
	delete[]u_before1;

	delete[]u_after4;
	delete[]u_after3;
	delete[]u_after2;
	delete[]u_after1;

	delete_2d(dpdz, NX);
	delete_2d(dpdx, NX);
	delete_2d(dvzdz, NX);
	delete_2d(dvxdx, NX);

	delete_2d(dpdzt, NX);
	delete_2d(dpdxt, NX);
	delete_2d(dvzdzt, NX);
	delete_2d(dvxdxt, NX);

	delete_2d(vz, NX);
	delete_2d(vx, NX);
	delete_2d(p, NX);
	delete_2d(px, NX);
	delete_2d(pz, NX);

	delete_2d(b_alpha1, NX);
	delete_2d(b_alpha2, NZ);
	delete_2d(b_alpha3, NX);
	delete_2d(b_alpha4, NZ);

	delete_2d(DT, NX);
	delete_2d(wavefront_p, NX);
	delete_2d(record_p, NX);
	delete[]trace;
	///////////////
	delete_2d(Parameters, v_area);
	delete[]v_parameters;
	delete[]Mv;
	delete[]duw;
	delete[]du;
	return 0;
}
