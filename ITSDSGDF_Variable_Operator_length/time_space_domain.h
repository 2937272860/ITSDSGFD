#pragma once
#ifndef _TIME_SPACE_DOMAIN_H
#define _TIME_SPACE_DOMAIN_H
#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>
#include<time.h>
using namespace std;
#define PI 3.1415926
double** area_2d(int NX, int NZ);//����NX��NZ�ж�ά����


double* linspace(int L);//����L����һά����

int* linspace_int(int L);

void delete_2d(double** b, int NX);//deleteNX�еĶ�ά����

double factorial(int N);


void LU(double** a, double* b, double* x, int  N);//���ϵ������A���������b��Nά���Է����飬���ؽ�x�������ʽΪһά����


void solve_fn(int M, double r, double* fn);

double q_coe(int m, int n, int u, int w);

double mnfmn(int m, int n, double* fn);


void fn_to_gn(double* fn, double* gn, double** duw, int M, int N);

double alpha(int n, int u);

double beta(int M, double* gn, int n);


void time_space_domain_coefficience(int M, int N, double r, double* x1, double* x2);


double ricker(double idt, double fm);//��Դ�����׿��Ӳ�������ʱ��ʱ��idt����Ƶfm��


void CatchUpForImplicit(int NX, double b, double* u_before, double* u_after);


void CatchUp(double** A, double* u_before, double* u_after, int N);//���Nά���ԽǾ���


double H(double beta, double r, int N, double** duw);



double error(double beta, int M, int N, double r, double* du, double b, double** duw);




double ErrorA(double beta, int M, int N, double r, double* du, double b, double** duw);
double ErrorR(double beta, int M, int N, double r, double* du, double b, double** duw);


double L2A(double** RombergT, int RomberOrder, double betamax, double r, int M, int N, double* du, double b, double** duw);

double L2R(double** RombergT, int RomberOrder, double betamax, double r, int M, int N, double* du, double b, double** duw);


double f1A(double beta, int u, int j, double r, int N, double** duw);

double f1R(double beta, int u, int j, double r, int N, double** duw);

double a1R(double** RombergT, int RomberOrder, double betamax, int u, int j, double r, int N, double** duw);

double a1A(double** RombergT, int RomberOrder, double betamax, int u, int j, double r, int N, double** duw);



double f2A(double beta, int j, double r, int N, double** duw);

double f2R(double beta, int j, double r, int N, double** duw);

double a2A(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw);
double a2R(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw);

double f3A(double beta, int u, double r, int N, double** duw);

double f3R(double beta, int u, double r, int N, double** duw);


double a3A(double** RombergT, int RomberOrder, double betamax, int u, double r, int N, double** duw);

double a3R(double** RombergT, int RomberOrder, double betamax, int u, double r, int N, double** duw);

double f4R(double beta, double r, int N, double** duw);

double f4A(double beta, double r, int N, double** duw);


double a4A(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw);

double a4R(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw);
double g1R(double beta, int j, double r, int N, double** duw);

double g1A(double beta, int j, double r, int N, double** duw);

double b1R(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw);

double b1A(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw);

double g2A(double beta, double r, int N, double** duw);
double g2R(double beta);

double b2A(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw);

double b2R(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw);

void LS_OptimizeR_eta(double eta, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB);

void LS_OptimizeR_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB);//����������ͽ����������Ӧ�����µ�L1���Ͳ��ϵ��

void Remez_Optimize_eta(double eta, int M, int N, double r, double* du, double b, double** duw, double* bak);
void Remez_Optimize_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* bak);

double dispersion4(double kh, double theta, double r, int M, int N, double* x1, double* x2);//ʱ����Ƶɢ����


void vtoMvRemez(int N, double dt, double dh, double fmax, double esp0, double eta0, int* Mv);
#endif // !_TIME_SPACE_DOMAIN_H

