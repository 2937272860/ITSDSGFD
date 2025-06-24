#pragma once
#ifndef _TIME_SPACE_DOMAIN_H
#define _TIME_SPACE_DOMAIN_H
#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>
#include<time.h>
using namespace std;
#define PI 3.1415926

/******************************************/
//basic function
double** area_2d(int NX, int NZ);//Generate a two-dimensional array with NX columns and NZ rows.

double* linspace(int L);//Define a one-dimensional double array with length L.

int* linspace_int(int L);//Define a one-dimensional int array with length L.

void delete_2d(double** b, int NX);//Release a two-dimensional array with NX columns.

double factorial(int N);//Factorial function

void LU(double** a, double* b, double* x, int  N);
//Solve the N-dimensional linear system of equations with coefficient matrix A and nonhomogeneous term b using LU decomposition,
// and return the solution x in the form of a one-dimensional array.

/******************************************/
//ITSDSGFD coefficients calculated by TE method
/**********/
void solve_fn(int M, double r, double* fn);
//Calculate the values of fn for ( n = 1, 2, ..., M, M+1 ) and store them in the array fn.

double q_coe(int m, int n, int u, int w);
double mnfmn(int m, int n, double* fn);
void fn_to_gn(double* fn, double* gn, double** duw, int M, int N);
double alpha(int n, int u);
double beta(int M, double* gn, int n);
/**********/

/*************************************************/
//optimization function
double H(double beta, double r, int N, double** duw);

double ErrorR(double beta, int M, int N, double r, double* du, double b, double** duw);//Relative error


double L2R(double** RombergT, int RomberOrder, double betamax, double r, int M, int N, double* du, double b, double** duw);//The L2 norm of relative error.

/*************************************************/
/*
When optimizing the difference coefficients using the least squares method,it is necessary to solve a system of linear equations.
At this time, both the coefficient matrix and the non - homogeneous term of the system of equations are integral expressions,
and numerical integration methods need to be used to solve them.
*/
//All numerical integrations in this program group are solved using the Romberg integration formula.

double f1R(double beta, int u, int j, double r, int N, double** duw);

double a1R(double** RombergT, int RomberOrder, double betamax, int u, int j, double r, int N, double** duw);

double f2R(double beta, int j, double r, int N, double** duw);

double a2R(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw);

double f3R(double beta, int u, double r, int N, double** duw);

double a3R(double** RombergT, int RomberOrder, double betamax, int u, double r, int N, double** duw);

double f4R(double beta, double r, int N, double** duw);

double a4R(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw);

double g1R(double beta, int j, double r, int N, double** duw);

double b1R(double** RombergT, int RomberOrder, double betamax, int j, double r, int N, double** duw);

double g2R(double beta);

double b2R(double** RombergT, int RomberOrder, double betamax, double r, int N, double** duw);
/*************************************************/


// Optimize the spatial implicit FD coefficients with the Least square method(LS) under the premise of a given maximum wavenumber range.
// When outputting the difference coefficients, also output the error values corresponding to the wavenumber range, FD operator length M,N, Courant number r and temporal FD coefficients duw.
void LS_OptimizeR_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB);


// Optimize the spatial implicit FD coefficients using the Least square method(LS) under the premise of a given error limitation eta. 
// While outputting the difference coefficients, output the maximum wavenumber range that meets the error limitation eta.
//Use the bisection method to find the maximum wavenumber that satisfies the error limitation.
void LS_OptimizeR_eta(double eta, int M, int N, double r, double* du, double b, double** duw, double* BetaMaxNewAndB);

double dispersion4(double kh, double theta, double r, int M, int N, double* x1, double* x2);//phase velocity dispersion
//x1 stores temporal FD coefficients, and x2 stores spatial FD coefficients



// Optimize the spatial implicit FD coefficients with the Remez exchange algorithm (RA) under the premise of a given maximum wavenumber range.
// When outputting the difference coefficients, also output the error values E corresponding to the wavenumber range, FD operator length M,N, Courant number r and temporal FD coefficients duw.
void Remez_Optimize_beta(double betamax, int M, int N, double r, double* du, double b, double** duw, double* bak);

// Optimize the spatial implicit FD coefficients using the Least square method(LS) under the premise of a given error limitation eta. 
// While outputting the difference coefficients, output the maximum wavenumber range that meets the error limitation eta.
//Use the bisection method to find the maximum wavenumber that satisfies the error limitation.
void Remez_Optimize_eta(double eta, int M, int N, double r, double* du, double b, double** duw, double* bak);



/************************************************/
//necessary function for modeling

double ricker(double idt, double fm);
//The source function is a Ricker wavelet, 
// which takes the time instant idt and the dominant frequency fm as inputs and outputs the wavefield value at the corresponding instant.



void CatchUp(double** A, double* u_before, double* u_after, int N);
//The CatchUp function is used to solve a system of linear equations with a tridiagonal matrix as the coefficient matrix. 
// Here, double** A is a two-dimensional array storing the coefficient matrix,
// double* u_before is the non-homogeneous term, 
// double* u_after is the solution to the equations,
// and int N is the order of the equations.

/*******************************************/

#endif // !_TIME_SPACE_DOMAIN_H

