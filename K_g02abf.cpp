/* The purpose of this C++ file is to mimic the NAG function g02abf from https://www.nag.co.uk/numeric/fl/nagdoc_fl26/html/g02/g02abf.html

In this file we will store matrices as arrays or vectors. For example, this 2x3 matrix:
   11 12 13
   21 22 23
would be stored as {11,12,13,21,22,23} where the user should remember the number of rows (2) and columns (3) */

#include <cmath>
#include "K_g02abf.h"
#include "K_jacobi.h"
#include <vector>
#include "K_matrix.h"
#include "K_borsdorf.h"


//The following 3 functions have been copied from https://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
//This website describes an algorithm for finding the nearest correlation matrix to a symmetric square matrix


//Function copied from https://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
//A must be symmetric square nxn matrix
//eigenvectors, eigenvalues, tempArray and result are all empty vectors of length n^2 (nxn matrices) which get filled by function
void K_g02abf::proj_spd(std::vector<double> & A, std::vector<double> & result, int n)
{
	//get eigenvectors and eigenvalues of A
	std::vector<double> eigenvectors;
	std::vector<double> eigenvalues;
	K_jacobi::JacobiEigenvectors(eigenvalues, eigenvectors, A, n);
	//set result to V * diag(max(diag(D), 0))
	for (int j = 0; j < n; j++)
	{
		//get max(0, jth eigenvalue)
		double ej = eigenvalues[j];
		if (ej < 0.0)
			ej = 0.0;
		//multiply jth eigenvector (column j of eigenvectors) by ej
		int in = 0; //wll be i*n
		for (int i = 0; i < n; i++)
		{
			//multiply eigenvectors[i j] by ej
			result[in + j] = ej * eigenvectors[in + j];
			in += n;
		}
	}

	//set tempArray to result*V'
	K_matrix::transposeSquareMatrix(eigenvectors, n);
	std::vector<double> tempArray = K_matrix::mult_ignoreZeros(result, eigenvectors, n, n, n, 0.0);
	//tempArray = K_matrix::MatrixMult_ignoreZeros(result, eigenvectors, n); //this takes ages unless result or eigenvectors is mostly zeros

	/*original matlab:
	[V, D] = eig(A)
	A = V * diag(max(diag(D), 0))*V'
	where D is the diagonal matrix such that Dii = ith eigenvalue and the ith column of V is the corressponding eigenvector
	All symmetric matrices are diagonalisable, so using eigendecomposition of matrices, we have A[i,j] = sum of (Dkk * Vik * Vjk) over k
	By setting the negative Dii to zero, we effectively do:
	A[i,j]  =  sum_k(Dkk * Vik * Vjk)  -  2 * sum_indicesChanged(Dkk * Vik * Vjk)
			=  original A              -  M

	Therefore a more efficient way to change A is to subtract M from it, however I tried this and it didn't work. TC*/

	//Ensure symmetry
	for (int i = 0; i < n; i++)
	{
		const int in = i * n;
		for (int j = 0; j < n; j++)
			result[in + j] = 0.5 * (tempArray[in + j] + tempArray[j * n + i]);
	}
}


//Function copied from https://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
//Replaces diagonal of nxn matrix A with ones
void K_g02abf::proj_unitdiag(std::vector<double> & A, int n)
{
	for (int i = 0; i < n * n; i = i + n + 1)
	{
		A[i] = 1.0;
	}
}


/*Function copied from https://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/ assuming FLAG = PRNT = 0
A is symmetric nxn matrix
X is empty vector before function called
At the end of the function X will be the nearest correllation matrix to A using tolerance tol.
nWeights is vector of n weights*/
std::vector<double> K_g02abf::nearcorr(std::vector<double> & A, int n, double tol, std::vector<double> & nWeights, bool useWeights, int maxits, int* iter)
{
	//   By N.J.Higham, 13 / 6 / 01, updated 30 / 1 / 13, 15 / 11 / 14, 07 / 06 / 15.
	//   Reference:  N.J.Higham, Computing the nearest correlation
	//   matrix-- - A problem from finance.IMA J.Numer.Anal., 22(3) : 329 - 343, 2002.

	//Declare n by n matrices used in the calculation
	std::vector<double> dS(n * n, 0.0);
	std::vector<double> R(n * n, 0.0);
	std::vector<double> R_wtd(n * n, 0.0);
	std::vector<double> Whalf(n * n, 0.0);
	std::vector<double> X(n * n, 0.0); //this will be the nearest correlation matrix
	std::vector<double> Xold(n * n, 0.0);
	std::vector<double> Y(n * n, 0.0);
	std::vector<double> Yold(n * n, 0.0);

	//set w to either nWeights or all 1s
	std::vector<double> w(n, 1.0);
	if (useWeights)
	{
		for (int i = 0; i < n; i++)
		{
			w[i] = nWeights[i];
		}
	}

	*iter = 0;
	double rel_diffX = 100.0 * tol;
	double rel_diffY = 100.0 * tol;
	double rel_diffXY = 100.0 * tol;

	std::vector<double> ww = K_matrix::mult(w, n, 1, w, n);

	//Set X=Y=A and Whalf
	for (int i = 0; i < n * n; i++)
	{
		X[i] = A[i];
		Y[i] = A[i];
		Whalf[i] = sqrt(ww[i]);
	}

	while ((rel_diffX > tol) || (rel_diffY > tol) || (rel_diffXY > tol)) //This while-loop takes over an hour with 700x700 matrices, unless they contain lots of zeros
	{
		//update nxn matrices
		for (int i = 0; i < n * n; i++)
		{
			Xold[i] = X[i];
			R[i] = Y[i] - dS[i];
			R_wtd[i] = Whalf[i] * R[i];
		}

		//Assume R_wtd not updated here
		proj_spd(R_wtd, X, n);

		for (int i = 0; i < n * n; i++)
		{
			X[i] = X[i] / Whalf[i];
			dS[i] = X[i] - R[i];
			Yold[i] = Y[i];
		}

		for (int i = 0; i < n * n; i++)
		{
			Y[i] = X[i];
		}
		proj_unitdiag(Y, n);
		//assume X not updated here

		//frobenius norms
		rel_diffX = K_matrix::FrobDiff(X, Xold, n * n) / K_matrix::FrobeniusNorm(X, n * n);
		rel_diffY = K_matrix::FrobDiff(Y, Yold, n * n) / K_matrix::FrobeniusNorm(Y, n * n);
		rel_diffXY = K_matrix::FrobDiff(Y, X, n * n) / K_matrix::FrobeniusNorm(Y, n * n);

		*iter = *iter + 1;
		if (*iter > maxits)
		{
			break;
		}
	}
	return X;
}

//Mimic NAG function g02abf from https://www.nag.co.uk/numeric/fl/nagdoc_fl26/html/g02/g02abf.html
//Condition: A must be symmetric square matrix
void K_g02abf::g02abf(double* G, const int* N, char OPT [[maybe_unused]], const double* W, const double* ERRTOL, const int* MAXITS, const int* MAXIT [[maybe_unused]], double* X, const int* ITER [[maybe_unused]], int* IFAIL [[maybe_unused]] )
{
	const int n = *N;
	//create vector versions of matrices G and X, and vector W
	std::vector<double> Gvec(n * n, 0.0);
	std::vector<double> Wvec(*N, 0.0);
	//set Wvec = W
	for (int i = 0; i < *N; i++)
		Wvec[i] = W[i];
	//set Gvec = G
	for (int i = 0; i < n * n; i++)
		Gvec[i] = G[i];
	//decide whether to use weights or not
	//If later we use nearcorr() instead of cornewton(), create a boolean called useWeights. If OPT=='A' set useWeights=false, otherwise set it to true

	//max no of iterations for gmres()
	int maxits = 2 * n;
	if (*MAXITS > 0)
		maxits = *MAXITS;
	//max no of iterations of newton method
	int maxit = 200;
	if (*MAXIT > 0)
		maxit = *MAXIT;
	//error tolerance
	double errTol = 0.00000001; //supposed to be *N * sqrt(machine precision)
	if (*ERRTOL > 0.0)
		errTol = *ERRTOL;
	//Get nearest correlation matrix, two ways of doing this (algorithms from sections A.2 and A.3) written in http://eprints.maths.manchester.ac.uk/1085/1/thesis.pdf)
	//To use A.3 type Xvec = K_g02abf::nearcorr(Gvec, n, errTol, Wvec, useWeights, maxits, ITER)
	//However, we use A.2 below:
	bool failed;
	std::vector<double> Xvec = K_borsdorf::cornewton(Gvec, *N, errTol, 0, maxits, maxit, &failed);
	if (failed)
		*IFAIL = 1;
	else
		*IFAIL = 0;
	//update X
	for (int i = 0; i < n * n; i++)
		X[i] = Xvec[i];
	//set G to (G+G^T)/2 with diagonal 1s
	K_matrix::transposeSquareMatrix(Gvec, *N);
	for (int i = 0; i < n * n; i++)
		G[i] = 0.5 * (G[i] + Gvec[i]);
	for (int i = 0; i < *N; i++)
		G[i * (*N) + i] = 1.0;
}
