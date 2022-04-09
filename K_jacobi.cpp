/*The purpose of this C++ file is to mimic the Jacobi Eigenvalue algorithm.
All four functions in this file have been copied from https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm

In this file we will store matrices as arrays or vectors.For example, this 2x3 matrix :
11 12 13
21 22 23
would be stored as { 11,12,13,21,22,23 } where the user should note the number of rows (2) and columns (3) */

#include <cmath>
#include "K_jacobi.h"
#include "KUtils.h"
#include "K_matrix.h"
#include <vector>

/*Function copied from https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm
 Change two elements of M as shown:
┌   ┐    ┌         ┐┌   ┐
│Mkl│    │cos  −sin││Mkl│
│   │ := │         ││   │
│Mij│    │sin   cos││Mij│
└   ┘    └         ┘└   ┘
Bear in mind that i,j,k,l will range from 0 to numColsM-1 */
void K_jacobi::rotate(int k, int l, int i, int j, std::vector<double>& M, int numColsM, double cosPhi, double sinPhi)
{
	const double Mkl = M[k * numColsM + l];
	const double Mij = M[i * numColsM + j];
	M[k * numColsM + l] = cosPhi * Mkl - sinPhi * Mij;
	M[i * numColsM + j] = sinPhi * Mkl + cosPhi * Mij;
}//end function


/*Function copied from https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm
Returns index (column) of largest right-to-diagonal element in row k of square nxn matrix S
Eg, if k=1 and S = 1 1 1 1 1 then returns 3 because row 1 has largest right-of-diagonal element in column 3
				   1 1 1 2 1
				   1 1 1 1 1
				   1 1 1 1 1
				   1 1 1 1 1
Function should never be called if k = n-1 */
int K_jacobi::maxind(std::vector<double>& S, int n, int k)
{
	int m = k + 1; //m will be maximum index
	for (int i = k + 2; i < n; i++) //from S[k,k+2] to S[k,numColsM-1]
	{
		const double Ski = S[k * n + i];
		const double Skm = S[k * n + m];
		if (abs(Ski) > abs(Skm))
			m = i;
	}//end for
	return m; //luckily function does not get called when k=n-1
}//end function


//Function copied from https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm
//Updates eigenvalues, changed and state
void K_jacobi::update(int k, double t, std::vector<double>& eigenvalues, std::vector<bool>& changed, int* state)
{
	const double y = eigenvalues[k];
	eigenvalues[k] = y + t;
	bool y_eq_ek = (abs(y - eigenvalues[k]) < 0.0000000001);
	if (changed[k] && y_eq_ek)
	{
		changed[k] = false;
		*state = (*state) - 1;
	}
	else if (!changed[k] && !y_eq_ek)
	{
		changed[k] = true;
		*state = (*state) + 1;
	}//end if
}//end function


//Function copied from https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm
//Function that looks at nxn symmetric matrix S and updates eigenvalues and eigenvectors, both of which start as empty vectors
//S is vector of length n^2, where the element S[i,j] is written in S[i*n + j]
//eigenvectors is a vector of length n^2, where the ith column is the eigenvector corresponding to eigenvalue[i]
//bear in mind indices span from 0 to n-1, and S changes during function
void K_jacobi::JacobiEigenvectors(std::vector<double>& eigenvalues, std::vector<double>& eigenvectors, std::vector<double>& S, int n)
{
	//vars
	int state;
	//set state, eigenvectors, eigenvalues, changed and ind with initial values
	state = n;
	std::vector<int> ind(n, 0);
	std::vector<bool> changed(n, true);
	for (int i = 0; i < n * n; i++)
		eigenvectors.push_back(0);
	for (int i = 0; i < n; i++)
	{
		eigenvectors[i * n + i] = 1.0; //identity nxn
		ind[i] = maxind(S, n, i);  //index of largest off-diagonal element in row i
		eigenvalues.push_back(S[i * n + i]); //Sii
	}//end for

	while (state != 0) //next rotation
	{
		int m = 0; //find index(k, l) of pivot p
		double SmMax = S[m * n + ind[m]];
		for (int i = 1; i < n - 1; i++)
		{
			const double SiMax = S[i * n + ind[i]];
			if (abs(SiMax) > abs(SmMax))
			{
				m = i;
				SmMax = SiMax;
			}//end if
		}//end for
		const int k = m;
		const int l = ind[m];
		const double p = S[k * n + l]; //Skl
		//calculate cos φ, s = sin φ
		const double y = 0.5 * (eigenvalues[l] - eigenvalues[k]);
		const double d = abs(y) + sqrt(p * p + y * y);
		const double r = sqrt(p * p + d * d);
		//if r=0, must be all zeros above diagonal (as required) so break now to avoid division by zero
		if (abs(r) < 0.0000000001)
			break;
		const double cosPhi = d / r;
		double sinPhi = p / r;
		//if d=0, must be all zeros above diagonal (as required) so break now to avoid division by zero
		if (abs(d) < 0.0000000001)
			break;
		double t = p * p / d;
		//if p=0, must be all zeros above diagonal (as required) so break now to avoid infinite loop
		if (abs(p) < 0.0000000001)
			break;
		if (y < 0)
		{
			sinPhi = -sinPhi;
			t = -t;
		}//end if
		S[k * n + l] = 0.0; //Skl
		update(k, -t, eigenvalues, changed, &state);
		update(l, t, eigenvalues, changed, &state);
		//rotate rows and columns k and l
		for (int i = 0; i < k; i++)
			rotate(i, k, i, l, S, n, cosPhi, sinPhi);
		for (int i = k + 1; i < l; i++)
			rotate(k, i, i, l, S, n, cosPhi, sinPhi);
		for (int i = l + 1; i < n; i++)
			rotate(k, i, l, i, S, n, cosPhi, sinPhi);
		//rotate eigenvectors
		for (int i = 0; i < n; i++)
			rotate(i, k, i, l, eigenvectors, n, cosPhi, sinPhi); //eigenvectors are columns, not rows
		//rows k, l have changed, update rows indk, indl
		ind[k] = maxind(S, n, k);
		ind[l] = maxind(S, n, l);
	}//end while
}//end function


//Function identical to K_jacobi::JacobiEigenvectors() except that eigenvalues and eigenvectors already have length n and n^2 respectively, and function does not modify S
void K_jacobi::JacobiEigenvectors2(std::vector<double>& eigenvalues, std::vector<double>& eigenvectors, std::vector<double>& S, int n)
{
	//copy of S
	std::vector<double> S_copy = S;
	//create empty vectors
	std::vector<double> tempEigVals;
	std::vector<double> tempEigVecs;
	//fill empty vectors
	JacobiEigenvectors(tempEigVals, tempEigVecs, S_copy, n);
	//modify eigenvalues and eigenvectors
	for (int i = 0; i < n; i++)
		eigenvalues[i] = tempEigVals[i];
	for (int i = 0; i < n * n; i++)
		eigenvectors[i] = tempEigVecs[i];
}


//Find inverse of nxn matrix A using eigendecomposition
//see https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix
//invA input as an empty vector
bool K_jacobi::EigInvSym(std::vector<double>& A, std::vector<double>& invA, int n)
{
	//decompose A into Q * D * invQ where Q is the square nxn matrix whose ith column is ith eigenvector of A, and D is a diagonal matrix where D[i,i] is ith eigenvalue of A.

	//get eigenvalues and eigenvectors
	std::vector<double> Q; //each column an eigenvector
	std::vector<double> eigenvalues; //diagonal of D
	JacobiEigenvectors(eigenvalues, Q, A, n);

	//get inverse of Q
	std::vector<double> invQ = K_matrix::transpose(Q, n, n);			//assume inverse of Q is transpose?

	//replace invQ with inverse(D) * invQ
	for (int i = 0; i < n; i++)
	{
		if ( abs(eigenvalues[i]) < 0.0000001) //if ith eigenvalue almost zero
			return false;
		const double invDii = 1 / eigenvalues[i];
		for (int j = 0; j < n; j++)
			invQ[i * n + j] *= invDii; //divide invQ[i,j] by ith eigenvalue
	}

	//set result as Q * invD * invQ
	invA = K_matrix::mult_ignoreZeros(Q, invQ, n, n, n, 0.0);

	return true;
}