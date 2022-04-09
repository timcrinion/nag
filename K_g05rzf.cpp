//The purpose of this file is to mimic the NAG function g05rzf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05rzf.html with n=1


#include <cmath>
#include "K_g05saf.h"
#include "K_g05rzf.h"
#include <vector>


//Function that does Cholesky decomposition according to https://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky%E2%80%93Banachiewicz_and_Cholesky%E2%80%93Crout_algorithms
//Given nxn matrix C, returns lower nxn triangular matrix L such that C=LL^T
//C must be a covariance matrix in the form of a vector
std::vector<double> K_g05rzf::CholeskyDecomposition_vec(std::vector<double>& C, int n, bool* fail)
{
	*fail = false;
	//create nxn matrix L full of zeros
	std::vector<double> L(n * n, 0.0);

	//fill L row by row, each time starting from left of row until we reach diagonal
	for (int i = 0; i < n; i++) //row i
	{
		for (int j = 0; j <= i; j++) //row j<=i
		{
			//want sum of all L_ik * L_jk from k = 0 to j-1
			double sum = 0.0;
			int ik = i * n; // will be i*n+k aka index of C such that C[ik] = C[i,k]
			int jk = j * n; // will be j*n+k aka index of C such that C[jk] = C[j,k]
			for (int k = 0; k < j; k++)
				sum += L[ik++] * L[jk++];
			//Set L[i,j]
			if (j < i) //if not on diagonal, set to (A_ij - sum)/L_jj
			{
				if (abs(L[j * n + j]) < 0.0000000001)
					* fail = true;
				else
					L[i * n + j] = (C[i * n + j] - sum) / L[j * n + j];
			}
			else //if on diagonal, set to sqrt(A_jj - sum)
			{
				if (abs(C[j * n + j] - sum) < 0.0)
					* fail = true; // NOT POSITIVE DEFINITE
				else
					L[i * n + j] = sqrt(C[j * n + j] - sum);
			}
		}
	}
	return L;
}


//Multiplies matrix M by column vector v
//v must have length numColsM
//Result Mv will have length numRowsM
std::vector<double> K_g05rzf::MatrixByVector(std::vector<double> M, std::vector<double> v, int numRowsM, int numColsM)
{
	std::vector<double> Mv(numRowsM, 0.0);
	int ij = 0; //index of M such that M[ij] = M[i,j]
	for (int i = 0; i < numRowsM; i++) //calculate Mv[i]
	{
		double Mvi = 0.0;
		for (int j = 0; j < numColsM; j++)
			Mvi += M[ij++] * v[j]; //increase result by M_ij * v_j (then increment ij by 1)
		Mv[i] = Mvi;
	}
	return Mv;
}


//Add two (double) vectors of length n
std::vector<double> K_g05rzf::addVectors(std::vector<double> & vec1, std::vector<double> & vec2, int n)
{
	std::vector<double> sum(n, 0.0);
	for (int i = 0; i < n; i++)
		sum[i] = vec1[i] + vec2[i];
	return sum;
}


//mimic NAG function g05rzf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05rzf.html with n=1
//Output random multivariate normal vector X
//meanX is mean of random vector X of length m
//covX is mxm covariance matrix of X
//state array used to generate random numbers
std::vector<double> K_g05rzf::qq_g05rzf(int* state, std::vector<double> & xMean, std::vector<double> & covX, int m) //maybe introduce new inputs for efficiency: vector lowTriang and boolean redoLowTriang
{
	bool fail;
	const std::vector<double> lowTriang = CholeskyDecomposition_vec(covX, m, &fail);
	//create a sample of m independent standard normal random numbers
	const std::vector<double> sample = K_g05saf::idtStdNormals(state, m);
	//X will be a normal random vector with covariance matrix covX
	//Use lower triangular matrix to create m variables whose covariance matrix is covX
	std::vector<double> X = MatrixByVector(lowTriang, sample, m, m);
	//adjust for mean
	X = addVectors(xMean, X, m);
	return X;
}


//Puts the elements from a double array into a vector
std::vector<double> K_g05rzf::doubleArrayToVector(const double* inArray, int length)
{
	std::vector<double> outVec(length, 0.0);
	for (int i = 0; i < length; i++)
	{
		outVec[i] = inArray[i];
	}
	return outVec;
}

//Check a n by n matrix M is symmetric. Returns true if symmetric, false otherwise.
bool K_g05rzf::isSymmetric(std::vector<double> & M, int n, double tolerance)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++) //for all j>i
		{
			if (abs(M[i * n + j] - M[j * n + i]) > tolerance) // if Mij not equal to Mji
				return false;
		}
	}
	return true;
}