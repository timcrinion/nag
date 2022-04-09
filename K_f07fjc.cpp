#include <cmath>
#include <vector>
#include "K_g05rzf.h"
#include "K_matrix.h"
#include "K_f07fjc.h"


//fills column_j with the values in jth column of matrix M (where 0 <= j < numCols)
void K_f07fjc::getColumn(std::vector<double>& M, std::vector<double>& column_j, int numRows, int numCols, int j)
{
	int ij = j; //this will run through all indices of M that are its jth column
	for (int i = 0; i < numRows; i++)
	{
		column_j[i] = M[ij];
		ij += numCols;
	}
}

//Sets jth column of matrix M equal to column_j (where 0 <= j < numCols)
void K_f07fjc::setColumn(std::vector<double>& M, std::vector<double>& column_j, int numRows, int numCols, int j)
{
	int ij = j; //this will run through all indices of M that are its jth column
	for (int i = 0; i < numRows; i++)
	{
		M[ij] = column_j[i];
		ij += numCols;
	}
}

//Returns inverse X of lower diagonal n by n matrix L
//Method 2 from section 2.1 of http://www.netlib.org/lapack/lawnspdf/lawn27.pdf
//Assume n>1
//fail set to true if L not invertible
std::vector<double> K_f07fjc::inverseLowerTriang(std::vector<double>& L, int n, bool *invertible)
{
	*invertible = true; //assume invertible until we find a diagonal element of L equal to zero

	//start with X n by n. Fill with zeros
	std::vector<double> X(n * n, 0.0);

	//fill in bottom right element of X. This completes its rightmost column
	if (abs(L[n * n - 1]) < 0.0000000001)
		* invertible = false;
	X[n * n - 1] = 1 / L[n * n - 1]; //X[n-1,n-1] = 1 / L[n-1,n-1]

	//initialise for next loop
	std::vector<double> LCol(n, 0.0); //make length n. will be jth column of L
	std::vector<double> XCol(n, 0.0); //make length n. will be jth column of X

	//fill up X column by column from the right
	for (int j = n - 2; j >= 0; j--)  // j = n-2, ..., 0
	{
		//set Xjj to 1/Ljj
		if ( abs(L[j * n + j]) < 0.0000000001) //if Ljj almost zero
			* invertible = false;
		const double Xjj = 1 / L[j * n + j];
		//set LCol to jth column of L
		getColumn(L, LCol, n, n, j);
		//set XCol to X * LCol
		XCol = K_g05rzf::MatrixByVector(X, LCol, n, n); //try new function that replaces input vector
		//scale XCol by -Xjj
		for (int i = 0; i < n; i++)
			XCol[i] *= -Xjj;
		//set jth column of X to XCol
		setColumn(X, XCol, n, n, j);
		//set jth diagonal element
		X[j * n + j] = Xjj;
	}

	return X;

}

//return inverse of symmetric n by n matrix S
std::vector<double> K_f07fjc::invSymmetric(std::vector<double>& S, int n, bool* fail)
{
	std::vector<double> L = K_g05rzf::CholeskyDecomposition_vec(S, n, fail);

	std::vector<double> invL = inverseLowerTriang(L, n, fail); //fail reset from Cholesky decomposition, change this.

	//get transpose of invL
	std::vector<double> invL_T = invL;
	K_matrix::transposeSquareMatrix(invL_T, n);

	//return invL_T * invL
	//return K_matrix::MatrixMult(invL_T, n, n, invL, n); //returns S^(-1)
	return K_matrix::mult_ignoreZeros(invL_T, invL, n, n, n, 0.0);
}

//Modifies matrix M by swapping row a with row b (0<= a,b < numRowsM)
void K_f07fjc::swapRows(std::vector<double>& M, int numCols, int a, int b)
{
	int aj = a * numCols; //such that M[aj]=M[a,j]
	int bj = b * numCols; //such that M[bj]=M[b,j]
	for (int j = 0; j < numCols; j++)
	{
		const double d = M[aj];
		M[aj++] = M[bj]; //update M[aj], then increment aj
		M[bj++] = d; //update M[bj], then increment bj
	}
}

//Modifies matrix M by multiplying row i by factor (0<= i < numRows)
void K_f07fjc::multRow(std::vector<double>& M, int numCols, int i, double factor)
{
	const int start = i * numCols; //index of M[i,0]
	const int finish = start + numCols; //index of M[i+1,0]
	for (int ij = start; ij < finish; ij++) //want M[ij]=M[i,j]
		M[ij] *= factor;
}

//Modifies matrix M by adding a factor times row b, to row a (0<= a,b < numRows)
void K_f07fjc::addMultRow(std::vector<double> & M, int numCols, int a, int b, double factor)
{
	const int start = a * numCols; //index of M[a,0]
	const int finish = start + numCols; //index of M[a+1,0]
	int bj = b * numCols; //want M[bj]=M[b,j]
	for (int aj = start; aj < finish; aj++) //for j = 0 to numCols, want M[aj] = M[a,j]
		M[aj] += factor * M[bj++]; //update M[aj], then increment bj
}


//Apply Gauss-Jordan elimination to a nxn matrix M
//X input as empty vector. Makes X inverse of M.
//If M invertible, returns true, otherwise false.
bool K_f07fjc::GaussJordanSquareInverse(std::vector<double> & M, std::vector<double> & X, int n)
{
	//make X = nxn identity matrix
	for (int i = 0; i < n * n; i++)
		X.push_back(0);
	for (int i = 0; i < n; i++)
		X[i * n + i] = 1; //Xii=1

	//for j=0 to n-1, set M[j,j]=1 and all other M[i,j] to zero
	for (int j = 0; j < n; j++) //j = 0 ,..., n-1
	{
		//find row>=j such that M[row,j] is closest to 1
		int bestRow = -1;
		for (int row = j; row < n; row++)
		{
			if ( abs(M[row * n + j]) > 0.00000001) //if M[row,j] not zero
			{
				if (bestRow == -1) //if this is the first nonzero row
					bestRow = row;
				else
				{
					if (abs(M[row * n + j]) - 1 < abs(M[bestRow * n + j]) - 1) //if M[row,j] closer to +-1 than M[bestRow,j] is
						bestRow = row;
				}
			}
		}
		if (bestRow == -1)
			return false;

		//Want M[j,j] nonzero. Therefore if necessary, swap row j with highestNonzero
		if (bestRow != j)
		{
			swapRows(M, n, j, bestRow); //swap rows so that M[j,j] nonzero
			swapRows(X, n, j, bestRow); //same for X
		}

		//scale row j so that M[j,j] = 1
		const double Mjj = M[j * n + j];
		multRow(M, n, j, 1 / Mjj); //now M[j,j]=1
		multRow(X, n, j, 1 / Mjj); //same for X

		//set rest of column j (except Mjj) to zero
		for (int i = 0; i < n; i++) //for each row i
		{
			if (i != j)
			{
				const double Mij = M[i * n + j];
				addMultRow(M, n, i, j, -Mij); //added multiple of row j to row i so now Mij=0
				addMultRow(X, n, i, j, -Mij); //same for X
			}
		}
	}
	return true;
}
