#include <cmath>
#include <vector>
#include "KUtils.h"
#include "K_matrix.h"

/*This file contains functions for dealing with matrices.
In this file we will store matrices as vectors by concatenating rows from top to bottom. For example, this 2x3 matrix:
   11 12 13
   21 22 23
would be stored as the vector {11,12,13,21,22,23} where the user should remember the number of rows (2) and columns (3)*/

//Function that multiplies matrix A by matrix B and stores answer in AB.
//A must have numRowsA rows and numColsA columns.
//B must have numColsA rows and numColsB columns.
//AB is matrix of dimensions numRowsA by numColsB.
std::vector<double> K_matrix::mult(std::vector<double>& A, int numRowsA, int numColsA, std::vector<double>& B, int numColsB)
{
	//Set size of AB, fill with zeros
	std::vector<double> AB(numRowsA * numColsB, 0.0);

	std::vector<double> BT = transpose(B, numColsA, numColsB); //note that BT has same dimensions as A

	int i0 = 0; //index of AB such that  A[i0] =  A[i,0]

	//For each row i of A, and column j of B, calculate entry [i,j] of AB:
	for (int i = 0; i < numRowsA; i++) //go through rows of A
	{
		for (int j = 0; j < numColsB; j++) //go through columns of B
		{
			//Want ABij = Ai0*B0j + Ai1*B1j + Ai2*B2j + etc
			double ABij = 0.0;
			int ik = i * numColsA;
			int jk = j * numColsA;
			for (int k = 0; k < numColsA; k++)
				ABij += A[ik++] * BT[jk++]; //increase ABij by A[i,k]*B[k,j]
			AB[i0 + j] = ABij;
		}//end for
		i0 += numColsB; //index of A such that A[i0] = A[i,0]
	}//end for
	return AB;
}//end function


//Function that transposes square nxn matrix M
void K_matrix::transposeSquareMatrix(std::vector<double>& M, int n) //can make more efficient using vector swap function?
{
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			const double Mij = M[i * n + j]; //pull out copy of Mij
			M[i * n + j] = M[j * n + i]; //set Mij = Mji
			M[j * n + i] = Mij; //set Mji = Mij
		}//end for
	}//end for
}//end function


//Function that transposes matrix M
//M[i,j] stored as M[i * numColsM + j]
std::vector<double> K_matrix::transpose(std::vector<double>& M, int numRowsM, int numColsM)
{
	//initialise result
	std::vector<double> result(numRowsM * numColsM, 0.0);

	//declare these now to make next for loops more efficient
	int ij = 0; //will be index i*numColsM+j such that M[ij] = M[i,j]

	for (int i = 0; i < numRowsM; i++)
	{
		int ji = i;
		for (int j = 0; j < numColsM; j++)
		{
			result[ji] = M[ij++];
			ji += numRowsM;
		}//end for
	}
	return result;
}//end function


//return Frobenius norm of array x of length n
double K_matrix::FrobeniusNorm(std::vector<double>& x, int n)
{
	double result = 0.0;
	for (int i = 0; i < n; i++)
		result += (x[i] * x[i]);
	return sqrt(result);
}

//return Frobenius difference between two arrays of length n
double K_matrix::FrobDiff(std::vector<double>& x, std::vector<double>& y, int n)
{
	double result = 0.0;
	for (int i = 0; i < n; i++)
	{
		const double diff = x[i] - y[i];
		result += diff * diff;
	}
	return sqrt(result);
}


//multiply A by B, output AB with dimensions numRowsA by numColsB
//B must have numColsA rows
//For efficiency, function ignores all entries <= tinyNumber
std::vector<double> K_matrix::mult_ignoreZeros(std::vector<double>& A, std::vector<double>& B, int numRowsA, int numColsA, int numColsB, double tinyNumber)
{
	const int numRowsB = numColsA; //must be the case for multiplication to be possible

	//Initialise result
	std::vector<double> AB(numRowsA * numColsB, 0.0); //numRowsA by numColsB

	//create vectors that will inform loop below of nonzero entries in B
	std::vector<int> numNonzerosB(numRowsB, 0); //length numRowsB
	std::vector<int> nonzerosB(numRowsB * numColsB, -1); //same dimensions as B
	//want numNonzerosB[i] = the number of nonzero values in ith row B[i,:]
	//want nonzerosB[i,j] = k to mean that B[i,k] is the jth nonzero in ith row B[i,:]
	//For example, numNonzerosB[i]=3 and nonzerosB[i,:] = { 0,2,5,-1,-1,-1 } means that only B[i,0], B[i,2] and B[i,5] are nonzero in B's ith row
	int i0 = 0; //such that B[i0] = B[i,0] = B[i * numColsB]
	int ij = 0; //such that B[ij] = B[i,j] = B[i * numColsB + j]
	for (int i = 0; i < numColsA; i++)//for row B[i,:]
	{
		for (int j = 0; j < numColsB; j++)//consider B[i,j]
		{
			if (abs(B[ij]) > tinyNumber) //if B[i,j] nonzero
			{
				nonzerosB[i0 + numNonzerosB[i]] = j;
				numNonzerosB[i]++;
			}
			ij++;
		}
		i0 += numColsB;
	}

	//want to set AB[i,j] = A[i,0]*B[0,j] + A[i,1]*B[1,j] + A[i,2]*B[2,j] + etc
	int ik = 0; //want A[ik] = A[i,k]
	i0 = 0; //want AB[i0] = AB[i,0]
	for (int i = 0; i < numRowsA; i++) //row i of A
	{
		int k0 = 0;
		for (int k = 0; k < numColsA; k++) //row k of A
		{
			const double Aik = A[ik++];
			if (abs(Aik) > tinyNumber) //if A[i,k] nonzero
			{
				for (int kk = 0; kk < numNonzerosB[k]; kk++) //for each nonzero in B[k,:]
				{
					const int j = nonzerosB[k0 + kk]; //nonzerosB[k,kk], so that Bkj is kkth nonzero on row B[k,:]
					AB[i0 + j] += Aik * B[k0 + j]; //increase ABij by Aik * Bkj
				}
			}
			k0 += numColsB;
		}
		i0 += numColsB;
	}

	return AB;
}


//Returns inner product u[1]v[1] + u[2]v[2] + ... u[n]v[n] of two vectors u and v of length n
double K_matrix::innerProduct(std::vector<double> & u, std::vector<double> & v, int n)
{
	double result = 0.0;
	for (int i = 0; i < n; i++)
		result += u[i] * v[i];
	return result;
}

//can do using  Gram–Schmidt process, Householder transformations, or Givens rotations.

//https://en.wikipedia.org/wiki/QR_decomposition#Using_Givens_rotations
//input:
//A, rowsA by colsA where rowsA >= colsA
//output:
//Q, unitary, rowsA by rowsA
//R, upper triangular, rowsA by colsA
void K_matrix::QR(std::vector<double> & A, int rowsA, int colsA, std::vector<double> & Q, std::vector<double> & R)
{
	//set Q to identity
	for (int i = 0; i < rowsA * rowsA; i++)
		Q[i] = 0.0;
	for (int i = 0; i < rowsA; i++)
		Q[i * rowsA + i] = 1.0;

	//set R to A
	for (int i = 0; i < rowsA * colsA; i++)
		R[i] = A[i];

	//apply givens rotations to Q and R
	for (int j = 0; j < colsA; j++) //for each column j of R
	{
		for (int i = j + 1; i < rowsA; i++)//for each i>j (Rij under diagonal)
		{
			double Rij = R[i * colsA + j];
			double Rjj = R[j * colsA + j];
			double div = sqrt(Rij * Rij + Rjj * Rjj);
			if (div > 0.0000000000000001)
			{
				//find cos and sin such that givens rotation sets Rij to zero
				double cos = Rjj / div;
				double sin = -Rij / div;
				//set Q to GQ
				K_matrix::givensLeft(cos, sin, Q, rowsA, i, j);
				//set R to GR
				K_matrix::givensLeft(cos, sin, R, colsA, i, j);
			}
		}
	}

	//transpose Q (inverses it)
	K_matrix::transposeSquareMatrix(Q, rowsA);

}


/*
function that replaces M with G*M, where G is the givens rotation below

	| 1            |
	|  cos  -sin   | <- row j
	|     1        |
G = |       1      |
	|  sin   cos   | <- row i
	|           1  |
	|             1|

G has dimensions rowsM by rowsM, square
M has dimensions rowsM by colsM
Cannot have i=j
*/
void K_matrix::givensLeft(double cos, double sin, std::vector<double> & M, int colsM, int i, int j)
{
	if (abs(cos - 1.0) + abs(sin) < 0.000000001) //if givens matrix is identity, return
		return;
	//Change row i and j of M:
	//	row j = cos*row j - sin*row i
	//	row i = sin*row j + cos*row i
	double Mik;
	double Mjk;
	int jk = j * colsM; //will be such that M[jk] = M[j,k] in loop below
	int ik = i * colsM; //will be such that M[ik] = M[i,k] in loop below
	for (int k = 0; k < colsM; k++)
	{
		Mik = M[ik];
		Mjk = M[jk];
		M[jk] = cos * Mjk - sin * Mik;
		M[ik] = sin * Mjk + cos * Mik;
		ik++;
		jk++;
	}
}




/*
function that replaces M with M*G, where G is the givens rotation below

	| 1            |
	|  cos  -sin   | <- row i
	|     1        |
G = |       1      |
	|  sin   cos   | <- row j
	|           1  |
	|             1|

G has dimensions rowsM by rowsM, square
M has dimensions rowsM by colsM
*/
//Same as givensLeft, only returns MG instead of GM
void K_matrix::givensRight(double cos, double sin, std::vector<double> & M, int colsM, int i, int j)
{
	if (abs(cos - 1.0) + abs(sin) < 0.000000001) //if givens matrix is identity, return
		return;
	//Change column i and j of M:
	//	column i =  cos*column i + sin*column j
	//	column j = -sin*column i + cos*column j
	double Mki;
	double Mkj;
	int ki = i; //will be such that M[ki] = M[k,i] in loop below
	int kj = j; //will be such that M[kj] = M[k,j] in loop below
	for (int k = 0; k < colsM; k++)
	{
		Mki = M[ki];
		Mkj = M[kj];
		M[ki] = cos * Mki + sin * Mkj;
		M[kj] = -sin * Mki + cos * Mkj;
		ki += colsM;
		kj += colsM;
	}
}





//use KUtils::min instead
double min2(double a, double b)
{
	if (a < b)
		return a;
	return b;
}


/*applies Singular Value Decomposition using algorithm 7 from https://www.cs.utexas.edu/~inderjit/public_papers/HLA_SVD.pdf using m=n
decompose square n by n matrix A into A = U * Sigma * V where
	U n by n with orthogonal vectors
	V n by n with orthogonal vectors
	Sigma is diagonal, n by n
*/
void K_matrix::SVD(std::vector<double> & A, int n, std::vector<double> & U, std::vector<double> & Sigma, std::vector<double> & V, double epsilon)
{
	//set B=A
	std::vector<double> B(n * n);
	for (int i = 0; i < n * n; i++)
		B[i] = A[i];

	//set U = I(n,n)
	for (int i = 0; i < n * n; i++)
		U[i] = 0.0;
	for (int i = 0; i < n; i++)
		U[i * n + i] = 1.0;

	//set V = I(n,n)
	for (int i = 0; i < n * n; i++)
		V[i] = 0.0;
	for (int i = 0; i < n; i++)
		V[i * n + i] = 1.0;

	//Next step begins "if m>n..." Ignore.

	//step 5
	double N2 = 0.0;
	for (int i = 0; i < n * n; i++)
		N2 += B[i] * B[i];
	double s = 0.0;

	//step 6
	int counter = 0;
	while (s > epsilon * epsilon * N2 || counter == 0)
	{
		//a
		s = 0.0;
		counter++;

		//b
		for (int i = 0; i < n - 1; i++) //i=0...n-2
		{
			for (int j = i + 1; j < n; j++) //j=i+1...n-1
			{
				double Bii = B[i * n + i];
				double Bij = B[i * n + j];
				double Bji = B[j * n + i];
				double Bjj = B[j * n + j];
				s += Bij * Bij + Bji * Bji;
				double sinLeft;
				double cosLeft;
				double sinRight;
				double cosRight;
				K_matrix::computeGivens(Bii, Bij, Bji, Bjj, sinLeft, cosLeft, sinRight, cosRight); //-sinLeft and -sinRight on row i
				//apply givens rotations to U, B and V
				K_matrix::givensLeft(cosLeft, -sinLeft, B, n, i, j); //row j has -sin
				K_matrix::givensRight(cosRight, sinRight, B, n, i, j); //row i has -sin
				K_matrix::givensLeft(cosLeft, -sinLeft, U, n, i, j);
				K_matrix::givensRight(cosRight, sinRight, V, n, i, j);
			}
		}

		if (counter > 100)
			break;
	}

	//step 7
	for (int i = 0; i < n; i++)
		Sigma[i * n + i] = B[i * n + i];

	//transpose U and V (inverses them)
	K_matrix::transposeSquareMatrix(U, n);
	K_matrix::transposeSquareMatrix(V, n);
}





/*
function that outputs sinL, cosL, sinR and cosR such that the matrix below has [1,2] = [2,1] = 0, as shown:

 |cosL -sinL|   |m11 m12|   |cosR -sinR|   |   0|
 |sinL  cosL| * |m21 m22| * |sinR  cosR| = |0   |

Want to satisfy equations 1 and 2
eq 1: (m12 cosL - m22 sinL) cosR = (m11 cosL - m21 sinL) sinR
eq 2: (m11 sinL + m21 cosL) cosR =-(m12 sinL + m22 cosL) sinR
*/

void K_matrix::computeGivens(double m11, double m12, double m21, double m22, double& sinL, double& cosL, double& sinR, double& cosR)
{
	//deal with special m12 = m21 = 0
	if (abs(m12) < 0.000000000001 && abs(m21) < 0.000000000001)
	{
		cosL = 1.0;
		sinL = 0.0;
		cosR = 1.0;
		sinR = 0.0;
		return;
	}

	//Deduce third equation from first two, regardless of sinR and cosR values:
	//	eq 3: (m12 cosL - m22 sinL)(m12 sinL + m22 cosL) = (m21 sinL - m11 cosL)(m11 sinL + m21 cosL)

	//Find sinL and cosL that satisfy eq 3:
	/*eq 3 true iff p(sinL^2 - cosL^2) = q sinL cosL, where
		  p = m12 m22 + m21 m11
		  q = m11^2 + m12^2 - m21^2 - m22^2
	*/
	double q = m11 * m11 + m12 * m12 - m21 * m21 - m22 * m22;
	double qq = q * q;

	//Thus eq 3 is satisfied by the following values for sinLeft and cosLeft
	if (qq < 0.0000000001) //q=0, so m11^2 + m12^2 = m21^2 + m22^2
	{
		sinL = sqrt(0.5);
		cosL = sinL; //setting both to -sqrt(0.5) works just as well
	}
	else //division by qq is safe
	{
		double p = m12 * m22 + m21 * m11;
		double r = 0.5 * sqrt(1.0 / (1.0 + (4.0 * p * p / qq)));
		if (p * q > 0)//if p and q have same sign
		{
			//make 0 <= cosL <= sinL
			sinL = sqrt(0.5 + r);
			cosL = sqrt(0.5 - r); //would also work if we multiplied both by -1, or if we swapped them and multiplied one by -1
		}
		else //if p and q have opposite signs
		{
			//make 0 <= sinL <= cosL
			sinL = sqrt(0.5 - r);
			cosL = sqrt(0.5 + r); //would also work if we multiplied both by -1, or if we swapped them and multiplied one by -1
		}
	}

	//Now eq 3 is true

	//Given sinL and cosL that satisfy eq3, we can find sinR and cosR that satisfy eq1 and eq2:
	cosR = m11 * cosL - m21 * sinL;
	sinR = m12 * cosL - m22 * sinL;
	double denominator = sqrt(cosR * cosR + sinR * sinR);
	if (denominator > 0.0000001) //division by denominator is safe
	{
		cosR = cosR / denominator;
		sinR = sinR / denominator; //these values satisfy eq1, and plugging them into eq2 then multiplying by denominator results eq3
	}
	else //in this case, m11 cosL = m21 sinL and m12 cosL = m22 sinL
	{//eq1 always true regardless of sinR and cosR, so satisfy eq2:
		cosR = -m12 * sinL - m22 * cosL;
		sinR = m11 * sinL + m21 * cosL;
		denominator = sqrt(cosR * cosR + sinR * sinR);
		if (denominator > 0.0000001) //division by denominator is safe
		{
			cosR = cosR / denominator;
			sinR = sinR / denominator; //these values satisfy eq2
		}
		else //in this case, m11 sinL = -m21 cosL and m12 sinL = -m22 cosL
		{//both eq1 and eq2 always true regardless of sinR and cosR
			cosR = 1.0;
			sinR = 0.0;
		}
	}
}


//swap rows i and j of M
void K_matrix::swapRows(std::vector<double> & M, int numColumns, int i, int j)
{
	double Mik;
	int ik = i * numColumns; //such that M[ik]=M[i,k] in next for loop
	int jk = j * numColumns; //such that M[jk]=M[j,k] in next for loop
	for (int k = 0; k < numColumns; k++) //column k
	{
		Mik = M[ik];
		M[ik++] = M[jk];
		M[jk++] = Mik;
	}
}


//swap columns i and j of M
void K_matrix::swapColumns(std::vector<double> & M, int numRows, int numColumns, int i, int j)
{
	double Mki;
	int ki = i; //such that M[ki]=M[k,i] in next for loop
	int kj = j; //such that M[kj]=M[k,j] in next for loop
	for (int k = 0; k < numRows; k++) //row k
	{
		Mki = M[ki];
		M[ki] = M[kj];
		M[kj] = Mki;
		ki += numColumns;
		kj += numColumns;
	}
}


//output x such that Ax = y
//A n by n, x and y both length n
//returns true if successful, false otherwise
bool K_matrix::uppTriSolution(std::vector<double> & A, std::vector<double> & x, std::vector<double> & y, int n)
{
	for (int i = n - 1; i >= 0; i--) // solve x[n-1] then x[n-2] ... then x[0]
	{
		x[i] = y[i]; //initial value of xi
		for (int j = i + 1; j < n; j++) //sum over all Aij where i<j
			x[i] -= A[i * n + j] * x[j]; //subtract xj Aij
		double Aii = A[i * n + i];
		if (abs(Aii) > 0.00000001) //division is safe
			x[i] /= Aii; //divide by Aii
		else //division impossible
			return false;
	}
	return true;
}


//output inverse V of n by n upper triangular matrix U
//output also upper triangular
//returns true if successful, false otherwise
bool K_matrix::uppTriInv(std::vector<double> & U, std::vector<double> & V, int n)
{
	for (int i = 0; i < n; i++) //row i
	{
		for (int j = 0; j < n; j++) //column j>=i
		{
			int ij = i * n + j; //such that U[i,j] = U[ij]
			if (i > j) //left of diagonal
				V[ij] = 0.0;
			else if (i == j) //on diagonal
			{
				if (abs(U[ij]) < 0.00000001)
					return false;
				V[ij] = 1.0 / U[ij];
			}
			else //i<j, right of diagonal
			{
				//by upper triangularity, VUij = (sum from k=i to j of Vik*Ukj) = 0
				V[ij] = 0.0;
				int ik = i * n + i; //want V[i,k] = V[ik], starting at k=i
				int kj = i * n + j; //want U[k,j] = U[kj], starting at k=i
				for (int k = i; k < j; k++) //for k=i to j-1
				{
					V[ij] -= V[ik++] * U[kj]; //subtract Uik*Vik
					kj += n;
				}
				if (abs(U[kj]) < 0.00000001) //if Ujj zero
					return false;
				V[ij] /= U[kj]; //divide by Ujj
			}
		}
	}
	return true;
}


//multiply 3 matrices
std::vector<double> K_matrix::mult3matrices(std::vector<double> & A, std::vector<double> & B, std::vector<double> & C, int rowsA, int rowsB, int rowsC, int colsC)
{
	// A times B
	std::vector<double> AB;
	AB = K_matrix::mult_ignoreZeros(A, B, rowsA, rowsB, rowsC, 0.0);
	// AB times C
	std::vector<double> ABC;
	ABC = K_matrix::mult_ignoreZeros(AB, C, rowsA, rowsC, colsC, 0.0);
	return ABC;
}

//multiply 4 matrices
std::vector<double> K_matrix::mult4matrices(std::vector<double> & A, std::vector<double> & B, std::vector<double> & C, std::vector<double> & D, int rowsA, int rowsB, int rowsC, int rowsD, int colsD)
{
	// A times B times C
	std::vector<double> ABC;
	ABC = K_matrix::mult3matrices(A, B, C, rowsA, rowsB, rowsC, rowsD);
	// ABC times D
	std::vector<double> ABCD;
	ABCD = K_matrix::mult_ignoreZeros(ABC, D, rowsA, rowsD, colsD, 0.0);
	return ABCD;
}