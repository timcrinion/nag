//This C++ file is used to reduce square matrices to upper Hessenberg form
//uses results from https://www.math.kth.se/na/SF2524/matber15/qrmethod.pdf

#include "K_matrix.h"
#include "K_UpperHessenberg.h"
#include "KUtils.h"
/*
A matrix is an upper Hessenberg matrix if all entries below the subdiagonal are zero, ie

x x x x x
x x x x x
0 x x x x
0 0 x x x
0 0 0 x x

*/


//sets n by n matrix M to identity
void K_UpperHessenberg::eye(std::vector<double>& M, int n)
{
	for (int i = 0; i < n*n; i++)
		M[i] = 0.0;
	for (int i = 0; i < n; i++)
		M[i*n + i] = 1.0;
}



/*
function that fills P (n by n) such that P is householder reflector I-2*v*v^T of a normalised vector v of length n
we have P = P^T = P^-1

	 | 1-2v1v1   -2v1v2   -2v1v3          |
	 |  -2v2v1  1-2v2v2   -2v2v3          |
P =  |  -2v3v1   -2v3v2  1-2v3v3          |
	 |                            1-2vnvn |

The idea is that if |u|=1 then Px reflects x in the plane defined by the vector u
*/
void K_UpperHessenberg::householderReflector(std::vector<double>& P, std::vector<double>& v, int n)
{
	//fill upper diagonal part
	for (int i = 0; i < n; i++) //row i
	{
		for (int j = i; j < n; j++) //column j >= i
			P[i*n + j] = -2.0 * v[i] * v[j]; // Pij = -2 vi vj
		P[i*n + i] += 1.0; // Pii = 1 - 2 vi vi
	}
	//fill lower diagonal part
	for (int i = 0; i < n; i++) //row i
	{
		for (int j = 0; j < i; j++) //column j < i
			P[i*n + j] = P[j*n + i]; // Pij = Pji
	}
}


//normalise a vector v of length n
//returns true if successful (ie original norm > 0)
bool K_UpperHessenberg::normalise(std::vector<double>& v, int n)
{
	//calculate norm
	double norm = K_matrix::FrobeniusNorm(v, n);
	//normalise v
	if (KUtils::AreDoublesEqual(norm, 0.0))
		return false;
	K_UpperHessenberg::multVecByScalar(v, n, 1.0 / norm);
	return true;
}


//multiply vector v of length n by scalar r
void K_UpperHessenberg::multVecByScalar(std::vector<double>& v, int n, double r)
{
	for (int i = 0; i < n; i++)
		v[i] *= r;
}


/*
Want P such that Px = y where |x| = |y|

Let z = x - y
Let u = z / |z|

Then
	z'z = |x|^2 + |y|^2 - 2x'y

	z'x = |x|^2 - x'y

	uu'x = z (z'x) / |z|^2 = (z'x) z / (z'z) = z/2

Let P = I - 2uu'

Then
	Px = (I-2uu')x = x - 2z/2 = y
*/


//Given n by n matrix M, returns subcolumn { M[j+1,j]  M[j+2,j] ...  M[n-1,j] } below diagonal in column j, length n-j-1
std::vector<double> K_UpperHessenberg::colBelowDiag(std::vector<double>& M, int n, int j)
{
	std::vector<double> result(n - j - 1);
	int ij = (j + 1) * n + j; //want M[ij] = M[i,j] in loop below
	for (int i = 0; i < n - j - 1; i++)
	{
		result[i] = M[ij];
		ij += n;
	}
	return result;
}




//Given n by n matrix A, fills unitary Q and upper hessenberg H such that H = Q * A * Q^T
void K_UpperHessenberg::hessenbergReduction(std::vector<double>& A, int n, std::vector<double>& Q, std::vector<double>& H)
{
	//set H to A
	for (int i = 0; i < n*n; i++)
		H[i] = A[i];
	//set Q to identity
	K_UpperHessenberg::eye(Q, n);
	//prepare n by n matrix P for use in loop below
	std::vector<double> P(n*n);

	//loop through, multiplying Q and H by P each time
	for (int j = 0; j < n - 2; j++) //go through all columns j, except n-1 and n-2
	{
		//set colj to A[j+1:n-1 , j] of length n-j-1
		std::vector<double> colj;
		colj = K_UpperHessenberg::colBelowDiag(H, n, j);
		//get sign of A[j+1 , j]
		int sign = 1;
		if (colj[0] < 0)
			sign = -1;
		//increase or decrease top entry of colj by |colj|, whichever gets it further from zero
		colj[0] += sign * K_matrix::FrobeniusNorm(colj, n - j - 1);
		//scale colj so that |colj|=1
		if (!K_UpperHessenberg::normalise(colj, n - j - 1))
			continue; //skip to next column if |colj| zero
		//fill P such that P * A * P has zeros below subdiagonal in column j
		K_UpperHessenberg::fillP(P, j, n, colj);
		//update H and Q
		K_UpperHessenberg::leftRightMult(P, H, P, n); //replace H by PHP
		K_UpperHessenberg::leftMult(P, Q, n); //replace Q by PQ
	}
}






/*    fill P for K_UpperHessenberg::hessenbergReduction
		   ___________________________
		  |             |             |
		  |  identity   |    zeros    | j+1
 want P = |_____________|_____________|
		  |             |             |
		  |    zeros    |    I-2uu'   | n-j-1
		  |_____________|_____________|
				j+1          n-j-1

 where u = input vector of length n-j-1       */
void K_UpperHessenberg::fillP(std::vector<double>& P, int j, int n, std::vector<double>& u)
{
	//set P to identity
	K_UpperHessenberg::eye(P, n);
	//set bottom right quadrant
	for (int row = j + 1; row < n; row++) //row of P
	{
		for (int col = j + 1; col < n; col++) //column of P
		{
			P[row*n + col] -= 2.0 * u[row - j - 1] * u[col - j - 1]; //subtract from P[row,col]                 could be more efficient!
		}
	}
}


//given two n by n matrices A and B, replace B by A*B
void K_UpperHessenberg::leftMult(std::vector<double>& A, std::vector<double>& B, int n)
{
	std::vector<double> AB;
	AB = K_matrix::mult_ignoreZeros(A, B, n, n, n, 0.0);
	for (int i = 0; i < n*n; i++)
		B[i] = AB[i];
}

//given two n by n matrices A and B, replace A by A*B
void K_UpperHessenberg::rightMult(std::vector<double>& A, std::vector<double>& B, int n)
{
	std::vector<double> AB;
	AB = K_matrix::mult_ignoreZeros(A, B, n, n, n, 0.0);
	for (int i = 0; i < n*n; i++)
		A[i] = AB[i];
}

//given three n by n matrices A,B,C, replace B by A*B*C
void K_UpperHessenberg::leftRightMult(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, int n)
{
	std::vector<double> ABC;
	ABC = K_matrix::mult3matrices(A, B, C, n, n, n, n);
	for (int i = 0; i < n*n; i++)
		B[i] = ABC[i];
}
