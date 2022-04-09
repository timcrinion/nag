//The purpose of this C++ file is to imitate the NAG function g02daf() from https://www.nag.com/numeric/fl/nagdoc_fl26.0/html/g02/g02daf.html

#include "K_g02daf.h"
#include "K_matrix.h"
#include <vector>
#include <algorithm>
#include "KUtils.h"


/*
Function that fills b such that y = Xb + epsilon where
	y column vector, length n, specified by user
	X matrix, n by p, specified by user
	b column vector, length p
	epsilon column vector, length n
	sum of squares of epsilon is minimised
method is described in https://www.nag.com/numeric/fl/nagdoc_fl26.0/html/g02/g02daf.html
We read X by row, ie X[i,j] = X[i*p + j] = jth variable during ith observation
*/
void K_g02daf::multLinReg(int n, std::vector<double> X, int p, std::vector<double> Y, std::vector<double>& B, double* se_array, double tol, int* ifail)
{
	*ifail = 0;

	//get QR decomposition of X, X = Q * RStar
	std::vector<double> RStar(n * p); //n by p, same dimensions as X, but upper triangular (bottom n-p rows zero)
	std::vector<double> Q(n * n); //n by n, square
	K_matrix::QR(X, n, p, Q, RStar);

	//shave off bottom rows of RStar to get p by p matrix R
	std::vector<double> R(p * p); //p by p
	for (int i = 0; i < p * p; i++)
		R[i] = RStar[i];

	//get rank of R
	//QR decomposition means Q is full rank, so rank(X) = rank(R) = number of nonzero diagonals of R
	int rank = 0;
	for (int i = 0; i < p; i++)
	{
		if (!KUtils::AreDoublesEqual(R[i * p + i], 0.0)) //if Rii nonzero
			rank++;
	}

	//create c1, length p
	std::vector<double> QT; //transpose of Q, n by n
	QT = K_matrix::transpose(Q, n, n);
	std::vector<double> c; // c = QT * Y, length n
	c = K_matrix::mult_ignoreZeros(QT, Y, n, n, 1, 0.0);
	std::vector<double> c1(p); //first p elements of c
	for (int i = 0; i < p; i++)
		c1[i] = c[i];


	//if R full rank
	if (rank == p)
	{
		//Calculate B
		if (!K_matrix::uppTriSolution(R, B, c1, p)) //fill B such that RB = c1
			* ifail = 1;

		if (n == p)
			* ifail = 5;
		else
		{//Get standard errors and put them in se_array
			double variance = K_g02daf::getVariance(Y, X, B, n, p, rank, ifail);
			//want diagonals of inverse(X^T * X)
			//Since X=QR, and inv(Q)=Q^T, inv(X^T X) = inv(RStar^T RStar) = inv(R^T R) = inv(R) inv(R)^T
			std::vector<double> invR(p * p); //inverse of R, p by p, upper triangular
			if (!K_matrix::uppTriInv(R, invR, p))
				* ifail = 1;
			std::vector<double> invRT; //transpose of inverse of R, p by p
			invRT = K_matrix::transpose(invR, p, p);
			std::vector<double> invRTR; //inverse of (R^T R)
			invRTR = K_matrix::mult_ignoreZeros(invR, invRT, p, p, p, 0.0);
			//get standard errors
			for (int i = 0; i < p; i++)
				se_array[i] = sqrt(variance * invRTR[i * p + i]);
		}
	}
	else //rank < p
	{
		//Calculate three p by p matrices such that R = QStar * DStar * PT where DStar diagonal with rank(R) nonzero elements
		std::vector<double> QStar(p * p); //p by p
		std::vector<double> DStar(p * p); //p by p
		std::vector<double> PT(p * p); //p by p
		K_matrix::SVD(R, p, QStar, DStar, PT, tol);

		//create vector (length p) of diagonal magnitudes of DStar
		std::vector<double> diagonals(p);
		for (int i = 0; i < p; i++)
			diagonals[i] = abs(DStar[p * i + i]);

		//create place such that diagonals[place[0]] <= diagonals[place[1]] <= ... <= diagonals[place[p-1]]
		std::vector<int> place(p);
		for (int i = 0; i < p; i++)
			place[i] = i;
		std::sort(std::begin(place), std::end(place), [&](int i1, int i2) { return diagonals[i1] < diagonals[i2]; });

		//create numSmaller such that diagonals[i] has numSmaller[i] diagonal magnitudes smaller than itself
		std::vector<int> numSmaller(p);
		for (int i = 0; i < p; i++)
			numSmaller[place[i]] = i;

		//define diagonal[i] as "zero" if numSmaller[i] < p-rank, "nonzero" if numSmaller[i] >= p-rank
		//want to permute DStar such that first rank values are "nonzero" and last p-rank are "zero"
		//set i=0...rank-1 and j=rank to p-1. Each time Dii is zero and Djj is nonzero, swap them
		int i = 0;
		int j = rank;
		while (i < rank && j < p)
		{
			//increase i until we hit a zero
			while (numSmaller[i] >= p - rank) //while diagonal[i] nonzero
			{
				i++;
				if (i == rank)
					break;
			}
			//increase j until we hit a nonzero
			while (numSmaller[j] < p - rank) //while diagonal[j] zero
			{
				j++;
				if (j == p)
					break;
			}
			//if i or j too big, break
			if (i >= rank || j >= p)
				break;
			//swap DStar_ii and DStar_jj
			K_matrix::swapColumns(DStar, p, p, i, j);
			K_matrix::swapRows(PT, p, i, j); //these two operations do not change DStar*PT
			K_matrix::swapColumns(QStar, p, p, i, j);
			K_matrix::swapRows(DStar, p, i, j); //these two operations do not change QStar*DStar
			//next i and next j
			i++;
			j++;
		} //Now the top-left rank diagonals of DStar should be nonzero

		//trim
		int k = rank;
		//QStar1 is first k columns of QStar
		K_matrix::transposeSquareMatrix(QStar, p); //now is k by p, simply ignore bottom rows later
		//P1 is the first k columns of P
		std::vector<double> P1; //p by k
		P1 = K_matrix::transpose(PT, k, p); //transpose top k rows of PT
		//D is the nonzero diagonals of DStar
		std::vector<double> Dinv(k * k, 0.0); //k by k
		for (int a = 0; a < k; a++)
			Dinv[a * k + a] = 1.0 / DStar[a * p + a];

		//set B = P1 * Dinv * QStar1^T * c1
		std::vector<double> PDQ; //p by p
		PDQ = K_matrix::mult3matrices(P1, Dinv, QStar, p, k, k, p); //we ignore bottom p-k rows of QStar
		std::vector<double> PDQc; //length p
		PDQc = K_matrix::mult_ignoreZeros(PDQ, c1, p, p, 1, 0.0);
		for (int a = 0; a < p; a++)
			B[a] = PDQc[a];

		if (n == p)
			* ifail = 5;
		else
		{//get standard errors and put them in se_array
			double variance = K_g02daf::getVariance(Y, X, B, n, p, rank, ifail);
			std::vector<double> PDQT; //transpose of PDQ, p by p
			PDQT = K_matrix::transpose(PDQ, p, p);
			std::vector<double> invRTR; //p by p, plays equivalent role of invRTR in rank==p case
			invRTR = K_matrix::mult_ignoreZeros(PDQ, PDQT, p, p, p, 0.0); //PDQ times PDQT
			//get standard errors
			for (int a = 0; a < p; a++)
				se_array[a] = sqrt(variance * invRTR[a * p + a]);
		}
	}

	//standard errors were calculated using
	//https://online.stat.psu.edu/stat501/lesson/5/5.3#:~:text=S%3D%E2%88%9AMSE,the%20simple%20linear%20regression%20setting
	//http://users.stat.umn.edu/~helwig/notes/mlr-Notes.pdf //this states that se_array proportional to diagonal of inv(XT*X)
}

//calculates variance for multLinReg
double K_g02daf::getVariance(std::vector<double> Y, std::vector<double> X, std::vector<double> B, int n, int p, int rank, int* ifail)
{
	//epsilon = XB - Y
	std::vector<double> epsilon; //length p
	epsilon = K_matrix::mult_ignoreZeros(X, B, n, p, 1, 0.0);
	for (int i = 0; i < n; i++)
		epsilon[i] -= Y[i];
	//variance
	double variance = 0.0;
	for (int i = 0; i < n; i++)
		variance += epsilon[i] * epsilon[i];
	double denominator = abs(static_cast<double>(n - rank)); //unsure why, but this is mentioned in the two links here:
	//https://online.stat.psu.edu/stat501/lesson/5/5.3#:~:text=S%3D%E2%88%9AMSE,the%20simple%20linear%20regression%20setting
	//http://users.stat.umn.edu/~helwig/notes/mlr-Notes.pdf
	if (KUtils::AreDoublesEqual(denominator, 0.0)) //if denominator zero
		* ifail = 5;
	else //division safe
		variance = variance / denominator;
	return variance;
}

//mimics NAG function g02daf
//n = number of observations
//m = number of independent variables
//X_array n by m
//ip = length of B. This will be the number of nonzero values in isx (plus 1 if mean=true). Cannot be more than n.
//http://mezeylab.cb.bscb.cornell.edu/labmembers/documents/supplement%205%20-%20multiple%20regression.pdf
//X_array read by column, ie X_array[i + j*n] = jth variable during ith observation
//setting mean to true is equivalent to 'M' in original NAG, false is equivalent to 'Z'
void K_g02daf::g02daf(bool mean, int n, const double* X_array, int m, const int* isx, int ip, double* Y_array, double* B_array, double* se_array, double tol, int* ifail)
{
	//check values ok
	if (n < 2 || m < 1 || tol<0.0 || ip <= 0 || ip > n)
	{
		*ifail = 1;
		return;
	}

	//check that ip (length of B_array) is correct
	int isxCount = 0;
	for (int i = 0; i < m; i++)
	{
		if (isx[i] > 0)
			isxCount++;
	}
	if (mean)
		isxCount++;
	if (ip != isxCount)
	{
		*ifail = 4;
		return;
	}


	//create new matrix X, dimensions n by ip, including only the columns of X_array desired by the user
	//unlike X_array, want X read by row, ie X[i*ip + j] = jth variable during ith observation
	std::vector<double> X(n * ip);
	int iX = 0; //next index of X to be filled
	for (int i = 0; i < n; i++) //row i
	{
		if (mean) //user wants y = b0 + sum_i(bi * xi)
			X[iX++] = 1.0;//leftmost column all 1s
		for (int j = 0; j < m; j++) //column j of X_array
		{
			if (isx[j] > 0) // user wants jth column of X_array to be included in X
				X[iX++] = X_array[i + j * n]; //observation i, variable j
		}
	}

	//vector versions of B and Y
	std::vector<double> Y(Y_array, Y_array + n); // length n
	std::vector<double> B(B_array, B_array + ip); // length ip

	//fill B and se_array such that Y = X*B + errors
	K_g02daf::multLinReg(n, X, ip, Y, B, se_array, tol, ifail);

	//fill B_array
	for (int i = 0; i < ip; i++)
		B_array[i] = B[i];
}