//The purpose of this C++ file is to imitate the NAG function e02daf() from https://www.nag.com/numeric/fl/nagdoc_fl26/pdf/e02/e02daf.pdf
//Matrices are stored by row, i.e. the entry in row i and column j is in index [i*numColumns + j]
/*also see
https://www.nag.com/numeric/fl/nagdoc_fl26/html/e02/e02dcf.html
https://www.nag.com/numeric/fl/nagdoc_fl26/html/e02/e02ddf.html
https://www.nag.com/numeric/fl/nagdoc_fl26/html/e02/e02def.html
https://www.nag.com/numeric/fl/nagdoc_fl26/html/e02/e02dff.html */

#include "K_e02bf.h"
#include "K_e02daf.h"
#include "K_matrix.h"
#include "K_minres.h"

#include "KUtils.h"
#include <algorithm>
#include <vector>


/*
User inputs m points (x0,y0,f0), (x1,y1,f1) ... (x[m-1],y[m-1],f[m-1]) with respective weights w[0] ... w[m-1]
User inputs px knots xmin <= xmin <= xmin <= xmin <= lambda[4] <= lambda[5] <= ... <= lambda[px-5] <= xmax <= xmax <= xmax <= xmax
User inputs py knots ymin <= ymin <= ymin <= ymin <= mu[4] <= mu[5] <= ... <= mu[py-5] <= ymax <= ymax <= ymax <= ymax

Output of the function are coefficients in (px-4) by (py-4) matrix c:
   c[0,0]    c[0,1]    ... c[0,py-5]
   c[1,0]    c[1,1]        c[1,py-5]
   :                       :
   c[px-5,0] c[px-5,1] ... c[px-5,py-5]       where c is stored as a vector so that c[i,j] = c[i*(py-4) + j]

Define s(x,y) = sum_ij(cij * Mix * Njy) from i=0 to px-5 and j=0 to py-5
where Mix = B(i,4,x,lambda), the cubic spline from lambda[i] to lambda[i+4]
	  Njy = B(j,4,y,mu),     the cubic spline from mu[j] to mu[j+4]

Use least squares method to find c that makes f very similar to s(x,y)
Define Sigma = sum_r(rho(r)^2) from r=0 to m-1 where rho(r) = wr * (s(xr,yr) - fr)
Want dSigma/dc[a,b] = 0 for each coefficient c[a,b]:
This boils down to sum_r(wr^2 Ma(xr) Nb(yr) (sum_ij[cij Mi(xr) Nj(yr)] - fr)) = 0
which results in GC=F below, as described by P. Dierckx in https://www.sciencedirect.com/science/article/pii/0771050X77900079

| w1 M1(x1) N1(y1)      w1 M1(x1) N2(y1)      w1 M1(x1) N3(y1)  ...   w1 Mi(x1) Nj(y1)   ...     w1 Mqx(x1) Nqy(y1) | |  c11   |     |w1 f1|
| w2 M1(x2) N1(y2)      w2 M1(x2) N2(y2)      w2 M1(x2) N3(y2)  ...   w2 Mi(x2) Nj(y2)   ...     w2 Mqx(x2) Nqy(y2) | |  c12   |     |w2 f2|
| w3 M1(x3) N1(y3)      w3 M1(x3) N2(y3)      w3 M1(x3) N3(y3)  ...   w3 Mi(x3) Nj(y3)   ...     w3 Mqx(x3) Nqy(y3) | |  c13   |     |w3 f3|
|                                                                                                                   | |   :    |  =  |  :  |
| wm M1(xm) N1(ym)      wm M1(xm) N2(ym)      wm M1(xm) N3(ym)  ...   wm Mi(xm) Nj(ym)   ...     wm Mqx(xm) Nqy(ym) | |c[qx,qy]|     |wm fm| where qx=px-4 and qy=py-4

Therefore the coefficients can be found by solving TGC=TF where T is the transpose of G
*/
std::vector<double> K_e02daf::leastSquares(int m, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<double>& w,
	int px, std::vector<double>& lambda, int py, std::vector<double>& mu, bool& possible, double* DL, double eps)
{
	//get span of f
	double fMin = f[0];
	double fMax = f[0];
	for (int i = 1; i < m; i++)
	{
		if (f[i] < fMin)
			fMin = f[i];
		if (f[i] > fMax)
			fMax = f[i];
	}
	double fSpan = fMax - fMin;

	//create vectors a and b of length m such that
	//     lambda[a[r]] <= x[r] < lambda[a[r+1]]
	//         mu[b[r]] <= y[r] <     mu[b[r+1]]
	//This means that B(i,k,x[r],lambda) will only be positive when a[r] <= i <= a[r]+k
	std::vector<int> a(m, 0);
	std::vector<int> b(m, 0);
	for (int r = 0; r < m; r++)
	{
		while (lambda[a[r]] <= x[r])
		{
			a[r]++;
			if (a[r] >= px)
				break;
		}
		a[r]--;
		while (mu[b[r]] <= y[r])
		{
			b[r]++;
			if (b[r] >= py)
				break;
		}
		b[r]--;
	}


	int qx = px - 4;
	//create qx by m matrix M such that M[i,r] = B(i,4,x[r],lambda)
	std::vector<double> M(qx * m, 0.0);
	for (int i = 0; i < qx; i++)//row i
	{
		for (int r = 0; r < m; r++)
			M[i * m + r] = K_e02bf::B(i, 4, x[r], lambda); //this means M[i,r] * M[j,r]=0 when |i-j|>3
	}


	int qy = py - 4;
	//create qy by m matrix N such that N[i,r] = B(i,4,y[r],mu)
	std::vector<double> N(qy * m, 0.0);
	for (int i = 0; i < qy; i++)//row i
	{
		for (int r = 0; r < m; r++)
			N[i * m + r] = K_e02bf::B(i, 4, y[r], mu); //this means N[i,r] * N[j,r]=0 when |i-j|>3
	}


	int q = qx * qy; //total number of coefficients
	//Let G be m by q matrix in function description
	std::vector<double> G(m * q, 0.0);
	//Let F be vector of length m in function description
	std::vector<double> F(m, 0.0);


	//Fill G and F
	//for each row r and each column ij, as per function description
	//row r
	for (int r = 0; r < m; r++)
	{
		//column ij of G

		//pick out columns i (i1 <= i <= i2) such that B(i,4,x[r],lambda) nonzero
		int i1 = std::max(0, a[r] - 4);
		int i2 = std::min(qx-1, a[r]);
		for (int i = i1; i <= i2; i++)
		{
			//pick out columns j (j1 <= j <= j2) such that B(j,4,y[r],mu) nonzero
			int j1 = std::max(0, b[r] - 4);
			int j2 = std::min(qy - 1, b[r]);
			for (int j = j1; j <= j2; j++) //could go from j=b-3 to b+3 instead for efficiency, to avoid all zero Nb*Nj
			{
				int ij = i * qy + j; //column number of G
				G[r * q + ij] += w[r] * M[i * m + r] * N[j * m + r]; //wr M[i,r] N[j,r]
			}
		}
		F[r] += w[r] * f[r];
	}

	//GC=F where G is m by q, C is a column vector of length q, and F is a column vector of length m
	//Let T = transpose(G), q by m
	//Then TGC=TF where TG is q by q and TF is a column vector of length q
	std::vector<double> T;
	T = K_matrix::transpose(G, m, q);

	std::vector<double> TG;
	TG = K_matrix::mult_ignoreZeros(T, G, q, m, q, 0.0); //q by q

	std::vector<double> TF;
	TF = K_matrix::mult_ignoreZeros(T, F, q, m, 1, 0.0); //length q

	//instead of gmres, another way to solve would be invSymmetric(TG, q)
	std::vector<double> C = TF;
	possible = K_minres::gmres(TG, TF, C, q, q, 0.00000000000001 * fSpan); //try max_iterations=q instead of 200. This worked for e02baf!

	//fill DL
	double meanww = 0.0; //mean squared weight
	for (int r = 0; r < m; r++)
		meanww += w[r] * w[r];
	meanww /= m;
	K_e02daf::reduce(TG, q, q, eps*meanww);
	for (int i = 0; i < q; i++)
	{
		DL[i] = TG[i * q + i];
		if (abs(DL[i]) < eps)
			DL[i] = 0.0;
	}

	return C;
}


/*
User inputs px knots xmin <= xmin <= xmin <= xmin <= lambda[4] <= lambda[5] <= ... <= lambda[px-5] <= xmax <= xmax <= xmax <= xmax
User inputs py knots ymin <= ymin <= ymin <= ymin <= mu[4] <= mu[5] <= ... <= mu[py-5] <= ymax <= ymax <= ymax <= ymax
User inputs coefficients in (px-4) by (py-4) matrix c:
   c[0,0]    c[1,0]    ... c[0,py-5]
   c[1,0]    c[1,1]        c[1,py-5]
   :                       :
   c[px-5,0] c[px-5,1] ... c[px-5,py-5]
   where c is stored as a vector so that c[i,j] = c[i*(py-4) + j]

Output s(x,y) = sum_ij(cij * Mix * Njy) from i=0 to px-5 and j=0 to py-5
where Mix = B(i,4,x,lambda), the cubic spline from lambda[i] to lambda[i+4]
	  Njy = B(j,4,y,mu),     the cubic spline from mu[i] to mu[i+4]
*/
double K_e02daf::s(double x, double y, int px, std::vector<double> & lambda, int py, std::vector<double> & mu, std::vector<double> & c)
{
	double result = 0.0;
	for (int i = 0; i < px - 4; i++)
	{
		for (int j = 0; j < py - 4; j++)
			result += (K_e02bf::B(i, 4, x, lambda) * K_e02bf::B(j, 4, y, mu)) * c[i * (py - 4) + j]; //this could be far far more efficient!
	}
	return result;
}


/*
px = surface.m_nuk
py = surface.m_nvk
lambda = surface.m_ukv
mu = surface.m_vk
*/
void K_e02daf::daf(int m, int px, int py, double* x, double* y, double* f, double* w, double* lambda, double* mu, const int* point [[maybe_unused]] , int npoint [[maybe_unused]] , double* dl,
	double* c, int nc [[maybe_unused]] , const double* ws [[maybe_unused]] , int nws [[maybe_unused]] , double eps , const double* sigma [[maybe_unused]] , const int* rank [[maybe_unused]] , int* ifail)
{
	//create inputs for leastSquares
	std::vector<double> xVec(x, x + m);
	std::vector<double> yVec(y, y + m);
	std::vector<double> fVec(f, f + m);
	std::vector<double> wVec(w, w + m);
	std::vector<double> lambdaVec(lambda, lambda + px); //assume this is correct?
	std::vector<double> muVec(mu, mu + py); //assume this is correct?
	bool isPossible;
	//get cVec from leastSquares, and fill dl
	std::vector<double> cVec = K_e02daf::leastSquares(m, xVec, yVec, fVec, wVec, px, lambdaVec, py, muVec, isPossible, dl, eps);
	//fill c from cVec
	int qx = px - 4;
	int qy = py - 4;
	for (int i = 0; i < qx * qy; i++)
		c[i] = cVec[i]; //hopefully it's as straightforward as this
	if (isPossible)
		*ifail = 0;
	else
		*ifail = 1;

	//could calculate other values, but they are not used
}



//Function that reduces M to reduced row echelon form
//copied from https://en.wikipedia.org/wiki/Row_echelon_form#Pseudocode_for_reduced_row_echelon_form
void K_e02daf::RREF(std::vector<double>& M, int rowCount, int columnCount)
{
	int lead = 0;
	for (int r = 0; r < rowCount; r++) //row r
	{
		if (columnCount <= lead)
			return;
		int i = r;
		while (KUtils::AreDoublesEqual(M[i * columnCount + lead], 0)) //while M[i, lead] = 0
		{
			i++;
			if (rowCount == i)
			{
				i = r;
				lead++;
				if (columnCount == lead)
					return;
			}
		}
		if (i != r)
		{
			//Swap rows i and r
			double dud;
			for (int j = 0; j < columnCount; j++)
			{
				dud = M[r * columnCount + j]; //dud = Mrj
				M[r * columnCount + j] = M[i * columnCount + j]; //Mrj = Mij
				M[i * columnCount + j] = dud; //Mij = dud
			}
		}
		//Divide row r by M[r, lead]
		double Mrlead = M[r * columnCount + lead];
		for (int j = 0; j < columnCount; j++)
			M[r * columnCount + j] /= Mrlead; //divide Mrj by Mrlead

		for (i = 0; i < rowCount; i++) // same i as above?
		{
			if (i != r)
			{
				//Subtract M[i, lead] multiplied by row r from row i
				double Milead = M[i * columnCount + lead];
				for (int j = 0; j < columnCount; j++)
					M[i * columnCount + j] -= Milead * M[r * columnCount + j]; //Mij -= Milead * Mrj
			}
		}
		lead++;
	}
}

//function that performs Gaussian elimination to make M upper-triangular
void K_e02daf::reduce(std::vector<double>& M, int numRows, int numCols, double eps)
{
	for (int j = 0; j < numRows; j++) //column j
	{
		//find least row i for which Mij nonzero and i>=j
		int i = j;
		while (abs(M[i * numCols + j]) < sqrt(eps)) //while Mij zero
		{
			i++;
			if (i == numRows)
				break;
		}

		//if all zero, move onto next row
		if (i == numRows)
			continue;

		//swap rows j and i to ensure Mjj nonzero
		if (i != j)
		{
			double dud;
			for (int k = 0; k < numCols; k++)
			{
				dud = M[i * numCols + k]; //dud = Mik
				M[i * numCols + k] = M[j * numCols + k]; //Mik = Mjk
				M[j * numCols + k] = dud; //Mjk = dud
			}
		}

		//ensure all Mij below Mjj are zero
		for (i = j + 1; i < numRows; i++)
		{
			double mult = M[i * numCols + j] / M[j * numCols + j]; // Mij/Mjj
			for (int k = 0; k < numCols; k++)
				M[i * numCols + k] -= mult * M[j * numCols + k]; //decrease Mik by Mjk*Mij/Mjj
		}
	}
}
