//The purpose of the C++ file is to mimic nag functions E02BAF and E02BBF
#include "K_e02bf.h"
#include "KUtils.h"
#include "K_minres.h"
#include <vector>

//B-splines:


/*
Given a sequence of knots t[0] < t[1] < ... < t[n]
define B-splines iteratively as follows:

w(i,k,x) = (x-t[i]) / (t[i+k]-t[i])       that is, the line passing through (0,t[i]) and (1,t[i+k])

B(i,1,x) = 1 if t[i] <= x < t[i+1]
		   0 otherwise

for k>0:
B(i,k+1,x) = w(i,k,x)B(i,k,x) + (1-w(i+1,k,x))B(i+1,k,x)
*/

//define w as above, with knots t
double K_e02bf::w(int i, int k, double x, std::vector<double>& t)
{
	const double denominator = t[i + k] - t[i];
	if (denominator < 0.000001)
		return 0.0;
	return (x - t[i]) / denominator;
}

//define B as above, with knots t
//positive iff t[i] <= x < t[i+k], zero otherwise
double K_e02bf::B(int i, int k, double x, std::vector<double>& t)
{
	if (x < t[i])
		return 0.0;
	if (x > t[i + k])
		return 0.0;
	//if we reach this point, must have t[i] <= x < t[i+k], output nonzero
	if (k == 1)
		return 1.0;
	else // k > 1
	{
		//Bik = w(i,k-1)B(i,k-1) + (1-w(i+1,k-1))B(i+1,k-1)
		double result = 0.0;
		//only bother adding w(i,k-1)B(i,k-1) if it's nonzero
		if (x < t[i + k - 1])
			result += K_e02bf::w(i, k - 1, x, t) * K_e02bf::B(i, k - 1, x, t);
		//only bother adding (1-w(i+1,k-1))B(i+1,k-1) if it's nonzero
		if (x >= t[i + 1])
			result += (1 - K_e02bf::w(i + 1, k - 1, x, t)) * K_e02bf::B(i + 1, k - 1, x, t);
		return result;
	}
}



/*
User inputs m points (x1,y1), (x2,y2) ... (xm,ym) with respective weights w1 ... wm (must have x1 <= x2 <= ... <= xm)
User inputs q+4 knots x1 <= x1 <= x1 <= x1 <= t5 <= t6 <= ... <= t[q] <= xm <= xm <= xm <= xm
Output of the function are coefficients c1,c2...cq
Want to define coefficients c such that s(x) = sum from i=1 to q of (ci * Ni(x)) where Ni(x) = B(i,4,x,t)
Use least squares method to make y similar to s(x)
Therefore define residue R = sum from r=1 to m of wr^2 (yr - s(xr))^2
Want dR/dci = 0 for each coefficient ci:
This boils down to sum_r(w^2 c1 N1 Ni) + sum_r(w^2 c2 N2 Ni) + ... + sum_r(w^2 cm Nm Ni) = sum_r(w^2 yr Ni)
This gives the matrix equation

| sum_r(wwN1N1)  sum_r(wwN2N1) ... sum_r(wwNqN1) |  | c1 |    | sum_r(wwyN1) |
| sum_r(wwN1N2)  sum_r(wwN2N2)          :        |  | c2 |    | sum_r(wwyN2) |
|     :                                          |  | :  |  = |      :       |
| sum_r(wwN1Nq)  ...               sum_r(wwNqNq) |  | cq |    | sum_r(wwyNq) |

Therefore c1...cq can be found by multiplying the inverse of the square symmetric matrix on the LHS by the vector on the RHS
This commented explanation uses 1-based numbering, whereas the C++ code below uses 0-based numbering
*/
std::vector<double> K_e02bf::leastSquares(int m, std::vector<double> & x, std::vector<double> & y, std::vector<double> & w, int q, std::vector<double> & t, bool& possible)
{
	//get span of y
	double yMin = y[0];
	double yMax = y[0];
	for (int i = 0; i < m; i++)
	{
		if (y[i] < yMin)
			yMin = y[i];
		if (y[i] > yMax)
			yMax = y[i];
	}
	double ySpan = yMax - yMin;

	//create vector of squared weights
	std::vector<double> ww(m);
	for (int r = 0; r < m; r++)
		ww[r] = w[r] * w[r];

	//create q by m matrix N such that N[i,r] = B(i,4,x[r],t)
	std::vector<double> N(q * m, 0.0);
	for (int i = 0; i < q; i++)
	{
		//get values r1 and r2 such that N[i,r] only nonzero when r1 <= r <= r2 ?
		for (int r = 0; r < m; r++)
			N[i * m + r] = B(i, 4, x[r], t); //this means N[i,r] * N[j,r]=0 when |i-j|>3
	}

	//Let M be q by q matrix on LHS
	std::vector<double> M(q * q, 0.0);
	//for each row
	for (int i = 0; i < q; i++)
	{
		//for each nonzero column
		for (int j = 0; j < q; j++)
		//for (int j = i; j < i + 4; j++)     <-     try this line instead after we are sure function works, much more efficient
		{
			if (j < q)
			{
				//calculate M[i,j]
				int ij = i * q + j;
				for (int r = 0; r < m; r++)
					M[ij] += ww[r] * N[i * m + r] * N[j * m + r]; //add wr * wr * B(i,4,xr,t) * B(j,4,xr,t)
			}
		}
	}

	//fill lower half of M
	for (int i = 0; i < q; i++)
	{
		for (int j = 0; j < i; j++)
			M[i * q + j] = M[j * q + i];
	}

	//Let v be vector on RHS
	std::vector<double> v(q, 0.0);
	for (int i = 0; i < q; i++)
	{
		for (int r = 0; r < m; r++)
			v[i] += ww[r] * y[r] * N[i * m + r]; //add wr * wr * yr * B(i,4,xr,t)
	}

	//want M^-1 * v:
	//std::vector<double> invM;
	//invM = f07fjc::invSymmetric(M, 5); //use this if convenient
	std::vector<double> c = v;
	possible = K_minres::gmres(M, v, c, q, q, 0.00000000000001 * ySpan); //setting max_iterations=q instead of 200 causes tina test 7376 to stop filling c with nans on its 9th and last execution of this line
	return c;
}

void K_e02bf::e02baf (const int m, const int ncap7, const double x[], const double y[], const double w[], double lamda[],
	const double /*work1*/[], const double /*work2*/[], double c[], const double* /*ss*/, int* ifail)
{
	std::vector<double> x_inner (m);
	std::vector<double> y_inner (m);
	std::vector<double> w_inner (m);
	for (size_t i = 0; i < m; i++)
	{
		x_inner[i] = x[i];
		y_inner[i] = y[i];
		w_inner[i] = w[i];
	}
	std::vector<double> t_inner (ncap7);
	for (size_t i = 0; i < ncap7; i++)
		t_inner[i] = lamda[i];
	std::vector<double> c_inner = e02baf_inner (m, x_inner, y_inner, w_inner, ncap7, t_inner, *ifail);
	for (size_t i = 0; i < ncap7; i++)
	{
		lamda[i] = t_inner[i];
		c[i] = (i < ncap7 - 4) ? c_inner[i] : 0.0;
	}
}


/*
x,y,w vectors of length m
User inputs m points(x0,y0),(x1,y1)... with respective weights w0,w1... (must have x0 <= x1 <= x2...)
t vector of length ncap7-8
t consists of knots such that x[0] < t[0] <= t[1] <= ... <= t[ncap7-9] < x[m-1]
function increases size of t to ncap7 by adding eight points
Output of the function is vector c of length q=ncap7-4, such that s(x) = sum from i = 0 to q-1 of(c[i] * Ni(x)) where Ni(x) is the cubic B-spline B(i, 4, x, t)
*/
std::vector<double> K_e02bf::e02baf_inner(int m, std::vector<double> & x, std::vector<double> & y, std::vector<double> & w, int ncap7, std::vector<double> & t, int& ifail)
{
	int q = ncap7 - 4;

	// first four values t[0] to t[3] are equal to minimum of x
	// last four values t[q] to t[q+3] are equal to maximum of x

	bool possible;
	std::vector<double> c = K_e02bf::leastSquares(m, x, y, w, q, t, possible);
	ifail = 1;
	if (possible)
		ifail = 0;

	return c;
}


//mimics nag function e02bbf https://www.nag.com/numeric/fl/nagdoc_fl25/pdf/e02/e02bbf.pdf
//lambda has length ncap7
//C has length q = ncap7 - 4
void K_e02bf::e02bbf_vec(int ncap7, std::vector<double> & lambda, std::vector<double> & C, double x, double& s, int& ifail)
{
	ifail = 0;
	const int q = ncap7 - 4;

	if (ncap7 < 8)
	{
		ifail = 2;
		return;
	}

	//if x outside allowed domain, quit
	if (x<lambda[3] || x>lambda[q])
	{
		ifail = 1;
		s = 0.0;
		return;
	}

	//special cases
	if (KUtils::AreDoublesEqual(x, lambda[3]))
	{
		s = C[0];
		return;
	}
	if (KUtils::AreDoublesEqual(x, lambda[q]))
	{
		s = C[q - 1];
		return;
	}

	//find k such that lambda[k-1] <= x < lambda[k]
	int k = 4;
	while (lambda[k] < x)
		k++;

	//set s
	s = 0.0;
	//Ni(x) can only be nonzero for i = k-4, k-3, k-2, k-1
	for (int i = k - 4; i < k; i++)
		s += C[i] * K_e02bf::B(i, 4, x, lambda); //increase s by Ci*Ni(x)
}


//mimics nag function e02bbf https://www.nag.com/numeric/fl/nagdoc_fl25/pdf/e02/e02bbf.pdf
//same as function above but inputs are double pointers, not double vectors
void K_e02bf::e02bbf(int ncap7, double* lambda, double* C, double x, double& s, int& ifail)
{
	//create vector versions of lambda and C
	const int q = ncap7 - 4;
	std::vector<double> lambdaVec(lambda, lambda + ncap7);
	std::vector<double> cVec(C, C + q);
	K_e02bf::e02bbf_vec(ncap7, lambdaVec, cVec, x, s, ifail);
}
