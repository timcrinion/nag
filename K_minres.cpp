//C++ implementation of MINRES method: https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Example_code
#include <cmath>
#include "K_minres.h"
#include "K_matrix.h"
#include "K_g05rzf.h"
#include "KUtils.h"


//A n by n matrix
//b vector of length n
//Want to solve Ax=b for x

//function[x, e] = gmres(A, b, x, max_iterations, threshold)
bool K_minres::gmres(std::vector<double>& A, std::vector<double>& b, std::vector<double>& x, int n, int max_iterations, double threshold)
{
	bool possible = true;

	//n = length(A)
	const int m = max_iterations;

	//r = K_matrix::MatrixMult(A, n, n, x, 1); //r=Ax
	std::vector<double> r = K_matrix::mult(A, n, n, x, 1);
	for (int i = 0; i < n; i++) //r=b-Ax
		r[i] = b[i] - r[i];

	//b_norm = norm(b)
	const double b_norm = K_matrix::FrobeniusNorm(b, n); //could return x as zero vector if b_norm very small?

	//initialize the 1D vectors
	//sn = zeros(m, 1)
	std::vector<double> sn(m, 0.0);
	//cs = zeros(m, 1)
	std::vector<double> cs(m, 0.0);
	//r_norm = norm(r)
	const double r_norm = K_matrix::FrobeniusNorm(r, n);
	//Q(:, 1) = r / r_norm
	std::vector<double> Q(n * (m + 1), 0.0); //Q n by m+1
	if (abs(r_norm) <= 0.0000000000001) //unlikely but do anyway to avoid dividing by zero
		return true;
	const double inv_r_norm = 1.0 / r_norm;
	for (int i = 0; i < n; i++)
		Q[i * (m + 1)] = r[i] * inv_r_norm; //set Q[i,0]
	//beta = r_norm * e1
	std::vector<double> beta(m + 1, 0.0); //length m+1
	beta[0] = r_norm;

	//for k = 1 : m
	std::vector<double> H((m + 1) * m, 0.0); //H m+1 by m
	int k; //want k to be available after for-loop
	for (k = 0; k < m; k++)
	{
		//run arnoldi
		//[H(1:k + 1, k) Q(:, k + 1)] = arnoldi(A, Q, k)
		arnoldi(H, A, Q, k, n, m);

		//eliminate the last element in H ith row and update the rotation matrix
		//[H(1:k + 1, k) cs(k) sn(k)] = apply_givens_rotation(H(1:k + 1, k), cs, sn, k)
		apply_givens_rotation(H, cs, sn, k, m);

		//update the residual vector
		beta[k + 1] = -sn[k] * beta[k];
		beta[k] = cs[k] * beta[k];
		const double error = abs(beta[k + 1]) / b_norm;

		if (error <= threshold)
		{
			k++; //since breaking, undo the k-- a few lines below
			break;
		}
	}

	k--; //number of iterations

	//calculate the result
	//y = H(1:k, 1 : k) \ beta(1:k)
	//find y such that H(1:k, 1 : k) * y = beta(1:k). H upper diagonal, so can work out yk...y2,y1,y0
	std::vector<double> y(k + 1, 0.0);
	for (int i = k; i >= 0; i--) //i = k to 0
	{
		const int i0 = i * m; //such that H[i0] = H[i,0]
		double sum = 0.0;
		for (int j = i + 1; j <= k; j++)
			sum += H[i0 + j] * y[j]; // Hij * yj
		if (abs(H[i0 + i]) < 0.0000000001) //if Hii close to 0
			possible = false;
		else //if Hii far from 0
			y[i] = (beta[i] - sum) / H[i0 + i]; // (beta[i] - sum) / Hii
	}

	//x = x + Q(:, 1 : k) * y
	for (int i = 0; i < n; i++) //Q m by n
	{
		const int i0 = i * (m + 1); //such that Q[i0] = Q[i,0]
		for (int j = 0; j <= k; j++)
			x[i] += Q[i0 + j] * y[j]; // Q[i, j] * y[j]
	}

	return possible;
} //end


//----------------------------------------------------
//                  Arnoldi Function
//----------------------------------------------------
//function[h, q] = arnoldi(A, Q, k)
//called using [H(1:k + 1, k) Q(:, k + 1)] = arnoldi(A, Q, k)
//Q is m by n, A is n by n, H m+1 by m
void K_minres::arnoldi(std::vector<double> & H, const std::vector<double> & A, std::vector<double> & Q, int k, int n, int m)
{
	const std::vector<double> Qk = getCol(Q, m + 1, n - 1, k);
	std::vector<double> q = K_g05rzf::MatrixByVector(A, Qk, n, n); //length n

	//for i = 1:k
	std::vector<double> h(k + 2, 0.0); //length k+2, will be imported into H[0:k+1,k]
	for (int i = 0; i <= k; i++)
	{
		//h(i) = q' * Q(:, i)
		for (int j = 0; j < n; j++)
			h[i] += q[j] * Q[j * (m + 1) + i]; //q[j] * Q[j,i]

		//q = q - h(i) * Q(:, i)
		for (int j = 0; j < n; j++)
			q[j] -= h[i] * Q[j * (m + 1) + i]; // - h[i] * Q[j,i]
	}
	//h(k + 1) = norm(q)
	h[k + 1] = K_matrix::FrobeniusNorm(q, n);
	//q = q / h(k + 1)
	for (int i = 0; i < n; i++)
		q[i] /= h[k + 1];
	//called using [H(1:k + 1, k) Q(:, k + 1)] = arnoldi(A, Q, k), therefore update H and Q:
	for (int i = 0; i <= k + 1; i++)
		H[i * m + k] = h[i]; //H[i,k] = h[i]
	for (int i = 0; i < n; i++)
		Q[i * (m + 1) + k + 1] = q[i]; //Q[i,k+1] = q[i]
}





//---------------------------------------------------------------------
//                  Applying Givens Rotation to H col
//---------------------------------------------------------------------
//function[h, cs_k, sn_k] = apply_givens_rotation(h, cs, sn, k)
//called using [H(1:k + 1, k) cs(k) sn(k)] = apply_givens_rotation(H(1:k + 1, k), cs, sn, k)
//cs and sn length m, H n+1 by n
void K_minres::apply_givens_rotation(std::vector<double> & H, std::vector<double> & cs, std::vector<double> & sn, int k, int m)
{
	std::vector<double> h = getCol(H, m, k + 1, k);

	//apply for ith column
	//for i = 1:k - 1
	for (int i = 0; i < k; i++)
	{
		const double temp = cs[i] * h[i] + sn[i] * h[i + 1];
		h[i + 1] = -sn[i] * h[i] + cs[i] * h[i + 1];
		h[i] = temp;
	}

	//update the next sin cos values for rotation
	//[cs_k sn_k] = givens_rotation(h(k), h(k + 1))
	double cs_k, sn_k;
	givens_rotation(&cs_k, &sn_k, h[k], h[k + 1]);

	//eliminate H(i + 1, i)
	h[k] = cs_k * h[k] + sn_k * h[k + 1];
	h[k + 1] = 0.0;

	//set columns
	setCol(H, h, m, k + 1, k);

	//outputs
	cs[k] = cs_k;
	sn[k] = sn_k;
}


//----Calculate the Given rotation matrix----
//function[cs, sn] = givens_rotation(v1, v2)
void K_minres::givens_rotation(double* cs, double* sn, double v1, double v2)
{
	if (KUtils::AreDoublesEqual(v1, 0.0))
	{
		*cs = 0.0;
		*sn = 1.0;
	}
	else
	{
		const double t = sqrt(v1 * v1 + v2 * v2);
		*cs = abs(v1) / t;
		*sn = (*cs) * v2 / v1;
	}
}

//Returns M[0:i,j] (where 0 <= j < numColsM)
std::vector<double> K_minres::getCol(std::vector<double> & M, int numColsM, int i, int j)
{
	std::vector<double> column(i + 1, 0.0); //will be column, length i+1
	int kj = j; //such that M[kj] = M[k,j] in loop below
	for (int k = 0; k <= i; k++)
	{
		column[k] = M[kj];
		kj += numColsM;
	}
	return column;
}

//Sets M[0:i,j] to column (where 0 <= j < numCols)
//column must have length of at least i+1
void K_minres::setCol(std::vector<double> & M, std::vector<double> & column, int numColsM, int i, int j)
{
	int kj = j; //such that M[kj] = M[k,j] in loop below
	for (int k = 0; k <= i; k++)
	{
		M[kj] = column[k];
		kj += numColsM;
	}
}