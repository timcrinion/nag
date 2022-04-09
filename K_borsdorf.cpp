#include <cmath>
#include "K_minres.h"
#include "K_jacobi.h"
#include "K_borsdorf.h"
#include "K_matrix.h"

//following from http://eprints.maths.manchester.ac.uk/1085/1/thesis.pdf

//returns n by n matrix that is all zeros except for the vector v of length n along its diagonal
std::vector<double> K_borsdorf::diagEnlarge(const std::vector<double>& v, int n)
{
	std::vector<double> result(n * n, 0.0);
	for (int i = 0; i < n; i++)
		result[i * n + i] = v[i]; //result[i,i] = v[i]
	return result;
}

//returns diagonal of n by n matrix M
std::vector<double> K_borsdorf::diagReduce(const std::vector<double> & M, int n)
{
	std::vector<double> v(n, 0.0);
	for (int i = 0; i < n; i++)
		v[i] = M[i * n + i]; //result[i,i] = v[i]
	return v;
}


//All function below this line have been translated into C++ from matlab algorithm 2 written in http://eprints.maths.manchester.ac.uk/1085/1/thesis.pdf (section A.2)

//PROGRAM FOR FINDING THE NEAREST CORRELATION MATRIX
//function[X, iter, f_eval, normgrad] = cornewton(G, error_tol, flag2, diaones, maxitmeth, maxit, pre, prnt)
//G n by n matrix. Function returns G's nearest correlation matrix.
std::vector<double> K_borsdorf::cornewton(std::vector<double> & G, int n, double error_tol, int diaones, int maxitmeth, int maxit, bool* ifail)
{
	double eps = 0.000000000000001; //from matlab
	double tempDouble = 0.0;
	*ifail = false;

	/* Original matlab shown:
	% [X, iter, f_eval, normgrad] = cornewton(G, error_tol, flag2, diaones, maxitmeth, maxit, pre, prnt)
	% finds the nearest correlation matrix to the symmetric matrix G.
	% error_tol is a termination tolerance
	% (default is length(G) * 10 ^ (-7)).


	% flag2 = 0: using for full eigendecomposition MATLAB function eig
	% (LAPACK - function DSYEV).
	% flag2 = 1 : using for full eigendecomposition MATLAB - NAG Toolbox
	% function f08fa
	13 % flag2 = 2 : using for full eigendecomposition MATLAB - NAG Toolbox
	% function f08fc(default)
	% flag2 = 3 : using for full eigendecomposition MATLAB - NAG Toolbox
	% function f08fd
	% diaones = 0 : the matrix X is not changed after
	18 % the computation(default)
	% diaones = 1 : the diagonal is set to 1 after computing X
	% diaones = 2 : the diagonal is rescaled to 1 and then set to 1
	% maxitmeth : is the number of the maximal iterations in the
	% iterative method, (default is 200)
	23 % maxit : is the number of maximal iterations in the Newton - method,
	% (default is 200)
	% pre = 1 : using a preconditioner in the iterative method
	% (default is 1)
	% prnt = 1 for display of intermediate output . (default is 0)
	*/
	//Test input values
	//Test the matrix

	//Ensure symmetry
	//line 38, G = (G + G ’) / 2
	for (int i = 0; i < n; i++)
	{
		int in = i * n;
		for (int j = 0; j < i; j++)
		{
			int ij = in + j;
			int ji = j * n + i;
			G[ij] = 0.5 * (G[ij] + G[ji]);
			G[ji] = G[ij];
		}
	}

	//line 42
	if (error_tol <= 0)
		error_tol = n * 0.0000001;


	//line 51, default diaones = 0

	if (maxitmeth <= 0)
		maxitmeth = 200;

	if (maxit <= 0)
		maxit = 200;

	//line 65, pre = true

	//line 73
	//Determination of constants
	int Iter_inner = 40; // Maximum number of Line Search in Newton method
	double tol = 0.000001; // constant for the angle condition
	double sigma_1 = 0.0001; // tolerance in the line search of the Newton method

	//b0 = ones(n, 1)
	std::vector<double> b0(n, 1.0); //length n

	//method = ’ pminres ’

	//method is a string which determinates the iterative method
		// method = ’CG ’ - Conjugate gradient method(NAG - routine)
		// method = ’ SYMMLQ ’ - SYMMLQ(NAG - routine)
		// method = ’ pminres ’ - MINRES(implemetation of David Silvester) (default)
		// method = ’ RGMRES ’ - restarted generalized minimal residual
		// method(NAG - routine)
		// method = ’ TFQMR ’ - transpose free quasi minimal residual method
		// (NAG - routine)
		// method = ’ BIGCSTAB ’ - bi - conjugate gradient stabilized method
		// (NAG - routine)
	//line 93
	//Initial values
	//y = b0 - diag(G)
	std::vector<double> y(b0); //length n
	for (int i = 0; i < n; i++)
		y[i] -= G[i * n + i]; //-Gii
	//x0=y
	std::vector<double> x0(y); //length n

	int f_eval = 0; //line 98

	//[P, lambda] = eigdecomp(G, y, flag2, n)
	std::vector<double> P(n * n, 0.0); //length n^2, eigenvectors are columns
	std::vector<double> lambda(n, 0.0); //length n, eigenvalues
	eigdecomp(P, lambda, G, y, n);

	// compute grad(theta(y))
	//[f0, Fy] = gradient(y, lambda, P, b0, n)
	double f0;
	std::vector<double> Fy; //length n
	gradient(&f0, Fy, y, lambda, P, b0, n);

	f_eval = f_eval + 1; //103

	//b = b0 - Fy
	std::vector<double> b(b0); //length n
	for (int i = 0; i < n; i++)
		b[i] -= Fy[i];

	//grad = -b
	std::vector<double> grad(b); //length n
	for (int i = 0; i < n; i++)
		grad[i] *= -1;

	//normgrad = norm(b)
	double normgrad = K_matrix::FrobeniusNorm(b, n);

	//Omega = omega_matrix(lambda, n)
	std::vector<double> Omega; //will have length n^2
	Omega = omega_matrix(lambda, n);

	//line 108
	int k = 0;
	//NEWTON - ITERATION STARTS
	int flag;
	std::vector<double> d(n, 0.0); //length n
	while (normgrad > error_tol && k < maxit)
	{
		//COMPUTE V_k *d = -grad
		//if pre (always true)
		//precond = precond_matrix(Omega, P)
		std::vector<double> precond; //will have length n
		precond = precond_matrix(Omega, P, n);
		//[d, flag, iterk] = solver(@( x) Jacobian_matrix(x, Omega, P, n, 0), n, b,	min(0.5, normgrad), maxitmeth, method, precond )
		for (int i = 0; i < n; i++)
			d[i] = precond[i];
		flag = gmres(Omega, P, b, d, n, maxitmeth, minTwoDoubles(0.5, normgrad));

		//TEST CONDITION 2 (angle condition, descent conditon)
		double slope = 0.0;
		if (flag == 0) //"minres converged to the desired tolerance tol within maxit iterations"
		{
			//norm_d = norm(d)
			double norm_d = K_matrix::FrobeniusNorm(d, n);
			if (abs(norm_d) <= eps)
			{
				//X = computeX(P, lambda, n, diaones)
				std::vector<double> X;
				X = computeX(P, lambda, n, diaones);
				//error_function(12)
				*ifail = true;
				return X;
			}
			//Test for descent direction and angle - condition
			double grad_times_d = 0;
			for (int i = 0; i < n; i++)
				grad_times_d += grad[i] * d[i];
			if (-grad_times_d / (norm_d * norm_d) < minTwoDoubles(tol, normgrad))
				flag = 1; //"minres iterated maxit times but did not converge"
			else
				slope = grad_times_d; //slope = (grad)’* d
		}
		//if iterative method was unsuccessful or the angle condition failed, use the negative gradient direction
		if (flag != 0)
		{
			//d = -grad
			for (int i = 0; i < n; i++)
				d[i] = -grad[i];
			//slope = -normgrad ^ 2
			slope = -normgrad * normgrad;
		}


		//Temporary update
		// y = x0 + d
		for (int i = 0; i < n; i++)
			y[i] = x0[i] + d[i];
		//line 148, [P, lambda] = eigdecomp(G, y, flag2, n)
		eigdecomp(P, lambda, G, y, n);
		//[f, Fy] = gradient(y, lambda, P, b0, n) //compute F(y)
		double f;
		gradient(&f, Fy, y, lambda, P, b0, n);
		//norm_x = norm(x0)
		double norm_x = K_matrix::FrobeniusNorm(x0, n);


		f_eval = f_eval + 1;

		//line 153, ARMIJO - CONDITION
		//test whether accuracy problem can occur
		if ((abs(static_cast<double>(f - f0) / (1 + abs(f) + abs(f0))) < 100 * eps) && (abs(0.01 * slope / (1 + abs(f0) + abs(f))) < eps))
		{
			//TEST convergence condition
			//line 158
			tempDouble = 0.0;
			for (int i = 0; i < n; i++)
			{
				double diff = Fy[i] - b0[i];
				tempDouble += diff * diff;
			}
			if (sqrt(tempDouble) / normgrad > 0.9) //if norm(Fy - b0) / normgrad > 0.9
			{
				//Take negative gradient direction
				//d = -grad
				for (int i = 0; i < n; i++)
					d[i] = -grad[i];
				if (K_matrix::FrobeniusNorm(d, n) / (10 + 10 * norm_x) < eps) //if norm(d) / (1 + norm_x) / 10 < eps
				{
					//error_function(1) (warning d too small)
					*ifail = true;
					break;
				}
				//y = x0 + d
				for (int i = 0; i < n; i++)
					y[i] = x0[i] + d[i];
				//[P, lambda] = eigdecomp(G, y, flag2, n)
				eigdecomp(P, lambda, G, y, n);
				//168
				//[f, Fy] = gradient(y, lambda, P, b0, n)
				gradient(&f, Fy, y, lambda, P, b0, n);

				f_eval = f_eval + 1;
			}
			//else
				//print_function(prnt, 14)
		}
		else
		{	//ARMIJO - BACK - TRACKING
			int k_inner = 0;
			while (k_inner <= Iter_inner && f > f0 + sigma_1 * pow(0.5, k_inner) * slope)
			{
				k_inner = k_inner + 1;
				//y = x0 + 0.5^ k_inner * d (backtracking)
				tempDouble = pow(0.5, k_inner);
				for (int i = 0; i < n; i++)
					y[i] = x0[i] + tempDouble * d[i];
				//[P, lambda] = eigdecomp(G, y, flag2, n)
				eigdecomp(P, lambda, G, y, n);
				//[f, Fy] = gradient(y, lambda, P, b0, n) (compute F(y))
				gradient(&f, Fy, y, lambda, P, b0, n);
				if ((abs(static_cast<double>(f - f0) / (1 + abs(f) + abs(f0))) < 100 * eps) && (abs(0.01 * slope / (1 + abs(f0) + abs(f))) < eps))
				{
					//TEST convergence condition
					//line 183
					//if norm(Fy - b0) / normgrad > 0.9
					tempDouble = 0.0;
					for (int i = 0; i < n; i++)
					{
						double diff = Fy[i] - b0[i];
						tempDouble += diff * diff;
					}
					if (sqrt(tempDouble) / normgrad > 0.9)
					{
						//d = 0 (next if - condition is satisfied)
						for (int i = 0; i < n; i++)
							d[i] = 0.0;
						break;
					}
					break;
				}
			}//loop for while
			//line 193, Check whether step is sufficiently large
			if (k_inner > 0 && 0.1 * pow(0.5, k_inner) * K_matrix::FrobeniusNorm(d, n) / (1 + norm_x) < eps)
			{
				//Take negative gradient direction
				//d = -grad
				for (int i = 0; i < n; i++)
					d[i] = -grad[i];
				if (0.1 * K_matrix::FrobeniusNorm(d, n) / (1 + norm_x) < eps)
				{
					//line 198, error_function(1) (warning d too small)
					*ifail = true; //unsure about this line
					break;
				}
				//y = x0 + d
				for (int i = 0; i < n; i++)
					y[i] = x0[i] + d[i];
				//line 203, [P, lambda] = eigdecomp(G, y, flag2, n)
				eigdecomp(P, lambda, G, y, n);
				//[f, Fy] = gradient(y, lambda, P, b0, n)
				gradient(&f, Fy, y, lambda, P, b0, n);
			}
			//print_function(prnt, 5, k, k_inner) (print the iterations number)
			f_eval = f_eval + k_inner;
		} //line 208, end

		//UPDATE
		f0 = f; //int
		//x0 = y
		for (int i = 0; i < n; i++)
			x0[i] = y[i];
		k = k + 1;
		//b = b0 - Fy
		for (int i = 0; i < n; i++)
			b[i] = b0[i] - Fy[i];
		//norm_b = norm(b)
		double norm_b = K_matrix::FrobeniusNorm(b, n);
		//TEST whether norm of the gradient has become smaller
		if (0.1 * abs(norm_b - normgrad) / (1 + norm_b + normgrad) < eps)
		{
			//line 218, error_function(1) (warning minimal value achieved)
			*ifail = true;
			break;
		}
		normgrad = norm_b;
		//grad = -b
		for (int i = 0; i < n; i++)
			grad[i] = -b[i];
		//Omega = omega_matrix(lambda, n)
		Omega = omega_matrix(lambda, n);
	} //end loop for while
	if (k > maxit)
	{
		//error_function(2) (warning Maximal iteration number achieved)
		*ifail = true;
	}
	//line 228
	int iter = k;
	//compute X
	double minEig = minVec(lambda, n);
	if (iter == 0 && minEig >= 0.0)
	{
		//Setting the diagonal to 1 by our chosen
		//starting value was
		//sufficient to satisfy the stopping criterion
		//and all eigenvalues are nonnegative
		//X = G - diag(diag(G)) + eye(n)
		std::vector<double> X(G);
		for (int i = 0; i < n; i++)
			X[i * (n + 1)] = 1.0; //Xii = 1
		return X;
	}
	//line 238
	std::vector<double> X;
	X = computeX(P, lambda, n, diaones);
	return X;
} //end cornewton


//function for computing the solution X
//function[X] = computeX(P, lambda, n, diaones) where default diaones=0
std::vector<double> K_borsdorf::computeX(std::vector<double> & P, std::vector<double> & lambda, int n, int diaones) //output n by n, P n by n, lambda length n
{
	std::vector<double> C = K_matrix::transpose(P, n, n);
	//for i = 1: n
	for (int i = 0; i < n; i++) //row i
	{
		//C(i, :) = max(0, lambda(i)) * C(i, :)
		const int in = i * n;
		if (lambda[i] <= 0.0)
		{
			for (int j = 0; j < n; j++)
				C[in + j] = 0.0; //Cij = 0
		}
		else //lambda[i] > 0
		{
			for (int j = 0; j < n; j++)
				C[in + j] *= lambda[i]; //multiply Cij by ith eigenvalue
		}
	}
	std::vector<double> X = K_matrix::mult_ignoreZeros(P, C, n, n, n, 0.0);
	if (diaones == 1) //set diagonal to 1
	{
		//X = X - diag(diag(X)) + eye(n)
		for (int i = 0; i < n; i++)
			X[i * (n + 1)] = 1.0; //Xii = 1
	}
	if (diaones == 2) //rescale diagonal to 1 and set diagonal to 1
	{
		std::vector<double> scale = diagReduce(X, n);
		//scale(find(scale <= 1.0 e - 6)) = 1.0 e - 6
		for (int i = 0; i < n; i++)
		{
			if (scale[i] < 0.000001)
				scale[i] = 0.000001;
		}
		//scale = scale . ^ (-1 / 2)
		for (int i = 0; i < n; i++)
			scale[i] = pow(scale[i], -0.5);
		//X = diag(scale) * X * diag(scale)
		for (int i = 0; i < n; i++)
		{
			const int in = i * n;
			for (int j = 0; j < n; j++)
				X[in + j] *= scale[i] * scale[j]; //multiply Xij by scale[i] and scale[j]
		}
		//X = X - diag(diag(X)) + eye(n)
		for (int i = 0; i < n; i++)
			X[i * (n + 1)] = 1.0; //Xii
	}
	//X = (X + X ’) / 2 (Ensuring symmetry_
	for (int i = 0; i < n; i++)//row i, fill in lower triangular section
	{
		const int in = i * n; //such that X[in] = X[i,0]
		int ji = i; //such that X[ji] = X[j,i] in loop below
		for (int j = 0; j < i; j++) //column j<i, ij below diagonal
		{
			const int ij = in + j;
			X[ij] += X[ji];
			X[ij] *= 0.5;
			ji += n;
		}
	}
	for (int i = 0; i < n; i++)//row i, fill in upper triangular section
	{
		const int in = i * n; //such that X[in] = X[i,0]
		int ji = (i + 1) * n + i; //such that X[ji] = X[j,i] in loop below
		for (int j = i + 1; j < n; j++) //column j>i, ij above diagonal
		{
			X[in + j] = X[ji]; //Xij = Xji
			ji += n;
		}
	}
	return X;
}


//function for computing the eigenvalue decomposition of G + Diag(y)
//function [P, lambda] = eigdecomp(G, y, flag2, n)
//P n by n, G n by n, lambda length n, y length n
void K_borsdorf::eigdecomp(std::vector<double> & P, std::vector<double> & lambda, const std::vector<double> & G, std::vector<double> & y, int n) //flag2 = 0, 1, 2 or 3. Assume 0.
{
	// C = G + diag(y)
	std::vector<double> C(G); //C = G, length n^2
	for (int i = 0; i < n; i++)
		C[i * n + i] += y[i]; //Cii = Gii + yi
	//assume case 0 for now
	//P has length n^2, its columns will be eigenvectors of C
	//lambda has length n, its elements are eigenvalues of C
	K_jacobi::JacobiEigenvectors2(lambda, P, C, n);
}


//function for generating J(y)
//all vectors length n except P which is n by n
void K_borsdorf::gradient(double* f, std::vector<double> & userFy, std::vector<double> & y, std::vector<double> & lambda, std::vector<double> & P, std::vector<double> & b0, int n)
{
	std::vector<double> H(P); //set H=P, length n^2

	int in;

	//for j=1:n
	//	H(:, j) = max(lambda(j), 0) * H(:, j)
	for (int i = 0; i < n; i++)//row i
	{
		in = i * n;
		for (int j = 0; j < n; j++)//column j
		{
			if (lambda[j] < 0.0)
				H[in + j] = 0.0; //Hij=0
			else
				H[in + j] *= lambda[j]; //multiply Hij by jth eigenvalue
		}
	}
	//Fy = sum(P.*H, 2), length n (ith element = sum of ith row of P.*H)
	std::vector<double> Fy(n, 0.0);
	for (int i = 0; i < n; i++)//row i
	{
		in = i * n;
		for (int j = 0; j < n; j++)//column j
			Fy[i] += P[in + j] * H[in + j]; //increase Fy[i] by Pij*Hij
	}
	userFy = Fy;
	//f = sum(max(lambda, 0) .^ 2)
	*f = 0.0;
	for (int i = 0; i < n; i++)
	{
		if (lambda[i] > 0)
			* f += lambda[i] * lambda[i];
	}
	//f = 0.5* f - b0 ’* y
	*f = 0.5 * (*f);
	for (int i = 0; i < n; i++)
		* f -= b0[i] * y[i];
}

//line 393 % function for generating the constructed matrix M_y
//lambda has length n, output is n by n matrix
std::vector<double> K_borsdorf::omega_matrix(std::vector<double> & lambda, int n)
{
	std::vector<double> omega(n * n, 1.0); //omega = ones(n,n)
	for (int i = 0; i < n; i++)
	{
		const int in = i * n;
		for (int j = 0; j < n; j++)
		{
			//set omega[i,j]
			if (abs(lambda[i] - lambda[j]) > 0.0000000001) //if abs(lambda(i) - lambda(j)) > 1.0e - 10
				omega[in + j] = (maxTwoDoubles(0.0, lambda[i]) - maxTwoDoubles(0.0, lambda[j])) / (lambda[i] - lambda[j]);//omega(i, j) = (max(0, lambda(i)) - max(0, lambda(j))) / (lambda(i)-lambda(j))
			else if (maxTwoDoubles(lambda[i], lambda[j]) < 0.000000000000001) // else if max(lambda(i), lambda(j)) <= 1.0 e - 15
				omega[in + j] = 0.0; //omega(i, j) = 0
		}
	}
	return omega;
}

double K_borsdorf::maxTwoDoubles(double a, double b)
{
	if (a < b)
		return b;
	return a;
}

double K_borsdorf::minTwoDoubles(double a, double b)
{
	if (a < b)
		return a;
	return b;
}

double K_borsdorf::minVec(std::vector<double> vec, int vecLen)
{
	double result = vec[0];
	for (int i = 1; i < vecLen; i++)
	{
		if (vec[i] < vec[0])
			result = vec[i];
	}
	return result;
}


//inputs both n by n matrices, output has length n
std::vector<double> K_borsdorf::precond_matrix(std::vector<double> & Omega, std::vector<double> & P, int n)
{

	//PS = P . ^ 2
	std::vector<double> PS(n * n, 0.0);
	for (int i = 0; i < n * n; i++)
		PS[i] = P[i] * P[i];
	std::vector<double> tempVec = K_matrix::mult_ignoreZeros(PS, Omega, n, n, n, 0.0);
	for (int i = 0; i < n * n; i++)
		tempVec[i] *= PS[i];
	std::vector<double> c(n, 0.0);
	for (int i = 0; i < n; i++)
	{
		const int in = i * n;
		for (int j = 0; j < n; j++)
			c[i] += tempVec[in + j]; //tempVec[i,j]
	}
	//c(find(c <= 1.0 e - 6)) = 1.0 e - 6
	for (int i = 0; i < n; i++)
	{
		if (c[i] <= 0.000001)
			c[i] = 0.000001;
	}
	return c;
}


//line 543
//[d, flag, iterk] = solver(@( x) Jacobian_matrix(x, Omega, P, n, 0), n, b,	min(0.5, normgrad), maxitmeth, method, precond )
//Omega and P n by n matrices, x and output length n
std::vector<double> K_borsdorf::Jacobian_matrix(std::vector<double> & x, std::vector<double> & Omega, std::vector<double> & P, int n)//assume ep=0
{
	//H = P
	std::vector<double> H(P);

	for (int i = 0; i < n; i++)
	{
		const int in = i * n;
		//H(i, :) = x(i) * H(i, :), H = diag(x)* P
		for (int j = 0; j < n; j++)
			H[in + j] *= x[i]; //multiply Hij by x
	}
	K_matrix::transposeSquareMatrix(H, n); //P' * diag(x)
	//HTP = K_matrix::MatrixMult_ignoreZeros(H, P, n); //this line needs to be efficient because cornewton's while-loop calls gmres, gmres's for-loop calls arnoldi, which calls Jacobian_matrix
	std::vector<double> HTP = K_matrix::mult_ignoreZeros(H, P, n, n, n, 0.0);
	for (int i = 0; i < n * n; i++)
		H[i] = Omega[i] * HTP[i];
	//H = P * H
	//H = K_matrix::MatrixMult_ignoreZeros(P, H, n); //this line needs to be efficient because cornewton's while calls gmres, gmres's for calls arnoldi, which calls Jacobian_matrix
	H = K_matrix::mult_ignoreZeros(P, H, n, n, n, 0.0);
	//Ax = sum(P.*H, 2) + ep * x
	//assume ep=0
	std::vector<double> Ax(n, 0.0);
	for (int i = 0; i < n; i++)
	{
		const int in = i * n;
		for (int j = 0; j < n; j++)
			Ax[i] += P[in + j] * H[in + j];
	}
	return Ax;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following two functions are copied and pasted from K_minres.cpp
//Same as originals except A*x is replaced by Jacobian_matrix(x)
//Also, an integer flag is returned, similar to matlab function minres()
//Therefore K_borsdorf::gmres(...x...) will modify x such that b=Jacobian_matrix(x,Omega,P,n)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//function[x, e] = gmres(A, b, x, max_iterations, threshold)
//same as K_minres::gmres() except multiplication by A is replaced by Jacobian_matrix()
//Omega and P n by n, b and x length n
int K_borsdorf::gmres(std::vector<double> & Omega, std::vector<double> & P, std::vector<double> & b, std::vector<double> & x, int n, int max_iterations, double threshold)
{
	int flag = 1; //"minres iterated maxit times but did not converge" (matlab minres)
	const int m = max_iterations;
	std::vector<double> r = Jacobian_matrix(x, Omega, P, n); //instead of r=A*x
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
	if (r_norm < 0.00000000001)
		return true; //unlikely, but do anyway
	//Q(:, 1) = r / r_norm
	std::vector<double> Q(n * (m + 1), 0.0); //let Q be n by m+1
	const double inv_r_norm = 1.0 / r_norm;
	for (int i = 0; i < n; i++)
		Q[i * (m + 1)] = r[i] * inv_r_norm; //set Q[i,0]
	//beta = r_norm * e1    (Note: this is not the beta scalar in section "The method" above but the beta scalar multiplied by e1)
	std::vector<double> beta(m + 1, 0.0); //length m+1
	beta[0] = r_norm;
	//for k = 1 : m
	std::vector<double> H((m + 1) * m, 0.0); //H m+1 by m
	int k;
	for (k = 0; k < m; k++) //maybe perform check to see if iteration has changed, if not set flag to 3?
	{
		//run arnoldi
		//[H(1:k + 1, k) Q(:, k + 1)] = arnoldi(A, Q, k)
		arnoldi(H, Omega, P, Q, k, n, m); //send Omega and P to function instead of A

		//eliminate the last element in H ith row and update the rotation matrix
		//[H(1:k + 1, k) cs(k) sn(k)] = apply_givens_rotation(H(1:k + 1, k), cs, sn, k)
		K_minres::apply_givens_rotation(H, cs, sn, k, m);

		//update the residual vector
		beta[k + 1] = -sn[k] * beta[k];
		beta[k] = cs[k] * beta[k];
		const double error = abs(beta[k + 1]) / b_norm;

		if (error <= threshold)
		{
			flag = 0; //"minres converged to the desired tolerance tol within maxit iterations" (matlab minres)
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
			flag = 4; //"One of the scalar quantities calculated during minres became too small or too large to continue computing" (matlab minres)
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
	return flag;
} //end

//same as minres::arnoldi except multiplication by A is replaced by Jacobian_matrix()
//Q is m by n, A is n by n, H m+1 by m
void K_borsdorf::arnoldi(std::vector<double> & H, std::vector<double> & Omega, std::vector<double> & P, std::vector<double> & Q, int k, int n, int m)
{
	std::vector<double> Qk = K_minres::getCol(Q, m + 1, n - 1, k);
	std::vector<double> q = Jacobian_matrix(Qk, Omega, P, n); //instead of r=A*Qk
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

