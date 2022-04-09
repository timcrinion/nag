#pragma once
//C++ implementation of MINRES method: https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Example_code
#include <vector>


class K_minres
{
private:

public:

	//A n by n matrix
	//b vector of length n
	//Want to solve Ax=b for x

	//function[x, e] = gmres(A, b, x, max_iterations, threshold)
	static bool gmres(std::vector<double>& A, std::vector<double>& b, std::vector<double>& x, int n, int max_iterations, double threshold);


	//function[h, q] = arnoldi(A, Q, k)
	static void arnoldi(std::vector<double>& H, const std::vector<double>& A, std::vector<double>& Q, int k, int n, int m);


	//function[h, cs_k, sn_k] = apply_givens_rotation(h, cs, sn, k)
	static void apply_givens_rotation(std::vector<double>& H, std::vector<double>& cs, std::vector<double>& sn, int k, int m);


	//%%----Calculate the Given rotation matrix----%%
	//function[cs, sn] = givens_rotation(v1, v2)
	static void givens_rotation(double* cs, double* sn, double v1, double v2);


	//Returns M[0:i,j] (where 0 <= j < numColsM)
	static std::vector<double> getCol(std::vector<double>& M, int numColsM, int i, int j);

	//Sets M[0:i,j] to column (where 0 <= j < numCols)
	//column must have length of at least i+1
	static void setCol(std::vector<double>& M, std::vector<double>& column, int numColsM, int i, int j);
};