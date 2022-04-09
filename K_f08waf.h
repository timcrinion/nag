#pragma once

#include <vector>

class K_f08waf
{
	// member functions
public:

	static void nagcal(int NDIMEN, int IMAX, double* DA, double* DB, double* R, double* V, int* LNAGOK);

private:

	//Function that outputs generalised right-eigenvectors and eigenvalues of A and B
	//We say k (number) is a right generalised eigenvalue of A and B, and v (n-vector) is a right generalised eigenvector of A and B, if A * v = k * B * v
	//n eigenvales placed in R (length n), n vectors placed in V (n by n), such that ith row of V corressponds to eigenvalue R[i]
	//Convention here is that the entry of the 2D matrix M (row i, column j) is stored in M[i * numColumnsM + j]
	static void nagcal_inner(int n, int IMAX, std::vector<double>& A, std::vector<double>& B, std::vector<double>& R, std::vector<double>& V, bool& LNAGOK);

	//calculate eigenvectors and eigenvalues of n by n matrix A
	//eigVals length n
	//eigVecs n by n, such that 
	//WARNING: This algorithm is fine for small matrices but inefficient for larger ones
	//returns true if eigenvalues real, false if complex
	static bool QR_algorithm(std::vector<double>& A, int n, std::vector<double>& eigVecs, std::vector<double>& eigVals);


	//Given three existing n by n matrices A,B,C, set A to B*C 
	static void setAtoBC(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, int n);


	//get eigenvectors of upper triangular n by n matrix U
	//eigenvalues are diagonal values Uii
	//place eigenvectors in eigVecs such that jth column corressonds to jth eigenvalue Ujj
	static void eigUppTri(std::vector<double>& U, int n, std::vector<double>& eigVecs);


	//function that calculates eigenvecter corressponding to jth eigenvalue of n by n matrix U
	//jth eigenvalue is Ujj
	//answer stored in eigVec
	static void eigUppTri_colj(std::vector<double>& U, int n, std::vector<double>& eigVec, int j);


	//returns true if n by n matrix M is upper triangular, false otherwise
	static bool isUppTri(std::vector<double>& M, int n);

	//scales each column of n by n matrix M such that each has norm 1
	//if a column is all zeros, it is left untouched
	static void normalizeColumns(std::vector<double>& M, int n);

	//returns true if n by n matrix M is symmetric, false otherwise
	static bool isSymmetric(std::vector<double>& M, int n);

	static bool genEigSym(std::vector<double>& A, std::vector<double>& B, int n, std::vector<double>& eigVals, std::vector<double>& eigVecs);

	static bool genEigSym2(std::vector<double>& A, std::vector<double>& B, int n, std::vector<double>& eigVals, std::vector<double>& eigVecs);

	static double vecDiff(std::vector<double>& A, std::vector<double>& B, int length);

};
