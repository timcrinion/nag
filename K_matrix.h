#pragma once
#include <cmath>
#include <vector>


class K_matrix
{
public:
	//Function that multiplies matrix A by matrix B and stores answer in AB.
	//A must have numRowsA rows and numColsA columns.
	//B must have numColsA rows and numRowsA columns.
	//AB is square matrix of dimensions numRowsA by numRowsA.
	static std::vector<double> mult(std::vector<double>& A, int numRowsA, int numColsA, std::vector<double>& B, int numColsB);


	//Function that transposes square nxn matrix M
	static void transposeSquareMatrix(std::vector<double>& M, int n);


	//Function that transposes matrix M
	//M[i,j] stored as M[i * numColsM + j]
	static std::vector<double> transpose(std::vector<double>& M, int numRowsM, int numColsM);


	//return Frobenius norm of array x of length n
	static double FrobeniusNorm(std::vector<double>& x, int n);

	//return Frobenius difference between two arrays of length n
	static double FrobDiff(std::vector<double>& x, std::vector<double>& y, int n);


	//multiply A by B, output AB with dimensions numRowsA by numColsB
	//B must have numColsA rows
	//For efficiency, function ignores all entries <= tinyNumber
	static std::vector<double> mult_ignoreZeros(std::vector<double>& A, std::vector<double>& B, int numRowsA, int numColsA, int numColsB, double tinyNumber);


	static double innerProduct(std::vector<double>& u, std::vector<double>& v, int n);
	static void QR(std::vector<double>& A, int rowsA, int colsA, std::vector<double>& Q, std::vector<double>& R);
	static void givensLeft(double cos, double sin, std::vector<double>& M, int colsM, int i, int j);
	static void givensRight(double cos, double sin, std::vector<double>& M, int colsM, int i, int j);
	static void computeGivens(double m11, double m12, double m21, double m22, double& sinL, double& cosL, double& sinR, double& cosR);
	static void SVD(std::vector<double>& A, int n, std::vector<double>& U, std::vector<double>& Sigma, std::vector<double>& V, double epsilon);
	static void swapRows(std::vector<double>& M, int numColumns, int i, int j);
	static void swapColumns(std::vector<double>& M, int numRows, int numColumns, int i, int j);
	static bool uppTriSolution(std::vector<double>& A, std::vector<double>& x, std::vector<double>& y, int n);
	static bool uppTriInv(std::vector<double>& U, std::vector<double>& V, int n);
	static std::vector<double> mult3matrices(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, int rowsA, int rowsB, int rowsC, int colsC);
	static std::vector<double> mult4matrices(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, std::vector<double>& D, int rowsA, int rowsB, int rowsC, int rowsD, int colsD);

};