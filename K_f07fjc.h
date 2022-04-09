//This is the header file for K_f07fjc.cpp
#pragma once
#include <vector>

class K_f07fjc
{
private:

public:
	//fills column_j with the values in jth column of matrix M (where 0 <= j < numCols)
	static void getColumn(std::vector<double>& M, std::vector<double>& column_j, int numRows, int numCols, int j);

	//Sets jth column of matrix M equal to column_j (where 0 <= j < numCols)
	static void setColumn(std::vector<double>& M, std::vector<double>& column_j, int numRows, int numCols, int j);

	//Returns inverse X of lower diagonal n by n matrix L
	//Method 2 from section 2.1 of http://www.netlib.org/lapack/lawnspdf/lawn27.pdf
	//Assume n>1
	static std::vector<double> inverseLowerTriang(std::vector<double>& L, int n, bool* invertible);

	//The following may never be used

	//return inverse of symmetric n by n matrix S
	static std::vector<double> invSymmetric(std::vector<double>& S, int n, bool* fail);

	//Modifies matrix M by swapping row a with row b (0<= a,b < numRows)
	static void swapRows(std::vector<double>& M, int numCols, int a, int b);

	//Modifies matrix M by multiplying row i by factor (0<= i < numRows)
	static void multRow(std::vector<double>& M, int numCols, int i, double factor);

	//Modifies matrix M by adding a factor times row b, to row a (0<= a,b < numRows)
	static void addMultRow(std::vector<double>& M, int numCols, int a, int b, double factor);

	//Apply Gauss-Jordan elimination to a nxn matrix M
	//Makes X inverse of M. If M invertible, returns true, otherwise false.
	static bool GaussJordanSquareInverse(std::vector<double>& M, std::vector<double>& X, int n);
};