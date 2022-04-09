#pragma once

#include <vector>

class K_UpperHessenberg
{
	// member functions
public:

	//sets n by n matrix M to identity
	static void eye(std::vector<double>& M, int n);

	//Given n by n matrix A, fills unitary Q and upper hessenberg H such that H = Q * A * Q^T
	static void hessenbergReduction(std::vector<double>& A, int n, std::vector<double>& Q, std::vector<double>& H);

	//given two n by n matrices A and B, replace B by A*B
	static void leftMult(std::vector<double>& A, std::vector<double>& B, int n);

	//given two n by n matrices A and B, replace A by A*B
	static void rightMult(std::vector<double>& A, std::vector<double>& B, int n);

private:

	/*
	function that fills P (n by n) such that P is householder reflector I-2*v*v^T of a normalised vector v of length n
	we have P = P^T = P^-1

		 | 1-2v1v1   -2v1v2   -2v1v3          |
		 |  -2v2v1  1-2v2v2   -2v2v3          |
	P =  |  -2v3v1   -2v3v2  1-2v3v3          |
		 |                            1-2vnvn |

	The idea is that if |u|=1 then Px reflects x in the plane defined by the vector u 	*/
	static void householderReflector(std::vector<double>& P, std::vector<double>& v, int n);

	//normalise a vector v of length n
	//returns true if successful (ie original norm > 0)
	static bool normalise(std::vector<double>& v, int n);

	//multiply vector v of length n by scalar r
	static void multVecByScalar(std::vector<double>& v, int n, double r);

	//Given n by n matrix M, returns subcolumn { M[j+1,j]  M[j+2,j] ...  M[n-1,j] } below diagonal in column j, length n-j-1
	static std::vector<double> colBelowDiag(std::vector<double>& M, int n, int j);

	/*    fill P for K_f08waf::hessenbergReduction
			   _______________________________
			  |               |               |
			  |   identity    |     zeros     | j+1
	 want P = |_______________|_______________|
			  |               |               |
			  |     zeros     |     I-2uu'    | n-j-1
			  |_______________|_______________|
					 j+1            n-j-1

	 where u = input vector of length n-j-1       */
	static void fillP(std::vector<double>& P, int j, int n, std::vector<double>& u);

	//given three n by n matrices A,B,C, replace B by A*B*C
	static void leftRightMult(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, int n);

};
