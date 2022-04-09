#pragma once
//This is the header file for K_jacobi.cpp
#include <vector>

class K_jacobi
{
//All of these functions (except EigInvSym) have been copied from https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm#Algorithm
public:
	static void JacobiEigenvectors(std::vector<double>& eigenvalues, std::vector<double>& eigenvectors, std::vector<double>& S, int n); //finds eigenvectors and eigenvalues of real symmetric matrix S
	static void JacobiEigenvectors2(std::vector<double>& eigenvalues, std::vector<double>& eigenvectors, std::vector<double>& S, int n); //same as function above except eigenvalues and eigenvectors already full
	//Find inverse of nxn matrix A using eigendecomposition. invA input as an empty vector
	static bool EigInvSym(std::vector<double>& A, std::vector<double>& invA, int n); //this function might never be used

private:
	static int maxind(std::vector<double>& S, int n, int k);
	static void rotate(int k, int l, int i, int j, std::vector<double>& M, int numColsM, double cosPhi, double sinPhi);
	static void update(int k, double t, std::vector<double>& eigenvalues, std::vector<bool>& changed, int* state);
};

//JacobiEigenvectors does not work with S = { 4 , -30 , 60 , -35 , -30 , 300 , -675 , 420 , 60 , -675 , 1620 , -1050 , -35 , 420 , -1050 , 700 }; //4x4 matrix