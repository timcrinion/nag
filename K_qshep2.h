#pragma once
#include <vector>

class K_qshep2
{
	//Following taken from https://jblevins.org/mirror/amiller/toms660.f90

private:
	static void qshep2(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nq, int nw, int nr,
		std::vector<int>& lcell, std::vector<int>& lnext, double* xmin, double* ymin, double* dx, double* dy, double* rmax, std::vector<double>& rsq, std::vector<double>& a, int* ier);

	static void getnp2(double px, double py, std::vector<double>& x, std::vector<double>& y, int nr, std::vector<int>& lcell, std::vector<int>& lnext,
		double xmin, double ymin, double dx, double dy, int* np, double* dsq, std::vector<bool>& marked);

	static void givens(std::vector<double>& B, int j, int irow, double* c, double* s); //CALL givens(b(j, j), b(j, irow), c, s) //not sure about first 2 arguments

	static void rotate(int n, double c, double s, std::vector<double>& b, int jp1, int j, int irow);

	static void setup2(double xk, double yk, double fk, double xi, double yi, double fi, double s1, double s2, double r, std::vector<double>& M, int irow);

	static void store2(int n, std::vector<double>& x, std::vector<double>& y, int nr, std::vector<int>& lcell, std::vector<int>& lnext, double* xmin, double* ymin, double* dx, double* dy, int* ier);

	static double qs2val(double px, double py, int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f,
		int nr, std::vector<int>& lcell, std::vector<int>& lnext, double xmin, double ymin, double dx, double dy, double rmax, std::vector<double>& rsq, std::vector<double>& a);

	static void qs2grd(double px, double py, int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nr, std::vector<int>& lcell, std::vector<int>& lnext,
		double xmin, double ymin, double dx, double dy, double rmax, std::vector<double>& rsq, std::vector<double>& a, double* q, double* qx, double* qy, int* ier);

public:
	static int removeDuplicates(std::vector<double>& X, std::vector<double>& Y, std::vector<double>& F, std::vector<double>& newX, std::vector<double>& newY,
		std::vector<double>& newF, int n, double dist);

	static void e01sgf(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nq, int nw, int nr, std::vector<int>& lcell, std::vector<int>& lnext,
		double& xmin, double& ymin, double& dx, double& dy, double& rmax, std::vector<double>& rsq, std::vector<double>& a, int& ifail);

	static void e01shf(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nr, std::vector<int>& lcell, std::vector<int>& lnext, double xmin, double ymin,
		double dx, double dy, double rmax, std::vector<double>& rsq, std::vector<double>& a, int lenUV, const double* U, const double* V, double* Q, double* QX, double* QY, int& IFAIL);
};
