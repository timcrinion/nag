//This is the header file for K_interpTri.cpp
#pragma once
#include <vector>

class K_interpTri
{
private:
	static std::vector<int> MergeIntVectors(std::vector<int>& A, int lenA, std::vector<int>& B, int lenB, int& lenC);
	static int Neighbours2(int n, int p, std::vector<bool>& isOutput, std::vector<int>& output, std::vector<int>& qNeighbours,
		std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts, int discSize);


public:

	/*algorithm from section 2.2 of https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf
	inputs:
	n, number of nodes (this means there can't be more than 2n triangles)
	x, vector such that x[i] is x coordinate of ith node (x horizontal, increases from left to right)
	y, vector such that y[i] is y coordinate of ith node (y vertical, increases from down to up)
	y,
	outputs:
	TriAdj, a 2n by 3 vector such that triAdj[i,:] = indices of triangles adjacent to ith triangle (counterclockwise, -1 for none, as per Fig 2.1)
	TriPts, a 2n by 3 vector such that TriPts[i,:] = indices of nodes, counterclockwise. TriPts[i,j] joins TriAdj[i,j] and TriAdj[i,j+1] as shown

			 2
			/ \
		  2/   \1
		  /     \
		 0_______1
			 0
	*/




	static bool WeightedLeastSquares(int n, std::vector<double>& W, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, double& a, double& b, double& c, double& d, double& e, double f);

	static int Neighbours(int p, int Tp, std::vector<int>& TriAdj, std::vector<int>& TriPts, std::vector<int>& neighboursOut, int n);
	static int MergeOrderedNoRepeats(int lenA, std::vector<int>& A, int lenB, std::vector<int>& B, std::vector<int>& ABmerged);
	static int MergeInto(int lenA, std::vector<int>& A, int lenB, std::vector<int>& B);
	static void Nearest(int p, int k, int n, std::vector<double>& x, std::vector<double>& y, std::vector<int>& TriAdj, std::vector<int>& TriPts, std::vector<int>& randTri,
		std::vector<int>& kNearest, std::vector<bool>& isNear, std::vector<int>& qNeighbours, std::vector<int>& newDisc);
	static bool GetPartials(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& DX, std::vector<double>& DY,
		std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts);
	static double Interpolate(double xx, double yy, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<double>& fx, std::vector<double>& fy, bool linear);
	static void E01saf(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<int>& TriAdj, std::vector<int>& TriPts,
		std::vector<double>& dfdx, std::vector<double>& dfdy, std::vector<int>& grid, int gridDim, double& dx, double& dy, double& xMin, double& yMin, int& ifail, bool linear);
	static double E01sbf(double Px, double Py, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<int>& TriAdj, std::vector<int>& TriPts,
		std::vector<double>& dfdx, std::vector<double>& dfdy, std::vector<int>& grid, int gridDim, double dx, double dy, double xMin, double yMin, int& ifail, bool linear);
};

