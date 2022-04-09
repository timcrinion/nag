//This is the header file for K_Cartesian2D.cpp
#pragma once

#include <vector>

class K_Cartesian2D
{
private:

public:

	static double DistSquared(double x1, double y1, double x2, double y2);
	static void LineEquidistant(double x1, double y1, double x2, double y2, double& x, double& y, double& dx, double& dy);
	static bool Intersection(double x1, double y1, double dx1, double dy1, double x2, double y2, double dx2, double dy2, double& x, double& y);
	static bool Circle3Pts(double x1, double y1, double x2, double y2, double x3, double y3, double& x0, double& y0, double& rSquared);
	static double PseudoAng(double x, double y, double xOrigin, double yOrigin);
	//function that tests where the point (x,y) is in relation to the line travelling from start to end
	//function returns -1 if point (x,y) right of the line travelling from start point to end point
	//function returns  0 if point (x,y) on top of the line travelling from start point to end point
	//function returns  1 if point (x,y) left of the line travelling from start point to end point
	static int Left(double xStart, double yStart, double xEnd, double yEnd, double x, double y);
	static bool TwoFiniteLinesMeet(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double& x, double& y);
	static bool PtFromTriangle(double x, double y, int& edge, double x0, double y0, double x1, double y1, double x2, double y2);
	static double PseudoAngFromP(double x, double y, double xOrigin, double yOrigin, double angleAtPointP);
	static bool PtOnFiniteLine(double Ax, double Ay, double Bx, double By, double Px, double Py);

	static void OrderPointsByDistance(double Px, double Py, std::vector<double>& x, std::vector<double>& y, std::vector<int>& indicesToSort, int numIndicesToSort);
	static std::vector<int> Nearest_k(double px, double py, int k, int n, int m, std::vector<double>& x, std::vector<double>& y, std::vector<int>& grid, std::vector<int>& next, double dx, double dy, double xMin, double yMin);
	/*create m by m grid that divides the area spanned by n nodes up into rectangles, then places each node into its rectangle
	in:
		n,m,x,y
		x and y are vectors of length n. The n input nodes' positions are represented by x and y. "node i" is at (x[i],y[i])
	out:
		grid,next,dx,dy,xMin,yMix
		On input, grid and next should be full of -1s
		grid has length m*m. Set grid[a,b] = grid [a*m+b] = i if:
			a <= (x[i]-xMin)/dx < a+1 and
			b <= (y[i]-yMin)/dy < b+1 and
		next has length n. next[i] = j means point j lies in same cell as point i (unless j=-1)
		xMin will be the top of the grid, xMax the bottom, yMin the left, yMax the right (not the min and max values of input vectors x and y)
		dx will be the height of a cell, dy the width of a cell*/
	static void SetupGrid(int n, int m, std::vector<double>& x, std::vector<double>& y, std::vector<int>& grid, std::vector<int>& next, double& dx, double& dy, double& xMin, double& yMin);
	static void RecordNodesInCell(std::vector<int>& record, int& numRec, std::vector<int>& grid, int m, int i, int j, std::vector<int>& next);
	static void RecordVertically(std::vector<int>& record, int& numRec, std::vector<int>& grid, int m, int i1, int i2, int j, std::vector<int>& next);
	static void RecordHorizontally(std::vector<int>& record, int& numRec, std::vector<int>& grid, int m, int i, int j1, int j2, std::vector<int>& next);

};

