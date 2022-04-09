//This C++ file contains cartesian 2D functions. As per tradition, x axis is positive to the right, y axis positive upwards
#include <algorithm>
#include <cmath>
#include <vector>
#include "K_Cartesian2D.h"
#include "KUtils.h"

//function that calculates SQUARED distance between two points (x1,y2) and (x2,y2)
double K_Cartesian2D::DistSquared(double x1, double y1, double x2, double y2)
{
	double dx = x1 - x2;
	double dy = y1 - y2;
	return dx * dx + dy * dy;
}

//function that calculates line equidistant from two points (x1,y1) and (x2,y2)
//returns line in the form (x,y) + r * (dx,dy) where r can be any real number
//In: x1, y1, x2, y2
//Out: x, y, dx, dy
void K_Cartesian2D::LineEquidistant(double x1, double y1, double x2, double y2, double& x, double& y, double& dx, double& dy)
{
	//let (x,y) be mid-point
	x = (x1 + x2) * 0.5;
	y = (y1 + y2) * 0.5;
	//swap dx and dy and multiply one of them by -1
	dx = y1 - y2;
	dy = x2 - x1;
}

//function that calculates intersection (x,y) of two lines
//line 1 is (x1, y1) + r * (dx1, dy1) where r can be any real number
//line 2 is (x2, y2) + r * (dx2, dy2) where r can be any real number
//returns true if unique intersection exists, false otherwise
bool K_Cartesian2D::Intersection(double x1, double y1, double dx1, double dy1, double x2, double y2, double dx2, double dy2, double& x, double& y)
{
	/*
	want r1 and r2 such that
	(x1,y1) + r1*(dx1,dy1) = (x2,y2) + r2*(dx2,dy2)
	or
	| dx1  -dx2 | |r1|   |x2-x1|
	| dy1  -dy2 | |r2| = |y2-y1|
	or
	|r1|   | dx1  -dx2 |^-1 |x2-x1|
	|r2| = | dy1  -dy2 |    |y2-y1|
	*/
	//calculate inverse of determinant of inverse 2x2 matrix above
	double detInv = -dx1 * dy2 + dy1 * dx2;
	if (abs(detInv) < 0.00000001) //lines parallel
		return false;
	//calculate r1
	double r1 = -dy2 * (x2 - x1) + dx2 * (y2 - y1);
	r1 /= detInv;
	//calculate intersection
	x = x1 + r1 * dx1;
	y = y1 + r1 * dy1;
	return true;
}

//finds centre and radius of circle intersecting 3 points (x1,y1), (x2,y2) and (x3,y3)
//output centre=(x0,y0) and squared radius rSquared
//returns true if circle finite, false if 3 input points colinear
bool K_Cartesian2D::Circle3Pts(double x1, double y1, double x2, double y2, double x3, double y3, double& x0, double& y0, double& rSquared)
{
	//find pair of input points furthest apart
	double dist12 = K_Cartesian2D::DistSquared(x1, y1, x2, y2);
	double dist13 = K_Cartesian2D::DistSquared(x1, y1, x3, y3);
	double dist23 = K_Cartesian2D::DistSquared(x2, y2, x3, y3);
	//relabel points 1,2,3 as points a,b,c such that a and c furthest apart
	//if 1 and 3 furthest apart
	double xa = x1;
	double ya = y1; //a=1
	double xb = x2;
	double yb = y2; //b=2
	double xc = x3;
	double yc = y3; //c=3
	if (dist23 > dist12 && dist23 > dist13) //2 and 3 furthest apart
	{
		xa = x2;
		ya = y2; //a=2
		xb = x1;
		yb = y1; //b=1
	}
	if (dist12 > dist13 && dist12 > dist23) //1 and 2 furthest apart
	{
		xb = x3;
		yb = y3; //b=3
		xc = x2;
		yc = y2; //c=2
	}

	//line (x4,y4) + r*(dx4,dy4) equidistant from a and b
	double x4, y4, dx4, dy4;
	K_Cartesian2D::LineEquidistant(xa, ya, xb, yb, x4, y4, dx4, dy4);
	// line (x5,y5) + r*(dx5,dy5) equidistant from b and c
	double x5, y5, dx5, dy5;
	K_Cartesian2D::LineEquidistant(xb, yb, xc, yc, x5, y5, dx5, dy5);

	//find intersection of line 4 and 5
	bool intersectionExists = K_Cartesian2D::Intersection(x4, y4, dx4, dy4, x5, y5, dx5, dy5, x0, y0);

	//radius squared
	if (intersectionExists)
		rSquared = K_Cartesian2D::DistSquared(x0, y0, x1, y1);
	return intersectionExists;
}

/*function that returns pseudoangle of (x,y) from origin where
x is positive from left to right
y is positive from down to up
Euclidean angle measured anticlockwise from negative x, i.e. atan(opp=yOrigin-y, adj=xOrigin-x)
However, instead of returning euclidean angle, function returns pseudoangle that
increases, stays the same, or decreases, at the same time as euclidean angle.
i.e., pseudoangle(x1,y1) < pseudoangle(x2,y2) if and only if euclidean angle(x1,y1) < euclidean angle(x2,y2)
(x,y) cannot be the origin*/
double K_Cartesian2D::PseudoAng(double x, double y, double xOrigin, double yOrigin)
{
	//(x,y) relative to origin
	double xRel = x - xOrigin;
	double yRel = y - yOrigin;
	//instead of returning atan(opp = -yRel, adj = -xRel), return this:
	double result = abs(yRel) / (abs(xRel) + abs(yRel)); //0 on x axis, increasing to 1 on y axis
	//At this point result is what we want only if point xRel,yRel in bottom left quadrant
	//But we want results to be 0 left of origin, 1 below origin, 2 right of origin, 3 above origin
	bool left = xRel <= 0.0;
	bool bottom = yRel <= 0.0;
	if (!left && bottom) //bottom right quadrant
		result = 2.0 - result; //1 below origin, increasing to 2 right of origin
	else if (!left && !bottom) //top right quadrant
		result = 2.0 + result; //2 right of origin, increasing to 3 above origin
	else if (left && !bottom) //top left quadrant
		result = 4.0 - result; //3 above origin, increasing to 4 right of origin
	return result;
}

//Similar to pseudoAng() except datum changed such that the angle of a point P is set to zero
double K_Cartesian2D::PseudoAngFromP(double x, double y, double xOrigin, double yOrigin, double angleAtPointP)
{
	double result = K_Cartesian2D::PseudoAng(x, y, xOrigin, yOrigin) - angleAtPointP; // between -4 and 4
	if (result < 0.0)
		result += 4.0; //between 0 and 4
	return result;
}

//function that tests where the point (x,y) is in relation to the line travelling from start to end
//function returns -1 if point x,y right of the line travelling from start point to end point
//function returns  0 if point x,y on top of the line travelling from start point to end point
//function returns  1 if point x,y left of the line travelling from start point to end point
int K_Cartesian2D::Left(double xStart, double yStart, double xEnd, double yEnd, double x, double y)
{
	//set (xStart,xEnd) as origin
	//All points (x,y) that lie on the line from Start to End have (x-xStart) * (yEnd-yStart) = (y-yStart) * (xEnd-xStart)
	//Consider four cases:
	//y increases and (xEnd-xStart) positive: Then RHS increases and (x,y) left of line
	//y decreases and (xEnd-xStart) positive: Then RHS decreases and (x,y) right of line
	//y increases and (xEnd-xStart) negative: Then RHS decreases and (x,y) right of line
	//y decreases and (xEnd-xStart) negative: Then RHS increases and (x,y) left of line
	//In all four cases, LHS stays the same. Therefore:
	double LHS = (x - xStart) * (yEnd - yStart);
	double RHS = (y - yStart) * (xEnd - xStart);
	if (LHS < RHS)
		return 1; //left
	if (LHS > RHS)
		return -1; //right
	if (KUtils::AreDoublesEqual(LHS, RHS))
		return 0; //colinear
	return 999; //something wrong!!!
}


//function that checks to see if two lines AB and CD of finite lengths intersect
//calculates their intersection (x,y)
//returns:
//true if finite lines intersect and are not parallel
//false if finite lines do not intersect or are parallel
bool K_Cartesian2D::TwoFiniteLinesMeet(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double& x, double& y)
{
	//min and max x and y coordinates of lines
	double xMinAB = std::min(Ax, Bx);
	double xMaxAB = std::max(Ax, Bx);
	double yMinAB = std::min(Ay, By);
	double yMaxAB = std::max(Ay, By);
	double xMinCD = std::min(Cx, Dx);
	double xMaxCD = std::max(Cx, Dx);
	double yMinCD = std::min(Cy, Dy);
	double yMaxCD = std::max(Cy, Dy);

	//if x-span or y-span of AB disjoint from that of CD, return false
	if ((xMaxAB < xMinCD) || (xMaxCD < xMinAB)) //disjoint in x axis
		return false;
	if ((yMaxAB < yMinCD) || (yMaxCD < yMinAB)) //disjoint in y axis
		return false;

	//get intersection of lines if they were to stretch to infinity
	bool notParallel = K_Cartesian2D::Intersection(Ax, Ay, Bx - Ax, By - Ay, Cx, Cy, Dx - Cx, Dy - Cy, x, y);
	if (!notParallel)
		return false;
	//return false if intersection outside AB
	if (abs(Ax - Bx) > abs(Ay - By)) //if AB fatter than tall
	{
		//if x outside AB
		if ((x < Ax && x < Bx) || (Ax < x && Bx < x))
			return false;
	}
	else //if AB taller than fat
	{
		//if y outside AB
		if ((y < Ay && y < By) || (Ay < y && By < y))
			return false;
	}
	//return false if intersection outside CD
	if (abs(Cx - Dx) > abs(Cy - Dy)) //if CD fatter than tall
	{
		//if x outside CD
		if ((x < Cx && x < Dx) || (Cx < x && Dx < x))
			return false;
	}
	else //if CD taller than fat
	{
		//if y outside CD
		if ((y < Cy && y < Dy) || (Cy < y && Dy < y))
			return false;
	}
	//intersection exists and is inside both AB and CD
	return true;
}


//returns true if point P inside triangle 012
//Let C be centroid of triangle 012, then 0 travels to 1 travels to 2 travels to 0, anticlockwise around C
//sets edge to 0 if C to P passes through 01
//sets edge to 1 if C to P passes through 12
//sets edge to 2 if C to P passes through 02
//function fine as long as P not too close to C, 0, 1 or 2
bool K_Cartesian2D::PtFromTriangle(double x, double y, int& edge, double x0, double y0, double x1, double y1, double x2, double y2)
{
	//centre of triangle
	double Cx = (x0 + x1 + x2) / 3;
	double Cy = (y0 + y1 + y2) / 3;
	//pseudoangles anticlockwise from CP
	double angP = K_Cartesian2D::PseudoAng(x, y, Cx, Cy);
	double ang0 = K_Cartesian2D::PseudoAngFromP(x0, y0, Cx, Cy, angP);
	double ang1 = K_Cartesian2D::PseudoAngFromP(x1, y1, Cx, Cy, angP);
	double ang2 = K_Cartesian2D::PseudoAngFromP(x2, y2, Cx, Cy, angP);
	//find edge that CP passes through
	if (ang0 < ang1)		//order 01
	{
		if (ang1 < ang2)		//order 012
			edge = 2; //02
		else if (ang2 < ang0)	//order 201
			edge = 1; //12
		else					//order 021
			edge = 0; //01
	}
	else					//order 10
	{
		if (ang0 < ang2)		//order 102
			edge = 1; //12
		else if (ang2 < ang1)	//order 210
			edge = 2; //02
		else					//order 120
			edge = 0; //01
	}
	//if edge intersects CP then P must be outside triangle
	double duff1, duff2;
	if (edge == 0)
	{
		if (K_Cartesian2D::PtOnFiniteLine(x0, y0, x1, y1, x, y)) //if P lies on edge 0
			return true;
		return !K_Cartesian2D::TwoFiniteLinesMeet(x0, y0, x1, y1, Cx, Cy, x, y, duff1, duff2); //returns true if P inside triangle, false otherwise
	}
	else if (edge == 1)
	{
		if (K_Cartesian2D::PtOnFiniteLine(x1, y1, x2, y2, x, y)) //if P lies on edge 1
			return true;
		return !K_Cartesian2D::TwoFiniteLinesMeet(x1, y1, x2, y2, Cx, Cy, x, y, duff1, duff2);
	}
	else //edge == 2
	{
		if (K_Cartesian2D::PtOnFiniteLine(x2, y2, x0, y0, x, y)) //if P lies on edge 2
			return true;
		return !K_Cartesian2D::TwoFiniteLinesMeet(x2, y2, x0, y0, Cx, Cy, x, y, duff1, duff2);
	}
}


//given three points A,B,P, returns true if P lies on the finite line AB
//Assumption: A cannot equal B
bool K_Cartesian2D::PtOnFiniteLine(double Ax, double Ay, double Bx, double By, double Px, double Py)
{
	if (K_Cartesian2D::Left(Ax, Ay, Bx, By, Px, Py) == 0) //if P lies on the infinite line joining A and B
	{
		if (abs(Ax - Bx) < abs(Ay - By)) //if AB taller than fat
		{
			if ((Ay <= Py && Py <= By) || (By <= Py && Py <= Ay)) //if Py between Ay and By
				return true;
		}
		else //if AB fatter than tall
		{
			if ((Ax <= Px && Px <= Bx) || (Bx <= Px && Px <= Ax)) //if Px between Ax and Bx
				return true;
		}
	}
	return false;
}


//Given a point P = (Px,Py), and coordinates (x[0],y[0]), (x[1],y[1]), (x[2],y[2]) etc, this function orders elements of indicesToSort by distance from P
//Eg if Q = indicesToSort[0] then (x[Q],y[Q]) is the closest coordinate to P
void K_Cartesian2D::OrderPointsByDistance(double Px, double Py, std::vector<double> & x, std::vector<double> & y, std::vector<int> & indicesToSort, int numIndicesToSort)
{
	//create vector of squared distances of indices from P
	std::vector<double> distances(numIndicesToSort);
	for (int i = 0; i < numIndicesToSort; i++)
		distances[i] = K_Cartesian2D::DistSquared(Px, Py, x[indicesToSort[i]], y[indicesToSort[i]]);
	//order indices of distances
	std::vector<int> distIndices(numIndicesToSort); //set to 0,1,2,etc, then order
	for (int i = 0; i < numIndicesToSort; i++)
		distIndices[i] = i;
	std::sort(std::begin(distIndices), std::end(distIndices), [&](int i1, int i2) { return distances[i1] < distances[i2]; });
	//order indicesToSort
	for (int i = 0; i < numIndicesToSort; i++)
		distIndices[i] = indicesToSort[distIndices[i]];
	for (int i = 0; i < numIndicesToSort; i++)
		indicesToSort[i] = distIndices[i];
}


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
void K_Cartesian2D::SetupGrid(int n, int m, std::vector<double> & x, std::vector<double> & y, std::vector<int> & grid, std::vector<int> & next, double& dx, double& dy, double& xMin, double& yMin)
{
	//min and max values
	xMin = x[0];
	yMin = y[0];
	double xMax = x[0];
	double yMax = y[0];
	for (int i = 0; i < n; i++)
	{
		if (x[i] < xMin)
			xMin = x[i];
		if (x[i] > xMax)
			xMax = x[i];
		if (y[i] < yMin)
			yMin = y[i];
		if (y[i] > yMax)
			yMax = y[i];
	}

	//make a bit wider, just in case
	xMin -= (xMax - xMin) / 100;
	xMax += (xMax - xMin) / 100;
	yMin -= (yMax - yMin) / 100;
	yMax += (yMax - yMin) / 100;
	double xSpan = xMax - xMin;
	double ySpan = yMax - yMin;

	//dimensions of cells of grid
	dx = xSpan / m; //cell height
	dy = ySpan / m; //cell width

	//place each node in its rightful cell
	int a, b; //place nodes in grid[a,b]
	for (int i = 0; i < n; i++) //for node i
	{
		a = static_cast<int>((x[i] - xMin) / dx);
		b = static_cast<int>((y[i] - yMin) / dy);
		next[i] = grid[a * m + b]; //grid[a,b]
		grid[a * m + b] = i; //grid[a,b]
	}
}


//records all nodes of grid[i,j] into record
//such that record[numRec], record[numRec+1], record[numRec+2] etc are the nodes
//numRec is updated accordingly
void K_Cartesian2D::RecordNodesInCell(std::vector<int> & record, int& numRec, std::vector<int> & grid, int m, int i, int j, std::vector<int> & next)
{
	//only proceed if i and j legal
	if (i < 0 || i >= m || j < 0 || j >= m)
		return;
	//get first node p
	int p = grid[i * m + j]; //grid[i,j]
	//while p is a node (nonnegative), add it to record
	while (p > -1)
	{
		record[numRec++] = p;
		p = next[p]; //set p to next node in grid[i,j]
	}
}

//travels vertically, recording all nodes of grid[i1,j], of grid[i2,j], and of all cells in-between
//updates record and numRec accordingly
void K_Cartesian2D::RecordVertically(std::vector<int> & record, int& numRec, std::vector<int> & grid, int m, int i1, int i2, int j, std::vector<int> & next)
{
	for (int i = std::min(i1, i2); i <= std::max(i1, i2); i++) //could make more efficient here by starting from max(0,iStart) and ending at min(m,iEnd)
		K_Cartesian2D::RecordNodesInCell(record, numRec, grid, m, i, j, next);
}

//travels horizontally, recording all nodes of grid[i,j1], of grid[i,j2], and of all cells in-between
//updates record and numRec accordingly
void K_Cartesian2D::RecordHorizontally(std::vector<int> & record, int& numRec, std::vector<int> & grid, int m, int i, int j1, int j2, std::vector<int> & next)
{
	for (int j = std::min(j1, j2); j <= std::max(j1, j2); j++)
		K_Cartesian2D::RecordNodesInCell(record, numRec, grid, m, i, j, next);
}

//finds nearest k points to point p=p(px,py), given the outputs of the function setupGrid()
std::vector<int> K_Cartesian2D::Nearest_k(double px, double py, int k, int n, int m, std::vector<double> & x, std::vector<double> & y, std::vector<int> & grid, std::vector<int> & next, double dx, double dy, double xMin, double yMin)
{
	//output vector
	std::vector<int> found(n);
	int numFound = 0;

	//find cell [a,b] of grid in which p lies
	auto a = static_cast<int>((px - xMin) / dx);
	auto b = static_cast<int>((py - yMin) / dy);

	//record points of grid[a,b] in output vector
	K_Cartesian2D::RecordNodesInCell(found, numFound, grid, m, a, b, next);

	//increase "searched square" outwards around [a,b] until we find k nodes
	//"searched square" = 1 cell, then 9 cells, then 25 cells etc
	int a0 = a;
	int a1 = a;
	int b0 = b;
	int b1 = b;
	int size = 1; //such that 2*size + 1 = length of side of searched square
	while (numFound < k)
	{
		//stretch searched square upward by recording row above it
		a0--;
		K_Cartesian2D::RecordHorizontally(found, numFound, grid, m, a0, b0, b1, next);
		//stretch searched square downward by recording row below it
		a1++;
		K_Cartesian2D::RecordHorizontally(found, numFound, grid, m, a1, b0, b1, next);
		//stretch searched square leftward by recording row to its left
		b0--;
		K_Cartesian2D::RecordVertically(found, numFound, grid, m, a0, a1, b0, next);
		//stretch searched square rightward by recording row to its right
		b1++;
		K_Cartesian2D::RecordVertically(found, numFound, grid, m, a0, a1, b1, next);
		//size of searched square increased
		size++;
	}

	//at this point found consists of at least k points which are close, but not necessarily the k closest

	//get max distance of points in found from p
	double maxDist = 0.0;
	double dist;
	for (int i = 0; i < numFound; i++)
	{
		dist = K_Cartesian2D::DistSquared(px, py, x[found[i]], y[found[i]]);
		if (maxDist < dist)
			maxDist = dist;
	}
	maxDist = sqrt(maxDist);

	//calculate number of cells required to cover circle of radius maxDist centred at p
	int numCellsHigh = static_cast<int>(maxDist / dx) + 1;
	int numCellsWide = static_cast<int>(maxDist / dy) + 1;
	int numCells = std::max(numCellsHigh, numCellsWide);

	//increase "searched square" until it covers circle of radius maxDist
	while (size < numCells)
	{
		//stretch searched square upward by recording row above it
		a0--;
		K_Cartesian2D::RecordHorizontally(found, numFound, grid, m, a0, b0, b1, next);
		//stretch searched square downward by recording row below it
		a1++;
		K_Cartesian2D::RecordHorizontally(found, numFound, grid, m, a1, b0, b1, next);
		//stretch searched square leftward by recording row to its left
		b0--;
		K_Cartesian2D::RecordVertically(found, numFound, grid, m, a0, a1, b0, next);
		//stretch searched square rightward by recording row to its right
		b1++;
		K_Cartesian2D::RecordVertically(found, numFound, grid, m, a0, a1, b1, next);
		//size of searched square increased
		size++;
	}

	//at this point found certainly contains k closest points, but possibly others

	//sort found by distance from p so that the k closest are at the front
	K_Cartesian2D::OrderPointsByDistance(px, py, x, y, found, numFound);
	return found;
}