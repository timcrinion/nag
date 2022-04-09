//This is the header file for K_Thiessen.cpp
#pragma once
#include <vector>

class K_Thiessen
{
private:

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
	static void TRIGRD(int n, std::vector<double>& x, std::vector<double>& y, std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts);
	//insert boundary point qk just after qi
	//boundaryNext and boundaryPrev are linked lists, for example
	//boundaryNext[11] = 17
	//boundaryPrev[17] = 11
	//means 17 is immediately after node 11 if you travel around the boundary anticlockwise
	static void InsertBoundaryPoint(int qi, int qk, std::vector<int>& boundaryNext, std::vector<int>& boundaryPrev);
	//delete boundary point qi
	//boundaryNext and boundaryPrev are linked lists, for example
	//boundaryNext[11] = 17
	//boundaryPrev[17] = 11
	//means 17 is immediately after node 11 if you travel around the boundary anticlockwise
	static void DeleteBoundaryPoint(int qi, std::vector<int>& boundaryNext, std::vector<int>& boundaryPrev);
	//returns index of edge of triangle joining qi and qj such that TriAdj[triangle,edge] = other triangle sharing both qi and qj
	static int GetEdge(int triangle, int qi, int qj, std::vector<int>& TriPts);
	//sets t to triangle containing p = (x,y)
	//sets edge to that of t facing p
	//this function will only work if boundary is always convex (or at least if p is in a convex subsection of the triangulation containing first guess t=0)
	//returns true if t contains p, false if no triangle does (however, edge of t will still be on boundary facing p)
	static bool TRFIND(int& edge, int& t, double x, double y, std::vector<double>& X, std::vector<double>& Y, std::vector<int>& TriAdj, std::vector<int>& TriPts);
	//function that returns a triangle t with edge connecting points a and b
	//If a and b are on the boundary of the triangulation, there is only one possible output
	//Ta is any triangle connected to node a
	static int TriWithEdge(int a, int b, int Ta, std::vector<int>& TriAdj, std::vector<int>& TriPts);
	static void Optimise(int T1, int T2, std::vector<double>& X, std::vector<double>& Y, std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts);
	static void GetABCD(int T1, int T2, int& A, int& B, int& C, int& D, std::vector<int>& TriAdj, std::vector<int>& TriPts);
	//Given Table (set to either TriAdj or TriPts) this function replaces the value oldVal in triangle t's row t with replacement value newVal
	static void ReplacePt(std::vector<int>& TriTable, int t, int oldVal, int newVal);
	static bool ShouldSwap(double xA, double yA, double xB, double yB, double xC, double yC, double xD, double yD);
	static int PtIndex(int p, int t, std::vector<int>& TriPts);

};
