//This C++ file mimics the NAG function f08waf from https://www.nag.com/numeric/fl/nagdoc_fl26.0/html/f08/f08waf.html

#include "K_f08waf.h"
#include "K_f07fjc.h"
#include "K_g05rzf.h"
#include "K_matrix.h"
#include "K_UpperHessenberg.h"
#include "KUtils.h"
#include "K_jacobi.h"

//for use in L:\Dev\Paddington\grclib\prshul\nagcal.f

int ntimes = 0;
int nproblem = 0;

void K_f08waf::nagcal(int NDIMEN, int IMAX, double* DA, double* DB, double* R, double* V, int* LNAGOK)
{
	ntimes++;
	std::vector<double> A_inner (NDIMEN * NDIMEN);
	std::vector<double> B_inner (NDIMEN * NDIMEN);
	std::vector<double> R_inner (NDIMEN);
	std::vector<double> V_inner (NDIMEN * NDIMEN);

	for (int i = 0; i < NDIMEN; i++)
	{
		for (int j = 0; j < NDIMEN; j++)
		{
			A_inner[i * NDIMEN + j] = DA[i * IMAX + j];
			B_inner[i * NDIMEN + j] = DB[i * IMAX + j];
			V_inner[i * NDIMEN + j] = V [i * IMAX + j];
		}
	}
	bool LNAGOK_inner = false;
	//FILE* fp = fopen ("L:\\nagcal.txt", "w");
	//fprintf (fp, "NDIMEN = %d\n", NDIMEN);
	//fprintf (fp, "IMAX = %d\n", IMAX);
	//for (int i = 0; i < NDIMEN * NDIMEN; i++)
	//	fprintf (fp, "A[%d] = %16.14e\n", i, A_inner[i]);
	//for (int i = 0; i < NDIMEN * NDIMEN; i++)
	//	fprintf (fp, "B[%d] = %16.14e\n", i, B_inner[i]);
	//fclose (fp);

	nagcal_inner (NDIMEN, IMAX, A_inner, B_inner, R_inner, V_inner, LNAGOK_inner);

	*LNAGOK = LNAGOK_inner ? 1 : 0;
	for (int i = 0; i < NDIMEN; i++)
	{
		R[i] = R_inner[i];
		for (int j = 0; j < NDIMEN; j++)
		{
			V[i * IMAX + j] = V_inner[i * NDIMEN + j];
		}
	}
}

//useful link:
//https://www.math.kth.se/na/SF2524/matber15/qrmethod.pdf
//http://pi.math.cornell.edu/~web6140/TopTenAlgorithms/QRalgorithm.html
//http://www.mosismath.com/Eigenvalues/EigenvalsQR.html

//The K_jacobi::JacobiEigenvectors function only works when the input is a symmetric matrix.
//K_f08waf works when the function is not symmetric, but not when eigenvalues are complex

//Function that outputs generalised right-eigenvectors and eigenvalues of A and B
//We say k (number) is a right generalised eigenvalue of A and B, and v (n-vector) is a right generalised eigenvector of A and B, if A * v = k * B * v
//n eigenvalues placed in R (length n), n vectors placed in V (n by n), such that ith row of V corressponds to eigenvalue R[i]
//Convention here is that the entry of the 2D matrix M (row i, column j) is stored in M[i * numColumnsM + j]
//LNAGOK set to true if function ok, set to false if eigenvalues are imaginary or if neither A nor B are invertible
void K_f08waf::nagcal_inner(int n, int IMAX [[maybe_unused]], std::vector<double>& A, std::vector<double>& B,
	std::vector<double>& R, std::vector<double>& V, bool& LNAGOK) //SUBROUTINE NAGCAL(NDIMEN, IMAX, DA, DB, R, V, LNAGOK)
{
	LNAGOK = true;

	//want right eigenvectors

	//try inverse of A, or inverse of B
	std::vector<double> copyA = A; //copy of A to be discarded
	std::vector<double> copyB = B; //copy of B to be discarded
	std::vector<double> invA; //will be inverse of A
	std::vector<double> invB; //will be inverse of B

	//scale copies of A and B in case their entries are huge or tiny
	//this does not change the output eigenvectors
	//however, the eigenvalues of (copyA,copyB) will be frobB/frobA times those of (A,B). Need to undo this before exiting function
	double frobA = K_matrix::FrobeniusNorm(A, n * n);
	double frobB = K_matrix::FrobeniusNorm(B, n * n);
	double invFrobA = 1.0 / frobA;
	double invFrobB = 1.0 / frobB;
	for (int i = 0; i < n * n; i++)
	{
		copyA[i] *= invFrobA;
		copyB[i] *= invFrobB;
	}
	std::vector<double> copyB2 = copyB; //2nd copy of B



	//If B is invertible then Av = kBv means that (invB * A)v = kv. This means calculating eigenvalues and eigenvectors of C = invB * A
	if (K_f07fjc::GaussJordanSquareInverse(copyB, invB, n))
	{
		//C = invB * A
		std::vector<double> C;
		C = K_matrix::mult_ignoreZeros(invB, copyA, n, n, n, 0.0);

		//fill R with eigenvalues and V with eigenvectors of C (eigenvectors are columns of V)
		LNAGOK = K_f08waf::QR_algorithm(C, n, V, R);
	}
	else if (K_f07fjc::GaussJordanSquareInverse(copyA, invA, n)) //if A invertible, then Av = kBv means that (invA * B) v = v/k
	{
		//C = invA * B
		std::vector<double> C;
		C = K_matrix::mult_ignoreZeros(invA, copyB2, n, n, n, 0.0);

		//fill R with eigenvalues and V with eigenvectors of C (eigenvectors are columns of V)
		LNAGOK = K_f08waf::QR_algorithm(C, n, V, R);

		//if that didn't work, return
		if (!LNAGOK)
			return;

		//if Av = kBv then v/k = Cv
		//Therefore invert eigenvalues of C to get generalised eigenvalues of A and B
		for (int i = 0; i < n; i++)
		{
				if (abs(R[i]) < 0.00000000001) // if R[i] zero
			{
				LNAGOK = false;
				return;
			}
			R[i] = 1.0 / R[i];
		}
	}
	else //if neither A nor B invertible, fail
		LNAGOK = false;

	//since we scaled A and B at the start of the function, need to correct eigenvalues
	for (int i = 0; i < n; i++)
		R[i] *= frobA / frobB;


	std::vector<double> AV = K_matrix::mult(A, n, n, V, n);
	std::vector<double> BV = K_matrix::mult(B, n, n, V, n);
	for (int i = 0; i < n * n; i++)
		AV[i] = AV[i] / BV[i];

	//want eigenvectors to be rows, not columns:
	K_matrix::transposeSquareMatrix(V, n);

	return;
}



//returns true if n by n matrix M is symmetric, false otherwise
bool K_f08waf::isSymmetric(std::vector<double>& M, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (!KUtils::AreDoublesEqual(M[i * n + j], M[j * n + i])) //if Mij not equal to Mji
				return false;
		}
	}
	return true;
}






//calculate eigenvectors and eigenvalues of n by n matrix A
//eigVals length n
//eigVecs n by n, such that 
//WARNING: This algorithm is fine for small matrices but inefficient for larger ones
//returns true if eigenvalues real, false if complex
bool K_f08waf::QR_algorithm(std::vector<double>& A, int n, std::vector<double>& eigVecs, std::vector<double>& eigVals)
{

	//want upper-hessenberg matrix H similar to A and unitary matrix U such that H = U * A * U'
	std::vector<double> U(n * n); // n by n
	std::vector<double> H(n * n); // n by n
	K_UpperHessenberg::hessenbergReduction(A, n, U, H);


	//matrix B will always be similar to U, which is similar to A
	std::vector<double> B = H;

	//will be QR decompositions
	std::vector<double> Q(n * n);
	std::vector<double> R(n * n);

	//create unitary matrix QQQ such that B = QQQ' * A * QQQ
	std::vector<double> QQQ(n * n);
	K_UpperHessenberg::eye(QQQ, n);


	//do B=QR, then B=RQ, again and again (how long?)
	std::vector<double> RQ(n * n);
	std::vector<double> prevB(n*n);
	for (int i = 0; i < 1000; i++)
	{
		//fill Q and R such that B = QR
		K_matrix::QR(B, n, n, Q, R);

		//set B to RQ (= Q' B Q), after backing up B
		for (int j = 0; j < n*n; j++)
			prevB[j] = B[j];
		K_f08waf::setAtoBC(B, R, Q, n);

		//update QQQ such that B = QQQ' * A * QQQ
		K_UpperHessenberg::rightMult(QQQ, Q, n);

		//If B not changed, break
		if (K_f08waf::vecDiff(B, prevB, n*n) < 0.000000000000001)
			break;
	}

	//if B is upper triangular, then no complex eigenvalues
	if (!K_f08waf::isUppTri(B, n))
		return false;

	//get eigenvalues of B
	for (int i = 0; i < n; i++)
		eigVals[i] = B[i * n + i]; //R not B?

	//get eigenvectors of B
	K_f08waf::eigUppTri(B, n, eigVecs);

	//since A, H and B are similar, their eigenvalues are the same
	//since B = QQQ' H QQQ it follows that if v is an eigenvector of B then QQQv is an eigenvector of H
	//since H = U A U' it follows that if v is an eigenvector of H then U'v is an eigenvector of A
	//therefore set eigVecs = U' * QQQ * eigVecs so that they are eigenvectors of A
	K_UpperHessenberg::leftMult(QQQ, eigVecs, n);
	K_matrix::transposeSquareMatrix(U, n);
	K_UpperHessenberg::leftMult(U, eigVecs, n);

	//normalise eigenvectors before output
	K_f08waf::normalizeColumns(eigVecs, n);

	//all eigenvalues real
	return true;
}

//returns sum of differences between A and B
double K_f08waf::vecDiff(std::vector<double>& A, std::vector<double> & B, int length)
{
	double diff = 0.0;
	for (int i = 0; i < length; i++)
		diff += abs(A[i] - B[i]);
	return diff;
}


//Given three existing n by n matrices A,B,C, set A to B*C 
void K_f08waf::setAtoBC(std::vector<double> & A, std::vector<double> & B, std::vector<double> & C, int n)
{
	std::vector<double> BC;
	BC = K_matrix::mult_ignoreZeros(B, C, n, n, n, 0.0);
	for (int i = 0; i < n * n; i++)
		A[i] = BC[i];
}



//get eigenvectors of upper triangular n by n matrix U
//eigenvalues are diagonal values Uii
//place eigenvectors in eigVecs such that jth column corressonds to jth eigenvalue Ujj
void K_f08waf::eigUppTri(std::vector<double> & U, int n, std::vector<double> & eigVecs)
{
	//initialise jth eigenvector, to be changed in each for loop
	std::vector<double> vecj(n);
	for (int j = n - 1; j >= 0; j--) // n-1, n-2 ... 0
	{
		//fill jth eigenvector
		K_f08waf::eigUppTri_colj(U, n, vecj, j);
		//set as column j of eigVecs
		K_f07fjc::setColumn(eigVecs, vecj, n, n, j);
	}
}


//function that calculates eigenvecter corressponding to jth eigenvalue of n by n matrix U
//jth eigenvalue is Ujj
//answer stored in eigVec
void K_f08waf::eigUppTri_colj(std::vector<double> & U, int n, std::vector<double> & eigVec, int j)
{
	double eigVal = U[j * n + j]; // Ujj
	//set eigVec[j] to 1 and all later entries to zero
	eigVec[j] = 1.0;
	for (int i = j + 1; i < n; i++)
		eigVec[i] = 0.0;
	//check that no earlier diagonals are set to Ujj. This is important to avoid division by zero later
	for (int i = 0; i < j; i++)
	{
		if (KUtils::AreDoublesEqual(U[i * n + i], eigVal)) // if Uii = Ujj
		{ // if this happens, jth eigenvector is the same as ith eigenvector
			K_f08waf::eigUppTri_colj(U, n, eigVec, i);
			return;
		}
	}
	//if we reach this point then no earlier eigenvalues are equal to Ujj
	//fill entry j-1, j-2 ... 0
	for (int i = j - 1; i >= 0; i--)
	{
		double sum = 0.0; //want sum = U[i,i+1]v[i+1] + U[i,i+2]v[i+2] + ... + U[i,j]v[j]
		for (int k = i + 1; k <= j; k++)
			sum += U[i * n + k] * eigVec[k]; //incremement by Uik * vk
		eigVec[i] = sum / (eigVal - U[i * n + i]); //divide sum by (Ujj-Uii)
	}
}


//returns true if n by n matrix M is upper triangular, false otherwise
bool K_f08waf::isUppTri(std::vector<double> & M, int n)
{
	for (int j = 0; j < n; j++) // column j = 0, 1 ... n-1
	{
		for (int i = j + 1; i < n; i++) // row i = j+1, j+2 ... n-1
		{
			if (!KUtils::AreDoublesEqual(M[i * n + j], 0.0)) // if Mij nonzero
				return false;
		}
	}
	return true;
}


//scales each column of n by n matrix M such that each has norm 1
//if a column is all zeros, it is left untouched
void K_f08waf::normalizeColumns(std::vector<double> & M, int n)
{
	for (int j = 0; j < n; j++) //column j
	{
		double norm = 0.0;
		double Mij;
		for (int i = 0; i < n; i++) //row i
		{
			Mij = M[i * n + j]; //increase norm by Mij * Mij
			norm += Mij * Mij;
		}
		norm = sqrt(norm);
		for (int i = 0; i < n; i++) //row i
			M[i * n + j] /= norm; //divide Mij by norm
	}
}


//Function that outputs generalised right-eigenvectors and eigenvalues of n by n symmetric matrices A and B
//Only works if all eigenvalues of B are strictly greater than zero. Returns true if this is the case, false otherwise
//We say k (number) is a right generalised eigenvalue of A and B, and v (n-vector) is a right generalised eigenvector of A and B, if A * v = k * B * v
//n eigenvalues placed in eigVals (length n), n vectors placed in eigVecs (n by n), such that ith column of eigVecs corressponds to eigenvalue eigVec[i]
bool K_f08waf::genEigSym(std::vector<double>& A, std::vector<double>& B, int n, std::vector<double>& eigVals, std::vector<double>& eigVecs)
{
	//copied from http://fourier.eng.hmc.edu/e161/lectures/algebra/node7.html

	//eigens of B
	std::vector<double> eigValsB(n);
	std::vector<double> eigVecsB(n * n);
	K_jacobi::JacobiEigenvectors2(eigValsB, eigVecsB, B, n);

	//diagonal of inverse squares of B's eigenvalues
	std::vector<double> invSqrLambdaB(n * n, 0.0); // Lambda_B ^ -1/2
	for (int i = 0; i < n; i++)
	{
		double s = sqrt(abs(eigValsB[i]));
		if (eigValsB[i] < 0.00000000001)
		{
			nproblem++;
			return false;
		}
		invSqrLambdaB[i * n + i] = 1.0 / s;
	}

	std::vector<double> PsiB_; //Psi subscript B apostrophe
	PsiB_ = K_matrix::mult_ignoreZeros(eigVecsB, invSqrLambdaB, n, n, n, 0.0);

	std::vector<double> PsiB_T = PsiB_; //transpose of PsiB_
	K_matrix::transposeSquareMatrix(PsiB_T, n);

	std::vector<double> A_; //A apostrophe
	A_ = K_matrix::mult3matrices(PsiB_T, A, PsiB_, n, n, n, n);

	//eigens of A_
	std::vector<double> eigValsA_(n);
	std::vector<double> PsiA(n * n);
	K_jacobi::JacobiEigenvectors2(eigValsA_, PsiA, A_, n);

	std::vector<double> PsiAT = PsiA;
	K_matrix::transposeSquareMatrix(PsiAT, n);

	std::vector<double> Lambda;
	Lambda = K_matrix::mult3matrices(PsiAT, A_, PsiA, n, n, n, n);
	std::vector<double> Psi;
	Psi = K_matrix::mult3matrices(eigVecsB, invSqrLambdaB, PsiA, n, n, n, n);

	for (int i = 0; i < n; i++)
		eigVals[i] = Lambda[i * n + i];


	for (int i = 0; i < n * n; i++)
		eigVecs[i] = Psi[i];

	return true;
}


// Function that outputs generalised right - eigenvectors and eigenvalues of n by n symmetric matrices A and B
//Only works if all eigenvalues of B are strictly greater than zero. Returns true if this is the case, false otherwise
//    https://www.netlib.org/lapack/lug/node54.html
bool K_f08waf::genEigSym2(std::vector<double>& A, std::vector<double>& B, int n, std::vector<double>& eigVals, std::vector<double>& eigVecs)
{
	//do cholesky composition to get n by n lower triangular matrix L such tha B = LL^T
	bool choleskySuccess;
	std::vector<double> L = K_g05rzf::CholeskyDecomposition_vec(B, n, &choleskySuccess);
	if (!choleskySuccess)
		return false;

	//invert L
	bool L_invertible;
	std::vector<double> invL;
	invL = K_f07fjc::inverseLowerTriang(L, n, &L_invertible);
	if (!L_invertible)
		return false;

	//transpose invL
	std::vector<double> invLT = invL;
	K_matrix::transposeSquareMatrix(invLT, n);

	//set C = invL * A * invLT
	//C will be symmetric
	std::vector<double> C;
	C = K_matrix::mult3matrices(invL, A, invLT, n, n, n, n);

	//get eigenvectors and eigenvalues of C
	std::vector<double> eigVecsC(n * n);
	K_jacobi::JacobiEigenvectors2(eigVals, eigVecsC, C, n);

	//set eigVecs = invLT * eigVecsC
	K_f08waf::setAtoBC(eigVecs, invLT, eigVecsC, n);

	return true;
}