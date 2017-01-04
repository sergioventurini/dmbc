// matrix_utils.cpp

#include <R.h>
#include "tnt/jama_lu.h"

// *** Not all of the following functions are used in the package ***

// Conversion of the dissimilarities from matrix to vector
Array1D<double> dissM2V(const Array2D<double>& d){
	int n = d.dim1();
	int m = n*(n - 1)/2;
	int count = 0;
	Array1D<double> dvec(m, 0.);
	
	for(int j = 0; j < n; j++){
		for(int i = (j + 1); i < n; i++){
			dvec[count] = d[i][j];
			count++;
		}
	}

	return dvec;
}

// Matrix to vector
Array1D<double> mat2vec(const Array2D<double>& A, int j){
	int nrows = A.dim1();
	Array1D<double> v(nrows, 0.);
	
	for(int i = 0; i < nrows; i++){
		v[i] = A[i][j];
	}

	return v;
}

// Vector to matrix
Array2D<double> vec2mat(const Array2D<double>& A, int j, const Array1D<double>& v){
	int nrows = A.dim1();
	Array2D<double> B = A.copy();

	for(int i = 0; i < nrows; i++){
		B[i][j] = v[i];
	}

	return B;
}

// Matrix transposition
Array2D<double> transpose(const Array2D<double>& A){
	int nrows = A.dim1();
	int ncols = A.dim2();

	Array2D<double> B(ncols, nrows);

	for(int i = 0; i < ncols; i++){
		for(int j = 0; j < nrows; j++){
			B[i][j] = A[j][i];
		}
	}

	return B;
}

// Computation of the dot product of two vectors
double dotprod(const Array1D<double>& a, const Array1D<double>& b){
	int nrows = a.dim();
	double dot = 0;

	for(int i = 0; i < nrows; i++){
		dot += a[i]*b[i];
	}

	return dot;
}

// Matrix column sums
void colsums(double* colsums, const double* A, int nrows, int ncols){
	for(int j = 0; j < ncols; j++){
		for(int i = 0; i < nrows; i++){
			colsums[j] += A[i + nrows*j];
		}
	}
}

// Matrix row sums
void rowsums(double* rowsums, const double* A, int nrows, int ncols){
	for(int i = 0; i < nrows; i++){
		for(int j = 0; j < ncols; j++){
			rowsums[i] += A[i + nrows*j];
		}
	}
}

// Matrix column sums (TNT)
Array1D<double> colsums_TNT(const Array2D<double>& A){
	int nrows = A.dim1();
	int ncols = A.dim2();
	Array1D<double> csums(ncols, 0.);

	for(int j = 0; j < ncols; j++){
		for(int i = 0; i < nrows; i++){
			csums[j] += A[i][j];
		}
	}

	return csums;
}

// Matrix row sums (TNT)
Array1D<double> rowsums_TNT(const Array2D<double>& A){
	int nrows = A.dim1();
	int ncols = A.dim2();
	Array1D<double> rsums(nrows, 0.);

	for(int i = 0; i < nrows; i++){
		for(int j = 0; j < ncols; j++){
			rsums[i] += A[i][j];
		}
	}

	return rsums;
}

// Back-solve an upper triangular system 
Array1D<double> backward(const Array2D<double>& A, const Array1D<double>& RHS){
	int nrows = A.dim1();
	double temp;

	Array1D<double> b(nrows, 0.);

	// initialize recursion on last row of all right-hand-sides
	if(A[nrows - 1][nrows - 1] != 0){
		b[nrows - 1] = RHS[nrows - 1]/A[nrows - 1][nrows - 1];
	}
	else{
		error("Singular system!\n");
	}

	// now work through the remaining rows
	for(int i = (nrows - 2); i >= 0 ; i--){
		if(A[i][i] != 0){
			temp = RHS[i];
			for(int k = (i + 1); k < nrows; k++){
				temp -= b[k]*A[i][k]; 
			}
			b[i] = temp/A[i][i];
		}
		else{
			error("Singular system!\n");
		}
	}

	return b;
}

// Computation of the matrix inverse
Array2D<double> inverse(const Array2D<double>& A){
	int nrows = A.dim1();
	int ncols = A.dim2();

	if(nrows != ncols){
		error("Matrix inverse is available only for square matrices!\n");
	}

	// Set up identity matrix
	Array2D<double> eye(nrows, nrows, 0.);
	for(int i = 0; i < nrows; i++){
		eye[i][i] = 1;
	}
	
	JAMA::LU<double> ALU(A);
	Array2D<double> Ainv = ALU.solve(eye);

	return Ainv;
}

// Computation of the Mahalanobis distance
Array1D<double> mahalanobis(const Array2D<double>& x, const Array1D<double>& center, const Array2D<double>& cov){
  int nrows = x.dim1();
  int ncols = x.dim2();
  Array2D<double> centered = x.copy();
  Array2D<double> covinv = inverse(cov);

  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      centered[i][j] -= center[j];
    }
  }

  Array1D<double> m = rowsums_TNT(centered * matmult(centered, covinv));
  
  return m;
}
