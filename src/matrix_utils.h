#ifndef DMBC_MATRIX_UTILS_H
#define DMBC_MATRIX_UTILS_H

#include "tnt/tnt.h"

TNT::Array1D<double> dissM2V(const TNT::Array2D<double>& d);
TNT::Array1D<double> mat2vec(const TNT::Array2D<double>& A, int j);
TNT::Array2D<double> vec2mat(const TNT::Array2D<double>& A, int j, const TNT::Array1D<double>& v);
TNT::Array2D<double> transpose(const TNT::Array2D<double>& A);
double dotprod(const TNT::Array1D<double>& a, const TNT::Array1D<double>& b);
void colsums(double* colsums, const double* A, int nrows, int ncols);
void rowsums(double* rowsums, const double* A, int nrows, int ncols);
TNT::Array1D<double> colsums_TNT(const TNT::Array2D<double>& A);
TNT::Array1D<double> rowsums_TNT(const TNT::Array2D<double>& A);
TNT::Array1D<double> backward(const TNT::Array2D<double>& A, const TNT::Array1D<double>& RHS);
TNT::Array2D<double> inverse(const TNT::Array2D<double>& A);
TNT::Array1D<double> mahalanobis(const TNT::Array2D<double>& x, const TNT::Array1D<double>& center, const TNT::Array2D<double>& cov);

#endif
