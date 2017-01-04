#ifndef DMBC_UTILS_H
#define DMBC_UTILS_H

void logit(double* res, const double* p, int n);
void expit(double* res, const double* x, int n);
double euclidean(const double *x, int nr, int nc, int i1, int i2);
void dist(double* d, const double* x, int nr, int nc);
bool any_na_nan(const double* x, int n);
void sample_no_rep(int n, double* p, int* perm, int nans, int* ans);
void tableC(int* counts, const int* x, int nelem, int ndistelem);
int factorial(int x);
void permutations(int* perm, int n, int nperm, int byrow);
void which_min(int* ans, const double* r, int n);

#endif
