#ifndef DMBC_DISTRIB_H
#define DMBC_DISTRIB_H

void dprodber(double* prob, const int* d, const double* pi, int m, int logscale);
void dmultinorm(double* dens, const double* x, const double* mean, const double* sigma, int n, int p, int logscale = 0);
void rmultinorm(double* dev, int n, const double* mean, const double* sigma, int p);
void dinvgamma(double* dens, const double* x, const double alpha, const double beta, int n, int logscale);
void rinvgamma(double* dev, int n, const double alpha, const double beta);
void ddirichlet(double* dens, const double* x, const double* par, int n, int p, int logscale);
void rdirichlet(double* dev, int n, const double* par, int p);
void loglik_rbmds(double* loglik, const int* d, const double* z, const double alpha, int nr, int nc, int S);
void loglik_dmbc(double* loglik, const int* d, const double* z, const double* alpha, const double* lambda, const int* x, int nr, int nc, int S, int G);

#endif
