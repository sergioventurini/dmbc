#ifndef DMBC_RELABEL_H
#define DMBC_RELABEL_H

void relabel_theta(double* theta, double* z_chain, double* alpha_chain, double* eta_chain, double* sigma2_chain, double* lambda_chain,
		double* prob_chain, int* x_ind_chain, int init, int n, int p, int S, int M, int R, int G);
void pack_par(double* theta, const double* z, const double* alpha, const double* lambda, int n, int p, int M, int G);

#endif
