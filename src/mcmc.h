#ifndef DMBC_MCMC_H
#define DMBC_MCMC_H

void dmbc_mcmc(double* z_chain, double* alpha_chain, double* eta_chain, double* sigma2_chain, double* lambda_chain, double* prob_chain, double* x_chain,
		double* x_ind_chain, double* accept, double* loglik, double* logprior, double* logpost, int* Dm, double* z, int* x, int* ng, double* alpha,
		double* eta, double* sigma2, double* lambda, const double* hyper_eta_a, const double* hyper_eta_b, const double* hyper_lambda, const double gamma_z,
		const double gamma_alpha, const double hyper_sigma2_a, const double hyper_sigma2_b, int totiter, int n, int p, int S, int G, int verbose);

#endif
