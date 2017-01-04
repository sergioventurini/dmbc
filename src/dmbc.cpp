// dmbc.cpp

#include <R.h>
#include <Rinternals.h>
#include "mcmc.h"
#include "relabel.h"

extern "C" {

SEXP DMBC_MCMC
(
	SEXP raiD,
	SEXP raix,
	SEXP raing,
	SEXP radalpha,
	SEXP rn,
	SEXP rp,
	SEXP rG,
	SEXP rS,
	SEXP rtotiter,
	SEXP radZ,
	SEXP rgamma_z,
	SEXP reta,
	SEXP rgamma_alpha,
	SEXP rsigma2,
	SEXP rlambda,
	SEXP rhyper_eta_a,
	SEXP rhyper_eta_b,
	SEXP rhyper_sigma2_a,
	SEXP rhyper_sigma2_b,
	SEXP rhyper_lambda,
	SEXP rverbose
)
{
	SEXP rAns = NULL;
	const int rAnsItems = 12;

	int n = asInteger(rn);
	int p = asInteger(rp);
	int G = asInteger(rG);
	int S = asInteger(rS);
	int totiter = asInteger(rtotiter);
	double gamma_z = asReal(rgamma_z);
	double gamma_alpha = asReal(rgamma_alpha);
	double hyper_sigma2_a = asReal(rhyper_sigma2_a);
	double hyper_sigma2_b = asReal(rhyper_sigma2_b);
	int verbose = INTEGER(rverbose)[0];

	SEXP rz_chain = PROTECT(allocVector(REALSXP, (totiter*n*p*G)));
	SEXP ralpha_chain = PROTECT(allocVector(REALSXP, (totiter*G)));
	SEXP reta_chain = PROTECT(allocVector(REALSXP, (totiter*G)));
	SEXP rsigma2_chain = PROTECT(allocVector(REALSXP, (totiter*G)));
	SEXP rlambda_chain = PROTECT(allocVector(REALSXP, (totiter*G)));
	SEXP rprob_chain = PROTECT(allocVector(REALSXP, (totiter*S*G)));
	SEXP rx_chain = PROTECT(allocVector(REALSXP, (totiter*S)));
	SEXP rx_ind_chain = PROTECT(allocVector(REALSXP, (totiter*S*G)));
	SEXP raccept = PROTECT(allocVector(REALSXP, (2*G)));
	SEXP rloglik = PROTECT(allocVector(REALSXP, totiter));
	SEXP rlogprior = PROTECT(allocVector(REALSXP, totiter));
	SEXP rlogpost = PROTECT(allocVector(REALSXP, totiter));

	dmbc_mcmc(REAL(rz_chain), REAL(ralpha_chain), REAL(reta_chain), REAL(rsigma2_chain), REAL(rlambda_chain), REAL(rprob_chain), REAL(rx_chain),
		REAL(rx_ind_chain), REAL(raccept), REAL(rloglik), REAL(rlogprior), REAL(rlogpost), INTEGER(raiD), REAL(radZ), INTEGER(raix), INTEGER(raing),
		REAL(radalpha), REAL(reta), REAL(rsigma2), REAL(rlambda), REAL(rhyper_eta_a), REAL(rhyper_eta_b), REAL(rhyper_lambda), gamma_z,
		gamma_alpha, hyper_sigma2_a, hyper_sigma2_b, totiter, n, p, S, G, verbose);

	// packing results
	PROTECT(rAns = allocVector(VECSXP, rAnsItems));

	SET_VECTOR_ELT(rAns, 0, rz_chain);
	SET_VECTOR_ELT(rAns, 1, ralpha_chain);
	SET_VECTOR_ELT(rAns, 2, reta_chain);
	SET_VECTOR_ELT(rAns, 3, rsigma2_chain);
	SET_VECTOR_ELT(rAns, 4, rlambda_chain);
	SET_VECTOR_ELT(rAns, 5, rprob_chain);
	SET_VECTOR_ELT(rAns, 6, rx_chain);
	SET_VECTOR_ELT(rAns, 7, rx_ind_chain);
	SET_VECTOR_ELT(rAns, 8, raccept);
	SET_VECTOR_ELT(rAns, 9, rloglik);
	SET_VECTOR_ELT(rAns, 10, rlogprior);
	SET_VECTOR_ELT(rAns, 11, rlogpost);

	// cleanup and return
	UNPROTECT(13);  // rAns, rz_chain, ralpha_chain, reta_chain, rsigma2_chain, rlambda_chain, rprob_chain, rx_ind_chain, rx_chain,
					// raccept, rloglik, rlogprior, rlogpost

	return rAns;
}

SEXP RELABEL
(
	SEXP radtheta,
	SEXP radz,
	SEXP radalpha,
	SEXP radeta,
	SEXP radsigma2,
	SEXP radlambda,
	SEXP radprob,
	SEXP raix_ind,
	SEXP rinit,
	SEXP rn,
	SEXP rp,
	SEXP rS,
	SEXP rM,
	SEXP rR,
	SEXP rG
)
{
	SEXP rAns = NULL;
	const int rAnsItems = 8;

	int init = asInteger(rinit);
	int n = asInteger(rn);
	int p = asInteger(rp);
	int S = asInteger(rS);
	int M = asInteger(rM);
	int R = asInteger(rR);
	int G = asInteger(rG);

	relabel_theta(REAL(radtheta), REAL(radz), REAL(radalpha), REAL(radeta), REAL(radsigma2), REAL(radlambda),
			REAL(radprob), INTEGER(raix_ind), init, n, p, S, M, R, G);

	// packing results
	PROTECT(rAns = allocVector(VECSXP, rAnsItems));

	SET_VECTOR_ELT(rAns, 0, radtheta);
	SET_VECTOR_ELT(rAns, 1, radz);
	SET_VECTOR_ELT(rAns, 2, radalpha);
	SET_VECTOR_ELT(rAns, 3, radeta);
	SET_VECTOR_ELT(rAns, 4, radsigma2);
	SET_VECTOR_ELT(rAns, 5, radlambda);
	SET_VECTOR_ELT(rAns, 6, radprob);
	SET_VECTOR_ELT(rAns, 7, raix_ind);

	// cleanup and return
	UNPROTECT(1);  // rAns

	return rAns;
}

SEXP PACK_PAR
(
	SEXP radz,
	SEXP radalpha,
	SEXP radlambda,
	SEXP rn,
	SEXP rp,
	SEXP rM,
	SEXP rG
)
{
	int n = asInteger(rn);
	int p = asInteger(rp);
	int M = asInteger(rM);
	int G = asInteger(rG);
	int r = n*(n - 1)/2;

	SEXP rAns = PROTECT(allocVector(REALSXP, M*(r + 1)*G));

	pack_par(REAL(rAns), REAL(radz), REAL(radalpha), REAL(radlambda), n, p, M, G);

	// cleanup and return
	UNPROTECT(1);  // rAns

	return rAns;
}

} // end extern "C"
