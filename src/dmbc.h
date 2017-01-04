#ifndef DMBC_H
#define DMBC_H

extern "C" {

SEXP DMBC_MCMC(SEXP raiD, SEXP raix, SEXP raing, SEXP radalpha, SEXP rn, SEXP rp, SEXP rG, SEXP rS, SEXP rtotiter, SEXP radZ, SEXP rgamma_z, SEXP reta,
		SEXP rgamma_alpha, SEXP rsigma2, SEXP rlambda, SEXP rhyper_eta_a, SEXP rhyper_eta_b, SEXP rhyper_sigma2_a, SEXP rhyper_sigma2_b, SEXP rhyper_lambda,
		SEXP rverbose);
SEXP RELABEL(SEXP radtheta, SEXP radz, SEXP radalpha, SEXP radeta, SEXP radsigma2, SEXP radlambda, SEXP radprob, SEXP raix_ind, SEXP rinit, SEXP rn, SEXP rp,
		SEXP rS, SEXP rM, SEXP rR, SEXP rG);
SEXP PACK_PAR(SEXP radz, SEXP radalpha, SEXP radlambda, SEXP rn, SEXP rp, SEXP rM, SEXP rG);

} // end extern "C"

#endif
