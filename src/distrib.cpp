// distrib.cpp

#include <climits>
#include <R.h>
#include <Rmath.h>
#include "matrix_utils.h"
#include "utils.h"
#include "tnt/jama_eig.h"
#include "tnt/jama_svd.h"

//  Probability function for a multivariate normal random variable
void dprodber(double* prob, const int* d, const double* pi, int m, int logscale){
	double prob_tmp = 0;
	if(logscale == 1){
		*prob = 0;
		for(int i = 0; i < m; i++){
			prob_tmp = d[i]*log(pi[i]) + (1 - d[i])*log(1 - pi[i]);
			if(ISNAN(prob_tmp)){
				prob_tmp = log(pow(pi[i], d[i])*pow(1 - pi[i], 1 - d[i]));
			}
			*prob += prob_tmp;
		}
	}
	else{
		*prob = 1;	
		for(int i = 0; i < m; i++){
			prob_tmp = pow(pi[i], d[i])*pow(1 - pi[i], 1 - d[i]);
			*prob *= prob_tmp;
		}
	}
}

//  Density function for a multivariate normal random variable
void dmultinorm(double* dens, const double* x, const double* mean, const double* sigma, int n, int p, int logscale){
	Array2D<double> xTNT(n, p);
	Array1D<double> meanTNT(p);
	Array2D<double> sigmaTNT(p, p);
	for(int j = 0; j < p; j++){
		for(int i = 0; i < n; i++){
			xTNT[i][j] = x[n*j + i];
		}
		meanTNT[j] = mean[j];
		for(int k = 0; k < p; k++){
			sigmaTNT[j][k] = sigma[p*k + j];
		}
	}

	double logdet = 0;

	JAMA::Eigenvalue<double> sigmaeig(sigmaTNT);
	Array1D<double> lambda(p);
	sigmaeig.getRealEigenvalues(lambda);
	for(int i = 0; i < p; i++){
		logdet += log(lambda[i]);
	}
	Array1D<double> distval = mahalanobis(xTNT, meanTNT, sigmaTNT);
	for(int i = 0; i < n; i++){
		dens[i] = -(p * log(2*M_PI) + logdet + distval[i])/2;
	}

	if(logscale == 0){
		for(int i = 0; i < n; i++){
			dens[i] = exp(dens[i]);
		}
	}
}

//  Multivariate normal random deviates
void rmultinorm(double* dev, int n, const double* mean, const double* sigma, int p){
	Array1D<double> meanTNT(p);
	Array2D<double> sigmaTNT(p, p);
	for(int j = 0; j < p; j++){
		meanTNT[j] = mean[j];
		for(int i = 0; i < p; i++){
			sigmaTNT[i][j] = sigma[i + p*j];
		}
	}

	JAMA::SVD<double> sigmaSVD(sigmaTNT);
	Array2D<double> V(p, p);
	Array2D<double> U(p, p);
	Array1D<double> d(p);
	sigmaSVD.getV(V);
	sigmaSVD.getU(U);
	sigmaSVD.getSingularValues(d);

	Array2D<double> tmp(p, p);
	for(int i = 0; i < p; i++){
		for(int j = 0; j < p; j++){
			tmp[i][j] = transpose(U)[i][j]*sqrt(d[i]);
		}
	}
	Array2D<double> sigmaSVDsqrt = transpose(matmult(V, tmp));
	Array2D<double> normdev(n, p);
	for(int j = 0; j < p; j++){
		for(int i = 0; i < n; i++){
			normdev[i][j] = rnorm(0, 1);
		}
	}
	Array2D<double> out_tmp = matmult(normdev, sigmaSVDsqrt);
	for(int j = 0; j < p; j++){
		for(int i = 0; i < n; i++){
			dev[n*j + i] = out_tmp[i][j] + meanTNT[j];
		}
	}
}

//  Density function for an inverse gamma random variable
void dinvgamma(double* dens, const double* x, const double alpha, const double beta, int n, int logscale){
	if((alpha <= 0) || (beta <= 0)){
		error("alpha (shape) and/or beta (scale) parameters in dinvgamma() need to be both strictly positive.\n");
	}

	double lbeta = log(beta);
	double lgalpha = lgammafn(alpha);
	for(int i = 0; i < n; i++){
		dens[i] = alpha * lbeta - lgalpha - (alpha + 1) * log(x[i]) - (beta/x[i]);
	}

	if(logscale == 0){
		for(int i = 0; i < n; i++){
			dens[i] = exp(dens[i]);
		}
	}
}

//  Inverse gamma random deviates
void rinvgamma(double* dev, int n, const double alpha, const double beta){
	if((alpha <= 0) || (beta <= 0)){
		error("alpha (shape) and/or beta (scale) parameters in rinvgamma() need to be both strictly positive.\n");
	}

	for(int i = 0; i < n; i++){
		dev[i] = 1.0/rgamma(alpha, 1.0/beta);
	}
}

//  Density function for a Dirichlet random variable
void ddirichlet(double* dens, const double* x, const double* par, int n, int p, int logscale){
	double tmp = 0;
	Array2D<double> xTNT(n, p);
	Array1D<double> parTNT(p);
	for(int j = 0; j < p; j++){
		for(int i = 0; i < n; i++){
			xTNT[i][j] = x[n*j + i];
		}
		parTNT[j] = par[j];
	}

	tmp = 0;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < p; j++){
			tmp += xTNT[i][j];
			if(xTNT[i][j] < 0 || xTNT[i][j] > 1){
				error("Some elements of x outside the [0, 1] range.\n");
			}
		}
		if(std::fabs(tmp - 1.0) >= sqrt(std::numeric_limits<double>::epsilon())){
			error("Some rows of x sum to a value different from 1.\n");
		}
	}

	double logD = 0, salpha = 0, s = 0;
	for(int j = 0; j < p; j++){
		logD += lgammafn(par[j]);
		salpha += par[j];
	}
	logD = logD - lgammafn(salpha);
	for(int i = 0; i < n; i++){
		s = 0;
		for(int j = 0; j < p; j++){
			s += (par[j] - 1)*log(xTNT[i][j]);
		}
		dens[i] = s - logD;
		if(logscale == 0){
			dens[i] = exp(dens[i]);
		}
	}
}

//  Dirichlet random deviates
void rdirichlet(double* dev, int n, const double* par, int p){
	Array2D<double> z(n, p, 0.);
	Array1D<double> s(n, 0.);

	for(int j = 0; j < p; j++){
		for(int i = 0; i < n; i++){
			z[i][j] = rgamma(par[j], 1.0);
			s[i] += z[i][j];
		}
	}
	for(int i = 0; i < n; i++){
		for(int j = 0; j < p; j++){
			dev[n*j + i] = z[i][j]/s[i];
		}
	}
}

//  Loglikelihood of the model (contribution to the loglikelihood from a single cluster)
//  (d is a vector which stacks the dissimilarities by individual)
void loglik_rbmds(double* loglik, const int* d, const double* z, const double alpha, int n, int p, int S){
	int m = n*(n - 1)/2;

	double* delta = new double[m];
	dist(delta, z, n, p);

	int* diss = new int[m];
	for(int k = 0; k < m; k++){
		diss[k] = 0;
		for(int s = 0; s < S; s++){
			diss[k] += d[k + m*s];
		}
		delta[k] = alpha + delta[k];
	}

	double* pi = new double[m];
	expit(pi, delta, m);

	*loglik = 0;
	for(int k = 0; k < m; k++){
		*loglik += log(pow(pi[k], diss[k])*pow(1 - pi[k], S - diss[k]));
	}

	delete[] delta;
	delete[] diss;
	delete[] pi;
}

//  Overall loglikelihood of the model
//  (d is a vector which stacks the dissimilarities by individual)
void loglik_dmbc(double* loglik, const int* d, const double* z, const double* alpha, const double* lambda, const int* x, int n, int p, int S, int G){
	int m = n*(n - 1)/2;
	int ind = 0;

	int* ng = new int[G];
	for(int g = 0; g < G; g++){
		ng[g] = 0;
	}
	tableC(ng, x, S, G);

	double* z_g = new double[n*p];
	int* D_g = new int[m*S];
	double* log_f_g = new double[G];
	*loglik = 0;
	for(int g = 0; g < G; g++){
		for(int j = 0; j < p; j++){
			for(int i = 0; i < n; i++){
				z_g[i + n*j] = z[i + n*j + n*p*g];
			}
		}
		for(int s = 0; s < S; s++){
			if(x[s] == (g + 1)){
				for(int k = 0; k < m; k++){
					D_g[ind] = d[k + m*s];
					ind++;
				}
			}
		}
		loglik_rbmds(&log_f_g[g], D_g, z_g, alpha[g], n, p, ng[g]);
		*loglik += ng[g]*log(lambda[g]) + log_f_g[g];
		ind = 0;
	}

	delete[] ng;
	delete[] z_g;
	delete[] D_g;
	delete[] log_f_g;
}
