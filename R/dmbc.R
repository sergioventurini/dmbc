dmbc <- function(D, p = 2, G = 3, burnin = 10000, nsim = 5000, prm.prop, prm.prior, random.start = TRUE, verbose = FALSE) {
	if (any(sapply(D, class) != "dist"))
		stop("D must be a list of 'dist' objects.")
	if (length(unique(sapply(D, length))) != 1)
		stop("the elements of D must have the same length.")
	if (any(is.na(D)))
		stop("NA values not allowed in the dissimilarity matrix D.")
	if (p < 1)
		stop("the number of latent dimensions p must be at least one.")
	if (G < 1)
		stop("the number of clusters/groups G must be at least one.")
	
	S <- length(D)
	n <- attr(D[[1]], "Size")
	m <- n*(n - 1)/2
	totiter <- burnin + nsim
	p <- as.integer(p)
	G <- as.integer(G)
	
	if (G >= S)
		stop("number of groups/clusters needs to be smaller than the number of subjects.")
	
	z.chain <- z.proc <- array(NA, dim = c(totiter, n, p, G))
	eta.chain <- alpha.chain <- sigma2.chain <- lambda.chain <- array(NA, dim = c(totiter, G))
	x.chain <- array(NA, dim = c(totiter, S))
	prob.chain <- x.ind.chain <- array(0, dim = c(totiter, S, G))
	loglik <- logprior <- logpost <- numeric(totiter)
	
	if (verbose) cat("Initialization of the algorithm...")
	
	dmbc.setup <- dmbc.init(D, p, G, random.start)
	z <- dmbc.setup$z
	x <- dmbc.setup$x
	ng <- dmbc.setup$ng
	alpha <- dmbc.setup$alpha
	eta <- dmbc.setup$eta
	sigma2 <- dmbc.setup$sigma2
	lambda <- dmbc.setup$lambda
	
	if (verbose) cat("done!\n")
	
	# set prior hyperparameters
	hyper.eta.a <- rep(1.5, G)
	if (p == 1) {
		hyper.eta.b <- .5*eta
	} else {
		hyper.eta.b <- .5*eta
	}
	hyper.sigma2.a <- prm.prior[["sigma2"]][1]
	hyper.sigma2.b <- prm.prior[["sigma2"]][2]
	hyper.lambda <- prm.prior[["lambda"]]
	
	# set proposal parameters
	gamma_z <- prm.prop[["z"]]
	gamma_alpha <- prm.prop[["alpha"]]

	# start iteration
	if (verbose) cat("Running the MCMC simulation...\n")
	
	res.mcmc <- .Call(CXX_DMBC_MCMC,
		raiD = as.integer(unlist(D)),
		raix = as.integer(x),
		raing = as.integer(ng),
		radalpha = as.double(alpha),
		rn = as.integer(n),
		rp = as.integer(p),
		rG = as.integer(G),
		rS = as.integer(S),
		rtotiter = as.integer(totiter),
		radZ = as.double(z),
		rgamma_z = as.double(gamma_z),
		reta = as.double(eta),
		rgamma_alpha = as.double(gamma_alpha),
		rsigma2 = as.double(sigma2),
		rlambda = as.double(lambda),
		rhyper_eta_a = as.double(hyper.eta.a),
		rhyper_eta_b = as.double(hyper.eta.b),
		rhyper_sigma2_a = as.double(hyper.sigma2.a),
		rhyper_sigma2_b = as.double(hyper.sigma2.b),
		rhyper_lambda = as.double(hyper.lambda),
		rverbose = as.integer(verbose)
	)

	z.chain <- array(res.mcmc[[1]], c(totiter, n, p, G))
	alpha.chain <- array(res.mcmc[[2]], c(totiter, G))
	eta.chain <- array(res.mcmc[[3]], c(totiter, G))
	sigma2.chain <- array(res.mcmc[[4]], c(totiter, G))
	lambda.chain <- array(res.mcmc[[5]], c(totiter, G))
	prob.chain <- array(res.mcmc[[6]], c(totiter, S, G))
	x.chain <- array(res.mcmc[[7]], c(totiter, S))
	x.ind.chain <- array(res.mcmc[[8]], c(totiter, S, G))
	accept <- t(array(res.mcmc[[9]], c(G, 2)))
	loglik <- as.numeric(res.mcmc[[10]])
	logprior <- as.numeric(res.mcmc[[11]])
	logpost <- as.numeric(res.mcmc[[12]])

	# post-processing:
	if (verbose) cat("Post-processing the chain:\n")

	## Procrustes transformation of Z_g
	if (verbose) cat("   - applying Procrustes transformation...")
	for (niter in 1:totiter) {
		for (g in 1:G) {
			if (p == 1) {
				z.proc[niter, , , g] <- as.numeric(procrustes(as.matrix(z.chain[niter, , , g]), as.matrix(z.chain[totiter, , , g]), 
								translation = TRUE, dilation = FALSE)$X.new)
			} else {
				z.proc[niter, , , g] <- procrustes(z.chain[niter, , , g], z.chain[totiter, , , g], translation = TRUE, dilation = FALSE)$X.new
			}
		}
	}
	if (verbose) cat("done!\n")

	# relabeling of the model parameters
	if (G > 1) {
		if (totiter > 10) {
			if (verbose) cat("   - relabeling the parameters chain...")
			init <- ifelse(totiter <= 100, 5, 100)
			
			theta <- .Call(CXX_PACK_PAR,
				radz = as.double(z.proc),
				radalpha = as.double(alpha.chain),
				radlambda = as.double(lambda.chain),
				rn = as.integer(n),
				rp = as.integer(p),
				rM = as.integer(totiter),
				rG = as.integer(G)
			)

			theta.relab <- .Call(CXX_RELABEL,
				radtheta = as.double(theta),
				radz = as.double(z.proc),
				radalpha = as.double(alpha.chain),
				radeta = as.double(eta.chain),
				radsigma2 = as.double(sigma2.chain),
				radlambda = as.double(lambda.chain),
				radprob = as.double(prob.chain),
				raix_ind = as.integer(x.ind.chain),
				rinit = as.integer(init),
				rn = as.integer(n),
				rp = as.integer(p),
				rS = as.integer(S),
				rM = as.integer(totiter),
				rR = as.integer(m + 1),
				rG = as.integer(G)
			)

			theta <- array(theta.relab[[1]], c(totiter, (m + 1), G))  # this is not needed elsewhere
			z.proc <- array(theta.relab[[2]], c(totiter, n, p, G))
			alpha.chain <- array(theta.relab[[3]], c(totiter, G))
			eta.chain <- array(theta.relab[[4]], c(totiter, G))
			sigma2.chain <- array(theta.relab[[5]], c(totiter, G))
			lambda.chain <- array(theta.relab[[6]], c(totiter, G))
			prob.chain <- array(theta.relab[[7]], c(totiter, S, G))
			x.ind.chain <- array(theta.relab[[8]], c(totiter, S, G))
			x.chain <- t(apply(x.ind.chain, 1, function(x) as.integer(x %*% 1:G)))

			if (verbose) cat("done!\n")
		} else {
			warning("the number of iterations is too small for relabeling; relabeling skipped.", call. = FALSE, immediate. = TRUE)
		}
	}
	
	out <- new("dmbc",
		z.chain = z.chain,
		z.chain.p = z.proc,
		alpha.chain = alpha.chain,
		eta.chain = eta.chain,
		sigma2.chain = sigma2.chain,
		lambda.chain = lambda.chain,
		prob.chain = prob.chain,
		x.ind.chain = x.ind.chain,
		x.chain = x.chain,
		accept = accept,
		obsdiss = D,
		dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
		control = list(burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = prm.prior, hyper.eta.a = hyper.eta.a, hyper.eta.b = hyper.eta.b),
		dim = list(n = n, p = p, G = G, S = S)
	)
	return(out)
}

dmbc.init <- function(D, p, G, random.start) {
	S <- length(D)
	n <- attr(D[[1]], "Size")

	# initialize x (cluster labels)
	if (random.start) {
		x <- sample(1:G, S, replace = TRUE)
		while (length(unique(x)) < G) {
			x <- sample(1:G, S, replace = TRUE)
		}
	} else {
		dmat <- list2matrix(D)
		d.clust <- hclust(dist(dmat, method = "manhattan"), method = "ward")
		x <- cutree(d.clust, k = G)
	}
	ng <- table(factor(x, levels = 1:G))
	
	Dm <- list2array(D)
	z <- array(NA, dim = c(n, p, G))
	alpha <- eta <- sigma2 <- numeric(G)

	for (g in 1:G) {
		# initialize Z_g
		Dg <- Dm[, , x == g]
		Dm_bar <- apply(Dg, c(1, 2), mean)
		D_bar <- as.dist(Dm_bar)
		d_ov.mean <- mean(D_bar)
		Dm_above <- as.matrix(as.dist(Dm_bar > d_ov.mean))
		
		Dg <- Dm[, , x == g]
		Dm_sum <- apply(Dg, c(1, 2), sum)
		z_mds <- cmdscale(d = Dm_sum, k = p)
		z[, , g] <- scale(z_mds)
		
		# initialize alpha_g
		alpha.glm <- glm(as.numeric(as.dist(Dm_above)) ~ 1, family = "binomial", offset = as.numeric(dist(z[, , g])))
		alpha[g] <- coef(alpha.glm)
		
		# initialize sigma2_g
		sigma2[g] <- vcov(alpha.glm)
	}
	
	# initialize eta
	if (p == 1) {
		eta <- apply(z, 3, var)
	} else {
		eta <- apply(z, 3, function(x) mean(diag(cov(x))))
	}
	
	# initialize lambda
	lambda <- ng/S

	return(list(z = z, x = x, ng = ng, alpha = alpha, eta = eta, sigma2 = sigma2, lambda = lambda))
}

dmbc.logLik.rbmds <- function(D, Z, alpha) {
	S <- length(D)
	diss <- list.sum(D)
	delta <- dist(Z)
	pi <- expit(alpha + delta)
	# out <- sum(diss*log(pi) + (S - diss)*log(1 - pi))
	out <- sum(log((pi^diss)*((1 - pi)^(S - diss))))
	
	return(out)
}

dmbc.logLik <- function(D, Z, alpha, lambda, x) {
	G <- length(alpha)
	ng <- as.numeric(table(factor(x, levels = 1:G)))
	logfg <- numeric(G)
	ll <- 0
	for (g in 1:G) {
		if (ng[g] > 0) {
			xg <- which(x == g)
			logfg[g] <- dmbc.logLik.rbmds(D[xg], Z[, , g], alpha[g])
			ll <- ll + ng[g]*log(lambda[g]) + logfg[g]
		}
	}

	return(ll)
}
