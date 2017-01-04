dprodber <- function(x, probs, log.p = FALSE) {
	if (any(probs <= 0) | any(probs >= 1))
		stop("elements of probs must be strictly positive.")
	if (any(!(x %in% c(0, 1))))
		stop("elements of x must be either 0 or 1.")
	if (length(x) != length(probs))
		stop("vectors x and probs must be of the same length.")
	return(ifelse(log.p, sum(x*log(probs) + (1 - x)*log(1 - probs)), prod((probs^x)*((1 - probs)^(1 - x)))))
}

dmultinorm <- function(x, mean, sigma, log = FALSE) {
	if (is.vector(x)) {
		x <- matrix(x, ncol = length(x))
	}
	if (missing(mean)) {
		mean <- rep(0, length = ncol(x))
	}
	if (missing(sigma)) {
		sigma <- diag(ncol(x))
	}
	if (NCOL(x) != NCOL(sigma)) {
		stop("x and sigma have non-conforming sizes.")
	}
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
		stop("sigma must be a symmetric matrix.")
	}
	if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming sizes.")
	}
	distval <- mahalanobis(x, center = mean, cov = sigma)
	logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
	logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
	if (log) {
		return(logretval)
	} else {
		return(exp(logretval))
	}
}

rmultinorm <- function(n, mean, sigma) {
	if (missing(mean) | missing(sigma)) {
		stop("both the mean vector and the covariance matrix must be specified.")
	}
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
		stop("sigma must be a symmetric matrix.")
	}
	if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming sizes.")
	}

	vs <- La.svd(sigma)
	vsqrt <- t(t(vs$vt) %*% (t(vs$u) * sqrt(vs$d)))
	
	p <- ncol(sigma)
	out <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
	out <- sweep(out, 2, mean, "+")
	dimnames(out) <- list(NULL, dimnames(sigma)[[2]])
	
	return(drop(out))
}

dinvgamma <- function(x, alpha, beta = 1, log = FALSE) {
	if ((alpha <= 0) | (beta <= 0)) {
		stop("alpha (shape) and/or beta (scale) parameter negative in dinvgamma().\n")
	}
	log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
	if (log) {
		return(log.density)
	} else return(exp(log.density))
}

rinvgamma <- function(n, alpha, beta = 1) {
	if ((alpha <= 0) | (beta <= 0)) {
		stop("alpha (shape) and/or beta (scale) parameter negative in rinvgamma().\n")
	}
	return(1/rgamma(n, shape = alpha, rate = beta))
}

ddirichlet <- function(x, alpha) {
	dirichlet1 <- function(x, alpha) {
		logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
		s <- sum((alpha - 1) * log(x))
		exp(sum(s) - logD)
	}
	if (!is.matrix(x)) 
		if (is.data.frame(x)) 
			x <- as.matrix(x)
		else x <- t(x)
	if (!is.matrix(alpha)) 
		alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), byrow = TRUE)

	if (any(alpha <= 0))
		stop("the elements of the alpha vector must be strictly positive.")
	if (any(dim(x) != dim(alpha))) 
		stop("mismatch between dimensions of x and alpha.")

	pd <- vector(length = nrow(x))
	for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, ])
	pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
	pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
	pd
}

rdirichlet <- function (n, alpha) {
	if (any(alpha <= 0))
		stop("the elements of the alpha vector must be strictly positive.")

	k = length(alpha)
	z = array(0, dim = c(n, k))
	s = array(0, dim = c(n, 1))
	for (i in 1:k) {
		z[, i] = rgamma(n, shape = alpha[i])
		s = s + z[, i]
	}
	for (i in 1:k) {
		z[, i] = z[, i]/s
	}
	return(z)
}

logit <- function(p) {
	return(log(p/(1 - p)))
}

expit <- function(x) {
	return(1/(1 + exp(-x)))
}

list2matrix <- function(D) {
	return(t(sapply(D, as.numeric)))
}

list2array <- function(D) {
	S <- length(D)
	n <- nrow(as.matrix(D[[1]]))
	p <- ncol(as.matrix(D[[1]]))
	
	out <- array(NA, dim = c(n, p, S))
	for (s in 1:S) {
		out[, , s] <- as.matrix(D[[s]])
	}
	
	return(out)
}

list.sum <- function(x) {
	m <- length(x)
	z <- x[[1]]
	if (m == 1) {
		return(z)
	}
	for (j in 2:m) {
		z <- z + x[[j]]
	}
	return(z)
}
