dmbc.get.postmean <- function(res) {
	burnin <- res@control$burnin
	nsim <- res@control$nsim
	totiter <- burnin + nsim
	
	z.postmean <- apply(res@z.chain.p[(burnin + 1):totiter, , , , drop = FALSE], c(2, 3, 4), mean)
	alpha.postmean <- colMeans(res@alpha.chain[(burnin + 1):totiter, , drop = FALSE])
	eta.postmean <- colMeans(res@eta.chain[(burnin + 1):totiter, , drop = FALSE])
	sigma2.postmean <- colMeans(res@sigma2.chain[(burnin + 1):totiter, , drop = FALSE])
	lambda.postmean <- colMeans(res@lambda.chain[(burnin + 1):totiter, , drop = FALSE])
	prob.postmean <- apply(res@prob.chain[(burnin + 1):totiter, , , drop = FALSE], c(2, 3), mean)
	class.postmean <- apply(prob.postmean, 1, which.max)
	out <- list(z.postmean = z.postmean, alpha.postmean = alpha.postmean, eta.postmean = eta.postmean, sigma2.postmean = sigma2.postmean,
		lambda.postmean = lambda.postmean, prob.postmean = prob.postmean, class.postmean = class.postmean)

	return(out)
}

dmbc.get.ml <- function(res) {
	burnin <- res@control$burnin
	nsim <- res@control$nsim
	totiter <- burnin + nsim
	ll <- res@dens$loglik

	i.ml <- which.max(ll[(burnin + 1):totiter])
	out <- list(z.ml = res@z.chain[(burnin + i.ml), , , ], alpha.ml = res@alpha.chain[(burnin + i.ml), ],
		eta.ml = res@eta.chain[(burnin + i.ml), ], sigma2.ml = res@sigma2.chain[(burnin + i.ml), ],
		lambda.ml = res@lambda.chain[(burnin + i.ml), ], loglik = ll[burnin + i.ml])

	return(out)
}

dmbc.get.map <- function(res) {
	burnin <- res@control$burnin
	nsim <- res@control$nsim
	totiter <- burnin + nsim
	lpost <- res@dens$logpost

	i.map <- which.max(lpost[(burnin + 1):totiter])
	out <- list(z.map = res@z.chain[(burnin + i.map), , , ], alpha.map = res@alpha.chain[(burnin + i.map), ], 
		eta.map = res@eta.chain[(burnin + i.map), ], sigma2.map = res@sigma2.chain[(burnin + i.map), ],
		lambda.map = res@lambda.chain[(burnin + i.map), ], logpost = lpost[burnin + i.map])

	return(out)
}
