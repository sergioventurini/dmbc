setClass(Class = "dmbc",
	representation(
		z.chain = "array",
		z.chain.p = "array",
		alpha.chain = "matrix",
		eta.chain = "matrix",
		sigma2.chain = "matrix",
		lambda.chain = "matrix",
		prob.chain = "array",
		x.ind.chain = "array",
		x.chain = "matrix",
		accept = "matrix",
		obsdiss = "list",
		dens = "list",
		control = "list",
		dim = "list"
	)
)

setMethod("initialize", "dmbc",
		function(
			.Object,
			z.chain = array(),
			z.chain.p = array(),
			alpha.chain = matrix(),
			eta.chain = matrix(),
			sigma2.chain = matrix(),
			lambda.chain = matrix(),
			prob.chain = array(),
			x.ind.chain = array(),
			x.chain = matrix(),
			accept = matrix(),
			obsdiss = list(),
			dens = list(),
			control = list(),
			dim = list()
		)
		{
			.Object@z.chain <- z.chain
			.Object@z.chain.p <- z.chain.p
			.Object@alpha.chain <- alpha.chain
			.Object@eta.chain <- eta.chain
			.Object@sigma2.chain <- sigma2.chain
			.Object@lambda.chain <- lambda.chain
			.Object@prob.chain <- prob.chain
			.Object@x.ind.chain <- x.ind.chain
			.Object@x.chain <- x.chain
			.Object@accept <- accept
			.Object@obsdiss <- obsdiss
			.Object@dens <- dens
			.Object@control <- control
			.Object@dim <- dim
			.Object
		}
)

setMethod("summary",
		"dmbc",
		function(object, summary.Z = FALSE, summary.alpha = TRUE, summary.eta = TRUE, summary.sigma2 = TRUE, summary.lambda = TRUE,
				quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
			n <- object@dim[["n"]]
			p <- object@dim[["p"]]
			G <- object@dim[["G"]]
			burnin <- object@control[["burnin"]]
			nsim <- object@control[["nsim"]]
			totiter <- burnin + nsim
			out <- list()
			if (summary.Z) {
				dat <- object@z.chain.p
				attr(dat, "dim") <- c(totiter, n*p*G)
				colnm <- character(0)
				for (g in 1:G) {
					for (j in 1:p) {
						for (i in 1:n) {
							colnm <- append(colnm, paste("z_", i, "_", j, "_", g, sep = ""))
						}
					}
				}
				colnames(dat) <- colnm
				z.coda <- as.mcmc(dat)
				cat("\n***  Z LATENT POSITIONS  ***\n")
				print(summary(z.coda, quantiles = quantiles))
				out[["z"]] <- print(summary(z.coda, quantiles = quantiles))
			}
			if (summary.alpha) {
				dat <- object@alpha.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("alpha_", g, sep = ""))
				}
				colnames(dat) <- colnm
				alpha.coda <- as.mcmc(dat)
				cat("\n***  ALPHA PARAMETERS  ***\n")
				out[["alpha"]] <- print(summary(alpha.coda, quantiles = quantiles))
			}
			if (summary.eta) {
				dat <- object@eta.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("eta_", g, sep = ""))
				}
				colnames(dat) <- colnm
				eta.coda <- as.mcmc(dat)
				cat("\n***  ETA PARAMETERS  ***\n")
				out[["eta"]] <- print(summary(eta.coda, quantiles = quantiles))
			}
			if (summary.sigma2) {
				dat <- object@sigma2.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("sigma2_", g, sep = ""))
				}
				colnames(dat) <- colnm
				sigma2.coda <- as.mcmc(dat)
				cat("\n***  SIGMA2 PARAMETERS  ***\n")
				out[["sigma2"]] <- print(summary(sigma2.coda, quantiles = quantiles))
			}
			if (summary.lambda) {
				dat <- object@lambda.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("lambda_", g, sep = ""))
				}
				colnames(dat) <- colnm
				lambda.coda <- as.mcmc(dat)
				cat("\n***  LAMBDA PARAMETERS  ***\n")
				out[["lambda"]] <- print(summary(lambda.coda, quantiles = quantiles))
			}
			
			return(out)
		}
)

setMethod("plot",
		signature(x = "dmbc", y = "missing"),
		function(x, plot.Z = TRUE, plot.alpha = TRUE, plot.eta = TRUE, plot.sigma2 = TRUE, plot.lambda = TRUE, plot.dens = TRUE, ...) {
			n <- x@dim[["n"]]
			p <- x@dim[["p"]]
			G <- x@dim[["G"]]
			burnin <- x@control[["burnin"]]
			nsim <- x@control[["nsim"]]
			totiter <- burnin + nsim
			if (plot.Z) {
				to.plot <- FALSE
				if(interactive()) {
					answer <- readline(paste("The number of elements in the Z latent positions matrix\n",
						"is usually large, so many graphs are going to be produced.\nProceed anyway? [yes, no]", sep = ""))
					if (substr(tolower(answer), 1, 1) == "y")
						to.plot <- TRUE
				}
				if (to.plot) {
					dat <- x@z.chain.p
					attr(dat, "dim") <- c(totiter, n*p*G)
					colnm <- character(0)
					for (g in 1:G) {
						for (j in 1:p) {
							for (i in 1:n) {
								colnm <- append(colnm, paste("z_", i, "_", j, "_", g, sep = ""))
							}
						}
					}
					colnames(dat) <- colnm
					z.coda <- as.mcmc(dat)
					dev.new()
					plot(z.coda, ask = TRUE)
				}
			}
			if (plot.alpha) {
				dat <- x@alpha.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("alpha_", g, sep = ""))
				}
				colnames(dat) <- colnm
				alpha.coda <- as.mcmc(dat)
				dev.new()
				plot(alpha.coda, ask = FALSE)
			}
			if (plot.eta) {
				dat <- x@eta.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("eta_", g, sep = ""))
				}
				colnames(dat) <- colnm
				eta.coda <- as.mcmc(dat)
				dev.new()
				plot(eta.coda, ask = FALSE)
			}
			if (plot.sigma2) {
				dat <- x@sigma2.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("sigma2_", g, sep = ""))
				}
				colnames(dat) <- colnm
				sigma2.coda <- as.mcmc(dat)
				dev.new()
				plot(sigma2.coda, ask = FALSE)
			}
			if (plot.lambda) {
				dat <- x@lambda.chain
				colnm <- character(0)
				for (g in 1:G) {
					colnm <- append(colnm, paste("lambda_", g, sep = ""))
				}
				colnames(dat) <- colnm
				lambda.coda <- as.mcmc(dat)
				dev.new()
				plot(lambda.coda, ask = FALSE)
			}
			if (plot.dens) {
				dat <- cbind(x@dens[["loglik"]], x@dens[["logprior"]], x@dens[["logpost"]])
				colnames(dat) <- c("loglikelihood", "logprior", "logposterior")
				dens.coda <- as.mcmc(dat)
				dev.new()
				plot(dens.coda, ask = FALSE)
			}
		}
)

setClass(Class = "dmbcIC",
	representation(
		logprior = "matrix",
		logmlik = "matrix",
		logcorrfact = "matrix",
		DCIC = "matrix",
		post.mean = "list",
		control = "list",
		res_last_p = "list"
	)
)

setMethod("initialize", "dmbcIC",
		function(
			.Object,
			logprior = matrix(),
			logmlik = matrix(),
			logcorrfact = matrix(),
			DCIC = matrix(),
			post.mean = list(),
			control = list(),
			res_last_p = list()
		)
		{
			.Object@logprior <- logprior
			.Object@logmlik <- logmlik
			.Object@logcorrfact <- logcorrfact
			.Object@DCIC <- DCIC
			.Object@post.mean <- post.mean
			.Object@control <- control
			.Object@res_last_p <- res_last_p
			.Object
		}
)

setMethod("plot",
		signature(x = "dmbcIC"),
		function(x, ...) {
			pmax <- nrow(x@DCIC)
			Gmax <- ncol(x@DCIC)
			clrs <- colors(distinct = TRUE)[c(21, 410, 23, 235, 62, 405, 28, 319, 360, 394)]
			clrs <- c(clrs, setdiff(colors(distinct = TRUE)[-1], clrs))
			lgd <- character(0)

			plot(1:Gmax, x@DCIC[1, ], type = "n", xlab = "G", ylab = "DCIC", xaxt = "n", ylim = c(min(x@DCIC[is.finite(x@DCIC)]), max(x@DCIC[is.finite(x@DCIC)])),
				main = "DCIC")
			for (p.i in 1:pmax) {
				lines(1:Gmax, x@DCIC[p.i, ], type = "b", col = clrs[p.i])
				lgd <- append(lgd, paste("p = ", p.i, sep = ""))
			}
			axis(1, at = 1:Gmax)
			dcic.min.G <- ifelse((which.min(t(x@DCIC)) %% Gmax) == 0, Gmax, (which.min(t(x@DCIC)) %% Gmax))
			points(dcic.min.G, min(x@DCIC), col = "green", pch = 10, cex = 1.75, lwd = 1.5)
			legend("topleft", lgd, lty = rep(1, pmax), col = clrs[1:pmax])
		}
)
