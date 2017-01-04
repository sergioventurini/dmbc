library(dmbc)

# load data
data(simdiss, package = "dmbc")

# model choice
pmax <- 2
Gmax <- 2
prm.prop <- list(z = 1.5, alpha = .75)
prm.prior <- c(1e-1, 1e-1)	# note that when the IC are computed it is assumed that the lambda hyperparameters are set to rep(1, G) for every model
							# (otherwise the whole (long!) list of values should be provided as an argument)
burnin <- 1000
nsim <- 2000
seed <- 1809

set.seed(seed)
sim.ic <- dmbc.IC(D = simdiss, pmax = pmax, Gmax = Gmax, burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = prm.prior, random.start = TRUE)

pmax <- pmax + 1
Gmax <- Gmax + 2
new.ic <- dmbc.IC.more(D = simdiss, pmax = pmax, Gmax = Gmax, sim.ic, verbose = TRUE)

# plot the results
plot(new.ic)
