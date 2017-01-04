library(dmbc)

# load data
data(simdiss, package = "dmbc")

G <- 3
p <- 2
prm.prop <- list(z = 1.5, alpha = .75)
prm.prior <- list(sigma2 = c(1e-1, 1e-1), lambda = rep(1, G))
burnin <- 20000
nsim <- 10000
seed <- 1406

set.seed(seed)
sim.dmbc <- dmbc(simdiss, p, G, burnin, nsim, prm.prop, prm.prior, random.start = TRUE, verbose = TRUE)

# summarize and plot the results
summary(sim.dmbc)
plot(sim.dmbc)
