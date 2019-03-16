#' Estimation of a DMBC model.
#'
#' \code{dmbc()}, the main function of the package, estimates a DMBC model
#'   for a given set of \emph{S} dissimilarity matrices.
#'
#' @param data An object of class \code{dmbc_data} containing the data
#'   to analyze.
#' @param p A length-one numeric vector indicating the number of dimensions of the
#'   latent space.
#' @param G A length-one numeric vector indicating the number of cluster to
#'   partition the \emph{S} subjects.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution See
#'   \code{\link{dmbc_control}} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{dmbc_prior}} for more details.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created for the duration of the \code{dmbc()} call.
#' @return A \code{dmbc_fit_list} object.
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#' @seealso \code{\link{bmds}} for Bayesian (metric) multidimensional scaling.
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 20000
#' nsim <- 10000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "multicore")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' summary(sim.dmbc, include.burnin = FALSE)
#'
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("teal")
#' plot(sim.dmbc, what = "trace", regex_pars = "eta")
#'
#' z <- dmbc_get_configuration(sim.dmbc, chain = 1, est = "mean",
#'   labels = 1:16)
#' summary(z)
#' color_scheme_set("mix-pink-blue")
#' graph <- plot(z, size = 2, size_lbl = 3, label_object = FALSE, adjust = .1)
#' graph + panel_bg(fill = "gray90", color = NA)
#' graph + geom_text(aes(label = lbl), nudge_x = .75, nudge_y = 0, size = 3)
#' }
#' @export
dmbc <- function(data, p = 2, G = 3, control = dmbc_control(), prior = NULL, cl = NULL) {
  D <- data@diss

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
  
  ### [for future developments] ###
  family <- "binomial"
  if (is.null(family))
    stop("the family argument is required.")
  if (!(family %in% .dmbcEnv$allowedfamilies))
    stop("'family' not recognized.")

  .dmbcEnv$current_p <- p
  .dmbcEnv$current_G <- G
  .dmbcEnv$current_family <- family

  control <- check_list_na(control, dmbc_control())
  if (!check_control(control))
    stop("the control list is not correct; see the documentation for more details.")

  nsim <- control[["nsim"]]
  burnin <- control[["burnin"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  threads <- control[["threads"]]
  seed <- control[["seed"]]
  parallel <- control[["parallel"]]
  random.start <- control[["random.start"]]
  store.burnin <- control[["store.burnin"]]
  verbose <- control[["verbose"]]

  have_mc <- have_snow <- FALSE
  if (parallel != "no" && threads > 1L) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) {
      warning("number of cores forced to 1 (i.e. no parallel computing used).")
      threads <- 1L
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }

  S <- length(D)
  n <- attr(D[[1]], "Size")
  m <- n*(n - 1)/2
  totiter <- burnin + nsim
  p <- as.integer(p)
  G <- as.integer(G)
  
  if (G >= S)
    stop("number of groups/clusters needs to be smaller than the number of subjects.")
  
  # save current random number generator kind
  old.rng <- RNGkind()[1L]
  RNGkind(kind = "L'Ecuyer-CMRG")
  
  # perform MCMC simulation
  if (nchains > 1L && (have_mc || have_snow)) {
    dmbc_fit_parallel <- function(c, D.c, p.c, G.c, family.c, control.c, prior.c, lib) {
      suppressPackageStartupMessages(require(dmbc, lib.loc = lib))
      control.c[["verbose"]] <- FALSE
      # cat("Starting cluster node", c, "on local machine\n")
      start.c <- dmbc_init(D = D.c, p = p.c, G = G.c, family = family.c, random.start = control.c[["random.start"]])
      if (is.null(prior.c)) {
        prior.c <- dmbc_prior()
      } else {
        prior.c <- check_list_na(prior.c, dmbc_prior())
      }
      if (!check_prior(prior.c))
        stop("the prior hyperparameter list is not correct; see the documentation for more details.")
      dmbc_fit(D = D.c, p = p.c, G = G.c, family = family.c, control = control.c, prior = prior.c, start = start.c)
    }
    environment(dmbc_fit_parallel) <- .GlobalEnv # this prevents passing objects other than those needed for
                                                 # evaluating the dmbc_fit_parallel function (maybe it is not strictly
                                                 # needed)

    if (verbose) {
      devout <- ""
      if (.Platform$OS.type != "windows" && !have_mc) {
        cat("--- STARTING PARALLEL SIMULATION OF", nchains, "CHAINS ---\n")
      } else {
        cat("Performing parallel simulation of", nchains, "chains...\n")
      }
    } else {
      if (.Platform$OS.type != "windows") {
        devout <- '/dev/null'
      } else {
        devout <- 'nul:'
      }
    }

    res <- if (have_mc) {
             if (!is.null(seed)) {
               set.seed(seed)
               parallel::mc.reset.stream()
             }
             parallel::mclapply(seq_len(nchains), dmbc_fit_parallel, mc.cores = threads, mc.set.seed = TRUE,
               D.c = D, p.c = p, G.c = G, family.c = family, control.c = control, prior.c = prior,
               lib = .dmbcEnv$path.to.me)
           } else if (have_snow) {
             if (is.null(cl)) {
               cl <- parallel::makePSOCKcluster(rep("localhost", threads), outfile = devout) # outfile doesn't work on 
                                                                                             # Windows
               parallel::clusterSetRNGStream(cl, seed)
               res <- parallel::parLapply(cl, seq_len(nchains), dmbc_fit_parallel, D.c = D, p.c = p, G.c = G, 
                family.c = family, control.c = control, prior.c = prior, lib = .dmbcEnv$path.to.me)
               parallel::stopCluster(cl)
               res
             } else parallel::parLapply(cl, seq_len(nchains), dmbc_fit_parallel, D.c = D, p.c = p, G.c = G,
               family.c = family, control.c = control, prior.c = prior, lib = .dmbcEnv$path.to.me)
           }

    if (verbose) {
      if (.Platform$OS.type != "windows" && !have_mc){
        cat("--- END OF PARALLEL SIMULATION OF", nchains, "CHAINS ---\n")
      } else {
        # cat("done!\n")
      }
    }
  } else {
    res <- list()
    for (ch in 1:nchains) {
      if (verbose && nchains > 1L) cat("--- STARTING SIMULATION OF CHAIN", ch, "OF", nchains, "---\n")

      if (verbose) cat("Initialization of the algorithm...\n")
  
      dmbc.start <- dmbc_init(D, p, G, family, random.start)
      if (is.null(prior)) {
        prior <- dmbc_prior()
      } else {
        prior <- check_list_na(prior, dmbc_prior())
      }
      if (!check_prior(prior))
        stop("the prior hyperparameter list is not correct; see the documentation for more details.")
    
      if (verbose) {
        # cat("done!\n")
      }

      res[[ch]] <- dmbc_fit(D = D, p = p, G = G, family = family, control = control, prior = prior, start = dmbc.start)

      if (verbose && nchains > 1L) cat("--- END OF CHAIN", ch, "OF", nchains, "---\n\n")
    }
  }

  # restore previous random number generator kind
  RNGkind(kind = old.rng)

  res <- new("dmbc_fit_list", results = res)

  return(res)
}
