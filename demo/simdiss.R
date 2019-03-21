# Script to simulate the simdiss object
make_simdata <- function(prob, n, seed) {
  set.seed(seed)

  sample1 <- sample(1:3, n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  sample2 <- sample(1:3, n, replace = TRUE, prob = c(0.5, 0.2, 0.3))
  sample3 <- sample(1:3, n, replace = TRUE, prob = c(0.3, 0.2, 0.5))

  diss1 <- as.matrix(dist(sample1))
  diss1[diss1 > 0] <- 1
  diss2 <- as.matrix(dist(sample2))
  diss2[diss2 > 0] <- 1
  diss3 <- as.matrix(dist(sample3))
  diss3[diss3 > 0] <- 1

  set.seed(seed)
  diss1_ran1 <- diss1_ran2 <- diss1_ran3 <- diss1
  diss2_ran1 <- diss2_ran2 <- diss2_ran3 <- diss2_ran4 <- diss2_ran5 <- diss2
  diss3_ran1 <- diss3_ran2 <- diss3

  for (i in 1:(nrow(diss1) - 1)) {
    for (j in (i + 1):nrow(diss1)) {
      if (diss1[i, j] == 1) {
        diss1_ran1[i, j] <- diss1_ran1[j, i] <- rbinom(1, 1, prob)
        diss1_ran2[i, j] <- diss1_ran2[j, i] <- rbinom(1, 1, prob)
        diss1_ran3[i, j] <- diss1_ran3[j, i] <- rbinom(1, 1, prob)
      }
      if (diss1[i, j] == 0) {
        diss1_ran1[i, j] <- diss1_ran1[j, i] <- rbinom(1, 1, 0.1)
        diss1_ran2[i, j] <- diss1_ran2[j, i] <- rbinom(1, 1, 0.1)
        diss1_ran3[i, j] <- diss1_ran3[j, i] <- rbinom(1, 1, 0.1)
      }
      if (diss2[i, j] == 1) {
        diss2_ran1[i, j] <- diss2_ran1[j, i] <- rbinom(1, 1, prob)
        diss2_ran2[i, j] <- diss2_ran2[j, i] <- rbinom(1, 1, prob)
        diss2_ran3[i, j] <- diss2_ran3[j, i] <- rbinom(1, 1, prob)
        diss2_ran4[i, j] <- diss2_ran4[j, i] <- rbinom(1, 1, prob)
        diss2_ran5[i, j] <- diss2_ran5[j, i] <- rbinom(1, 1, prob)
      }
      if (diss2[i, j] == 0) {
        diss2_ran1[i, j] <- diss2_ran1[j, i] <- rbinom(1, 1, 0.1)
        diss2_ran2[i, j] <- diss2_ran2[j, i] <- rbinom(1, 1, 0.1)
        diss2_ran3[i, j] <- diss2_ran3[j, i] <- rbinom(1, 1, 0.1)
        diss2_ran4[i, j] <- diss2_ran4[j, i] <- rbinom(1, 1, 0.1)
        diss2_ran5[i, j] <- diss2_ran5[j, i] <- rbinom(1, 1, 0.1)
      }
      if (diss3[i, j] == 1) {
        diss3_ran1[i, j] <- diss3_ran1[j, i] <- rbinom(1, 1, prob)
        diss3_ran2[i, j] <- diss3_ran2[j, i] <- rbinom(1, 1, prob)
      }
      if (diss3[i, j] == 0) {
        diss3_ran1[i, j] <- diss3_ran1[j, i] <- rbinom(1, 1, 0.1)
        diss3_ran2[i, j] <- diss3_ran2[j, i] <- rbinom(1, 1, 0.1)
      }
    }
  }

  subjects <- c("diss1_1", "diss1_2", "diss1_3",
    "diss2_1", "diss2_2", "diss2_3", "diss2_4", "diss2_5",
    "diss3_1", "diss3_2")
   
  all_diss <- list()
  all_diss[[1]] <- as.dist(diss1_ran1)
  all_diss[[2]] <- as.dist(diss1_ran2)
  all_diss[[3]] <- as.dist(diss1_ran3)
  all_diss[[4]] <- as.dist(diss2_ran1)
  all_diss[[5]] <- as.dist(diss2_ran2)
  all_diss[[6]] <- as.dist(diss2_ran3)
  all_diss[[7]] <- as.dist(diss2_ran4)
  all_diss[[8]] <- as.dist(diss2_ran5)
  all_diss[[9]] <- as.dist(diss3_ran1)
  all_diss[[10]] <- as.dist(diss3_ran2)
   
  names(all_diss) <- subjects

  return(all_diss)
}

library(dmbc)

n <- 16
S <- 10

all_diss_9 <- make_simdata(prob = 0.9, n = n, seed = 120)
all_diss_8 <- make_simdata(prob = 0.8, n = n, seed = 120)
all_diss_7 <- make_simdata(prob = 0.7, n = n, seed = 120)
all_diss_6 <- make_simdata(prob = 0.6, n = n, seed = 120)
all_diss_5 <- make_simdata(prob = 0.5, n = n, seed = 120)

simdata <- new("dmbc_data", diss = all_diss_9, n = n, S = S, 
  family = "binomial")
