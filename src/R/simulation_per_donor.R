args <- commandArgs(TRUE)
if (length(args) >= 1) {
  donor <- args[1]
}

library(ggpubr)
library(tidyverse)
library(cardelino)
library(MASS)

in_dir = "data/cell_assignment/"
out_dir <- "data/simulations/"
setwd(out_dir)

carde <- readRDS(paste0(in_dir, "cardelino_results.",
                        donor, ".filt_lenient.cell_coverage_sites.rds"))

set.seed(0)

## re-assignment
carde_binom <- cell_assign_EM(carde$A, carde$D, carde$tree$Z,
                              Psi = NULL, model = "binomial")
carde[["prob_binom"]] <- carde_binom$prob
carde[["theta_binom"]] <- carde_binom$theta

carde_Bern <- cell_assign_EM(carde$A, carde$D, carde$tree$Z,
                             Psi = NULL, model = "Bernoulli")
carde[["prob_Bern"]] <- carde_Bern$prob
carde[["theta_Bern"]] <- carde_Bern$theta

carde_Gibbs <- cell_assign_Gibbs(carde$A, carde$D,
                                 carde$tree$Z, Psi = NULL,
                                 min_iter = 1000, wise = "variant")
carde[["prob_Gibbs"]] <- carde_Gibbs$prob
carde[["theta0_Gibbs"]] <- carde_Gibbs$theta0
carde[["theta1_Gibbs"]] <- carde_Gibbs$theta1


## save estimate of real data
saveRDS(carde, paste0(out_dir, donor,
                      ".filt_lenient.cell_coverage_sites.mult.rds"))


## simulation
demuxlet <- function(A, D, Config, theta1=0.5, theta0=0.01) {
  P0_mat <- dbinom(A, D, theta0, log = TRUE)
  P1_mat <- dbinom(A, D, theta1, log = TRUE)
  P0_mat[which(is.na(P0_mat))] <- 0
  P1_mat[which(is.na(P1_mat))] <- 0
  logLik_mat <- t(P0_mat) %*% (1 - Config) + t(P1_mat) %*% Config
  prob_mat <- exp(logLik_mat) / rowSums(exp(logLik_mat))
  prob_mat
}

fit.beta <- fitdistr(carde_Gibbs$theta1[, 1], "beta",
                     start = list(shape1 = 2, shape2 = 5))
beta_var <- c(fit.beta$estimate[[1]] + fit.beta$estimate[[2]])
beta_mean <- fit.beta$estimate[[1]] / beta_var

Config <- carde$tree$Z
Psi <- carde$tree$P[,1]
sim_dat <- sim_read_count(Config, carde$D, Psi = Psi, cell_num = 500,
                          wise0 = "global", wise1 = "variant",
                          means = c(carde_Gibbs$theta0, beta_mean),
                          vars = c(20, beta_var))

I_sim <- sim_dat$I_sim
A_sim <- sim_dat$A_sim
D_sim <- sim_dat$D_sim
A_germ_sim <- sim_dat$A_germ_sim
D_germ_sim <- sim_dat$D_germ_sim

sim_Bern <- cell_assign_EM(sim_dat$A_sim, sim_dat$D_sim, Config, Psi = NULL)
sim_binom <- cell_assign_EM(sim_dat$A_sim, sim_dat$D_sim, Config, Psi = NULL,
                            model = "binomial")
sim_Gibbs <- cell_assign_Gibbs(sim_dat$A_sim, sim_dat$D_sim, Config, Psi = NULL,
                               wise = "variant")
sim_dat[["prob_demux"]] <- demuxlet(sim_dat$A_sim, sim_dat$D_sim, Config,
                                    theta0 = carde_Gibbs$theta0)
sim_dat[["prob_Gibbs"]] <- sim_Gibbs$prob
sim_dat[["theta0_Gibbs"]] <- sim_Gibbs$theta0
sim_dat[["theta1_Gibbs"]] <- sim_Gibbs$theta1

sim_dat[["prob_Bern"]] <- sim_Bern$prob
sim_dat[["theta_Bern"]] <- sim_Bern$theta

sim_dat[["prob_binom"]] <- sim_binom$prob
sim_dat[["theta_binom"]] <- sim_binom$theta


## save samulation data
saveRDS(sim_dat, paste0(out_dir, donor, ".simulate.rds"))

