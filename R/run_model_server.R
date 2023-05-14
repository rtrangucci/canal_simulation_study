args <- commandArgs(trailingOnly = TRUE)
data_fl <- args[1]
data_idx <- as.numeric(args[2])
mod_nm <-  args[3]
seed <- data_idx
set.seed(seed)
print(args)
data_fit <- readRDS(data_fl)[[data_idx]]
inf_model <- readRDS(mod_nm)
prep <- basename(data_fl)
prep <- strsplit(prep, ".RDS")[[1]][1]
tot_sv_dir <- file.path("/scratch/stats_dept_root/stats_dept1/trangucc/",prep)
if (!dir.exists(tot_sv_dir)) {
  dir.create(tot_sv_dir)
}

options(CMDSTANR_NO_VER_CHECK = TRUE)
print("loading cmdstanr")
library(cmdstanr, quietly = T)
print("loaded cmdstanr")
source("R/summarise_simulations.R")
library(dplyr)
library(posterior)

print("loading fit_mod")
fit_mod <- function(model, data, prep, data_idx) {
  fit <- model$sample(
    data = data$data,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 1000,
    chains = 4, parallel_chains = 4,
    adapt_delta = 0.90,
    max_treedepth = 12
  )
  score <- score_model(fit, data, prep, data_idx)
  return(
    score
  )
}

print("loaded fit_mod")

print("fitting-model")
fit <- fit_mod(inf_model, data_fit, prep, data_idx)

saveRDS(fit, file.path(tot_sv_dir,paste0(prep,"_idx_",data_idx,".RDS")))
