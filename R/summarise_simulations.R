library(dplyr)
cover <- function(x, interval) {
  if (x < interval[1]  || x > interval[2])
    return(0)
  return (1)
}
imse <- function(draws, true, dx) {
  mean_error_draws <- sweep(draws, 2, true, '-')
  mean_error <- colMeans(mean_error_draws)
  mean_mse <- mean_error^2
  var_f <- apply(draws,2,var)
  classical_imse <- sum(mean_mse * dx)
  classical_ivar <- sum(var_f * dx)
  classical_ibias <- sum(mean_error * dx)
  metrics <- c(
      'classical_imse' = classical_imse,
      'classical_ivar' = classical_ivar,
      'classical_ibias' = classical_ibias
    )
  return(
    metrics
   )
}
gen_mu <- function(mod, dat) {
  analytical_risk_x <- function(x, M) {
    return(0.15 + x^2 / M^2)
  }
  
  analytical_risk_y <- function(x, M, C) {
    return(analytical_risk_x(M/2,M) + x^2 / C^2)
  }
  
  analytical_risk_x2 <- function(x, M, C) {
    return(analytical_risk_y(2*C/3, M, C) - 1/4  + x^2 / M^2)
  }
  post_r_fun <- mod$draws("f",format = "draws_matrix")
  col_nms_post_r <- colnames(post_r_fun)
  post_r_fun_x <- cbind(post_r_fun[, grepl("f\\[1,", col_nms_post_r)], post_r_fun[, grepl("f\\[2,", col_nms_post_r)])
  pts_x <- c(dat$data$pts[1,],dat$data$pts[2,])
  dx <- c(rep(dat$data$delta_t[1],length(dat$data$pts[1,])),rep(dat$data$delta_t[2],length(dat$data$pts[2,])))
  true_x <- sapply(pts_x, function(x) analytical_risk_x(x, dat$M))
  
  post_r_fun_x2 <- cbind(post_r_fun[, grepl("f\\[3,", col_nms_post_r)], post_r_fun[, grepl("f\\[4,", col_nms_post_r)])
  pts_x2 <- c(dat$data$pts[3,],dat$data$pts[4,])
  dx2 <- c(rep(dat$data$delta_t[3],length(dat$data$pts[3,])),rep(dat$data$delta_t[4],length(dat$data$pts[4,])))
  true_x2 <- sapply(pts_x2, function(x) analytical_risk_x2(x, dat$M, dat$C))
  
  post_r_fun_y <- cbind(post_r_fun[, grepl("f\\[5,", col_nms_post_r)], post_r_fun[, grepl("f\\[6,", col_nms_post_r)])
  pts_y <- c(dat$data$pts[5,],dat$data$pts[6,])
  dy <- c(rep(dat$data$delta_t[5],length(dat$data$pts[5,])),rep(dat$data$delta_t[6],length(dat$data$pts[6,])))
  true_y <- sapply(pts_y, function(x) analytical_risk_y(x, dat$M, dat$C))
  
  return(
    list(
      mod_x = post_r_fun_x,
      true_x = true_x,
      pts_x = pts_x,
      dx = dx,
      mod_x2 = post_r_fun_x2,
      true_x2 = true_x2,
      pts_x2 = pts_x2,
      dx2 = dx2,
      mod_y = post_r_fun_y,
      true_y = true_y,
      pts_y = pts_y,
      dy = dy
    )
  )
}
gen_imse <- function(mod, dat) {
  mu_draws <- gen_mu(mod,dat)
  dx <- mu_draws$dx
  dx2 <- mu_draws$dx2
  dy <- mu_draws$dy
  
  imse_x <- imse(mu_draws$mod_x, mu_draws$true_x, dx)
  imse_x2 <- imse(mu_draws$mod_x2, mu_draws$true_x2, dx2)
  imse_y <- imse(mu_draws$mod_y, mu_draws$true_y, dy)
  imse_df <- list(
    x = imse_x,
    x2 = imse_x2,
    y = imse_y
  )
  return(
    imse_df
  )
}
suffix <- function(vec, suff) {
  names(vec) <- paste(names(vec), suff,sep='_')
  return(data.frame(t(vec)))
}
score_model <- function(mod_i, dat_i, data_nm, idx) {
  if(grepl('clustered',data_nm))
    spat_string <- 'c'
  else
    spat_string <- 'u'
  post_mean <- data.frame(
    alpha = rep(NA_real_,1),
    sd_alpha = rep(NA_real_,1),
    rho = rep(NA_real_,1),
    sd_rho = rep(NA_real_,1),
    cover_alpha = rep(NA_real_,1),
    cover_rho = rep(NA_real_,1)
  ) 
  max_treedepth <- mod_i$metadata()$max_treedepth
  sum_mod <- mod_i$summary()
  diags <- mod_i$sampler_diagnostics()
  mu_j <- mod_i$draws("mu",format = "draws_matrix")
  tot_risk_j <- mod_i$draws("total_risk",format = "draws_matrix")
  alpha <- mod_i$draws('alpha', format = "draws_matrix")
  beta <- mod_i$draws("beta", format = "draws_matrix")
  rho <- mod_i$draws("rho", format = "draws_matrix")
  interval_alpha <- quantile(alpha, c(0.1,0.9))
  interval_beta <- quantile(beta, c(0.1,0.9))
  interval_rho <- quantile(rho, c(0.1,0.9))
  intervals_p_50 <- t(apply(mu_j,2,quantile, c(0.25,0.75)))
  intervals_risk_50 <- t(apply(tot_risk_j,2,quantile, c(0.25,0.75)))
  intervals_p_80 <- t(apply(mu_j,2,quantile, c(0.10,0.90)))
  intervals_risk_80 <- t(apply(tot_risk_j,2,quantile, c(0.10,0.90)))
  imse_mod <- gen_imse(mod_i, dat_i)
  post_mean[1,'alpha'] <- mean(alpha) 
  post_mean[1,'beta'] <- mean(beta) 
  post_mean[1,'rho'] <- mean(rho) 
  post_mean[1,'bias_alpha'] <- post_mean[1,'alpha'] - dat_i$int
  post_mean[1,'bias_beta'] <- post_mean[1,'beta'] - dat_i$slope
  post_mean[1,'bias_rho'] <-  post_mean[1,'rho'] - dat_i$rho
  post_mean[1,'sd_beta'] <- sd(beta) 
  post_mean[1,'sd_alpha'] <- sd(alpha) 
  post_mean[1,'sd_rho'] <- sd(rho) 
  post_mean[1,'cor_rho_alpha'] <- cor(rho,alpha) 
  post_mean[1,'cor_rho_beta'] <- cor(rho,beta) 
  post_mean[1,'cor_alpha_beta'] <- cor(alpha,beta) 
  post_mean[1,'cover_alpha'] <- cover(dat_i$int,interval_alpha)
  post_mean[1,'cover_beta'] <- cover(dat_i$slope,interval_beta)
  post_mean[1,'cover_rho'] <- cover(dat_i$rho,interval_rho)
  post_mean[1,'pred_bias'] <- mean(colMeans(mu_j) - dat_i$canal_risks)
  post_mean[1,'int_pred_bias'] <- mean(colMeans(mu_j) + post_mean[1,'alpha'] - (dat_i$canal_risks + dat_i$int))
  post_mean[1,'risk_bias'] <- mean(colMeans(tot_risk_j) - dat_i$risks)
  post_mean[1,'pred_width_50'] <- mean(intervals_p_50[,2] - intervals_p_50[,1])
  post_mean[1,'pred_width_80'] <- mean(intervals_p_80[,2] - intervals_p_80[,1])
  post_mean[1,'risk_width_50'] <- mean(intervals_risk_50[,2] - intervals_risk_50[,1])
  post_mean[1,'risk_width_80'] <- mean(intervals_risk_80[,2] - intervals_risk_80[,1])
  post_mean[1,'pred_cover_50'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$canal_risks[x], intervals_p_50[x,])))
  post_mean[1,'pred_cover_80'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$canal_risks[x], intervals_p_80[x,])))
  post_mean[1,'risk_cover_50'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$risks[x], intervals_risk_50[x,])))
  post_mean[1,'risk_cover_80'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$risks[x], intervals_risk_80[x,])))
  post_mean[1,'pred_cover_50_01'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$canal_risks[x], intervals_p_50[x,]))[dat_i$canal_risks > 0.01])
  post_mean[1,'pred_cover_80_01'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$canal_risks[x], intervals_p_80[x,]))[dat_i$canal_risks > 0.01])
  post_mean[1,'risk_cover_50_01'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$risks[x], intervals_risk_50[x,]))[dat_i$canal_risks > 0.01])
  post_mean[1,'risk_cover_80_01'] <- mean(sapply(1:dat_i$data$J, function(x) cover(dat_i$risks[x], intervals_risk_80[x,]))[dat_i$canal_risks > 0.01])
  post_mean[1,'max_rhat'] <- max(sum_mod[,'rhat'])
  post_mean[1,'min_neff_bulk'] <- min(sum_mod[,'ess_bulk'])
  post_mean[1, "min_neff_tail"] <- min(sum_mod[, "ess_tail"])
  post_mean[1,'n_diverge'] <- sum(diags[,,'divergent__'])
  post_mean[1, "max_treedepth"] <- sum(diags[, , "treedepth__"] == max_treedepth)
  post_mean[1, "idx"] <- idx
  post_mean <- cbind(post_mean,suffix(imse_mod$x,'x'))
  post_mean <- cbind(post_mean,suffix(imse_mod$y,'y'))
  post_mean <- cbind(post_mean,suffix(imse_mod$x2,'x2'))
  post_mean <- post_mean %>%
    mutate(
      n = dat_i$data$J,
      G = dat_i$data$G,
      spat = spat_string,
      tot_n = n * G
    )
  return(post_mean)
}

read_files <- function(dir) {
  fls <- list.files(dir, recursive = TRUE,full.names = TRUE)
  res <- list()
  for (fl_i in seq_along(fls)) {
    res[[fl_i]] <- readRDS(fls[fl_i])
  }
  res <- do.call(rbind,res)
  return(res)
}
