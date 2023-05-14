l <- 0.3
M <- 10
C <- 4

cover <- function(x, interval) {
  if (x < interval[1]  || x > interval[2])
    return(0)
  return (1)
}

## Generate houshold coordinates from PRM 
# kern <- function(d, l, nu) {
#   (1 + d^2 / (nu * l^2)) ^ (-(nu + 1) / 2)
# }

kern <- function(d, l, nu) {
  exp(-d/l)
}

kern_ab <- function(a, b, l, nu) {
  exp(-sqrt(a^2 + b^2)/l)
}

analytical_risk_x <- function(x, M) {
  return(0.15 + x^2 / M^2)
}

analytical_risk_y <- function(x, M, C) {
  return(analytical_risk_x(M/2,M) + x^2 / C^2)
}

analytical_risk_x2 <- function(x, M, C) {
  return(analytical_risk_y(2*C/3, M, C) - 1/4  + x^2 / M^2)
}

kernel_x <- function(x, a, b, l, nu, M) {
  kern_ab(x - a, b, l, nu) * analytical_risk_x(x, M)
}

kernel_y <- function(y, a, b, l, nu, M, C) {
  kern_ab(M/2 - a, b - y, l, nu) * analytical_risk_y(y, M, C)
}

kernel_x2 <- function(x, a, b, l, nu, M, C) {
  kern_ab(x - a, b - 2*C/3, l, nu) * analytical_risk_x2(x, M, C)
}

b_fun <- function(u, M) {
  return(M - M * sqrt(1 - u))
}

trunc_exp <- function(u,lambda,t) {
  u_star <- u * pexp(t,lambda)
  return(-1/lambda * log(1 - u_star))
}

l_trunc_exp <- function(u,lambda,t) {
  return(t + -1/lambda * log(1 - u))
}

lt_trunc_exp <- function(u,lambda,l,t) {
  u_star <- u * (exp(-lambda * l) - exp(-lambda * t))
  return(-1/lambda * log(exp(-l*lambda) - u_star))
}

pts_extract <- function(pts) {
  len <- length(pts)
  endpts <- pts[c(1,len)]
  intervals <- pts
  centroids <- sapply(2:length(intervals),
                      function(i)
                        mean(intervals[(i-1):i])
  )
  dt <- diff(intervals)[1]
  return(
    list(
      endpts = endpts,
      centroids = centroids,
      dt = dt
    )
  )
}

gen_house_locs <- function(n = 564,M = 10,C = 4,house_dist = "uniform") {
   if (house_dist == 'uniform') {
      xy <- randtoolbox::sobol(n, dim = 2)
      a <- xy[,1] * M
      b <- xy[,2] * C
   } else if (house_dist == 'clustered_uniform') {
     community_sizes <- as.numeric(table(sample(c(1,2,3,4),n,replace=T)))
     n_div_4 <- community_sizes[1]
     a <- runif(n_div_4) * M
     b <- sapply(runif(n_div_4), function(u) lt_trunc_exp(u,4, 0.005, 1/3*C) + 2/3 * C)
     n_div_4 <- community_sizes[2]
     a <- c(a,runif(n_div_4) * M)
     b_new <- sapply(runif(n_div_4), function(u) -lt_trunc_exp(u,4, 0.005, 1/3*C) + 2/3 * C)
     b <- c(b,b_new)
     n_div_4 <- community_sizes[3]
     a_new <- runif(n_div_4) * M
     b_new <- sapply(runif(n_div_4), function(u) lt_trunc_exp(u,4, 0.005, 2/3*C))
     b <- c(b,b_new)
     a <- c(a,a_new)
     n_div_4 <- community_sizes[4]
     small_communities <- as.numeric(table(sample(c(0,1),n_div_4,replace=T)))
     n_div_8 <- small_communities[1]
     b_new <- runif(n_div_8) * C
     a_new <- M/2 + sapply(runif(n_div_8), function(u) lt_trunc_exp(u,1, 0.005, M/4))
     b <- c(b,b_new)
     a <- c(a,a_new)
     n_div_8 <- small_communities[2]
     b_new <- runif(n_div_8) * C
     a_new <- M/2 - sapply(runif(n_div_8), function(u) lt_trunc_exp(u,1, 0.005, M/4))
     a <- c(a,a_new)
     b <- c(b,b_new)
   } else if (house_dist == 'clustered') {
     community_sizes <- as.numeric(table(sample(c(1,2,3,4),n,replace=T)))
     n_div_4 <- community_sizes[1]
     a <- runif(n_div_4) * M
     b <- sapply(runif(n_div_4), function(u) lt_trunc_exp(u,4, 0.005, 1/3*C) + 2/3 * C)
     n_div_4 <- community_sizes[2]
     a <- c(a,runif(n_div_4) * M)
     b_new <- sapply(runif(n_div_4), function(u) -lt_trunc_exp(u,4, 0.005, 1/3*C) + 2/3 * C)
     b <- c(b,b_new)
     n_div_4 <- community_sizes[3]
     a_new <- runif(n_div_4) * M
     b_new <- sapply(runif(n_div_4), function(u) lt_trunc_exp(u,4, 0.005, 2/3*C))
     b <- c(b,b_new)
     a <- c(a,a_new)
     n_div_4 <- community_sizes[4]
     a_cent <- M/2 + 0.1
     b_cent <- C/3
     theta <- runif(n_div_4, min = 0, max = 2*pi)
     r <- runif(n_div_4)
     a_new <- a_cent + r * cos(theta)
     b_new <- b_cent + r * sin(theta)
     a <- c(a,a_new)
     b <- c(b,b_new)
   }
   return(
     list(
       a = a,
       b = b,
       n = n,
       M = M,
       C = C,
       house_dist = house_dist
     )
   )
}

gen_comp_dataset <- function(M = 10, C = 4, l = 0.3, G = 3, 
                             int = 0.05, slope = -0.15, a, b) {
  n <- length(a)
  risk_bind_x <- function(a, b) {
    return(stats::integrate(kernel_x, 0, M, a, b, l, 5, M))
  }
  risk_bind_y <- function(a, b) {
    return(stats::integrate(kernel_y, 0, C, a, b, l, 5, M, C))
  }
  risk_bind_x2 <- function(a, b) {
    return(stats::integrate(kernel_x2, 0, M, a, b, l, 5, M, C))
  }
  total_risk_bind <- function(a, b) {
    return(risk_bind_x(a, b)$value + risk_bind_y(a, b)$value + risk_bind_x2(a, b)$value)
  }
  canal_risk <- mapply(function(x, y) total_risk_bind(x, y), a, b)
  x <- rnorm(n)
  x <- scale(x)[, 1]
  risks <- (canal_risk + int) * exp(x * slope) * 1
  z <- lapply(risks, function(x) rbinom(G, 1, 1 - exp(-x)))
  y_mat <- t(do.call(what = cbind, args = z))
  dat <- list(
    y = y_mat,
    coords = cbind(a, b),
    J = nrow(y_mat),
    M = M,
    G = G,
    C = C,
    x = x,
    canal_risk = canal_risk,
    risks = risks,
    int = int,
    slope = slope,
    l = l
  )
  return(dat)
}

gen_data <- function(J_sample, full_dat, S) {
  idx_sample <- sample(full_dat$J,J_sample,replace = FALSE)
  dat <- full_dat
  dat$y <- dat$y[idx_sample,]
  dat$coords <- dat$coords[idx_sample,]
  dat$J <- J_sample
  dat$x <- dat$x[idx_sample]
  dat$canal_risk <- dat$canal_risk[idx_sample]
  dat$risks <- dat$risks[idx_sample]
  M <- dat$M
  C <- dat$C
  pts_x1 <- pts_extract(seq(0,M/2,length.out = S + 1))
  pts_x2 <- pts_extract(seq(M/2,M,length.out = S + 1))
  pts_x3 <- pts_extract(seq(0,M/2,length.out = S + 1))
  pts_x4 <- pts_extract(seq(M/2,M,length.out = S + 1))
  pts_y1 <- pts_extract(seq(0,2*C/3,length.out = S+1))
  pts_y2 <- pts_extract(seq(2*C/3,C,length.out = S+1))

  x1_dist <- lapply(pts_x1$centroids,
                    function(s)
                      sqrt(rowSums(sweep(dat$coords,MARGIN = 2, STATS = c(s,0), FUN = '-')^2)
                            )
                    )
  x2_dist <- lapply(pts_x2$centroids,
                    function(s)
                      sqrt(rowSums(sweep(dat$coords,MARGIN = 2, STATS = c(s,0), FUN = '-')^2)
                            )
                    )
  x3_dist <- lapply(pts_x3$centroids,
                    function(s)
                      sqrt(rowSums(sweep(dat$coords,MARGIN = 2, STATS = c(s,2*C/3), FUN = '-')^2)
                            )
                    )
  x4_dist <- lapply(pts_x4$centroids,
                    function(s)
                      sqrt(rowSums(sweep(dat$coords,MARGIN = 2, STATS = c(s,2*C/3), FUN = '-')^2)
                            )
                    )
  y1_dist <- lapply(pts_y1$centroids,
                    function(s)
                      sqrt(rowSums(sweep(dat$coords,MARGIN = 2, STATS = c(M/2,s), FUN = '-')^2)
                          )
                    )
  y2_dist <- lapply(pts_y2$centroids,
                    function(s)
                      sqrt(rowSums(sweep(dat$coords,MARGIN = 2, STATS = c(M/2,s), FUN = '-')^2)
                          )
                    )
  dat_mod <- dat
  dat_mod$dists <- list(do.call(rbind,x1_dist),
                        do.call(rbind,x2_dist),
                        do.call(rbind,x3_dist),
                        do.call(rbind,x4_dist),
                        do.call(rbind,y1_dist),
                        do.call(rbind,y2_dist)
  )
  dat_mod$delta_t <- c(pts_x1$dt,pts_x2$dt,pts_x3$dt,pts_x4$dt,pts_y1$dt,pts_y2$dt)
  dat_mod$pts <- do.call(rbind,list(pts_x1$centroids,pts_x2$centroids,
                                    pts_x3$centroids,pts_x4$centroids,
                                    pts_y1$centroids,pts_y2$centroids))
  dat_mod$endpt_coords <- do.call(rbind,list(pts_x1$endpts, pts_x2$endpts,
                                              pts_x3$endpts, pts_x4$endpts,
                                              pts_y1$endpts, pts_y2$endpts))
  dat_mod$endpt_inds <- do.call(rbind,list(c(1,2),c(2,3),
                                            c(4,5),c(5,6),
                                            c(2,5),c(5,7)))
  dat_mod$N_seg <- 6
  dat_mod$N_endpts <- dat_mod$N_seg + 1

  dat_mod$S <- S
  dat_mod$alpha_loc <- 0
  dat_mod$alpha_scale <- 0.3
  dat_mod$beta_scale <- 0.3
  dat_mod$beta_loc <- 0
  dat_mod$ell_shape <- 4
  dat_mod$ell_rate <- 1
  dat_mod$rho_scale <- 0.5

  return(list(
    data = dat_mod,
    M = M,
    C = C,
    risks = dat$risks,
    canal_risks = dat$canal_risk,
    int = dat$int,
    slope = dat$slope,
    rho = dat$l
    )
  )
}
gen_full_data_sets <- function(S, l, L, J_sched, par_dir, G = G, int = 5e-2, house_locs) {
  data_sets <- list()
  full_data <- gen_comp_dataset(
    M = house_locs$M,
    C = house_locs$C,
    l = l,
    G = G,
    int = int,
    a = house_locs$a,
    b = house_locs$b
  )
  house_dist <- house_locs$house_dist
  for (J in J_sched) {
    for (i in 1:L) {
      data_sets[[i]] <- gen_data(S = S, J_sample = J, full_dat = full_data)
    }
    saveRDS(data_sets,
            paste0(par_dir, "/S_", S, "_n_", J, "_G_", G, "_", house_dist, "_int_0_05.RDS"))
  }
  
}

