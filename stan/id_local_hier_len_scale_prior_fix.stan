functions {
  real exp_quad(real d, real ell) {
     return exp(-0.5/ell^2 * d^2);    
  }
  real poisson_binomial_lpmf(int y, vector theta) {
    int N = rows(theta);
    matrix[N + 1, N + 1] alpha;
    vector[N] log_theta = log(theta);
    vector[N] log1m_theta = log1m(theta);

    if (y < 0 || y > N)
      reject("poisson binomial variate y out of range, y = ", y,
             " N = ", N);
    for (n in 1:N)
      if (theta[n] < 0 || theta[n] > 1)
        reject("poisson binomial parameter out of range,",
               " theta[", n, "] =", theta[n]);

    if (N == 0)
      return y == 0 ? 0 : negative_infinity();

    // dynamic programming with cells
    // alpha[n + 1, tot + 1] = log prob of tot successes in first n trials
    alpha[1, 1] = 0;
    for (n in 1:N) {
      // tot = 0
      alpha[n + 1, 1] = alpha[n, 1] + log1m_theta[n];

      // 0 < tot < n
      for (tot in 1:min(y, n - 1))
        alpha[n + 1, tot + 1]
            = log_sum_exp(alpha[n, tot] + log_theta[n],
                          alpha[n, tot  + 1] + log1m_theta[n]);

      // tot = n
      if (n > y) continue;
      alpha[n + 1, n + 1] = alpha[n, n] + log_theta[n];
    }
    return alpha[N + 1, y + 1];
  }  
  real kernel(real d, real rho) {
    return exp(-d/rho);
  }
  matrix inv_22(real off_diag, real diag) {
    matrix[2, 2] K_inv;
    K_inv[1,1] = diag;
    K_inv[2,2] = diag;
    K_inv[1,2] = -off_diag;
    K_inv[2,1] = K_inv[1,2];
    return K_inv * inv(square(diag) - off_diag^2);
  }
  vector kernel_vec(real[] d, real rho) {
    int dim_d = size(d);
    vector[dim_d] kernels;
    for (i in 1:dim_d)
      kernels[i] = kernel(d[i], rho);
    return kernels; /// sum(kernels);
  }
//  y_arr[j, 1] = N_obs[j];
//  y_arr[j, 2] = N_canal[j];
//  y_arr[j,3] = D;
//  y_arr[j,4] = S;
//  y_arr[j,5:(N_obs[j]+ 4)] = y[j,1:N_obs[j]];
//  y_arr[j,N_obs[j]+ 4 + 1] = ll[j];
//  X_aug[j] = rep_matrix(0, max_obs, D + max_canal_pts + 1);
//  X_aug[j, 1:(D*N_obs[j]) = to_array_1d(X[j][1:N_obs[j],]);
//  X_aug[j, (D*N_obs[j] + 1):(N_canal[j] + D*N_obs[j])] = dists[1:N_canal[j],j];
//  X_aug[j, D + 1] = delta_t;

//  shared[1:S] = f;
//  shared[S+1] = alpha;
//  shared[S+2] = rho_scale;
//  shared[(S + 3):(S + 2 + D)] = beta;
// shared[(S + 2 + D + 1):(S + 2 + D + L)] = eta_rhos

  vector lp_fun(vector shared, vector hh_spec, real[] x, int[] y) {
    real canal = x[y[3] * y[1] + y[2] + 1] * dot_product(shared[1:y[4]], //  y_arr[j,(N_obs[j] + 5):(N_obs[j] + N_canal[j] + 4)] = canal_pts[j,1:N_canal[j]];
                                 kernel_vec(x[(y[3] * y[1] + 1):(y[2] + y[3] * y[1])], //  X_aug[j][1, (D + 1):(N_canal[j] + D)] = dists[1:N_canal[j],j];
                                            // y[5 + y[1]]
                   // shared[S+2] = rho_scale; shared[(S+2 + D + 1):(S+2 + D + L)] = eta_rhos;
                                            shared[y[4] + 2] * shared[y[4] + y[3] + 2 + y[5 + y[1]]]));
    return [bernoulli_lpmf(y[5:(y[1] + 4)] | 1 - exp(-(canal + shared[y[4] + 1]  //  shared[S+1] = alpha;
                                                           + hh_spec[1]) *
                                                             exp(to_matrix(x[1:(y[3]*y[1])], y[1], y[3]) * shared[(y[4] + 3):(y[4] + 2 + y[3])])))]';
  }
  row_vector kernel_rvec(vector d, real rho) {
    int dim_d = rows(d);
    row_vector[dim_d] kernels;
    for (i in 1:dim_d)
      kernels[i] = kernel(d[i], rho);
    return kernels;// / sum(kernels);
  }
  matrix kernel_mat(matrix d, real rho) {
    int rows_d = rows(d);
    int n_units = cols(d);
    matrix[n_units,rows_d] kernels;
    for (j in 1:n_units) {
      kernels[j,] = kernel_rvec(d[,j], rho);
    }
    return kernels;
  }
  vector integrate(matrix K, vector f, real dx) {
    int row_K = rows(K);
    int col_K = cols(K);
    vector[row_K] integral = K[,1] * f[1] + K[,col_K] * f[col_K];
    integral += 2 * block(K,1,2,row_K,col_K-2) * f[2:(col_K-1)];
    return integral * dx / 2;
  }
  vector integrate_simple(matrix K, vector f, real dx) {
    int row_K = rows(K);
    int col_K = cols(K);
    vector[row_K] integral = K * f;
    return integral * dx;
  }
  vector gen_mu(real[] pts, real ell, vector cond_vals, vector eta, real mag, int n_end_pts) {
    int S = size(pts);
    matrix[S, S] K = cov_exp_quad(pts, mag, ell) + diag_matrix(rep_vector(1e-6,S));
    int S_m_pts = S - n_end_pts;
    matrix[n_end_pts, n_end_pts] K_inv; 
    matrix[S_m_pts,n_end_pts] K_12;
    matrix[S_m_pts, S_m_pts] K_cond;
    if (n_end_pts == 2) {
      K_inv = inv_22(K[S,S-1],1+1e-6); 
    } else if (n_end_pts > 2) {
      matrix[n_end_pts,n_end_pts] K_block = block(K,S_m_pts + 1,S_m_pts + 1,n_end_pts, n_end_pts);
      matrix[n_end_pts,n_end_pts] K_block_chol = cholesky_decompose(K_block);
      matrix[n_end_pts,n_end_pts] K_L_inv = mdivide_left_tri_low(K_block_chol, diag_matrix(rep_vector(1,n_end_pts)));
      K_inv = K_L_inv' * K_L_inv;
    }
    K_12 = block(K,1,S_m_pts + 1, S_m_pts, n_end_pts);
    if (n_end_pts > 1) {
      K_cond = block(K,1,1,S_m_pts,S_m_pts) - quad_form_sym(K_inv,K_12');
      for (s in 1:S_m_pts) {
        K_cond[s,s] += 1e-6;
      }
      for (i in 1:(S_m_pts - 1))
        for (j in (i + 1):S_m_pts)
          K_cond[i,j] = K_cond[j,i];
      return exp(K_12 * K_inv * cond_vals + cholesky_decompose(K_cond) * eta); // * inv_mag_sq
    } else {
      K_cond = block(K,1,1,S_m_pts,S_m_pts) - K_12 * K_12';
      for (s in 1:S_m_pts) {
        K_cond[s,s] += 1e-6;
      }
      for (i in 1:(S_m_pts - 1))
        for (j in (i + 1):S_m_pts)
          K_cond[i,j] = K_cond[j,i];
      return exp(K_12 * cond_vals + cholesky_decompose(K_cond) * eta); // * inv_mag_sq
    }
  }
}
data {
  int<lower=1> J; // Number of households
  int<lower=1> S; // Number of canal segments
  int<lower=1> D; // Dimension of covariates
  int<lower=1> N_canal[J]; // Number of canal points per household
  int<lower=1> N_obs[J]; // Number of observations per household
  int<lower=1> D_Z;
  int<lower=1> max_obs; // max(N_obs) 
  int<lower=1> max_canal_pts; // max(N_canal) 
  int<lower=0, upper=1> y[J, max_obs]; // Observations
  int canal_pts[J, max_canal_pts];
  int<lower=1> K; // Number of localidad comb
  int<lower=1, upper = K> kk[J]; // Map from house to localidad comb
  int<lower=1> L; // Number of localidad
  int<lower=1, upper = L> ll[J]; // Map from house to localidad
  matrix[max_obs, D] X[J]; // Covariates
  matrix[max_canal_pts, J] dists; // HH dist to integration points int canal_pts[J, max_canal_pts]; // Canal index for integration points
  vector[J] dist;
  matrix[J, D_Z] Z;
//  real pts[S];
  real delta_t;
  matrix[D, D] R_inv;
  // canal ppc
  real delta_distance;
  int n_dist_intervals;
  int n_seg;
  int max_seg_pts;
  int num_pts[n_seg];
  matrix[max_seg_pts,n_seg] canal_pos;
  int canal_idx_mat[max_seg_pts,n_seg];
  // end_pts 
  int max_end_pts;
  int n_end_pts[n_seg];
  int end_pts_idx[n_seg,max_end_pts];
  matrix[n_seg,max_end_pts] end_pts_pos;
  int tot_end_pts;
  matrix[tot_end_pts,tot_end_pts] pt_dist_mat;
  int max_dep_pts;
  int dep_pts[tot_end_pts,max_dep_pts];
  int comp_order[tot_end_pts];
  int n_dep_pts[tot_end_pts];
  real<lower=0> gam_shape;
  real<lower=0> gam_rate;
}
transformed data {
  int y_arr[J,max_obs + 4 + 1] = rep_array(0, J, max_obs + 4 + 1);
  real X_aug[J,D*max_obs + max_canal_pts + 1];
  matrix[max_seg_pts + max_end_pts,n_seg] aug_pts;
  real<lower=0> delta = 1e-10;
  int lower_bounds[n_dist_intervals];
  matrix[J, D_Z] scale_Z;
  real weighting = sqrt(0.5);
  for (d in 1:n_dist_intervals)
    lower_bounds[d] = d - 1;
  for (j in 1:J) {
    y_arr[j,1] = N_obs[j];
    y_arr[j,2] = N_canal[j];
    y_arr[j,3] = D;
    y_arr[j,4] = S;
    y_arr[j,5:(N_obs[j]+ 4)] = y[j,1:N_obs[j]];
    y_arr[j,N_obs[j]+ 4 + 1] = ll[j];
    X_aug[j] = rep_array(0, D*max_obs + max_canal_pts + 1);
    X_aug[j, 1:(D * N_obs[j])] = to_array_1d(X[j][1:N_obs[j],]);
    X_aug[j, (D * N_obs[j]+ 1):(S + D * N_obs[j])] = to_array_1d(dists[,j]);
    X_aug[j, S + D * N_obs[j] + 1] = delta_t;
  }
  for (d in 1:D_Z)
    scale_Z[,d] = Z[,d] / sd(Z[,d]);
  for (seg in 1:n_seg) {
    aug_pts[1:num_pts[seg],seg] = canal_pos[1:num_pts[seg],seg];
    if (n_end_pts[seg] > 0) {
      aug_pts[(num_pts[seg] + 1):(num_pts[seg] + n_end_pts[seg]),seg] =
         end_pts_pos[seg,1:n_end_pts[seg]]';
    }
  }
}
parameters {
  real<lower=0> rho_scale;
  vector<lower=0>[L] eta_rhos;
  vector[S] eta_f;
  real<lower=0> len_scale;
  real<lower=0> alpha;
  vector[D] theta;
  vector[K-1] eta_loc;
//  vector[J-1] eta_hh;
//  real<lower=0> sigma_hh;
  vector[tot_end_pts] eta_end_pts;
//  vector[D_Z] gamma;
}
transformed parameters {
  vector[S + 2 + D + L] shared;
  vector[1] hh[J];
  vector[tot_end_pts] end_pts;
  vector[S] f;
  {
    vector[tot_end_pts] end_pts_mean;
    for (pt in comp_order) {
      if (dep_pts[pt,1] == 0) {
        end_pts[pt] = eta_end_pts[pt];
        end_pts_mean[pt] = 0;
      } else {
        if (n_dep_pts[pt] == 1) {
          int idx_dep = dep_pts[pt,1];
          real cov_exp = exp_quad(pt_dist_mat[pt,idx_dep],len_scale);
          end_pts_mean[pt] = cov_exp * (end_pts[idx_dep] - end_pts_mean[idx_dep]);
          end_pts[pt] = end_pts_mean[pt] + sqrt(1 - cov_exp^2) * eta_end_pts[pt];  
        } else if (n_dep_pts[pt] == 2) {
          int idx_dep[2] = dep_pts[pt,];
          row_vector[2] cov_exp = weighting * [exp_quad(pt_dist_mat[pt,idx_dep[1]], len_scale),
                                               exp_quad(pt_dist_mat[pt,idx_dep[2]], len_scale)];
          end_pts_mean[pt] = cov_exp * (end_pts[idx_dep] - end_pts_mean[idx_dep]);
          end_pts[pt] = end_pts_mean[pt] + sqrt(1 - dot_self(cov_exp)) * eta_end_pts[pt];
        }
      }
    }
    for (seg in 1:n_seg) {
      int num_pts_seg = num_pts[seg];
      int n_end_pts_seg = n_end_pts[seg];
      int idx_f[num_pts_seg] = canal_idx_mat[1:num_pts_seg,seg];
      if (n_end_pts_seg == 0) {
        matrix[num_pts_seg,num_pts_seg] cov_mat = cov_exp_quad(to_array_1d(canal_pos[1:num_pts_seg,seg]),
                                                 1.0, len_scale);
        for (n in 1:num_pts_seg)
          cov_mat[n,n] += 1e-6;
        f[idx_f] = exp(cholesky_decompose(cov_mat) * eta_f[idx_f]);
      }
      else {
        // Make aug_pts with end_pt_pos appended to end of 
        f[idx_f] = gen_mu(to_array_1d(aug_pts[1:(num_pts_seg + n_end_pts_seg),seg]),
                          len_scale, end_pts[end_pts_idx[seg,1:n_end_pts_seg]] - end_pts_mean[end_pts_idx[seg,1:n_end_pts_seg]],
                          eta_f[idx_f], 1.0, n_end_pts_seg);
      }
    }
  }

  shared[1:S] = f[canal_pts[1,]];
  shared[S+1] = alpha;
  shared[S+2] = rho_scale;
  shared[(S + 3):(S + 2 + D)] = theta;
  shared[(S + 2 + D + 1):(S + 2 + D + L)] = eta_rhos;
  for (j in 1:J) {
    if (kk[j] < K)
      hh[j][1] = exp(eta_loc[kk[j]]);
    else 
      hh[j][1] = 0;
  }
}
model {
  rho_scale ~ normal(0, 0.1);
  eta_rhos ~ normal(0, 1);
  eta_f ~ normal(0, 1);
  alpha ~ normal(0, 0.3);
  theta ~ normal(0, 0.3);
//  gamma ~ normal(0, 1);
  eta_loc ~ normal(0, 1);
//  eta_hh ~ normal(0, 1);
//  sigma_hh ~ normal(0, 0.5);
  len_scale ~ gamma(gam_shape, gam_rate);
  eta_end_pts ~ normal(0, 1);
  target += sum(map_rect(lp_fun, shared, hh, X_aug, y_arr));
}
generated quantities {
  vector[J] mu;
  vector[D] beta = R_inv * theta;
  vector[J] log_lik;
  int y_rep[J,max_obs] = rep_array(0,J,max_obs);
  matrix[J,max_obs] logit_p = rep_matrix(0,J,max_obs);
  int case_count_dist[n_dist_intervals] = rep_array(0, n_dist_intervals);
  vector[L] rhos = eta_rhos * rho_scale;
  for (j in 1:J) {
    real min_dist = dist[j];
    real mu_j = kernel_rvec(dists[,j], rhos[kk[j]]) * f[canal_pts[1,]] * delta_t;
    vector[N_obs[j]] log_pred = (hh[j][1] + alpha + mu_j) * exp(X[j][1:N_obs[j],] * theta);
    mu[j] = mu_j;
    logit_p[j,1:N_obs[j]] = log_pred';
    log_lik[j] = poisson_binomial_lpmf(sum(y[j,1:N_obs[j]]) | 1 - exp(- log_pred));
    y_rep[j,1:N_obs[j]] = bernoulli_rng(1 - exp(- log_pred));
    for (d in 1:n_dist_intervals) {
      if (min_dist >= lower_bounds[d] * delta_distance && min_dist < (lower_bounds[d] + 1) * delta_distance) {
        case_count_dist[d] += sum(y_rep[j,1:N_obs[j]]);
        break;
      }
    }
  }
}

