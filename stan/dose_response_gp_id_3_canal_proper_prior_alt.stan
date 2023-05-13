functions {
  real kernel(real d, real rho) {
    return exp(-d / rho);
  }
  real exp_quad(real d, real ell) {
    return exp(-0.5 / ell ^ 2 * d ^ 2);
  }
  row_vector kernel_vec(vector d, real rho) {
    int dim_d = rows(d);
    row_vector[dim_d] kernels;
    for (i in 1 : dim_d) 
      kernels[i] = kernel(d[i], rho);
    return kernels;
  }
  matrix kernel_mat(matrix d, real rho) {
    int rows_d = rows(d);
    int cols_d = cols(d);
    matrix[cols_d, rows_d] kernels;
    for (j in 1 : cols_d) {
      kernels[j,  : ] = kernel_vec(d[ : , j], rho);
    }
    return kernels;
  }
  array[] matrix kernel_list(array[] matrix d, real rho) {
    int N_seg = size(d);
    array[N_seg] matrix[cols(d[1]), rows(d[1])] kernels;
    for (seg in 1 : N_seg) 
      kernels[seg] = kernel_mat(d[seg], rho);
    return kernels;
  }
  vector integrate(matrix K, vector f, real dx) {
    return K * f * dx;
  }
  matrix inv_22(real off_diag) {
    matrix[2, 2] K_inv;
    K_inv[1, 1] = 1;
    K_inv[2, 2] = 1;
    K_inv[1, 2] = -off_diag;
    K_inv[2, 1] = K_inv[1, 2];
    return K_inv * inv(1 - off_diag ^ 2);
  }
  vector gen_mu(array[] real pts, real ell, vector cond_vals, vector eta,
                real mag) {
    int S = size(pts);
    int Sm2 = S - 2;
    matrix[S, S] K = add_diag(gp_exp_quad_cov(pts, mag, ell), 1e-6);
    matrix[2, 2] K_22_inv = inv_22(K[S, S - 1]);
    matrix[Sm2, 2] K_12 = block(K, 1, S - 1, Sm2, 2);
    matrix[Sm2, Sm2] K_cond = block(K, 1, 1, Sm2, Sm2)
                              - quad_form(K_22_inv, K_12'); // * inv_mag_sq;
    for (s in 1 : Sm2) {
      K_cond[s, s] += 1e-6;
    }
    return exp(K_12 * K_22_inv * cond_vals + cholesky_decompose(K_cond) * eta); // * inv_mag_sq
  }
}
data {
  int<lower=1> J; // Number of households
  int<lower=1> S; // Number of integration points
  int<lower=1> G; // Number of observations per household
  int<lower=1> N_seg;
  int<lower=1> N_endpts;
  array[J, G] int<lower=0, upper=1> y; // Observations
  array[N_seg] matrix[S, J] dists; // HH dist to integration points 
  real M;
  real C;
  array[N_seg, S] real pts;
  array[N_seg] real delta_t;
  matrix[N_seg, 2] endpt_coords;
  array[N_seg, 2] int<lower=1, upper=N_endpts> endpt_inds;
  real<lower=0> rho_scale;
  vector[J] x;
  real alpha_loc;
  real<lower=0> alpha_scale;
  real beta_loc;
  real<lower=0> beta_scale;
  real<lower=0> ell_shape;
  real<lower=0> ell_rate;
}
transformed data {
  array[J] int<lower=0, upper=G> Y;
  // int ell_idx[N_seg];
  real scale = 1.0;
  array[N_seg, S + 2] real pts_aug;
  real weighting = sqrt(0.5);
  for (j in 1 : J) 
    Y[j] = sum(y[j,  : ]);
  for (seg in 1 : N_seg) {
    for (s in 1 : S) {
      pts_aug[seg, s] = pts[seg, s];
    }
    pts_aug[seg, (S + 1) : (S + 2)] = to_array_1d(endpt_coords[seg,  : ]);
  }
  // ell_idx[1:2] = rep_array(1,2);
  // ell_idx[3:4] = rep_array(2,2);
  // ell_idx[5:6] = rep_array(3,2);
}
parameters {
  //  real<lower=0> ell[N_seg];
  real<lower=0> ell;
  real<lower=0> rho;
  matrix[S, N_seg] eta_f;
  real<lower=0> alpha;
  real beta;
  vector[N_endpts] eta_end_pts;
}
transformed parameters {
  vector[J] mu = rep_vector(0, J);
  array[N_seg] vector[S] f;
  vector[N_endpts] end_pts;
  {
    vector[N_endpts] end_pts_mean;
    row_vector[2] c_5_67 = weighting
                           * [exp_quad(M / 2., ell), exp_quad(C / 3, ell)];
    row_vector[2] c_2_35 = weighting
                           * [exp_quad(M / 2., ell),
                              exp_quad(2. * C / 3, ell)];
    real c_1_2 = exp_quad(M / 2., ell);
    real c_4_5 = exp_quad(M / 2., ell);
    end_pts_mean[6] = 0;
    end_pts_mean[7] = 0;
    end_pts_mean[3] = 0;
    end_pts[6 : 7] = eta_end_pts[6 : 7] + end_pts_mean[6 : 7];
    end_pts[3] = eta_end_pts[3] + end_pts_mean[3];
    end_pts_mean[5] = c_5_67 * (end_pts[6 : 7] - end_pts_mean[6 : 7]);
    end_pts[5] = end_pts_mean[5] +
                 + sqrt(1 - dot_self(c_5_67)) * eta_end_pts[5];
    end_pts_mean[2] = c_2_35 * [end_pts[3] - end_pts_mean[3], end_pts[5] - end_pts_mean[5]]';
    end_pts[2] = end_pts_mean[2]
                 + sqrt(1 - dot_self(c_2_35)) * eta_end_pts[2];
    end_pts[1] = c_1_2 * (end_pts[2] - end_pts_mean[2])  + sqrt(1 - c_1_2 ^ 2) * eta_end_pts[1];
    end_pts[4] = c_4_5 * (end_pts[5] - end_pts_mean[5]) + sqrt(1 - c_4_5 ^ 2) * eta_end_pts[4];
  }
  {
    array[N_seg] matrix[J, S] k_dist = kernel_list(dists, rho);
    for (seg in 1 : N_seg) {
      //  real ell_use = ell[ell_idx[seg]];
      f[seg] = gen_mu(pts_aug[seg], ell, end_pts[endpt_inds[seg,  : ]],
                      eta_f[ : , seg], 1.0);
      mu += integrate(k_dist[seg], f[seg], delta_t[seg]);
    }
  }
}
model {
  for (j in 1 : J) {
    real risk = exp(x[j] * beta) * (alpha + mu[j]) * scale;
    target += Y[j] * log1m_exp(-risk) + (G - Y[j]) * -risk;
  }
  to_vector(eta_f) ~ std_normal();
  rho ~ normal(0, rho_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  beta ~ normal(beta_loc, beta_scale);
  ell ~ gamma(ell_shape, ell_rate);
  eta_end_pts ~ std_normal();
}
generated quantities {
  vector[J] total_risk;
  for (j in 1 : J) 
    total_risk[j] = exp(x[j] * beta) * (alpha + mu[j]) * scale;
}
