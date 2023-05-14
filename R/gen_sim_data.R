source("R/fun_sim_data.R")
S <- 20
G <- 3
n <- 100
l <- 0.1
L <- 200
set.seed(123)
house_locs_unif <- gen_house_locs(n = 2000)
house_locs_clustered <- gen_house_locs(n = 2000, house_dist = "clustered_uniform")
gen_full_data_sets(S = S, G = 3, J_sched = c(100,500,1000,2000), l = l, L = L, house_locs = house_locs_unif, par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 10, J_sched = c(100, 500, 1000, 2000), l = l, L = L, house_locs = house_locs_unif, par_dir = "fake-data-sets")
gen_full_data_sets(S = S, G = 100, J_sched = c(100, 500, 1000, 2000), l = l, L = L, house_locs = house_locs_unif, par_dir = "fake-data-sets")


gen_full_data_sets(S = S, G = 3, J_sched = c(100, 500, 1000, 2000), l = l, L = L, house_locs = house_locs_clustered, par_dir = "fake-data-sets")
gen_full_data_sets(S = S, G = 10, J_sched = c(100, 500, 1000, 2000), l = l, L = L, house_locs = house_locs_clustered, par_dir = "fake-data-sets")
gen_full_data_sets(S = S, G = 100, J_sched = c(100, 500, 1000, 2000), l = l, L = L, house_locs = house_locs_clustered, par_dir = "fake-data-sets")
