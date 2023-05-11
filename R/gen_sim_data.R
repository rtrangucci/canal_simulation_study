source("fun_sim_data.R")
S <- 100
G <- 3
n <- 100
l <- 0.1
house_dist <- 'uniform'
gen_full_data_sets(S = S, G = G, n = n, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = G, n = 500, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = G, n = 1000, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 10, n = n, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 10, n = 500, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = 10, n = 1000, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = 100, n = n, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 100, n = 500, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 100, n = 1000, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets')

gen_full_data_sets(S = S, G = G, n = n, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = G, n = 500, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = G, n = 1000, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 10, n = n, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 10, n = 500, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = 10, n = 1000, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = 100, n = n, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 100, n = 500, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets')
gen_full_data_sets(S = S, G = 100, n = 1000, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets')

gen_full_data_sets(S = S, G = 10, n = 150, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = 10, n = 150, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets') #

gen_full_data_sets(S = S, G = 6, n = 500, l = l, L = 50, house_dist = house_dist, par_dir = 'fake-data-sets') #
gen_full_data_sets(S = S, G = 6, n = 500, l = l, L = 50, house_dist = 'clustered_uniform', par_dir = 'fake-data-sets') #

dat <- gen_data(S=100, G=10, seed = 123, int = 0.05, n=100, l=0.1,house_dist = 'uniform')
