## set up runs to evolve different robustness initial sub populations with 8 degrees of robustness

## 8.8.17
## randinit, no recombination
## note: robust_flag set to 0. (so that each genotype represents one scale, which case mutation is not additive)
bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-5 0 0 8  diff_init_rob_evol_8_degree

## if robust_flag set to 1 (so that when mutate, each robustness genotype can only +1 or -1 at one mutation)
## the correct additive case (using 7 bit)
bash set_init_diff_rob_evol.sh  8000000  1 1 1 7 0 1e-5 0 0 7  diff_init_rob_evol_8_degree

## 8.9.17
## set up for 16 degrees, either using 4 bits or using 15 bits additively
bash set_init_diff_rob_evol.sh  8000000  1 1 0 4 0 1e-5 0 0 16  diff_init_rob_evol_16_degree

## additive
bash set_init_diff_rob_evol.sh  8000000  1 1 1 15 0 1e-5 0 0 15  diff_init_rob_evol_16_degree


## will use 3 bits for 8 states
## and change mutation rates (higher or lower)

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-4 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-3 0 0 8  diff_init_rob_evol_8_degree



bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-2 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 0.1 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-6 0 0 8  diff_init_rob_evol_8_degree


bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-7 0 0 8  diff_init_rob_evol_8_degree

## just set up on 8.13.17
bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0 1e-8 0 0 8  diff_init_rob_evol_8_degree



## for recombination (note: when rate is 0.3, it will take longer to run)
bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.1 1e-5 0 0 8  diff_init_rob_evol_8_degree


bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.2 1e-5 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.3 1e-5 0 0 8  diff_init_rob_evol_8_degree


## to do
## the next two are waiting to be run



## change for both recomb and mutation rate (not yet set up)
bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.1 1e-6 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.1 1e-7 0 0 8  diff_init_rob_evol_8_degree

## set up on 9.8.17 (set up)
bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.2 1e-6 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.2 1e-7 0 0 8  diff_init_rob_evol_8_degree

bash set_init_diff_rob_evol.sh  8000000  1 1 0 3 0.2 1e-8 0 0 8  diff_init_rob_evol_8_degree

## for stabil
bash set_init_diff_rob_evol.sh  8000000  0 1 0 3 0 1e-5 0 0 8  diff_init_rob_evol_8_degree
