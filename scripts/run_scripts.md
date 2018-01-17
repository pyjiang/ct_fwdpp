This document shows how to run simulations using these scripts (set shared parameters for a set of runs with different random number seeds). The parameters pass through these scripts will pass to the c++ executable as explained in the README.md file, other non-specified parameters will use default values as in the main text.



### Basic simulations with gamma
```bash
bash  basic_runs.sh  gamma stype
```
This will generate 20 runs for the basic simulation with specified gamma value and stype (with either optimal initial environment stype=0, or with a random initial environment stype=1). Other parameters set with default values as specified in the paper. The random seed being used will be output into the log file.



### Simulations with fixed alpha
#### First need to generate initial populations
```bash
bash set_init_pop_for_scale.sh  recomb_rate  stype ga index_begin index_end
```
  - index_begin: the beginning index of  simulated populations
  - index_end: the ending index of  simulated populations; for generating 20 simulations, index_begin=1, index_end=20

#### Then generate populations with different alpha values from each initial populations
```bash
bash set_simu_diff_alpha.sh   ngen   stype index_begin index_end  init_ga recomb_rate alpha_list
```
  - init_ga: the gamma value used for the initial population
  - alpha_list: 
