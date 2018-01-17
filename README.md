
Note: to display math symbols, please download chrome plug-in [GitHub with MathJax](https://chrome.google.com/webstore/detail/github-with-mathjax/ioemnmodlmafdkllaclgeombjnmnbima?hl=en)

## Build requirement
### packages: [fwdpp](https://github.com/molpopgen/fwdpp), [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- I use fwdpp version 0.5.1 (with fixing a small bug in the #ifndef statement in fwdpp/tags/gamete_tags.hpp; alternatively, a version with this bug fixed is 0.5.3), eigen 3.3.1.

## Installation
- Compiler should support C++11
- Download the required packages (only C++ header files are needed)
- Edit Makefile in this directory (change CPPFLAGS -I to the package directory)
- `make all` to make all the targets
- or optionally, only make specific target
  +  `make haploid_save_new` to make the target for basic simulations with a fixed level of robustness with  $\alpha$
  + `make haploid_diff_init_rob_evol` to make the target for simulations which can evolve robustness with $\gamma$

## Run
- For basic simulations with a fixed level of robustness with  $\alpha$

  usage:

```bash
  haploid_save_new  N  theta  ngens seed nz Lv  c pv1 ga_flag ga stype out_common_name
```
  + N: population size (as if the individuals are diploids)
  + theta: 2Nu (u is the mutation rate per individual per generation)
  + ngens: how many generations of forward simulations to run
  + seed: the random number seed used to initialize parameters
  + nz: length of z vector (phenotype vector)
  + Lv: length of v vector (genotype vector)
  + c: the percentage of non-zero in the T matrix
  + pv1: percent of 1 in the v vector
  + ga_flag: flag for $\gamma$. if gamma_flag =0 : pass $\gamma$ (next argument) as a parameter;
    if gamma_flag =1 :  pass index of $\gamma$ in the list (the list contains 10 $\gamma$ which are equally separated by their robustness value)
  + ga: $\gamma$ value
  + stype: which type of simulation to run. if type=0: with optimal initial environment; if type=1, with a random initial environment
  + out_common_name: the prefix for the output files.

- For simulations which can evolve robustness with $\gamma$

  usage:

```bash
  haploid_diff_init_rob_evol N theta rmb ngens  seed mu1_base nz Lv  robust_bits  c pv1 ga stype out_common_name   robust_flag  low_bound high_bound  scale_flag
```

  + N: population size (as if the individuals are diploids)
  + theta: 2Nu (u is the mutation rate per individual per generation)
  + rmb: a value between 0 and 0.5, where 0.5 means random shuffling of robustness locus and the regular genotype locus, and 0 means complete linkage
  + ngens: how many generations of forward simulations to run
  + seed: the random number seed used to initialize parameters
  + mu1_base: mutation rate per bit for the robustness locus
  + nz: length of z vector (phenotype vector)
  + Lv: length of v vector (genotype vector)
  + robust_bits: the number of bits that encodes robustness locus
  + c: the percentage of non-zero in the T matrix
  + pv1: percent of 1 in the v vector
  + ga: $\gamma$ value
  + stype: which type of simulation to run. if type=0: with optimal initial environment; if type=1, with a random initial environment
  + out_common_name: the prefix for the output files
  + robust_flag: if robust_flag =1: evolve robustness with the robustness locus, if robust_flag=2, do not evolve robustness locus
  + low_bound: a number that is the log2 lower bound of the $\alpha$ value
  + high_bound:  a number that is the log2 higher bound of the $\alpha$ value
  + scale_flag: specify how robustness genotype maps to robustness set by $\alpha$ value
    if scale_flag = 8, it is exactly as described in the paper.

### Run with scripts
Alternatively, you can run with scripts that are in ./log directory.
scripts: generates the same simulation results as in the paper.
or you can run with different seed values with
