#ifndef _READ_SAVED_STATE_H_
#define _READ_SAVED_STATE_H_

#include <vector>
#include "ct_multiloc.hpp"
#include "ct_IO.tcc" // include IO
#include "ct_eigen.hpp"
#include <iostream>
#include <fstream> //ofstream
#include <sstream>
#include "ct_IO_high_level.hpp"
#include "ct_serialization.hpp"
#include "ct_singlepop.hpp"

#include "ct_mutation.hpp"

// temporarily change to mutation_base for compatible with basic type (hack: in order to run basic type)
// using mtype = KTfwd::mutation_base;
using mtype = CT::ct_mutation ;

// for two loc
// using twolocus_pop_t = CT::multiloc_t<mtype>;

// for current verion of pybind (v2.0), it needs to specify base class as well
// using popbase_2loc_pop_t = CT::popbase_2loc_t<mtype>;


using popbase_pop_two_t = CT::popbase_two_t <mtype>;


// create a new struct that has two pop and also b vector. [for interface output]
// using twolocus_pop_w_b_t = CT::twolocus_w_b_t<mtype>;

// twolocus_pop_w_b_t read_saved_state_two_loc(std::string out_common_name, std::string prev_ngen_s);


// for one linked loc
// using popbase_pop_t = CT::popbase_t<mtype>;

// using singlepop_t = CT::ct_singlepop_with_h<mtype>;

// using single_pop_h_w_b_t = CT::singlepop_h_w_b_t<mtype>;
// create a new struct that has two pop and also b vector. [for interface output]

template <typename pop_t,
          typename Read_pop_func>
pop_t read_saved_state_one_linked_loc(std::string out_common_name, std::string prev_ngen_s, gsl_rng * r , double & mu, std::size_t & prev_ngen,  int flag, const Read_pop_func & read_pop_func);

// template <typename pop_t,
//           typename Read_pop_func>
// pop_t read_saved_state_one_linked_loc_2env(std::string out_common_name, std::string prev_ngen_s,  gsl_rng * r, double & mu, std::size_t & prev_ngen, uint32_t & robust_bits, int & robust_flag,
//                                             const Read_pop_func & read_pop_func);


template <typename pop_t,
          typename Read_pop_func>
pop_t read_saved_state_two_loc(std::string out_common_name, std::string prev_ngen_s, gsl_rng * r, double & mu1, double & mu2, double & recomb,
                               std::size_t & prev_ngen,
                               uint32_t & robust_bits, int & robust_flag, const Read_pop_func & read_pop_func);

#endif
#include "read_saved_state.tcc"
