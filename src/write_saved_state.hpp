// header file for write_saved_state
#ifndef _WRITE_SAVED_STATE_H_
#define _WRITE_SAVED_STATE_H_

#include "ct_serialization.hpp"
#include "ct_IO.tcc"
#include <fstream>
#include "ct_IO_high_level.hpp"
#include <iomanip>

template<typename pop_t,
         typename Write_pop_func>
void write_saved_state_one_linked_loc(std::string out_common_name, uint ngens, const pop_t & pop,
                                               double mu, double littler, const gsl_rng * r, const Write_pop_func & write_pop_func);

template<typename pop_t,
        typename Write_pop_func>
void write_saved_state_two_loc(std::string out_common_name, uint ngens, const pop_t & pop,
                              double mu1, double mu2, double recomb, const gsl_rng * r, const uint32_t robust_bits,  const int robust_flag, // flag of robustness
                              const Write_pop_func & write_pop_func);

template<typename pop_t>
void write_fix(std::string out_common_name, const pop_t & pop);

#endif
#include "write_saved_state.tcc"
