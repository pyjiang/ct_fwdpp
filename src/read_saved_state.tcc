// modified : read saved state also for one big locus
#ifndef _READ_SAVED_STATE_TCC_
#define _READ_SAVED_STATE_TCC_

// #include "read_saved_state.hpp"
#include "ct_IO_high_level.hpp"
#include "ct_serialization.hpp"
#include "ct_IO.tcc"


// for basic type
// modifed on 3.27.17
// add keep_gamete_flag
template <typename pop_t,
          typename Read_pop_func>
pop_t read_saved_state_one_linked_loc(std::string out_common_name, std::string prev_ngen_s,  gsl_rng * r, double & mu,
                                      std::size_t & prev_ngen, int keep_gamete_flag, const Read_pop_func & read_pop_func)
{
  std::string in_state_file= out_common_name + "_"+ prev_ngen_s +".state";
  std::string rd_state= out_common_name + "_"+ prev_ngen_s + ".rd";


  FILE * frd_state= fopen(rd_state.c_str(),"rb");

  gsl_rng_fread(frd_state,r); // read in random number state
  fclose(frd_state);

  //create a string to hold data
  std::string data;

  CT::gz_to_string(in_state_file,data);

  std::istringstream in(data); // initialize istringstream obj


  // create empty gametes, mutations, diploids
  typename pop_t::gcont_t gametes;
  typename pop_t::mcont_t mutations;
  typename pop_t::dipvector_t haploids;
  std::vector<uint>  mcounts;
  Cell_type ct;
  // std::size_t prev_ngen;
  double  recomb;
  ArrayXb b_vec;

  std::string old_out_common_name; // read in old out common name, but not use it !


  read_pop_func(gametes, mutations, haploids, std::bind(CT::mutation_reader<mtype>(),std::placeholders::_1),
                mcounts, ct, &prev_ngen, &mu, &recomb, old_out_common_name, b_vec, in, keep_gamete_flag);
    // need to calculate z_val for cell type (because it only recovers v, w, h)
    // ct.z_val = calc_cell_type_z_val(ct.w, ct.v, ct.h);
    // ct.wv= calc_cell_type_wv(ct.w,ct.v,Lv-robust_bits).array();


    // need to initialize pop with gametes, mutations etc.
    pop_t pop(gametes,mutations,mcounts,haploids,b_vec,ct);

    return pop;

}



// template <typename pop_t,
//           typename Read_pop_func>
// pop_t read_saved_state_one_linked_loc_2env(std::string out_common_name, std::string prev_ngen_s,  gsl_rng * r, double & mu, std::size_t & prev_ngen, const Read_pop_func & read_pop_func)
// {
//   std::string in_state_file= out_common_name + "_"+ prev_ngen_s +".state";
//   std::string rd_state= out_common_name + "_"+ prev_ngen_s + ".rd";
//
//
//   FILE * frd_state= fopen(rd_state.c_str(),"rb");
//
//   gsl_rng_fread(frd_state,r); // read in random number state
//   fclose(frd_state);
//
//   //create a string to hold data
//   std::string data;
//
//   CT::gz_to_string(in_state_file,data);
//
//   std::istringstream in(data); // initialize istringstream obj
//
//
//   // create empty gametes, mutations, diploids
//   typename pop_t::gcont_t gametes;
//   typename pop_t::mcont_t mutations;
//   typename pop_t::dipvector_t  haploids;
//   std::vector<uint>  mcounts;
//   Cell_type ct;
//   // std::size_t prev_ngen;
//   double  recomb;
//   ArrayXb b_vec;
//
//   std::string old_out_common_name; // read in old out common name, but not use it !
//
//
//   ArrayXb b2_vec;
//
//   read_pop_func(gametes, mutations, haploids, std::bind(CT::mutation_reader<mtype>(),std::placeholders::_1),
//                 mcounts, ct, &prev_ngen, &mu, &recomb, old_out_common_name, b_vec, b2_vec, in, 1);
//
//
//   // need to calculate z_val for cell type (because it only recovers v, w, h)
//   // ct.z_val = calc_cell_type_z_val(ct.w, ct.v, ct.h);
//   // ct.wv= calc_cell_type_wv(ct.w,ct.v,Lv-robust_bits).array();
//
//
//   // need to initialize pop with gametes, mutations etc.
//   pop_t pop(gametes, mutations, mcounts, haploids, b_vec, b2_vec, ct);
//
//   return pop;
//
//
//
// }


// read saved state for all two locus cases
// modified on 2.6.17
// add an extra parameter of recomb, indicates the rate of recombination between the two locus
// modified on 1.11.17
// pass robust bits and robust_flag as parameter
template <typename pop_t,
          typename Read_pop_func>
pop_t read_saved_state_two_loc(std::string out_common_name, std::string prev_ngen_s, gsl_rng * r, double & mu1, double & mu2, double & recomb,
                               std::size_t & prev_ngen,
                               uint32_t & robust_bits, int & robust_flag, const Read_pop_func & read_pop_func)
{
  std::string in_state_file= out_common_name + "_"+ prev_ngen_s +".state";
  std::string rd_state= out_common_name + "_"+ prev_ngen_s + ".rd";


  FILE * frd_state= fopen(rd_state.c_str(),"rb");

  // gsl_rng * rtest = (gsl_rng *)malloc(fsize);
  gsl_rng_fread(frd_state,r); // read in random number state
  fclose(frd_state);

  //create a string to hold data
  std::string data;

  CT::gz_to_string(in_state_file,data);

  std::istringstream in(data); // initialize istringstream obj


  // create empty gametes, mutations, diploids
  typename pop_t::gcont1_t gametes1; // gametes for locus 1
  typename pop_t::gcont2_t gametes2; // gametes for locus 2
  typename pop_t::mcont_t mutations;
  typename pop_t::dipvector_t haploids;
  std::vector<uint>  mcounts;
  Cell_type ct;
  ArrayXb b_vec;
  // uint32_t robust_bits;
  // int robust_flag; // added on 12.7.16


  std::string old_out_common_name; // read in old out common name, but not use it !
  read_pop_func(gametes1, gametes2, mutations, haploids, std::bind(CT::mutation_reader<mtype>(),std::placeholders::_1),
                mcounts, ct, &prev_ngen, &mu1, &mu2, &recomb, old_out_common_name, b_vec, &robust_bits, &robust_flag, in);



  // need to initialize pop with gametes, mutations etc.
  pop_t pop(gametes1, gametes2, mutations, mcounts, haploids, b_vec, ct);

  return pop;
  // return ct;
}

#endif


// int main()
// {
//   // read_saved_state_two_loc("two_locus/two_locus_rand_nz1000_lv10000_ga1_500_1","8799999");
//   read_saved_state_one_linked_loc("test_hap_robust_bits_v21/with_h_hap_rand_nz_1000_rvz_10_ga_1_500_7","8799999");
// }
