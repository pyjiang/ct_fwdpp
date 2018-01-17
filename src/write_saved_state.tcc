// created on 1.5.17
// write saved state
#ifndef _WRITE_SAVED_STATE_TCC_
#define _WRITE_SAVED_STATE_TCC_



// for basic type
// write saved state of the current population state
// note: this new pop type includes Cell_type b_vec
template<typename pop_t,
         typename Write_pop_func>
void write_saved_state_one_linked_loc(std::string out_common_name, uint ngens, const pop_t & pop,
                                      double mu, double littler, const gsl_rng * r,  const Write_pop_func & write_pop_func)
{

  // write gsl state into a separate file
  std::string out_rd_state = out_common_name + "_"+ std::to_string(ngens) +".rd";
  FILE * frd_state= fopen(out_rd_state.c_str(),"wb"); // because gsl_rng_fwrite takes only C type file stream!
  gsl_rng_fwrite(frd_state,r); // write gsl state to a separate file
  fclose(frd_state);

  //write it to a buffer:
  std::ostringstream buffer;
  // first searlize to buffer !
  // will save gsl state into a separate file!!

  // this function will take general parameters

  write_pop_func(pop.gametes, pop.mutations, pop.haploids, std::bind(CT::mutation_writer(),std::placeholders::_1,std::placeholders::_2),
                 pop.mcounts, pop.ct, ngens, mu, littler, out_common_name, pop.b_vec, buffer);



  std::string s = buffer.str();

  std::string out_state_file= out_common_name + "_"+ std::to_string(ngens) + ".state";
  CT::stringstream_to_gz(out_state_file,s);

}

// for all two locus case
// modified on 1.9.17, the state file will have the index of how many generations have been simulated
template<typename pop_t,
         typename Write_pop_func>
void write_saved_state_two_loc(std::string out_common_name, uint ngens, const pop_t & pop,
                               double mu1, double mu2, double recomb, const gsl_rng * r, const uint32_t robust_bits,  const int robust_flag, // flag of robustness
                               const Write_pop_func & write_pop_func)
{

  // write gsl state into a separate file
  std::string out_rd_state = out_common_name + "_"+ std::to_string(ngens) +".rd";
  FILE * frd_state= fopen(out_rd_state.c_str(),"wb"); // because gsl_rng_fwrite takes only C type file stream!
  gsl_rng_fwrite(frd_state,r); // write gsl state to a separate file
  fclose(frd_state);

  //write it to a buffer:
  std::ostringstream buffer;
  // first searlize to buffer !
  // will save gsl state into a separate file!!

  // this function will take general parameters
  write_pop_func(pop.gametes1, pop.gametes, pop.mutations, pop.haploids, std::bind(CT::mutation_writer(),std::placeholders::_1,std::placeholders::_2),
                       pop.mcounts, pop.ct, ngens, mu1, mu2,  recomb, out_common_name, pop.b_vec, robust_bits, robust_flag, buffer);


  std::string s = buffer.str();

  std::string out_state_file= out_common_name + "_"+ std::to_string(ngens) + ".state";
  CT::stringstream_to_gz(out_state_file,s);

}


// modified on 2.21.17
// will also write robustness associated with every fixed mutations (in order to see distribution)
// note: now the mutation type needs to have scale (it is not back compatible for now)
// modified on 2.3.17
// this version should work for both one locus and two locus
// this will output the position in double (rather than round) : will be converted outside of file
template<typename pop_t>
void write_fix(std::string out_common_name, const pop_t & pop)
{
  // output fixation
  std::ofstream out_fix;
  std::string fixation_file= out_common_name  + "_fix.txt";
  out_fix.open(fixation_file,std::ofstream::out | std::ofstream::app);

  // uint Lv = pop.ct.v.size();

  for(uint i=0; i< pop.fixations.size();i++)
  {
    out_fix << pop.fixation_times[i] << "," << std::setprecision(12) << pop.fixations[i].pos << ","  << pop.fixations[i].scale << std::endl ;
  }
  out_fix.close();
}


#endif
