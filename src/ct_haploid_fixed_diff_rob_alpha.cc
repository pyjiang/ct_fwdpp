// created on 5.23.17
// will read in the initial state generated from ct_haploid_fixed_rob_alpha_init, so that initial population has the same genotype, phenotype
// then this file will apply different scales (arbitrary) for the same genotype to phenotype settings (because want to test a larger range of robustness values)


// modifed on 2.3.17
// take parameter scale_flag instead of previous scale_func_flag
// update on 2.1.17
// should be compatible with the new version (no initial scale), and initial w matrix is given by gamma
// this version will work with shorter bits (10 bits ) that encodes robustness
// and increased mutation rates on these bits
// this version will use wij
// in stead of with h
// created on 1.23.17
// to initialize the initial population having 3 different degrees of robustnes
// 1/3 have 1/2 gamma, 1/3 have gamma, and 1/3 have 2*gamma
// modifed on 1.20.17
// add gamma flag
// if gamma_flag =0 : pass gamma as a parameter
// if gamma_flag =1 :  pass index of gamma (because list of gamma is saved )
// created on 1.12.17
// a new version for haploid simulation, without robust bits

#include <vector>
#include <iterator>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/debug.hpp>
// #include <Sequence/SimData.hpp>



#include "sample_haploid.tcc"
#include "ct_mut_model.hpp"
#include "ct_util.hpp"
// #include "ct_fitness_models.hpp" (does not use this)
#include "ct_fitness_policy.tcc"
#include "ct_mutation_internal.hpp"
#include "ct_mutation.tcc"

#include "ct_eigen.hpp"


#include "initialize_pop.tcc"
#include "write_saved_state.hpp"

#include "robust_dist.tcc"
#include "ct_multiloc.hpp"
#include "read_saved_state.hpp"

// #include "static_member.hpp"
#include "ct_mutation.hpp"
#include <iomanip>

using mtype = CT::ct_mutation ;
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
// using single_pop_h_w_b_t = CT::singlepop_h_w_b_t<mtype>; // use singlepop_h_w_b_t type (because need h vector), but set robust bits =0
// using single_pop_wij_b_t = CT::singlepop_wij_b_t<mtype>; // without h in gamete

using twolocus_pop_no_h_t = CT::twolocus_no_h_t<mtype>;

extern double LOW_SCALE_BOUND;
extern double HIGH_SCALE_BOUND;


int main(int argc, char ** argv)
{

  if (argc != 6)
    {
      // change parameter from rvz to Lv (total length of v); and also add a parameter of robust_bits
      // do not use gamma to initialze
      // take parameters to set low and high scale boundary
      // calc_scale_flag ==0: additive calc_scale_func; calc_scale_func ==1: threshold (only for 10 bits)
      // modified on 7.11.17
      // if seed = 0, it means keep the previous state and keep running
      // modified on 7.12.17
      // add another argument (input both prev_ngens and ngens)
      std::cerr << "Too few or too more arguments\n"
    << "Usage: diploid_ind  prev_ngens  ngens seed   saved_state_common_name      scale    \n"; // type=0: stabil, type=1: randinit
      exit(10);
    }

    // read in arguments
    int argument=1;

    const unsigned prev_ngens = unsigned(atoi(argv[argument++]));       //Number of generations that have simulated

    const unsigned ngens = unsigned(atoi(argv[argument++]));       //Number of generations to simulate
    // int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate

    const unsigned seed = unsigned(atoi(argv[argument++]));        //Random number seed

    // const double mu = theta/double(4*N);                 //per-gamete mutation rate

    std::string in_common_name= argv[argument++];

    std::string scale_str = argv[argument++];
    double scale =  stod(scale_str); // scale value

    // if seed ==0 (continue running, will output to the same mean fit file)
    // else, will add scale to output file
    std::string out_common_name ;
    if (seed ==0)
    {
      out_common_name = in_common_name ;
    }
    else
    {
      out_common_name = in_common_name + "_" + scale_str ;
    }

    //Write the command line to stderr
    std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
    std::cerr << '\n';


    //Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)

    double mu1, mu2, rho;
    uint32_t robust_bits;
    int robust_flag;

    // for test, use old r for now
    gsl_rng * r_old = gsl_rng_alloc (gsl_rng_mt19937); // will read from state, but not be used in the following

    GSLrng r(seed);

    std::size_t prev_ngen;

    // modified on 7.12.17, read from previous saved state (not necessarily gen=0)
    // read from previous state ngen=0
    twolocus_pop_no_h_t pop= read_saved_state_two_loc<twolocus_pop_no_h_t>(in_common_name, std::to_string(prev_ngens), r_old, mu1, mu2, rho, prev_ngen, robust_bits, robust_flag,
                            std::bind(CT::read_binary_pop_hap(),std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                            std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                            std::placeholders::_9,std::placeholders::_10,  std::placeholders::_11, std::placeholders::_12, std::placeholders::_13,
                            std::placeholders::_14, std::placeholders::_15, std::placeholders::_16 ));

    // give the option to use the same seed or change seed
    const gsl_rng * r_use = (seed==0) ? r_old : r.get(); // the pointer of r in the model

    unsigned N = pop.haploids.size()/2;

    std::size_t Lv = pop.ct.v.size() - robust_bits;

    // then change scale factor in the gamete
    pop.gametes1[0].scale = scale;

    // write file handle for mean fitness and mean robustness
    std::string out_mean_fit_file = out_common_name +"_mean_fit.txt";

    std::ofstream out_mean_fit(out_mean_fit_file,std::ofstream::out | std::ofstream::app);

    unsigned generation;
    double wbar;
    // double mean_robust;

    unsigned twoN = 2*N;


    // store alpha vector
    // std::vector<double> alpha_vector = {pop.gametes[0].alpha, pop.gametes[1].alpha, pop.gametes[2].alpha };

    // output alpha list as the first line!
    // out_rob_dist << alpha_vector[0] << "\t" << alpha_vector[1] << "\t" << alpha_vector[2] << "\t" << std::endl;

    VectorXi v1 ; //dummy
    // VectorXi v2= pop.ct.v;

    // simulations starts here!
    for( generation = prev_ngens; generation < prev_ngens + ngens; ++generation )
    {

        wbar= CT::sample_haploid(r_use,
                          				pop.gametes1,
                                  pop.gametes,
                          				pop.haploids,
                          				pop.mutations,
                          				pop.mcounts,
                          				N,
                                  N,
                          				0, // no mutation on robustness genotype
                                  mu2,
                                  std::bind(CT::ct_mut_model(),std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup),
                                [r_use](){return ((double) gsl_rng_uniform_int(r_use,1)) / 1 ;}), //dummy (placeholder)
                                  std::bind(CT::ct_mut_model(),std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup),
                                [r_use,Lv](){return ((double) gsl_rng_uniform_int(r_use,Lv)) / Lv +1  ;}),
                          			  std::bind(CT::ct_fit_simple(),std::placeholders::_1, pop.ct.h, pop.b_vec, 1),
                          				// pop.neutral,
                          				// pop.selected,
                                  pop.ct.w,
                                  v1,
                                  pop.ct.v,
                                  rho,
                                  0,
                                  0,
                                  0);

        CT::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,twoN,v1,pop.ct.v);

        assert(KTfwd::check_sum(pop.gametes,twoN));

        out_mean_fit << generation << "," << wbar << std::endl;

      }
      // close mean fit file
      out_mean_fit.close();
      // out_rob_dist.close();

      // // now needs to copy v1 and v2 back to ct.v (so that it will be easier to write to file)
      // copy_back( pop.ct.v,v2);

      write_fix(out_common_name,pop);

      write_saved_state_two_loc(out_common_name, prev_ngens + ngens, pop, 0, mu2, rho, r_use, 0, 0,
                                      std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                      std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                      std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12, std::placeholders::_13,
                                      std::placeholders::_14, std::placeholders::_15, std::placeholders::_16  ) );

      return 0;


}
