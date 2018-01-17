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

  if (argc != 13)
    {
      // change parameter from rvz to Lv (total length of v); and also add a parameter of robust_bits
      // do not use gamma to initialze
      // take parameters to set low and high scale boundary
      // calc_scale_flag ==0: additive calc_scale_func; calc_scale_func ==1: threshold (only for 10 bits)
      std::cerr << "Too few or too more arguments\n"
    << "Usage: diploid_ind N theta rho   seed  nz Lv  c pv1  ga type out_common_name  scale    \n"; // type=0: stabil, type=1: randinit
      exit(10);
    }

    // read in arguments
    int argument=1;
    const unsigned N = unsigned(atoi(argv[argument++]));           //Number of diploids
    const double theta = atof(argv[argument++]);         //4*n*mutation rate.  Note: mutation rate is per REGION, not SITE!!
    const double rho = atof(argv[argument++]);           // use this parameter to pass in proportion of shuffled gametes
    // const unsigned ngens = unsigned(atoi(argv[argument++]));       //Number of generations to simulate
    // int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
    const unsigned seed = unsigned(atoi(argv[argument++]));        //Random number seed

    // const double mu = theta/double(4*N);                 //per-gamete mutation rate


    uint32_t nz= unsigned(atoi(argv[argument++])); // number of cell types
    // double rvz=  atof(argv[argument++]);    //ratio of length of v to length of z
    uint32_t Lv= unsigned(atoi(argv[argument++]));


    double c= atof(argv[argument++]); // ratio of non-zeros in w matrix
    double pv1= atof(argv[argument++]); // ratio of 1s in v

    double  ga =  atof(argv[argument++]);

    const double mu2 = theta/double(4*N); // mutation rate per locus for regulation v


    int type= atoi(argv[argument++]); // which type of simulation (randint or stabil)
    std::string out_common_name= argv[argument++]; // out common name

    // int fixation_flag=atoi(argv[argument++]);// whether output fixation or not, fixation_flag =1: output fixation; fixation_flag = 0: do not output

    double scale =  atof(argv[argument++]); // scale value

    // const double littler = rho/double(4*N); // haploid does not have recombination, so this is a dummy variable


    //Write the command line to stderr
    std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
    std::cerr << '\n';

    // calculate robust_degrees, based on robust_flag
    // uint robust_degrees;

    //Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);



    twolocus_pop_no_h_t pop= initialize_pop_w_rob_no_h<twolocus_pop_no_h_t>(N, r.get(), nz, Lv,    c, pv1, ga, type, theta,
                                                                            scale,  &eval_fitness);


    // write the initial state of the population to file
    write_saved_state_two_loc(out_common_name, 0, pop, 0, mu2, rho, r.get(), 0, 0,
                                    std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                    std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                    std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12, std::placeholders::_13,
                                    std::placeholders::_14, std::placeholders::_15, std::placeholders::_16  ) );



    return 0;

}
