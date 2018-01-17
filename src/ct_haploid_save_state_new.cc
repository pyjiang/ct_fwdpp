
#include <vector>
#include <iterator>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/debug.hpp>
// #include <Sequence/SimData.hpp>


#include "sample_haploid.tcc"
#include "ct_mut_model.hpp"
#include "ct_util.hpp"

#include "ct_fitness_policy.tcc"
#include "ct_mutation_internal.hpp"
#include "ct_mutation.tcc"

#include "ct_eigen.hpp"

#include "initialize_pop.tcc"
#include "write_saved_state.hpp"

#include "ct_mutation.hpp"

#include <iomanip>


using mtype = CT::ct_mutation ;

using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
using singlepop_basic_pop_t = CT::singlepop_basic_t<mtype>; // for basic type



int main(int argc, char ** argv)
{

  if (argc != 13)
    {
      // always output fixations
      std::cerr << "Too few or too more arguments\n"
    << "Usage: diploid_ind N theta  ngens  seed nz Lv  c pv1 ga_flag ga type out_common_name \n"; //
      exit(10);
    }

    // read in arguments
    int argument=1;
    const unsigned N = unsigned(atoi(argv[argument++]));           //Number of diploids
    const double theta = atof(argv[argument++]);         //4*n*mutation rate.  Note: mutation rate is per REGION, not SITE!!

    const unsigned ngens = unsigned(atoi(argv[argument++]));       //Number of generations to simulate

    const unsigned seed = unsigned(atoi(argv[argument++]));        //Random number seed

    const double mu = theta/double(4*N);                 //per-gamete mutation rate

    uint32_t nz= unsigned(atoi(argv[argument++])); // number of cell types
    uint32_t Lv= unsigned(atoi(argv[argument++]));


    double c= atof(argv[argument++]); // ratio of non-zeros in w matrix
    double pv1= atof(argv[argument++]); // ratio of 1 in v
    int ga_flag = atoi(argv[argument++]); // ga_flag ==0: ga is passed as a parameter; ga_flag==1: ga is passed as index to the gamma value in the ga_list

    // stored gamma value, which are equally separated in robustness by gamma
    // 10 digit accuracy (in string)
    double  ga_list[10]={0.00100000, 0.11247546, 0.22826615, 0.35353067, 0.49564361, 0.66695182, 0.89196037, 1.23201245, 1.92302210, 15.7945323 };
    uint ga_index;
    double ga;
    // if flag == 0: gamma is passed as a parameter
    assert (ga_flag ==0 || ga_flag ==1);

    if (ga_flag==0)
    {
      ga =  atof(argv[argument++]);
    }
    else  // if flag==1: gamma is passed as index to the saved gamma values
    {
      ga_index = atoi(argv[argument++]);
      if( ga_index >= 10)
      {
        std::cerr << "gamma index not in range!" << std::endl;
        exit(10);
      }
      else
      ga = ga_list[ga_index];
    }

    int type= atoi(argv[argument++]); // which type of simulation (randint or stabil), type=0: stabil, type=1: randinit
    std::string out_common_name= argv[argument++]; // out common name

    const double littler = 0; // this haploid model does not have recombination, so this is a dummy variable


    //Write the command line to stderr
    std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
    std::cerr << '\n';


    //Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    singlepop_basic_pop_t pop= initialize_pop<singlepop_basic_pop_t>(N, r.get(), nz, Lv, c, pv1, ga, type, theta,
                                                      &eval_fitness);

    // write the initial state of the population to file
    write_saved_state_one_linked_loc(out_common_name, 0, pop, mu, littler, r.get(),
                                    std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                    std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                    std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12  ) );

    // write file handle for mean fitness and mean robustness
    std::string out_mean_fit_file = out_common_name +"_mean_fit.txt";

    std::ofstream out_mean_fit(out_mean_fit_file,std::ofstream::out | std::ofstream::app);

    unsigned generation;
    double wbar;


    unsigned twoN = 2*N;

    // simulations starts here!
    for( generation = 0; generation < ngens; ++generation )
    {
        //Iterate the population through 1 generation
        wbar = CT::sample_haploid(r.get(),
                                  pop.gametes,  //non-const reference to gametes
                                  pop.haploids, //non-const reference to haploids
                                  pop.mutations, //non-const reference to mutations
                                  pop.mcounts, // count of each mutation in the mutation list
                                  N,     //current pop size, remains constant
                                  N,
                                  mu,    //mutation rate per gamete
                                   /*
                                   The mutation model (KTfwd::infsites) will be applied by
                             sample_diploid in order to add mutations to gametes each generation.
                                       */
                                 std::bind(CT::ct_mut_model(),std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup),
                               [&r,Lv](){return ((double) gsl_rng_uniform_int(r.get(),Lv)) / Lv ;}),
                                 //The function to generation recombination positions:
                                //  [](){return 0.;}, // no recombination so far
                                 std::bind(CT::ct_fit_simple(), std::placeholders::_1, pop.ct.h , pop.b_vec),
                                //  pop.neutral,
                                //  pop.selected,
                                 pop.ct.w,
                                 pop.ct.v);


        CT::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation, twoN, pop.ct.v); //for updating fixed or lost mutations

        assert(KTfwd::check_sum(pop.gametes,twoN));

        out_mean_fit << generation << "," << wbar << std::endl;


      }
      // close mean fit file
      out_mean_fit.close();

      write_fix(out_common_name,pop);

      write_saved_state_one_linked_loc(out_common_name, ngens, pop, mu, littler, r.get(),
                                      std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                      std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                      std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12  ) );

      return 0;
}
