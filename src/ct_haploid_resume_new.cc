// created on 1.12.17
// new main function to resume running from saved state using
// new interface functions!

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
#include "read_saved_state.hpp"

#include <iomanip>

// using mtype = KTfwd::mutation_base;
using mtype = CT::ct_mutation ;
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
using uint_t = KTfwd::uint_t ;
using singlepop_basic_pop_t = CT::singlepop_basic_t<mtype>; // for basic type


int main(int argc, char ** argv)
{
  // add 6 more arguments specific to my problem
  if (argc != 4)
    {
      std::cerr << "Too few or too more arguments\n"
		<< "Usage: diploid_ind saved_state_common_name prev_ngen  ngens \n";  //add fixation flag (in order to output fixation as well )
      exit(10);
    }

    int argument=1;
    std::string out_common_name= argv[argument++]; // out file common name
    std::string prev_ngen_s= argv[argument++]; // the stopping generation for the previous run
    const unsigned ngens= unsigned(atoi(argv[argument++])); // number of generations to run

    // int fixation_flag=atoi(argv[argument++]);// whether output fixation or not, fixation_flag =1: output fixation; fixation_flag = 0: do not output


    // gsl_rng * rtest = (gsl_rng *)malloc(fsize);
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    double mu;
    std::size_t prev_ngen;


    singlepop_basic_pop_t pop= read_saved_state_one_linked_loc<singlepop_basic_pop_t>(out_common_name,prev_ngen_s, r, mu, prev_ngen, 1,
                            std::bind(CT::read_binary_pop_hap(),std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                            std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                            std::placeholders::_9,std::placeholders::_10,  std::placeholders::_11, std::placeholders::_12, std::placeholders::_13));



    // write file handle for mean fitness and mean robustness
    std::string out_mean_fit_file = out_common_name +"_mean_fit.txt";

    std::ofstream out_mean_fit(out_mean_fit_file,std::ofstream::out | std::ofstream::app);

    unsigned N = pop.haploids.size()/2;

    std::size_t Lv = pop.ct.v.size();

    double wbar;
    unsigned twoN = 2*N;
    double littler =0; // set recombination to be 0

    // simulations starts here!
    for( uint generation = prev_ngen; generation < (prev_ngen+ ngens); ++generation  )
    {
      //Iterate the population through 1 generation
      wbar = CT::sample_haploid(r,
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
                             [&r,Lv](){return ((double) gsl_rng_uniform_int(r,Lv)) / Lv ;}),
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

      write_saved_state_one_linked_loc(out_common_name, prev_ngen+ ngens, pop, mu, littler, r,
                                      std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                      std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                      std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12  ) );

      return 0;
}
