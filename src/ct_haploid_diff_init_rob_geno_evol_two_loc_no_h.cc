// this initialize pop with diff robustness genotype


// when scale_flag=8, the first 3 bits in genotype represents robust bits
// 0: represents the lowest digit and vice versa. For example, 4 (100) will have mutation in bit 2.
// modified on 8.8.17
// added scale_flag =8, 8 degrees of robustness
// updated on 2.2.17
// remove calc_scale_func in the argument  (scale_flag is sufficient)
// modified on 2.1.17
// different initial robustness, and also have the ability to evolve robustness
// modified initial robustness set up; the initial geno will reflect the robust bits (when evolve)
// so now the initial gamma (alpha) is specified by parameter gamma, not by default calculation from robust bits (because that was a different funtion for initlize and further evolution)


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

// #include "scale_vec.hpp"
#include "ct_mutation.hpp"
#include <iomanip>

using mtype = CT::ct_mutation ;
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;


using twolocus_pop_no_h_t = CT::twolocus_no_h_t<mtype>;

extern double LOW_SCALE_BOUND;
extern double HIGH_SCALE_BOUND;


int main(int argc, char ** argv)
{

  if (argc != 19)
    {
      // updated on 3.6.17
      // remove robust_val in the argument
      // change parameter from rvz to Lv (total length of v); and also add a parameter of robust_bits
      // do not use gamma to initialze
      // take parameters to set low and high scale boundary
      // calc_scale_flag ==0: additive calc_scale_func; calc_scale_func ==1: threshold (only for 10 bits)
      // scale_flag inidicates how many different initial differet robustness , scale_flag =3, 5, 10.
      // robust_flag indicates whether the robust geno is considered using all 2^n (=0) or considered additively (=1)
      std::cerr << "Too few or too more arguments\n"
    << "Usage: diploid_ind N theta rho ngens  seed mu1_base nz Lv  robust_bits  c pv1 ga type out_common_name   robust_flag  low_bound high_bound  scale_flag\n"; // type=0: stabil, type=1: randinit
      exit(10);
    }

    // read in arguments
    int argument=1;
    const unsigned N = unsigned(atoi(argv[argument++]));           //Number of diploids
    const double theta = atof(argv[argument++]);         //4*n*mutation rate.  Note: mutation rate is per REGION, not SITE!!
    const double rho = atof(argv[argument++]);           // use this parameter to pass in proportion of shuffled gametes
    const unsigned ngens = unsigned(atoi(argv[argument++]));       //Number of generations to simulate
    // const unsigned samplesize1 = unsigned(atoi(argv[argument++])); //Sample size to draw from the population
    // int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
    const unsigned seed = unsigned(atoi(argv[argument++]));        //Random number seed

    // const double mu = theta/double(4*N);                 //per-gamete mutation rate

    // mutation rate of robust bits
    double mu1_base = atof(argv[argument++]);

    uint32_t nz= unsigned(atoi(argv[argument++])); // number of cell types
    // double rvz=  atof(argv[argument++]);    //ratio of length of v to length of z
    uint32_t Lv= unsigned(atoi(argv[argument++]));

    uint robust_bits = unsigned(atoi(argv[argument++]));


    double c= atof(argv[argument++]); // ratio of non-zeros in w matrix
    double pv1= atof(argv[argument++]); // ratio of 1s in v

    double  ga =  atof(argv[argument++]);

    const double mu1= mu1_base*robust_bits; // mutation rate per locus for robust bits
    const double mu2 = theta/double(4*N); // mutation rate per locus for regulation v


    int type= atoi(argv[argument++]); // which type of simulation (randint or stabil)
    std::string out_common_name= argv[argument++]; // out common name

    // int fixation_flag=atoi(argv[argument++]);// whether output fixation or not, fixation_flag =1: output fixation; fixation_flag = 0: do not output

    int robust_flag = atoi(argv[argument++]);


    // low and high scale boundary
    LOW_SCALE_BOUND = atof(argv[argument++]);
    HIGH_SCALE_BOUND = atof(argv[argument++]);

    // add flag to specify which function to use to calculate scale for robustness
    // int scale_func_flag = atoi(argv[argument++]);

    int scale_flag = atoi(argv[argument++]);

    //Write the command line to stderr
    std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
    std::cerr << '\n';

    // calculate robust_degrees, based on robust_flag
    uint robust_degrees;

    //Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    // function pointer to the current calculate scale function
    double (*update_scale_func)(uint32_t,uint32_t);

    // make sure shuffle_ratio is [0,0.5]
    assert(rho >=0 || rho <=0.5 ) ;

    // added flag=10 on 2.21.17 (for 2 based, 10 digits)
    // add flag = 7 , for 8 degree, additive case
    // flag=15 is the additive case for flag=16
    assert(scale_flag ==0  || scale_flag ==8 || scale_flag ==7 || scale_flag ==16 || scale_flag ==15);


    std::vector<int> robust_geno_vec;

    // added on 2.20.17
    // add an option to use just 2 bits, 3 degrees to encode robustness
    // note: robust flag can only take 0,7,8,15,16 now
    if (scale_flag ==0)
    {
       robust_degrees =3;
       update_scale_func = &calc_scale_2base; // modified on 2.26.17, 2-based three degrees of robustness
       robust_geno_vec = {0,1,2};
    }
    else if(scale_flag ==8 || scale_flag ==7) // scale_flag=7 is the additive case, using 7 bits, the robust_geno represents how many 1s in the genotype ; while sclae_flag = 8 is using 3 bit to represent 8 state
    {
        robust_degrees = 8;
        update_scale_func = &calc_scale_8_degree;
        robust_geno_vec = {0,1,2,3,4,5,6,7};
    }
    else // scale_flag == 15 or ==16  new, added 16 degrees of alpha, spanning robustness
    {
      robust_degrees = 16;
      update_scale_func = &calc_scale_16_degree;
      robust_geno_vec = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    }

    // scale of robustness
    // the total range is from 0.5 to 2.
    // std::vector<double> scale_vec = {0.5000736 ,  0.58100263,  0.67502875,  0.78427151,  0.9111935 ,
    //                               1.05865581,  1.22998257,  1.42903587,  1.66030281,  1.92899666};


    // // initialize static function pointer by the scale_func_flag
    // initialize_calc_function_ptr( scale_func_flag);

    // std::vector<double> scale_vec = {0.5, 0.6, 0.7, 1, 1.2, 1.5, 1.7, 2};

    // single_pop_wij_b_t pop= initialize_pop_diff_rob<single_pop_wij_b_t>(N, r.get(), nz, Lv, c, pv1, ga, type, theta, 2, 500, 1000, scale_vec,
    //                                                   &eval_fitness);

    // intialize scale from the degrees
    // output initial scale as the first line !

    // initial robust_geno vector
    // std::vector<int> robust_geno_vec={1,5,9};

    twolocus_pop_no_h_t pop= initialize_pop_diff_rob<twolocus_pop_no_h_t>(N, r.get(), nz, Lv +robust_bits, c, pv1, ga, type, theta, robust_degrees, robust_flag, robust_bits, robust_bits,
                                                                          robust_geno_vec, scale_flag);

    // write the initial state of the population to file
    write_saved_state_two_loc(out_common_name, 0, pop, mu1, mu2, rho, r.get(), robust_bits, robust_flag,
                                    std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                    std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                    std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12, std::placeholders::_13,
                                    std::placeholders::_14, std::placeholders::_15, std::placeholders::_16  ) );



    // write file handle for mean fitness and mean robustness
    std::string out_mean_fit_file = out_common_name +"_mean_fit.txt";
    std::string out_rob_dist_file = out_common_name +"_robust_dist.txt";

    std::ofstream out_mean_fit(out_mean_fit_file,std::ofstream::out | std::ofstream::app);
    std::ofstream out_rob_dist(out_rob_dist_file,std::ofstream::out | std::ofstream::app);

    unsigned generation;
    double wbar;
    // double mean_robust;

    unsigned twoN = 2*N;

    // store alpha vector
    // std::vector<double> alpha_vector = {pop.gametes[0].alpha, pop.gametes[1].alpha, pop.gametes[2].alpha };

    // output alpha list as the first line!
    // out_rob_dist << alpha_vector[0] << "\t" << alpha_vector[1] << "\t" << alpha_vector[2] << "\t" << std::endl;

    // initial scale vector
    std::vector<double> scale_vec;
    for (auto i: robust_geno_vec)
    {
      double scale = update_scale_func(robust_degrees,i);
      scale_vec.push_back(scale);
      out_rob_dist << scale << "\t";
    }
    out_rob_dist << std:: endl;

    // for (int i=0; i< scale_vec.size() ;i++)
    // {
    //   out_rob_dist << scale_vec[i] << "\t";
    // }
    // out_rob_dist << std:: endl;


    VectorXi v1= pop.ct.v.head(robust_bits);
    VectorXi v2= pop.ct.v.tail(Lv);


    // simulations starts here!
    for( generation = 0; generation < ngens; ++generation )
    {

        wbar= CT::sample_haploid(r.get(),
                          				pop.gametes1,
                                  pop.gametes,
                          				pop.haploids,
                          				pop.mutations,
                          				pop.mcounts,
                          				N,
                                  N,
                          				mu1,
                                  mu2,
                                  std::bind(CT::ct_mut_model(),std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup),
                                [&r,robust_bits](){return ((double) gsl_rng_uniform_int(r.get(),robust_bits)) / robust_bits ;}) ,
                                  std::bind(CT::ct_mut_model(),std::placeholders::_1,std::placeholders::_2,std::ref(pop.mut_lookup),
                                [&r,Lv](){return ((double) gsl_rng_uniform_int(r.get(),Lv)) / Lv +1  ;}),
                          			  std::bind(CT::ct_fit_simple(),std::placeholders::_1, pop.ct.h, pop.b_vec, 1),
                          				// pop.neutral,
                          				// pop.selected,
                                  pop.ct.w,
                                  v1,
                                  v2,
                                  rho,
                                  robust_bits,
                                  robust_flag,
                                  robust_degrees);

        CT::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,twoN,v1,v2);

        assert(KTfwd::check_sum(pop.gametes,twoN));

        out_mean_fit << generation << "," << wbar << std::endl;

        calc_and_output_scale_dist_two_loc(out_rob_dist, pop.gametes1, pop.haploids, N, scale_vec );
      }
      // close mean fit file
      out_mean_fit.close();
      out_rob_dist.close();

      // now needs to copy v1 and v2 back to ct.v (so that it will be easier to write to file)
      copy_back( pop.ct.v,v1,v2);

      write_fix(out_common_name,pop);

      write_saved_state_two_loc(out_common_name, ngens, pop, mu1, mu2, rho, r.get(), robust_bits, robust_flag,
                                      std::bind(CT::write_binary_pop_hap(), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                                      std::placeholders::_5,std::placeholders::_6, std::placeholders::_7, std::placeholders::_8,
                                      std::placeholders::_9,std::placeholders::_10, std::placeholders::_11, std::placeholders::_12, std::placeholders::_13,
                                      std::placeholders::_14, std::placeholders::_15, std::placeholders::_16  ) );

      return 0;
}
