#ifndef _SAMPLE_HAPLOID_
#define _SAMPLE_HAPLOID_


//this is a template for haploid simulation!
// adapted from sample_diploid.tcc


#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/diploid_fitness_dispatch.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include "ct_gamete.hpp"
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
#include "ct_haploid_fitness_dispatch.hpp"

#include "math.h"
#include "ct_debug.hpp"
#include "ct_mutation.tcc"

//#include "new_haploids.hpp" // not used now (this new haploids have h saved in haploids)
#include "ct_sample_haploid_helpers.hpp"
#include <map>

//#include "comp_two_method.hpp"
#include "ct_multiloc.hpp"


namespace CT
{

      // haploid model, basic simulations for gamma
      // single deme, N changing
      // for basic gemete type
      // modified on 1.12.17
      // modified on 2.15.17
      // for wv, add scaling factor (set as 1)
      template< typename gamete_type,
      typename gamete_cont_type_allocator,
      typename mutation_type,
      typename mutation_cont_type_allocator,
      typename haploid_geno_t,
      typename haploid_vector_type_allocator,
      typename haploid_fitness_function,
      typename mutation_model,
      // typename recombination_policy,
      template<typename,typename> class gamete_cont_type,
      template<typename,typename> class mutation_cont_type,
      template<typename,typename> class haploid_vector_type,
      typename Mat,
      typename Vec,
      typename mutation_removal_policy = std::true_type,
      typename gamete_insertion_policy = emplace_back>
      double
      sample_haploid(const gsl_rng * r,
      gamete_cont_type<gamete_type,gamete_cont_type_allocator> & gametes,
      haploid_vector_type<haploid_geno_t,haploid_vector_type_allocator> & haploids,
      mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
      std::vector<uint_t> & mcounts,
      const uint_t & N_curr,
      const uint_t & N_next,
      const double & mu,
      const mutation_model & mmodel,
      // const recombination_policy & rec_pol,
      const haploid_fitness_function & ff, // fitness function not need to be const!
      // typename gamete_type::mutation_container & neutral,
      // typename gamete_type::mutation_container & selected,
      const Mat & ref_w,
      const Vec & ref_v,
      const mutation_removal_policy mp = mutation_removal_policy(),
      const gamete_insertion_policy & gpolicy_mut = gamete_insertion_policy())
      {
       /*
         The mutation and gamete containers contain both extinct and extant objects.
         The former are useful b/c the represent already-allocated memory.  The library
         uses these extinct objects to 'recycle' them into new objects.  The function calls
         below create FIFO queues of where extinct objects are.  These queues are passed to
         mutation and recombination functions and used to decide if recyling is possible or
         if a new object needs to be 'emplace-back'-ed into a container.

         The type of the FIFO queue is abstracted with the name KTfwd::fwdpp_internal::recycling_bin_t,
         which is a C++11 template alias.

         The details of recycling are implemented in fwdpp/internal/recycling.hpp
       */
       // make the zero count mutations and gametes into the recyling bin
       auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
       auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);

       // mutation first
       // mutation on each haploid!

       for (uint_t i=0; i<2*N_curr;i++)
       {
         haploids[i]= mutate_gamete_recycle (mut_recycling_bin, gam_recycling_bin, r, mu, gametes, mutations, haploids[i], ref_w, ref_v, mmodel,
                                            std::bind(CT_internal::add_N_mutations_recycle(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5, 1, std::placeholders::_6, std::placeholders::_7, 0),
                                            gpolicy_mut);

         // then update fitness! (now)
          ff(gametes[haploids[i]]); // modified on 1.12.17
       }
       // here will add assertion to make sure all the gamete flags are all 0.
       assert(check_gamete_flags(gametes)); // comment out for debug

       // then recalculate fitness
       std::vector<double> fitnesses(haploids.size());

       double wsum = 0. , wbar; //pop'n mean fitness
       for (uint_t i = 0 ; i < 2*N_curr ; ++i)
       {
         fitnesses[i] = gametes[haploids[i]].fitness;
         wsum += fitnesses[i];
       }


       wbar= wsum/ (2*N_curr); // mean fitness


       // then sample gamete of the next gen from their fitness
       std::vector<int> count(2*N_curr,0);// count of gamete in the next gen
       std::vector<int> lost_haploids; //lost haploid indecies
       std::vector<int> gained_haploids; //gained haploid indecies

       //set up a random generator pointer (for fixed size)
       fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(2*N_curr,fitnesses.data()));

       // sample haploid index for the next generation based on individual fitness
       for (uint_t i = 0 ; i < 2*N_curr ; ++i)
       {
         int index= gsl_ran_discrete(r,lookup.get());
         count[index]++;
       }
       // figure out which gamete are lost and which are gained
       for (uint_t i = 0 ; i < 2*N_curr ; ++i)
       {
         if (count[i]==0)
         {
           gametes[haploids[i]].n--;
           lost_haploids.emplace_back(i);
         }
         else if (count[i]>=1)
         {
           gametes[haploids[i]].n+= count[i]-1;
           gained_haploids.insert(gained_haploids.end(),count[i]-1,i);
         }
       }
       // update haploids
       for (uint_t i =0; i< lost_haploids.size(); i++)
       {
         haploids[lost_haploids[i]]= haploids[gained_haploids[i]];
       }

       // need to updat haploid and gamete structure!
       /*
         At the end of the above loop, we have a bunch of new diploids
         that are all recombined and mutated sampling of the parental generation.

         Our problem is that we no longer know how many times each mutation is present, which
         is corrected by the following call.

         Although the implementation of process_gametes is super-trivial, it is actually the
         most computationally-expensive part of a simulation once mutation rates are large.

         Further, the function is hard to optimize. Recall that gametes store mutations in order
         according to position.  Thus, when we go from a gamete to a position in mcounts, we are
         accessing the latter container out of order with respect to location in memory.  process_gametes
         is thus the "scatter" part of a "scatter-gather" idiom.  Modern x86 CPU have little available
         for vectorizing such cases.  I've experimented with CPU intrinsics to attempt memory prefetches,
         but never saw any performance improvement, and the code got complex, and possibly less portable.

         The implementation is in fwdpp/internal/sample_diploid_helpers.hpp
        */
       fwdpp_internal::process_gametes(gametes,mutations,mcounts);// update mcounts with new count of gametes
       // then update reference
       fwdpp_internal::gamete_cleaner(gametes,mutations,mcounts,2*N_next,mp);// here already remove mutation in the gamete mutation list!!
       return wbar;
      }


        // for all two locus cases, encoding robustness as a locus to evolve along with the phenotype
        // in this case h is not stored in ct_gamete1 (it only affects fitness calculation)
        // with recombination between the two locus
        template<typename ct_gamete2_t,
                 typename haploid_t,
                 typename gamete1_cont_type_allocator,
                 typename gamete2_cont_type_allocator,
                 typename mutation_type,
                 typename mutation_cont_type_allocator,
                 typename haploid_vector_type_allocator,
                 typename haploid_fitness_function,
                 typename mutation_model1,
                 typename mutation_model2,
                 template<typename,typename> class gamete_cont_type,
                 template<typename,typename> class mutation_cont_type,
                 template<typename,typename> class haploid_vector_type,
                 typename mutation_removal_policy = std::true_type,
                 typename gamete_insertion_policy = emplace_back
           >
        double
        sample_haploid(const gsl_rng * r,
          gamete_cont_type<ct_gamete_geno_no_h_t,gamete1_cont_type_allocator> & gametes1, // gamete container for the locus robust bits
          gamete_cont_type<ct_gamete2_t,gamete2_cont_type_allocator> & gametes2, // gamete container for genotype
          haploid_vector_type<haploid_t,haploid_vector_type_allocator> & haploids,
          mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
          std::vector<uint_t> & mcounts,
          const uint_t & N_curr,
          const uint_t & N_next,
          const double & mu1,
          const double & mu2,
          const mutation_model1 & mmodel1,
          const mutation_model2 & mmodel2,
          const haploid_fitness_function & ff, // fitness function not need to be const!
          // typename ct_gamete_geno_no_h_t::mutation_container & neutral,
          // typename ct_gamete_geno_no_h_t::mutation_container & selected,
          const MatrixXd & ref_w,
          const VectorXi & ref_v1,
          const VectorXi & ref_v2,
          double shuffle_ratio, // proportion of shuffling the allele (an aproximate to recombination)
          const int & robust_bits,
          const int & flag, //robustness flag
          const int & robust_degrees,
          const mutation_removal_policy & mp = mutation_removal_policy(),
          const gamete_insertion_policy & gpolicy_mut = gamete_insertion_policy())
        {

          auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
          auto gam1_recycling_bin = fwdpp_internal::make_gamete_queue(gametes1);
          auto gam2_recycling_bin = fwdpp_internal::make_gamete_queue(gametes2);

          // mutation first
          // mutation on each haploid!


          // genotype mutation step (happening in the gamete )
          // haploids will be defined as pair of index
          std::vector<int> gametes1_idx;
         //  std::vector<int> gametes2;
          int index1;

          // here only mutate loc 1
          for (uint_t i=0; i<2*N_curr;i++)
          {
             // index to the gamete of locus 1
             index1 = mutate_gamete_recycle (mut_recycling_bin,gam1_recycling_bin, r,mu1,gametes1,mutations,haploids[i].first, ref_v1,mmodel1,robust_bits,flag,robust_degrees,
                                             std::bind(CT_internal::add_N_mutations_recycle(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
                                             std::placeholders::_5,std::placeholders::_6,std::placeholders::_7,std::placeholders::_8,std::placeholders::_9),
                                             gpolicy_mut);



             gametes1_idx.push_back(index1);
          }

        // random segregation step
        // only need to shuffle the order of one of the vector (here only shuffle gamete1)
        // for test, comment out this line
        // only shuffle the first k gamete order.
        // if shuffle_ratio =0, do not shuffle; if ratio >= 0.5, random shuffle
        // if 0< shuffle_ratio < 0.5 , based on the probability, choose the number of random pairs and exchange index
        if (shuffle_ratio >= 0.5 )
        {
          gsl_ran_shuffle(r,gametes1_idx.data(),2*N_curr,sizeof(int));
        }
        else if (shuffle_ratio>0 )
        {
          int n_exchange = gsl_ran_poisson(r,shuffle_ratio*N_curr);
          int temp_idx[2];
          int temp_val ;
          std::vector<int> index_of_gametes1_idx(2*N_curr);
          for (int i=0; i<2*N_curr ; i++)
          {
            index_of_gametes1_idx[i]=i;
          }
          // then random choose by index
          for (int i=0; i< n_exchange; i++)
          {
             gsl_ran_choose(r,temp_idx, 2, index_of_gametes1_idx.data(),2*N_curr, sizeof(int) );
             // only exchange when the two index values are different
             if(gametes1_idx[temp_idx[0]] != gametes1_idx[temp_idx[1]])
             {
               temp_val = gametes1_idx[temp_idx[0]];
               gametes1_idx[temp_idx[0]] = gametes1_idx[temp_idx[1]];
               gametes1_idx[temp_idx[1]] = temp_val;
             }
          }
        }

        // then recalculate fitness
        std::vector<double> fitnesses(haploids.size());

        // create a std::map to store pair of index
        typedef std::pair<int, int> Key;
        std::map< Key , double> fitness_map;
        std::map< Key , double>::iterator it;


        // two things being done in the fitness function
        // 1) resample h for haploids
        // 2) recalculate fitness based on wv in g, and h in haploid
        double wsum = 0. , wbar; //pop'n mean fitness

        double ratio ; // the scaling factor of the robust genotype

        // then calcualte fitness
        for (uint_t i=0; i<2*N_curr;i++)
        {
          haploids[i].first=gametes1_idx[i];

          // use the robustness in robust genotype
          ratio = gametes1[haploids[i].first].get_robust_alpha_ratio();

          // index to the gamete of locus 2
          haploids[i].second = mutate_gamete_recycle (mut_recycling_bin,gam2_recycling_bin, r,mu2,gametes2,mutations,haploids[i].second,ref_w, ref_v2,mmodel2,
                                                      std::bind(CT_internal::add_N_mutations_recycle(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5, ratio, std::placeholders::_6, std::placeholders::_7, 1),
                                                      gpolicy_mut);


          it = fitness_map.find(haploids[i]);
          if (it == fitness_map.end()) // if it does not in the fitness_map, recalculate fitness and add to it
          {
            fitnesses[i]= ff(gametes2[haploids[i].second].wv);
            fitness_map[haploids[i]]= fitnesses[i];
          }
          else
          {
            fitnesses[i]= it->second; // if exist, just take fitness value stored in the fitness_map
          }

          wsum += fitnesses[i];
        }

          wbar= wsum/ (2*N_curr); // mean fitness


          // then sample gamete of the next gen from their fitness
          std::vector<int> count(2*N_curr,0);// count of gamete in the next gen
          std::vector<int> lost_haploids; //lost haploid indecies
          std::vector<int> gained_haploids; //gained haploid indecies

          //set up a random generator pointer (for fixed size)
          fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(2*N_curr,fitnesses.data()));

          // sample haploid index for the next generation based on individual fitness
          for (uint_t i = 0 ; i < 2*N_curr ; ++i)
          {
            int index= gsl_ran_discrete(r,lookup.get());
            count[index]++;
          }
          // figure out which gamete are lost and which are gained
          for (uint_t i = 0 ; i < 2*N_curr ; ++i)
          {
            if (count[i]==0)
            {
              gametes1[haploids[i].first].n--;
              gametes2[haploids[i].second].n--;
              lost_haploids.emplace_back(i);
            }
            else if (count[i]>=1)
            {
              gametes1[haploids[i].first].n+= count[i]-1;
              gametes2[haploids[i].second].n+= count[i]-1;
              gained_haploids.insert(gained_haploids.end(),count[i]-1,i);
            }
          }
          // update haploids
          for (uint_t i =0; i< lost_haploids.size(); i++)
          {
            haploids[lost_haploids[i]]= haploids[gained_haploids[i]];
          }


          // need to updat haploid and gamete structure!
          /*
            At the end of the above loop, we have a bunch of new diploids
            that are all recombined and mutated sampling of the parental generation.

            Our problem is that we no longer know how many times each mutation is present, which
            is corrected by the following call.

            Although the implementation of process_gametes is super-trivial, it is actually the
            most computationally-expensive part of a simulation once mutation rates are large.

            Further, the function is hard to optimize. Recall that gametes store mutations in order
            according to position.  Thus, when we go from a gamete to a position in mcounts, we are
            accessing the latter container out of order with respect to location in memory.  process_gametes
            is thus the "scatter" part of a "scatter-gather" idiom.  Modern x86 CPU have little available
            for vectorizing such cases.  I've experimented with CPU intrinsics to attempt memory prefetches,
            but never saw any performance improvement, and the code got complex, and possibly less portable.

            The implementation is in fwdpp/internal/sample_diploid_helpers.hpp
           */
          CT_internal::process_gametes(gametes1,gametes2,mutations,mcounts);// update mcounts with new count of gametes

          // caution!!  the below is tempted modification (not sure if it works fine or not)
          // then update reference
          fwdpp_internal::gamete_cleaner(gametes1,mutations,mcounts,2*N_next,mp);// here already remove mutation in the gamete mutation list!!
          fwdpp_internal::gamete_cleaner(gametes2,mutations,mcounts,2*N_next,mp);
          return wbar;

        }

}
#endif
