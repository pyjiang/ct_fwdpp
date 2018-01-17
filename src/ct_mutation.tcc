#ifndef _CT_MUTATION_TCC_
#define _CT_MUTATION_TCC_


//modified from mutations.tcc

#include <iostream>
#include <type_traits>
#include <numeric>
#include <gsl/gsl_randist.h>
// #include <fwdpp/internal/mutation_internal.hpp>
#include "ct_mutation_internal.hpp"
// #include "ct_gamete.hpp"

//for debug LD (if not include this, the compilation wouldn't pass)
// #include "ct_eigen.tcc"

#include "ct_eigen.hpp"


// PJ modified on 11.1.16
// will only allocate a new space when gamete count is larger than 1!
// add another parameter to this function, which takes in compiler time true/fase
namespace CT
{
  // define ct_gamete_with_h_t
  // using ct_gamete_with_h_t= ct_gamete_with_h<void,VectorXd>;

  // modified on 12.17.16
  // write some generalized version of mutate_gamete_recycle
  // so that it can take different functions of add_N_mutations_recycle, given different parameters outside of this function
  // add one more parameter of add_N_mutation_func.

  // updated on 12.7.16
  // do not think about back-compatibility for now (because add_N_mutations_recycle has been modified already)
  // now this version works with ct_gamete_robust_geno type (but without robust bits)

  // modified on 12.21.16
  // write two overloading function one with w, one without

  template< typename queue_type,
	    typename queue_type2,
	    typename mutation_model,
	    typename gamete_insertion_policy,
	    typename gcont_t,
	    typename mcont_t,
      typename Mat,
      typename Vec,
      typename add_N_mutation_func>
  std::size_t mutate_gamete_recycle( queue_type & recycling_bin,
				     queue_type2 & gamete_recycling_bin,
				     const gsl_rng * r,
				     const double & mu,
				     gcont_t & gametes,
				     mcont_t & mutations,
				     const std::size_t g, // the haploid referece (to gamete), if changed, will return the new referece!
             const Mat & ref_w,
             const Vec & ref_v,
				     const mutation_model & mmodel,
             const add_N_mutation_func & add_N_mutations_f,
				     const gamete_insertion_policy & gpolicy) // this is compiling time true/false value (false means gametes does not have h)
  {
    // static_assert( traits::valid_mutation_model<mutation_model,mcont_t,gcont_t>::value,
		//    "error: type mutation_model is not a dispatchable mutation model type!" );
    assert(g<gametes.size());
    assert( gamete_is_sorted_n(gametes[g],mutations) );
    assert( gamete_is_sorted_s(gametes[g],mutations));
    // if(! gamete_is_sorted_s(gametes[g],mutations) ) //for debug
    // {
    //   std::cout << "error";
    // }

    // assert( gamete_is_sorted_s(gametes[g],mutations) );
    unsigned nm = gsl_ran_poisson(r,mu); // number of mutations on this gamete
    if(!nm) return g;

    assert(gametes[g].n);

    if (gametes[g].n==1) // if it is the only count, just change itself!!
    {
      gametes[g].flag=1;
      // add N mutations to the current gamete (use from recycle is possible)
      add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[g],ref_w,ref_v);
      // add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[g],ref_v);
      return g;
    }
    else
    {
      gametes[g].n--;
      //Recycle an already-allocated gamete, if possible
      if (!gamete_recycling_bin.empty())
      {
        	auto idx = gamete_recycling_bin.front();
        	gamete_recycling_bin.pop();
        	assert(idx!=g);
        	assert(!gametes[idx].n);
        	// gametes[idx].mutations=gametes[g].mutations;
        	// gametes[idx].smutations=gametes[g].smutations;
          gametes[idx].copy(gametes[g]);

        	gametes[idx].n=1;

          // PJ added, set flag to 1 (that needs recalculate fitness!)
          gametes[idx].flag=1;

          // add N mutations to the current gamete (use from recycle is possible)
        	add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[idx],ref_w,ref_v);
          // add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[idx],ref_v);
        	return idx;
      }
      typename gcont_t::value_type ng( 1, gametes[g]);
      add_N_mutations_f(recycling_bin,mmodel,nm,mutations,ng,ref_w,ref_v);
      // add_N_mutations_f(recycling_bin,mmodel,nm,mutations,ng,ref_v);
      return gpolicy(std::move(ng),gametes);
    }
  }
  //
  // // this is to compare result of two gametes calcualtion
  // template< typename queue_type,
	//     typename queue_type2,
	//     typename mutation_model,
	//     typename gamete_insertion_policy,
	//     typename gcont1_t,
  //     typename gcont2_t,
	//     typename mcont_t,
  //     typename Mat,
  //     typename Vec,
  //     typename add_N_mutation_func>
  // std::size_t mutate_gamete_recycle( queue_type & recycling_bin,
	// 			     queue_type2 & gamete_recycling_bin,
	// 			     const gsl_rng * r,
	// 			     const double & mu,
	// 			     gcont1_t & gametes1,
  //            gcont2_t & gametes2,
	// 			     mcont_t & mutations,
	// 			     const std::size_t g, // the haploid referece (to gamete), if changed, will return the new referece!
  //            const Mat & ref_w,
  //            const Vec & ref_v,
	// 			     const mutation_model & mmodel,
  //            const add_N_mutation_func & add_N_mutations_f,
	// 			     const gamete_insertion_policy & gpolicy) // this is compiling time true/false value (false means gametes does not have h)
  // {
  //   // static_assert( traits::valid_mutation_model<mutation_model,mcont_t,gcont_t>::value,
	// 	//    "error: type mutation_model is not a dispatchable mutation model type!" );
  //   assert( g<gametes1.size());
  //   assert( gamete_is_sorted_n(gametes1[g],mutations) );
  //   assert( gamete_is_sorted_s(gametes1[g],mutations));
  //   // if(! gamete_is_sorted_s(gametes[g],mutations) ) //for debug
  //   // {
  //   //   std::cout << "error";
  //   // }
  //
  //   // assert( gamete_is_sorted_s(gametes[g],mutations) );
  //   unsigned nm = gsl_ran_poisson(r,mu); // number of mutations on this gamete
  //   if(!nm) return g;
  //
  //   assert(gametes1[g].n);
  //
  //   // make sure gamete1 and gamete2 has the same n!
  //   assert (gametes1[g].n == gametes2[g].n);
  //
  //
  //   if (gametes1[g].n==1) // if it is the only count, just change itself!!
  //   {
  //     gametes1[g].flag=1;
  //     gametes2[g].flag=1;
  //     // add N mutations to the current gamete (use from recycle is possible)
  //     add_N_mutations_f(recycling_bin,mmodel,nm,mutations, gametes1[g], gametes2[g], ref_w,ref_v);
  //     // add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[g],ref_v);
  //     return g;
  //   }
  //   else
  //   {
  //     gametes1[g].n--;
  //     gametes2[g].n--;
  //     //Recycle an already-allocated gamete, if possible
  //     if (!gamete_recycling_bin.empty())
  //     {
  //       	auto idx = gamete_recycling_bin.front();
  //       	gamete_recycling_bin.pop();
  //       	assert(idx!=g);
  //       	assert(!gametes1[idx].n);
  //       	// gametes[idx].mutations=gametes[g].mutations;
  //       	// gametes[idx].smutations=gametes[g].smutations;
  //         gametes1[idx].copy(gametes1[g]);
  //         gametes2[idx].copy(gametes2[g]);
  //
  //       	gametes1[idx].n=1;
  //         gametes2[idx].n=1;
  //
  //         // PJ added, set flag to 1 (that needs recalculate fitness!)
  //         gametes1[idx].flag=1;
  //         gametes2[idx].flag=1;
  //
  //         // add N mutations to the current gamete (use from recycle is possible)
  //       	add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes1[idx],gametes2[idx],ref_w,ref_v);
  //         // add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[idx],ref_v);
  //       	return idx;
  //     }
  //     typename gcont1_t::value_type ng1( 1, gametes1[g]);
  //     typename gcont2_t::value_type ng2( 1, gametes2[g]);
  //
  //
  //
  //     add_N_mutations_f(recycling_bin,mmodel,nm,mutations,ng1,ng2,ref_w,ref_v);
  //     // add_N_mutations_f(recycling_bin,mmodel,nm,mutations,ng,ref_v);
  //
  //     gpolicy(std::move(ng1),gametes1);
  //     auto y = gpolicy(std::move(ng2),gametes2);
  //
  //     //
  //     // // for debug only
  //     // if (y==232)
  //     // {
  //     //   std::cout << "here!";
  //     // }
  //
  //     return y;
  //
  //   }
  // }







  // this version does not have Mat
  // this function is specific to v that don't need to recalculate wv, so these v affect robust bits
  template< typename queue_type,
      typename queue_type2,
      typename mutation_model,
      typename gamete_insertion_policy,
      typename gcont_t,
      typename mcont_t,
      typename Vec,
      typename add_N_mutation_func>
  std::size_t mutate_gamete_recycle( queue_type & recycling_bin,
             queue_type2 & gamete_recycling_bin,
             const gsl_rng * r,
             const double & mu,
             gcont_t & gametes,
             mcont_t & mutations,
             const std::size_t g, // the haploid referece (to gamete), if changed, will return the new referece!
             const Vec & ref_v,
             const mutation_model & mmodel,
             const int & robust_bits,
             const int & flag,
             const int & robust_degrees,
             const add_N_mutation_func & add_N_mutations_f,
             const gamete_insertion_policy & gpolicy) // this is compiling time true/false value (false means gametes does not have h)
  {
    // static_assert( traits::valid_mutation_model<mutation_model,mcont_t,gcont_t>::value,
    //    "error: type mutation_model is not a dispatchable mutation model type!" );
    assert(g<gametes.size());
    assert( gamete_is_sorted_n(gametes[g],mutations) );
    assert( gamete_is_sorted_s(gametes[g],mutations));
    // if(! gamete_is_sorted_s(gametes[g],mutations) ) //for debug
    // {
    //   std::cout << "error";
    // }

    // assert( gamete_is_sorted_s(gametes[g],mutations) );
    unsigned nm = gsl_ran_poisson(r,mu); // number of mutations on this gamete
    if(!nm) return g;

    assert(gametes[g].n);

    if (gametes[g].n==1) // if it is the only count, just change itself!!
    {
      gametes[g].flag=1;
      // add N mutations to the current gamete (use from recycle is possible)
      add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[g],ref_v,robust_bits,flag,robust_degrees);
      return g;
    }
    else
    {
      gametes[g].n--;
      //Recycle an already-allocated gamete, if possible
      if (!gamete_recycling_bin.empty())
      {
          auto idx = gamete_recycling_bin.front();
          gamete_recycling_bin.pop();
          assert(idx!=g);
          assert(!gametes[idx].n);
          // gametes[idx].mutations=gametes[g].mutations;
          // gametes[idx].smutations=gametes[g].smutations;
          gametes[idx].copy(gametes[g]);

          gametes[idx].n=1;

          // PJ added, set flag to 1 (that needs recalculate fitness!)
          gametes[idx].flag=1;

          // add N mutations to the current gamete (use from recycle is possible)
          add_N_mutations_f(recycling_bin,mmodel,nm,mutations,gametes[idx],ref_v,robust_bits,flag,robust_degrees);
          return idx;
      }
      typename gcont_t::value_type ng( 1, gametes[g]);
      add_N_mutations_f(recycling_bin,mmodel,nm,mutations,ng,ref_v,robust_bits,flag,robust_degrees);
      return gpolicy(std::move(ng),gametes);
    }
  }







}
#endif /* _FWDPP_MUTATION_TCC_ */
