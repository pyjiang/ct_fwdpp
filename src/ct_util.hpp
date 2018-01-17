/// \file util.hpp
#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <set>
#include <map>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "ct_gamete.hpp"

namespace CT
{
  using namespace KTfwd;
  /*!
    Label all extinct and fixed variants for recycling

    \note: lookup must be compatible with lookup->erase(lookup->find(double))
  */
  template<typename mcont_t,
	   typename mutation_lookup_table>
  void update_mutations( mcont_t & mutations,
			 mutation_lookup_table & lookup,
			 std::vector<uint_t> & mcounts,
			 const unsigned twoN)
  {
    static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    assert(mcounts.size()==mutations.size());
    for(std::size_t i = 0 ; i < mcounts.size() ; ++i)
      {
	       assert(mcounts[i] <= twoN);
	        if(mcounts[i]==twoN || !mcounts[i] )
	         {
	            lookup.erase(mutations[i].pos);
	            mcounts[i]=0;
	         }
      }
  }

  /*!
    Label all extinct  variants for recycling

    \Note: lookup must be compatible with lookup->erase(lookup->find(double))
  */
  template<typename mcont_t,
	   typename mutation_lookup_table>
  void update_mutations( const mcont_t & mutations,
			 mutation_lookup_table & lookup,
			 std::vector<uint_t> & mcounts)

  {
    static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    for(std::size_t i = 0 ; i < mcounts.size() ; ++i)
      {
	if( !mcounts[i] )
	  {
	    lookup.erase(mutations[i].pos);
	    mcounts[i]=0;
	  }
      }
  }

  /*!
    Label all fixed and all extinct variants for recycling. Copy fixations and fixation times
    into containers.

    \note: lookup must be compatible with lookup->erase(lookup->find(double))
  */

  // added gametes as parameter 10.17.16
  template<typename mcont_t,
	         typename fixation_container_t,
	         typename fixation_time_container_t,
	         typename mutation_lookup_table,
           typename Vecd >
          //  typename ct_gamete>
  void update_mutations( mcont_t & mutations,
			                   fixation_container_t & fixations,
			                   fixation_time_container_t & fixation_times,
			                   mutation_lookup_table & lookup,
			                   std::vector<uint_t> & mcounts,
			                   const unsigned & generation,
                         const unsigned & twoN,
                        //  ct_gamete & gametes,
                         VectorXi & ref_v,
                         const MatrixXd & w,
                         const VectorXd & h,
                         Vecd  & val_ref )
  {
    static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    assert(mcounts.size()==mutations.size());

    int fix_flag=0; //flag to see if there is any fixation event

    uint_t Lv= ref_v.size();

    for(unsigned i=0;i<mcounts.size();++i)
    {
	      assert(mcounts[i] <= twoN);
        if(mcounts[i]==twoN)
         {
            fixations.push_back(mutations[i]);
            fixation_times.push_back(generation);
            // change ref_v
            // for test
            // std::cout.precision(8);
            // std::cout << i << " "; // print out fixed mutation position (in the mutation list)
            // std::cout << std::fixed << mutations[i].pos << "\t";
            // std::cout << "round" << std::round(mutations[i].pos*Lv) << "\t";
            // std::cout << "int " << int (mutations[i].pos*Lv) << "\n";
            // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " " <<std::flush;
            ref_v(std::round(mutations[i].pos*Lv))= 1-ref_v(std::round(mutations[i].pos*Lv));
            // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " "<<std::flush;
            fix_flag =1;
            mcounts[i]=0; //set count to zero to mark mutation as "recyclable"
            // lookup.erase(mutations[i].pos);
            // will remove all index in gamete list, added on 10.17.16
            // for (auto & gamete : gametes)
            // {
            //   if (gamete.n!=0)//remove mutation in smutation list!
            //   {
            //     std::remove(gamete.smutations.begin(),gamete.smutations.end(),i);
            //     // auto pos= std::find(gamete.smutations.begin(),gamete.smutations.end(),i); // find the pos in smutations list in gamete
            //     //
            //     // gamete.smutations.erase(pos);
            //   }
            // }
         }
  	    if(!mcounts[i] && std::abs(mutations[i].pos-1) > std::numeric_limits<double>::epsilon()) // avoid removing twice!! if mutation[i].pos not equal to 1 (not yet marked as removed), mark it and remove
        {
           lookup.erase(mutations[i].pos);
           mutations[i].pos= 1; // it needs to reset mutation position (a value that does not exist in the population for active mutations)
                                // so that it won't be searched!! when it has been removed already!
        }
    }

    // update val_ref if ref_v is changed
    if (fix_flag==1)
    {
      val_ref= calc_cell_type_z_val(w,ref_v,h);
    }
  }


  // modified for gamete with h, also have robust_bits
  // for bits that represent robustness: will also treat them as normal genotype
  // this actually is very similar to previous version, but instead of updating val_ref,
  // it will update wv in the population (because individual genotype will have their own h)
  template<typename mcont_t,
	         typename fixation_container_t,
	         typename fixation_time_container_t,
	         typename mutation_lookup_table,
           typename Vecd>
          //  typename ct_gamete>
  void update_mutations( mcont_t & mutations,
			                   fixation_container_t & fixations,
			                   fixation_time_container_t & fixation_times,
			                   mutation_lookup_table & lookup,
			                   std::vector<uint_t> & mcounts,
			                   const unsigned & generation,const unsigned & twoN,
                        //  ct_gamete & gametes,
                         int robust_bits,
                         VectorXi & ref_v,
                         const MatrixXd & w,
                         Vecd  & wv )
  {
    static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    assert(mcounts.size()==mutations.size());

    int fix_flag=0; //flag to see if there is any fixation event

    uint_t Lv= ref_v.size();

    for(unsigned i=0;i<mcounts.size();++i)
    {
	      assert(mcounts[i] <= twoN);
        if(mcounts[i]==twoN)
         {
            fixations.push_back(mutations[i]);
            fixation_times.push_back(generation);
            // change ref_v
            // for test
            // std::cout.precision(8);
            // std::cout << i << " "; // print out fixed mutation position (in the mutation list)
            // std::cout << std::fixed << mutations[i].pos << "\t";
            // std::cout << "round" << std::round(mutations[i].pos*Lv) << "\t";
            // std::cout << "int " << int (mutations[i].pos*Lv) << "\n";
            // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " " <<std::flush;
            ref_v(std::round(mutations[i].pos*Lv))= 1-ref_v(std::round(mutations[i].pos*Lv)); // change the bit in the ref_v
            // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " "<<std::flush;
            fix_flag =1;
            mcounts[i]=0; //set count to zero to mark mutation as "recyclable"
            // lookup.erase(mutations[i].pos);
            // will remove all index in gamete list, added on 10.17.16
            // for (auto & gamete : gametes)
            // {
            //   if (gamete.n!=0)//remove mutation in smutation list!
            //   {
            //     std::remove(gamete.smutations.begin(),gamete.smutations.end(),i);
            //     // auto pos= std::find(gamete.smutations.begin(),gamete.smutations.end(),i); // find the pos in smutations list in gamete
            //     //
            //     // gamete.smutations.erase(pos);
            //   }
            // }
         }
  	    if(!mcounts[i] && std::abs(mutations[i].pos-1) > std::numeric_limits<double>::epsilon()) // avoid removing twice!! if mutation[i].pos not equal to 1 (not yet marked as removed), mark it and remove
        {
           lookup.erase(mutations[i].pos);
           mutations[i].pos= 1; // it needs to reset mutation position (a value that does not exist in the population for active mutations)
                                // so that it won't be searched!! when it has been removed already!
        }
    }

    // update val_ref if ref_v is changed
    if (fix_flag==1)
    {
      wv= calc_cell_type_wv(w,ref_v,ref_v.size()-robust_bits);
    }
  }


  // used for basic type
  // Note: here needs to change function in previous version!!
  // modified on 12.19.16
  // remove the parameter w and robust_bits , which is not used in the function
  // new version, deal with gamete with h, with wv, so do not need to update common wv
  template<typename mcont_t,
           typename fixation_container_t,
           typename fixation_time_container_t,
           typename mutation_lookup_table>
          //  typename ct_gamete>
  void update_mutations( mcont_t & mutations,
                         fixation_container_t & fixations,
                         fixation_time_container_t & fixation_times,
                         mutation_lookup_table & lookup,
                         std::vector<uint_t> & mcounts,
                         const unsigned & generation,
                         const unsigned & twoN,
                        //  ct_gamete & gametes,
                        //  int robust_bits,
                         VectorXi & ref_v)
  {
    static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
       "mutation_type must be derived from KTfwd::mutation_base" );
    assert(mcounts.size()==mutations.size());

    uint_t Lv= ref_v.size();

    for(unsigned i=0;i<mcounts.size();++i)
    {
        assert(mcounts[i] <= twoN);
        if(mcounts[i]==twoN)
         {
            fixations.push_back(mutations[i]);
            fixation_times.push_back(generation);
            // change ref_v
            // for test
            // std::cout.precision(8);
            // std::cout << i << " "; // print out fixed mutation position (in the mutation list)
            // std::cout << std::fixed << mutations[i].pos << "\t";
            // std::cout << "round" << std::round(mutations[i].pos*Lv) << "\t";
            // std::cout << "int " << int (mutations[i].pos*Lv) << "\n";
            // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " " <<std::flush;
            ref_v(std::round(mutations[i].pos*Lv))= 1-ref_v(std::round(mutations[i].pos*Lv)); // change the bit in the ref_v
            // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " "<<std::flush;
            mcounts[i]=0; //set count to zero to mark mutation as "recyclable"
            // lookup.erase(mutations[i].pos);
            // will remove all index in gamete list, added on 10.17.16
            // for (auto & gamete : gametes)
            // {
            //   if (gamete.n!=0)//remove mutation in smutation list!
            //   {
            //     std::remove(gamete.smutations.begin(),gamete.smutations.end(),i);
            //     // auto pos= std::find(gamete.smutations.begin(),gamete.smutations.end(),i); // find the pos in smutations list in gamete
            //     //
            //     // gamete.smutations.erase(pos);
            //   }
            // }
         }
        if(!mcounts[i] && std::abs(mutations[i].pos-1) > std::numeric_limits<double>::epsilon()) // avoid removing twice!! if mutation[i].pos not equal to 1 (not yet marked as removed), mark it and remove
        {
           lookup.erase(mutations[i].pos);
           mutations[i].pos= 1; // it needs to reset mutation position (a value that does not exist in the population for active mutations)
                                // so that it won't be searched!! when it has been removed already!
        }
    }


  }

  // for all two locus cases
  // the mutations in the first locus will be [0,1), and the mutations in the second locus will be [1,2)
  template<typename mcont_t,
           typename fixation_container_t,
           typename fixation_time_container_t,
           typename mutation_lookup_table>
          //  typename ct_gamete>
  void update_mutations( mcont_t & mutations,
                         fixation_container_t & fixations,
                         fixation_time_container_t & fixation_times,
                         mutation_lookup_table & lookup,
                         std::vector<uint_t> & mcounts,
                         const unsigned & generation,const unsigned & twoN,
                        //  ct_gamete & gametes,
                        //  int robust_bits,
                         VectorXi & ref_v1,
                         VectorXi & ref_v2)
   {
     static_assert( typename traits::is_mutation_t<typename mcont_t::value_type>::type(),
        "mutation_type must be derived from KTfwd::mutation_base" );
     assert(mcounts.size()==mutations.size());

     uint_t len1= ref_v1.size();
     uint_t len2= ref_v2.size();

     for(unsigned i=0;i<mcounts.size();++i)
     {
         assert(mcounts[i] <= twoN);
         if(mcounts[i]==twoN)
          {
             fixations.push_back(mutations[i]);
             fixation_times.push_back(generation);
             // change ref_v
             // for test
             // std::cout.precision(8);
             // std::cout << i << " "; // print out fixed mutation position (in the mutation list)
             // std::cout << std::fixed << mutations[i].pos << "\t";
             // std::cout << "round" << std::round(mutations[i].pos*Lv) << "\t";
             // std::cout << "int " << int (mutations[i].pos*Lv) << "\n";
             // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " " <<std::flush;
             if(mutations[i].pos<1)
                ref_v1(std::round(mutations[i].pos*len1))= 1-ref_v1(std::round(mutations[i].pos*len1)); // change the bit in the ref_v
             else if (mutations[i].pos<2) // this is from [1,2)
                ref_v2(std::round((mutations[i].pos-1)*len2))= 1-ref_v2(std::round((mutations[i].pos-1)*len2));
             // std::cout << ref_v(std::round(mutations[i].pos*Lv)) << " "<<std::flush;
             mcounts[i]=0; //set count to zero to mark mutation as "recyclable"
             // lookup.erase(mutations[i].pos);
             // will remove all index in gamete list, added on 10.17.16
             // for (auto & gamete : gametes)
             // {
             //   if (gamete.n!=0)//remove mutation in smutation list!
             //   {
             //     std::remove(gamete.smutations.begin(),gamete.smutations.end(),i);
             //     // auto pos= std::find(gamete.smutations.begin(),gamete.smutations.end(),i); // find the pos in smutations list in gamete
             //     //
             //     // gamete.smutations.erase(pos);
             //   }
             // }
          }
         if(!mcounts[i] && std::abs(mutations[i].pos-2) > std::numeric_limits<double>::epsilon()) // avoid removing twice!! if mutation[i].pos not equal to 2 (not yet marked as removed), mark it and remove
         {
            lookup.erase(mutations[i].pos);
            mutations[i].pos= 2; // it needs to reset mutation position (a value that does not exist in the population for active mutations)
                                 // so that it won't be searched!! when it has been removed already!
         }
     }
   }




}
#endif /* _UTIL_HPP_ */
