#ifndef _CT_MUT_MODEL_
#define _CT_MUT_MODEL_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <type_traits>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/internal/recycling.hpp>
// this is a mutation model, modified from original infsites.hpp



namespace CT
{
   using uint_t = KTfwd::uint_t ;
  // store mutation index and a flag (whether it exists or not)
  struct mut_pos_indx_t
  {
    uint_t index; //index in the mutation list
    int flag;// if flag==1, exists in the mutation container; otherwise does not exist
  };

  // define a new struct, with not only mutation position, but also a flag of it whether exists in the mutation container or not
  template <typename position_result_t>
  struct mut_pos_t
  {
    position_result_t pos;// this position is the mutation position on genotype
    mut_pos_indx_t indx_obj;
    // uint_t index; //index in the mutation list
    // int flag;// if flag==1, exists in the mutation container; otherwise does not exist
    // mut_pos_t(){};
  };

  struct ct_mut_model
  {
    // generate a new mutation
    // it can be generated of any position and update the lookup table w/that position.
    // also update mutation counts
    // change lookup table from unordered_set to unordered_map (which has both key and value)

    template<typename position_t,
            typename lookup_table_t>
    inline mut_pos_t<typename std::result_of<position_t()>::type> generate_mut_pos( const position_t & posmaker,
                                                                                    lookup_table_t & lookup) const
    {
      static_assert( std::is_convertible<typename std::result_of<position_t()>::type,double>::value,
         "The return type of posmaker must be convertible to double." );
      // auto pos = posmaker();
      mut_pos_t<typename std::result_of<position_t()>::type>  pos_struct;
      pos_struct.pos =  posmaker();
      auto search=lookup.find(pos_struct.pos);

      if(search!= lookup.end())
      {
        //if find position
        // mutations[*search].xtra+=1; //use .xtra as mutation count
        pos_struct.indx_obj.index= (*search).second ; // record index
        pos_struct.indx_obj.flag=1;
      }
      else
      {
        lookup.insert({pos_struct.pos,0}); // the index in mutation list is unknown for now, will assign value out of this function!
        pos_struct.indx_obj.flag=0;
      }
      return pos_struct;
    }

    //copied from infsites.hpp
    /*!
      \param r A const gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
    */

    // modified on 2.6.17
    // remove unused parameter r
    // remove unused parameter generation
    // remove unused parameter neutral_mutation_rate
    // remove unused parameter selected_mutation_rate
    template<typename queue_t,
             typename mlist_t,
             typename lookup_table_t,
             typename position_t>
    inline mut_pos_indx_t  operator()(queue_t & recycling_bin,
                                      mlist_t & mutations,
                                      lookup_table_t & lookup,
                                      // const uint_t & generation,
                                      // const double & neutral_mutation_rate,
                                      // const double & selected_mutation_rate,
                                      const position_t & posmaker) const
    {
      // static_assert(std::is_same<typename mlist_t::value_type,KTfwd::popgenmut>::value,
      //   "mlist_t::value_type must be KTfwd::popgenmut");
      //Establish position of new mutation
      auto pos_struct = this->generate_mut_pos(posmaker,lookup);

      if (pos_struct.indx_obj.flag==0)// if need to put in a new mutation
      {
        pos_struct.indx_obj.index = KTfwd::fwdpp_internal::recycle_mutation_helper(recycling_bin,mutations,
                                                                                   pos_struct.pos,0
                                                                                 ); // from the third parameter and on, are parameters passed to mutation type constructor
                                                                                 // this function returns the mutation index for the newly arised mutation (to recycle as possible)
        // update mutation index in lookup table
        lookup[pos_struct.pos]= pos_struct.indx_obj.index ;
      }

      return pos_struct.indx_obj;
    }



  };
}

#endif
