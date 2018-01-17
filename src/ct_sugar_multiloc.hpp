#ifndef __CT_SUGAR_MULTILOC_MULTILOC_HPP__
#define __CT_SUGAR_MULTILOC_MULTILOC_HPP__

#include <type_traits>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
// #include <fwdpp/sugar/popbase.hpp>

#include "ct_sugar_popbase.hpp"


namespace CT {
  namespace sugar {
    /*!
      \brief Abstraction of what is needed to simulate a multilocus
      simulation using an individual-based sampler from fwdpp.

      All that is missing is the mutation_type and the container types.

      See @ref md_md_sugar for rationale, etc.

      \ingroup sugar
    */
    template<typename mutation_type,
	     typename mcont,
	     typename gcont1,
       typename gcont2,
	     typename dipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    struct multiloc : public popbase_2loc<mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type>
    {
      //! Dispatch tags for other parts of sugar layer
      using popmodel_t =  KTfwd::sugar::MULTILOCPOP_TAG;
      //! Typedef for base class
      using popbase_t = popbase_2loc<mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type>;

      using mcount_t = std::vector<uint_t>;

      using gcont1_t = gcont1;

      using gcont2_t = gcont2;

      //! Population size
      uint_t N;

      //! Container of individuals
      typename popbase_t::dipvector_t haploids; // change to haploid here

      // template <typename Vec > // gamete type
      // multiloc( const uint_t & popsize,
      //            const Vec & h_vec,
      //           //  const fitness_policy_type & ff,
      //            const uint & robust_geno,
      //            const int & robust_degrees,
      //            const Vec & wv,
      //            typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
      //            popbase_t(popsize,h_vec,robust_geno,robust_degrees, wv,reserve_size),
      //            N(popsize),
      //            haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0,0)))
      // {
      //
      // }

      // modified on 5.16.17
      // add back scale, add fitness, remove scale_func_flag, robust_geno, robust_degrees
      // modified on 2.1.17
      // remove scale in the parameter
      // modified on 1.31.16
      // gametes without h, for robustness don't evolve, but have any arbitrary scale (alpha) value
      template <typename Vec >
      multiloc( const uint_t & popsize,
                double fitness,
                //  const fitness_policy_type & ff,
                //  const uint & robust_geno,
                //  const int & robust_degrees,
                 double scale,
                //  int scale_func_flag,
                 const Vec & wv):
                 popbase_t(popsize, fitness, scale, wv),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0,0)))
      {

      }

      // addede on 2.1.17
      // gametes no h, with different initial robustness
      template <typename Vec,
                typename Cell_type_count_vec,
                typename Robust_Geno_Vec >
      multiloc( const uint_t & popsize,
                //  const fitness_policy_type & ff,
                const uint & robust_bits,
                 const uint & robust_degrees,
                 const Robust_Geno_Vec & robust_geno_vec,
                 const Cell_type_count_vec & ct_count_vec,
                 int scale_func_flag,
                 int robust_flag,
                 const Vec & wv):
                 popbase_t(popsize, robust_bits, robust_degrees, robust_geno_vec, ct_count_vec, scale_func_flag, robust_flag, wv),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0,0)))
      {
        // change haploid first index with different rob
        // change the initialize haploids that have the same count for the gamete
        uint prev_count =0;
        for (uint i =0 ; i< ct_count_vec.size();i++)
        {
          uint count = ct_count_vec[i];
          for (uint j= 0; j< count ;j++)
          {
            haploids[prev_count + j].first = i;
          }
          prev_count += count ;
        }

      }



      multiloc(const gcont1 & gametes1,
               const gcont2 & gametes2,
               const mcont & mutations,
               const mcount_t & mcounts,
               const dipvector & _haploids):
               popbase_t(gametes1,gametes2,mutations,mcounts),
               N(_haploids.size()),
               haploids(_haploids)
      {

      }

      // constructor from an existing obj
      multiloc(const multiloc<mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type> & pop):
               popbase_t(pop.gametes1,pop.gametes,pop.mutations,pop.mcounts),
               N(pop.haploids.size()),
               haploids(pop.haploids)
      {

      }

      bool operator==( const multiloc & rhs ) const
      {
	       return this->diploids == rhs.diploids && popbase_t::is_equal(rhs);
      }

      //! Empty all containers
      void clear()
      {
	       haploids.clear();
	       popbase_t::clear_containers();
      }
    };

  }
}

#endif
