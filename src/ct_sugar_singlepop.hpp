#ifndef __CT_SUGAR_SINGLEPOP_SINGLEPOP_HPP__
#define __CT_SUGAR_SINGLEPOP_SINGLEPOP_HPP__

/*
  A structure representing a single Wright-Fisher population.
  The user initizializes it with a population size, N
*/

#include <fwdpp/sugar/poptypes/tags.hpp>
// #include <fwdpp/sugar/popbase.hpp>
#include "ct_sugar_popbase.hpp"
#include <fwdpp/sugar/singlepop.hpp>
#include "ct_eigen.hpp"

namespace CT {
  namespace sugar {

    // using ct_gamete_robust_geno_t= ct_gamete_robust_geno<void,VectorXd>;

    // using namespace KTfwd::sugar;
    /*!
      \brief Abstraction of what is needed to simulate a single population
      using an individual-based sampler from fwdpp

      All that is missing is the mutation_type and the container types.

      See @ref md_md_sugar for rationale, etc.

      \ingroup sugar
    */
    // for haploid!
    template<typename mutation_type,
	           typename mcont,
	           typename gcont,
	           typename dipvector,
	           typename mvector,
	           typename ftvector,
             typename lookup_table_type>
    class singlepop : public popbase<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>
    {
    public:
      //! Population size
      uint_t N;

      //! Typedef for base class
      using popbase_t = popbase<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>;
      //! Dispatch tag for other parts of sugar layer
      using popmodel_t = KTfwd::sugar::SINGLEPOP_TAG;
      //! Fitness function signature compatible with this type
      using fitness_t = KTfwd::traits::fitness_fxn_t<typename popbase_t::dipvector_t,
						     typename popbase_t::gcont_t,
						     typename popbase_t::mcont_t>;

      //! Container of diploids
      typename popbase_t::dipvector_t haploids;


      // for basic type
      template <typename Vec>
      singlepop( const uint_t & popsize,
                 double fitness,
                 const Vec & wv,
		             typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100) :
	               popbase_t(popsize,fitness,wv,reserve_size),
	               N(popsize),
	               //All N diploids contain the only gamete in the pop
	               haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
      {
      }

      // added on 1.24.17
      // for initizialize empty popbase with just popsize
      singlepop( const uint_t & popsize,
		             typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100) :
	               popbase_t(popsize),
	               N(popsize),
	               //All N diploids contain the only gamete in the pop
	               haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
      {
      }


      // added on 1.23.17
      // for population initialize with cell types with different degrees of gamma
      template <typename Cell_type_vec,
                typename Cell_type_count_vec>
      singlepop( const uint_t & popsize,
                 double fitness,
                 const Cell_type_vec & ct_vec,
                 const Cell_type_count_vec & ct_count_vec,
                 double initial_alpha,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100) :
                 popbase_t(popsize,fitness,ct_vec,ct_count_vec,initial_alpha, reserve_size),
                 N(popsize),
                 //All N diploids contain the only gamete in the pop
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
      {
        // change the initialize haploids that have the same count for the gamete
        uint prev_count =0;
        for (uint i =0 ; i< ct_count_vec.size();i++)
        {
          uint count = ct_count_vec[i];
          for (uint j= 0; j< count ;j++)
          {
            haploids[prev_count + j] = i;
          }
          prev_count += count ;
        }

      }




      // inialize singlepop with gametes, mutations, etc
      template <typename mcount_t>
      singlepop( const gcont & gametes, const mcont & mutations, const mcount_t & mcounts,
                 const dipvector & _diploids,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100) :
                 popbase_t(gametes, mutations,mcounts,reserve_size),
                 N(_diploids.size()),
                 // constructor of diploid to copy from the paramter!
                 haploids(_diploids)
      {
      }




      // updated on 11.22.16
      //! Constructor for gametes with h
      // with robust_geno
      // enable it if gamete type is ct_gamete_with_h_t
      template <typename fitness_policy_type,
                typename Vec > // gamete type
      singlepop( const uint_t & popsize,
                 const Vec & h_vec,
                 const fitness_policy_type & ff,
                 const int & robust_geno,
                 const int & robust_degrees,
                 const Vec & wv,
                 const int dummy, // dummy variable, just used to distinguish this function from the one below
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
                 popbase_t(popsize,h_vec,ff,robust_geno,robust_degrees, wv,reserve_size),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
                 //All N diploids contain the only gamete in the pop
      {
      }

      // modified on 1.28.17
      // add scale
      // added on 1.23.17
      // for gamete type : ct_gamete_scale_wij
      template < typename Vec > // gamete type
      singlepop( const uint_t & popsize,
                 double fit,
                 const int & robust_geno,
                 const int & robust_degrees,
                 double scale,
                 const Vec & wv,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
                 popbase_t(popsize,fit,robust_geno,robust_degrees, scale, wv,reserve_size),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
                 //All N diploids contain the only gamete in the pop
      {
      }

      // created on 1.30.17
      // for gamete type : ct_gamete_scale_wij
      // but with initial different degrees of rob
      template < typename Vec,
                 typename Cell_type_count_vec,
                 typename Scale_Vec  >
      singlepop( const uint_t & popsize,
                 double fit,
                 const uint_t & robust_geno,
                 const Scale_Vec & scale_vec,
                 const Cell_type_count_vec  & ct_count_vec,
                 const Vec & wv,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
                 popbase_t(popsize,fit,robust_geno, scale_vec, ct_count_vec, wv, reserve_size),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
                 //All N diploids contain the only gamete in the pop
      {
        // change the initialize haploids that have the same count for the gamete
        uint prev_count =0;
        for (uint i =0 ; i< ct_count_vec.size();i++)
        {
          uint count = ct_count_vec[i];
          for (uint j= 0; j< count ;j++)
          {
            haploids[prev_count + j] = i;
          }
          prev_count += count ;
        }

      }

      // created on 2.2.17
      // pass robust_geno rather than scale_vec
      template < typename Vec,
                 typename Cell_type_count_vec,
                 typename Robust_Geno_Vec  >
      singlepop( const uint_t & popsize,
                 double fit,
                 const uint & Lv,
                 const uint & robust_degrees,
                 const Robust_Geno_Vec & robust_geno_vec,
                 const Cell_type_count_vec  & ct_count_vec,
                 int scale_function_flag,
                 const Vec & wv,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
                 popbase_t(popsize,fit,Lv,robust_degrees, robust_geno_vec, ct_count_vec, scale_function_flag, wv, reserve_size),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
                 //All N diploids contain the only gamete in the pop
      {
        // change the initialize haploids that have the same count for the gamete
        uint prev_count =0;
        for (uint i =0 ; i< ct_count_vec.size();i++)
        {
          uint count = ct_count_vec[i];
          for (uint j= 0; j< count ;j++)
          {
            haploids[prev_count + j] = i;
          }
          prev_count += count ;
        }

      }



      // added on 1.5.17
      // do not have fitness function in the constructor
      template < typename Vec > // gamete type
      singlepop( const uint_t & popsize,
                 const Vec & h_vec,
                 const int & robust_geno,
                 const int & robust_degrees,
                 const Vec & wv,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
                 popbase_t(popsize,h_vec,robust_geno,robust_degrees, wv,reserve_size),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0)))
                 //All N diploids contain the only gamete in the pop
      {
      }


      template <typename fitness_policy_type,
                typename Vec > // gamete type
      singlepop( const uint_t & popsize,
                 const Vec & h_vec,
                 const fitness_policy_type & ff,
                 const int & robust_geno,
                 const Vec & wv,
                 typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100):
                 popbase_t(popsize,wv,robust_geno,reserve_size),
                 N(popsize),
                 haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0,h_vec,ff)))
                 //All N diploids contain the only gamete in the pop
      {
        // initialize_singlepop<typename gcont::value_type>(popsize,h_vec,ff,robust_geno);
      }



      //  new constructor for haploids that have h_vec (to initialize them having the same initial h_vec )
      // enable if gamete type is ct_gamete_robust_geno





      // template <typename fitness_policy_type,
      //           typename Vec>
      // singlepop( const uint_t & popsize,
      //            const Vec & h_vec,
      //            const fitness_policy_type & ff,
      //            const int & robust_geno,
      //            typename popbase_t::gamete_t::mutation_container::size_type reserve_size = 100) :
      //            popbase_t(popsize,h_vec,ff,robust_geno,reserve_size),
      //            N(popsize),
      //            //All N diploids contain the only gamete in the pop
      //            haploids(typename popbase_t::dipvector_t(2*popsize,typename popbase_t::diploid_t(0,h_vec)))
      // {
      // }

      bool operator==( const singlepop & rhs ) const
      {
	       return this->haploids == rhs.haploids && popbase_t::is_equal(rhs);
      }

      //! Empty all the containers
      void clear()
      {
	       haploids.clear();
	       popbase_t::clear_containers();
      }
    };






  }
}
#endif
