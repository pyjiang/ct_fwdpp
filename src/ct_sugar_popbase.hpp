#ifndef CT_SUGAR_POPBASE_HPP
#define CT_SUGAR_POPBASE_HPP

#include <type_traits>
#include <vector>
#include <fwdpp/type_traits.hpp>

#include "ct_eigen.hpp"

// modified on 1.11.17
// reserve for gamete size change to popsize (not 4*popsize)
// this is modified from popbase.hpp (need to change the constructor for gamete!)
namespace CT
{
  namespace sugar
  {
    template<typename mutation_type,
	     typename mcont,
	     typename gcont,
	     typename dipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class popbase
    /*!
      \ingroup sugar
      \brief Base class for population objects
      \note Added in fwdpp 0.5.0
     */
    {
      static_assert(typename KTfwd::traits::is_gamete_t<typename gcont::value_type>::type(),
		    "gcont::value_type must be a gamete type" );
      static_assert(typename KTfwd::traits::is_mutation_t<typename mcont::value_type>::type(),
		    "mcont::value_type must be a mutation type" );

    public:
      //! Mutation type
      using mutation_t = mutation_type;
      //! Gamete type
      using gamete_t = typename gcont::value_type;
      //! Diploid vector type
      using dipvector_t = dipvector;
      //! Diploid type
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation vec type
      using mcont_t = mcont;
      //! Mutation count vector type
      using mcount_t = std::vector<uint_t>;
      //! Gamete vec type
      using gcont_t = gcont;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;
      //! container type for fixations
      using mvector_t = mvector;
      //! container type for fixation times
      using ftvector_t = ftvector;

      mcont_t mutations;
      /*!
	Used to keep track of mutation frequencies.

	Should have memory reserved externally,
	based on some good guess.
      */

      //! Contains number of times each mutation exists
      mcount_t mcounts;


      //! Container of gametes
      gcont gametes;

      /*!
	Vectors for holding copies of pointers to mutations during recombination.
	The requirement to declare these was introduced in fwdpp 0.3.3.

	In previous versions of the library, vectors like this had to be allocated
	for every crossover event for every generation.  The result was an excessive
	number of requests for memory allocation.

	Now, we create the vector once per replicate.  Further, we will reserve memory
	here, to minimize reallocs, etc., within fwdpp.

	Internally, fwdpp's job is to make sure that this vector is appropriately
	and efficiently cleared, but only when needed.

	\note: if not using the sugar features, you can create these vectors
	only once per simulation...
      */
      typename gamete_t::mutation_container neutral,selected;

      /*!
	\brief Can be used to track positions of segregating mutations.
	\note Must have interface like std::map or std::unordered_set
      */
      lookup_table_type mut_lookup;
      //! Vector of mutation_t to track fixations
      mvector fixations;
      /*! \brief vector<uint_t> records times when mutation_ts
	were added to mut_lookup
      */
      ftvector fixation_times;

      //! Constructor
      // for the basic type, add wv into gametes as well
      template<typename Vec>
      popbase( const uint_t & popsize,
               double fitness,
               const Vec & wv,
	             typename gamete_t::mutation_container::size_type reserve_size = 100) :
            	 mutations(mcont_t()),
            	 mcounts(mcount_t()),
            	 //The population contains a single gamete in 2N copies
            	 gametes(gcont(1,gamete_t(2*popsize,fitness,wv))),
            	 neutral(typename gamete_t::mutation_container()),
            	 selected(typename gamete_t::mutation_container()),
            	 mut_lookup(lookup_table_type()),
            	 fixations(mvector()),
            	 fixation_times(ftvector())
      {
      	//This is a good number for reserving,
      	//allowing for extra allocations when recycling is doing its thing
      	gametes.reserve(popsize);
      	//Reserve memory
      	neutral.reserve(reserve_size);
      	selected.reserve(reserve_size);
      }


      // // used for two locus initialization, this is for initilizing robustness
      // // created on 1.24.17
      // // initialize with only popsize
      // popbase( const uint_t & popsize,
      //          typename gamete_t::mutation_container::size_type reserve_size = 100) :
      //         //No muts in the population
      //         mutations(mcont_t()),
      //         mcounts(mcount_t()),
      //         //The population contains a single gamete in 2N copies
      //         gametes(gcont(1,gamete_t(2*popsize))),
      //         neutral(typename gamete_t::mutation_container()),
      //         selected(typename gamete_t::mutation_container()),
      //         mut_lookup(lookup_table_type()),
      //         fixations(mvector()),
      //         fixation_times(ftvector())
      // {
      //   //This is a good number for reserving,
      //   //allowing for extra allocations when recycling is doing its thing
      //   gametes.reserve(popsize);
      //   //Reserve memory
      //   neutral.reserve(reserve_size);
      //   selected.reserve(reserve_size);
      // }






      // Constructor 2
      // constructing popbase from mutations, gametes, as parameters
      popbase(const gcont_t & _gametes, const mcont_t & _mutations, const mcount_t & _mcounts, typename gamete_t::mutation_container::size_type reserve_size = 100):
              mutations(_mutations),
              mcounts(_mcounts),
              gametes(_gametes),
              mut_lookup(lookup_table_type()),
              fixations(mvector()),
              fixation_times(ftvector())
      {
        //Reserve memory
      	neutral.reserve(reserve_size);
      	selected.reserve(reserve_size);
        // need to initialize mut_lookup table using mutations and mcounts!
        for (uint_t i=0;i< _mcounts.size();i++)
        {
          if (_mcounts[i]!=0) // if it is non zero, need to insert into lookup table
          {
            mut_lookup.insert({_mutations[i].pos,i});
          }
        }
      }

        bool is_equal( const popbase & rhs ) const
        {
        	  return this->mutations == rhs.mutations &&
        	  this->mcounts == rhs.mcounts &&
        	  this->gametes == rhs.gametes &&
        	  this->fixations == rhs.fixations &&
        	  this->fixation_times == rhs.fixation_times;
      }


      //
      // // modified on 1.5.17
      // // pass initial fitness (rather than calculate it)
      // template<typename Vec>
      // popbase( const uint_t & popsize,
      //          const Vec & h_vec,
      //          double fitness,
      //          const int robust_geno,
      //          const int robust_degrees,
      //          const Vec & wv,
      //          typename gamete_t::mutation_container::size_type reserve_size = 100) :
      //         //No muts in the population
      //         mutations(mcont_t()),
      //         mcounts(mcount_t()),
      //         //The population contains a single gamete in 2N copies
      //         gametes(gcont(1,gamete_t(2*popsize,h_vec,fitness,robust_geno,robust_degrees,wv))),
      //         neutral(typename gamete_t::mutation_container()),
      //         selected(typename gamete_t::mutation_container()),
      //         mut_lookup(lookup_table_type()),
      //         fixations(mvector()),
      //         fixation_times(ftvector())
      // {
      //     //This is a good number for reserving,
      //     //allowing for extra allocations when recycling is doing its thing
      //     gametes.reserve(popsize);
      //     //Reserve memory
      //     neutral.reserve(reserve_size);
      //     selected.reserve(reserve_size);
      // }


      // this is for fixed level of alpha
      // modified on 1.31.17
      // add extra parameter, inidcating which scale function to use
      // modified on 1.28.17
      // add scale
      // created on 1.23.17
      // for gamete type ct_gamete_scale_wij, population with fixed scale value
      template<typename Vec>
      popbase( const uint_t & popsize,
               double fitness,
              //  const int robust_geno,
              //  const int robust_degrees,
               double scale,
              //  int scale_func_flag,
               const Vec & wv,
               typename gamete_t::mutation_container::size_type reserve_size = 100) :
              //No muts in the population
              mutations(mcont_t()),
              mcounts(mcount_t()),
              //The population contains a single gamete in 2N copies
              gametes(gcont(1,gamete_t(2*popsize,fitness,wv,scale))),
              neutral(typename gamete_t::mutation_container()),
              selected(typename gamete_t::mutation_container()),
              mut_lookup(lookup_table_type()),
              fixations(mvector()),
              fixation_times(ftvector())
      {
          //This is a good number for reserving,
          //allowing for extra allocations when recycling is doing its thing
          gametes.reserve(popsize);
          //Reserve memory
          neutral.reserve(reserve_size);
          selected.reserve(reserve_size);
      }


      // for evolving robustness, given initial different levels of robustness with different counts
      // modified on 2.2.17
      // here should pass Lv (total length, including robust bits), rather than robust bits here! (because it is one locus!, the distance is normalized by the total length)
      // created on 2.2.17
      // pass robust_geno_vec
      template<typename Vec,
               typename Cell_type_count_vec,
               typename Robust_Geno_Vec>
      popbase( const uint_t & popsize,
               double fitness,
               const uint & Lv,
               const uint & robust_degrees,
               const Robust_Geno_Vec & robust_geno_vec,
               const Cell_type_count_vec & ct_count_vec,
               int scale_function_flag,
               const Vec & wv,
               typename gamete_t::mutation_container::size_type reserve_size = 100) :
              //No muts in the population
              mutations(mcont_t()),
              mcounts(mcount_t()),
              //The population contains a single gamete in 2N copies
              gametes(gcont(robust_geno_vec.size(),gamete_t(2*popsize,fitness,scale_function_flag,wv))),
              neutral(typename gamete_t::mutation_container()),
              selected(typename gamete_t::mutation_container()),
              mut_lookup(lookup_table_type()),
              fixations(mvector()),
              fixation_times(ftvector())
      {
          //This is a good number for reserving,
          //allowing for extra allocations when recycling is doing its thing
          gametes.reserve(popsize);

          // also need to initialize initial mutations for the robust locus
          // get the biggest element in the robust geno vector
          auto iter = std::max_element(robust_geno_vec.begin(),robust_geno_vec.end());
          // then add these mutations to the mutation list
          for (uint i=0; i< *iter; i++)
          {
            // set the mutations to be not neutral
            double pos = (double) i / Lv;
            typename mcont::value_type mut( pos, false);
            this -> mutations.push_back(mut);
            this -> mcounts.push_back(0);
            // also need to initialize lookup tables! for the current mutations
            this -> mut_lookup.insert({pos,i});
          }


          // initialize the initial population consisiting of the different gametes with different robustness
          for (uint i =0; i< robust_geno_vec.size() ; i++)
          {

            gametes[i].n = ct_count_vec[i];
            gametes[i].fitness = fitness;
            // will use alpha in the parameter as gamma
            gametes[i].initialize_robustness(robust_geno_vec[i],robust_degrees);
            // gametes1[i].initialize_scale(scale_vec[i]);
            // std::cout << ct_vec[i].gamma << std::flush;

            // add counts to mcount
            for (uint j=0;j < robust_geno_vec[i]; j++)
            {
              this -> mcounts[j] += ct_count_vec[i];
            }
          }


          //Reserve memory
          neutral.reserve(reserve_size);
          selected.reserve(reserve_size);
      }



      // // do not currently use
      // // added on 1.5.17
      // // do not call fitness function when initialzing the pop (will call it later)
      // template<typename Vec >
      // popbase( const uint_t & popsize,
      //          const Vec & h_vec,
      //          const int robust_geno,
      //          const int robust_degrees,
      //          const Vec & wv,
      //          typename gamete_t::mutation_container::size_type reserve_size = 100) :
      //         //No muts in the population
      //         mutations(mcont_t()),
      //         mcounts(mcount_t()),
      //         //The population contains a single gamete in 2N copies
      //         gametes(gcont(1,gamete_t(2*popsize,h_vec,robust_geno,robust_degrees,wv))),
      //         neutral(typename gamete_t::mutation_container()),
      //         selected(typename gamete_t::mutation_container()),
      //         mut_lookup(lookup_table_type()),
      //         fixations(mvector()),
      //         fixation_times(ftvector())
      // {
      //     //This is a good number for reserving,
      //     //allowing for extra allocations when recycling is doing its thing
      //     gametes.reserve(popsize);
      //     //Reserve memory
      //     neutral.reserve(reserve_size);
      //     selected.reserve(reserve_size);
      // }



      // this is used to initialize for two locus
      // added on 12.19.16
      // added a constructor without fitness function
      // for gamete type ct_gamete_wv
      popbase( const uint_t & popsize,
               const VectorXd & wv,
               typename gamete_t::mutation_container::size_type reserve_size = 100) :
              //No muts in the population
              mutations(mcont_t()),
              mcounts(mcount_t()),
              //The population contains a single gamete in 2N copies
              gametes(gcont(1,gamete_t(2*popsize,wv))),

              neutral(typename gamete_t::mutation_container()),
              selected(typename gamete_t::mutation_container()),
              mut_lookup(lookup_table_type()),
              fixations(mvector()),
              fixation_times(ftvector())
      {
          //This is a good number for reserving,
          //allowing for extra allocations when recycling is doing its thing
          gametes.reserve(popsize);
          //Reserve memory
          neutral.reserve(reserve_size);
          selected.reserve(reserve_size);
      }



      //! Empty all the containers
      void clear_containers()
      {
        	mutations.clear();
        	mcounts.clear();
        	gametes.clear();
        	mut_lookup.clear();
        	fixations.clear();
        	fixation_times.clear();
      }
    };



    // created on 12.18.16
    // create a new type derived from popbase, for two locus case
    // in which case, will just add an extra member of gamete
    template<typename mutation_type,
	     typename mcont,
	     typename gcont1, // type of gamete2
       typename gcont2, // type of gamete2
	     typename dipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    struct popbase_2loc: public popbase<mutation_type,mcont,gcont2,dipvector,mvector,ftvector,lookup_table_type>
    {
      gcont1 gametes1; // locus 1 gamete  (will be the gamete storing mutations in robust bits)

      using gamete_t = typename gcont1::value_type;

      using popbase_t = popbase<mutation_type,mcont,gcont2,dipvector,mvector,ftvector,lookup_table_type>;

      using mcount_t = std::vector<uint_t>;

      // // constructor
      // // no need to calculate fitness before simulation
      // template<typename Vec>
      //         // typename fitness_policy_type>
      // popbase_2loc( const uint_t & popsize,
      //          const Vec & h_vec,
      //         //  const fitness_policy_type & ff,
      //          const uint robust_geno,
      //          const int robust_degrees,
      //          const Vec & wv,
      //          typename gamete_t::mutation_container::size_type reserve_size = 100)
      //             : popbase_t(popsize,wv),gametes1(gcont1(1,gamete_t(2*popsize,h_vec, robust_geno,robust_degrees)))
      // {
      //   // ff(this->gametes[0].wv, gametes1[0].h_vec); // calculate fitness
      //   gametes1.reserve(popsize);
      // }

      // modified on 5.16.17
      // add back scale
      // modified on 2.1.17
      // remove scale
      // modified on 1.31.17
      // gamete without h
      template<typename Vec>
              // typename fitness_policy_type>
      popbase_2loc( const uint_t & popsize,
              //  const fitness_policy_type & ff,
              //  const uint robust_geno,
              //  const int robust_degrees,
               double fitness,
               double scale,
              //  int scale_func_flag,
               const Vec & wv)
                  : popbase_t(popsize,fitness,wv),gametes1(gcont1(1,gamete_t(scale, 2*popsize)))
      {
        // ff(this->gametes[0].wv, gametes1[0].h_vec); // calculate fitness
        gametes1.reserve(popsize);
      }


      // currently being used
      // modified on 8.8.17
      // add robust_flag
      // modified on 2.6.17
      // modified the constructor in the gamete_t for robust bits no h.
      // added on 2.1.17
      // for initial pop with different rob
      template<typename Vec,
               typename Cell_type_count_vec,
               typename Robust_Geno_Vec >
              // typename fitness_policy_type>
      popbase_2loc( const uint_t & popsize,
              //  const fitness_policy_type & ff,
               const uint & robust_bits,
               const uint & robust_degrees,
               const Robust_Geno_Vec & robust_geno_vec,
               const Cell_type_count_vec & ct_count_vec,
               int scale_func_flag,
               int robust_flag,
               const Vec & wv)
                  : popbase_t(popsize,wv),gametes1(gcont1(robust_geno_vec.size(),gamete_t(2*popsize)))
      {
        // ff(this->gametes[0].wv, gametes1[0].h_vec); // calculate fitness
        gametes1.reserve(popsize);

        uint mut_count =0;
        // also need to initialize initial mutations for the robust locus
        // get the biggest element in the robust geno vector
        auto iter = std::max_element(robust_geno_vec.begin(),robust_geno_vec.end());
        // then add these mutations to the mutation list

        if(robust_flag==0) //using 2^n
          mut_count = robust_bits;
        else if(robust_flag==1)
          mut_count = *iter;
        // will use robust_flag here
        // if additive, use the element in the robust_geno vector
        // if using 2^n, will use robust bits
        for (uint i=0; i< mut_count ; i++)
        {
          // set the mutations to be not neutral
          double pos = (double) i / robust_bits;
          typename mcont::value_type mut( pos, false);
          this -> mutations.push_back(mut);
          this -> mcounts.push_back(0);
          // also need to initialize lookup tables! for the current mutations
          this -> mut_lookup.insert({pos,i});
        }


        // initialize the initial population consisiting of the different gametes with different robustness

        for (uint i =0; i< robust_geno_vec.size() ; i++)
        {

          gametes1[i].n = ct_count_vec[i];
          // will use alpha in the parameter as gamma
          gametes1[i].initialize_robustness_all(scale_func_flag, robust_geno_vec[i],robust_bits, robust_degrees);
          // gametes1[i].initialize_scale(scale_vec[i]);
          // std::cout << ct_vec[i].gamma << std::flush;

          // add counts to mcount
          if(robust_flag ==1)
          {
            for (uint j=0;j < robust_geno_vec[i]; j++)
            {
              this -> mcounts[j] += ct_count_vec[i];
            }
          }
          else if(robust_flag ==0)
          {
            int mask =1;
            int val;
            for (uint j=0; j< robust_bits; j++)
            {
              val = robust_geno_vec[i] & mask ;
              if(val)
                this -> mcounts[j] += ct_count_vec[i];
              mask <<= 1;
            }
          }
        }

      }



      popbase_2loc( const gcont1 _gametes1,
               const gcont2 _gametes2,
               const mcont & mutations,
               const mcount_t & mcounts)
                  : popbase_t(_gametes2,mutations,mcounts),gametes1(_gametes1)
      {
      }



    };

  }
}

#endif
