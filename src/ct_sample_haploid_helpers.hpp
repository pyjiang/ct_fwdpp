// write modified helper function for two lists of gametes
#ifndef CT_INTERNAL_SAMPLE_HAPLOID_HELPERS
#define CT_INTERNAL_SAMPLE_HAPLOID_HELPERS

namespace CT{
  namespace CT_internal
  {
    template<typename gcont1_t,
             typename gcont2_t,
	           typename mcont_t>
        inline void process_gametes( const gcont1_t & gametes1,
                                     const gcont2_t & gametes2,
				                             const mcont_t & mutations,
				                             std::vector<uint_t> & mcounts)
    /*!
      For every non-extinct gamete, increment the count of its mutations
      using the frequency of the gamete.

      This is usually the most expensive function call in a simulation.
    */
    {
      //zero out mcounts
      for(auto & mc : mcounts) mc=0;
      if(mutations.size()>mcounts.size())
    	{
    	  mcounts.resize(mutations.size(),0);
    	}
      // update for gametes1
      //update mutation counts
      for(const auto & g : gametes1)
    	{
    	  const auto n = g.n;
    	  if(n) //only do this for extant gametes
    	    {
    	      for(const auto & m : g.mutations) mcounts[m]+=n;
    	      for(const auto & m : g.smutations) mcounts[m]+=n;
    	    }
    	}
      // update for gametes2
      for(const auto & g : gametes2)
    	{
    	  const auto n = g.n;
    	  if(n) //only do this for extant gametes
    	    {
    	      for(const auto & m : g.mutations) mcounts[m]+=n;
    	      for(const auto & m : g.smutations) mcounts[m]+=n;
    	    }
    	}

    }
  }

}

#endif
