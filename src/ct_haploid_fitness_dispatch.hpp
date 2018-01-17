#ifndef __CT_INTERNAL_HAPLOID_FITNESS_DISPATCH_HPP__
#define __CT_INTERNAL_HAPLOID_FITNESS_DISPATCH_HPP__

#include <type_traits>
#include <fwdpp/tags/diploid_tags.hpp>

namespace CT {
  namespace CT_internal{
    //! Custom haploid
    template<typename fitness_policy_type,
	           typename haploid_t,
	           typename gcont_t>
    inline double haploid_fitness_dispatch( const fitness_policy_type & fp, const haploid_t & hap,
					     gcont_t & gametes) {
      return fp(gametes[hap]);
    }
  }
}


#endif
