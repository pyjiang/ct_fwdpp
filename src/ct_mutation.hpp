// added on 2.21.17
// custom type mutation type from mutation_base type
// added an extra value of scale
#ifndef _CT_MUTATION_
#define _CT_MUTATION_

#include <fwdpp/forward_types.hpp>

namespace CT
{
  struct ct_mutation : public KTfwd:: mutation_base {
    double scale;

    ct_mutation(const double & position, const double &_scale, const bool & isneutral = false, const std::uint16_t x = 0):
            KTfwd:: mutation_base(position, isneutral,x), scale(_scale)
    {
    }

    ct_mutation(const ct_mutation & ct_mut):
           KTfwd:: mutation_base(ct_mut), scale(ct_mut.scale)
    {
    }
  };
}

#endif
