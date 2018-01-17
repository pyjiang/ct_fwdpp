// write a few gsl functions that it will be easier to call from interface
#include "gsl_funcs.hpp"

// initialize gsl
// void initialize_gsl(gsl_rng *r,unsigned seed)
// {
//   gsl_rng_set(r,seed);
// }

// create gsl
GSLrng create_gsl(unsigned seed)
{
  GSLrng r(seed);
  return r;
}
