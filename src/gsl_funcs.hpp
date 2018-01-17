#ifndef _GSL_FUNC_
#define _GSL_FUNC_

#include <fwdpp/sugar/GSLrng_t.hpp>

using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

GSLrng create_gsl(unsigned seed);

#endif
