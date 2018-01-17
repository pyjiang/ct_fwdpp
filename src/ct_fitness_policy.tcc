// use Eigen package for matrix manipulation
// also use libigl which has wrapper for indexing in Eigen
// #include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

// add the header of ct_eigen.hpp (for calculating cell types)
#include "ct_eigen.tcc"
// #include "ct_gamete.hpp"
#include "ct_singlepop.hpp"
//#include "trace_robust_gametes.hpp"

namespace CT
{
  struct ct_fit_simple // this struct is to recalculate fitness
  {
        
        // calcualte fitness from basic gamete type (with wv in gamete type)
        // vector double type
        template <typename ct_gamete_t>
        double operator() ( ct_gamete_t & g,
                            const ArrayXd & h_vec,
                            const ArrayXb & b_vec) const
         {
               if(g.flag==0)
               {
                 return g.fitness;
               }
               else // recalculate fintess
               {
                 double fitness = calc_fitness_from_wv_h(g.wv, h_vec, b_vec );
                 g.set_fitness(fitness); // update fintess
                 g.set_flag(0); // set flag, inidicating that fitness has been calculated
                 return fitness;
               }
         }

        // two locus (with alpha )
        // PJ added a dummy variable for this case!
        // PJ modified this function, that it does not take template variable, but rather be specific
       double operator() (
                          // uint_t & flag, // flag in the gamete, indicating whether it needs recalculate or not (if it is non-zero, need to recalculate)
                          // double cur_g_fitness, // current fitness value of the gamete
                          const VectorXd & wv,
                          const VectorXd & h_vec,
                          const ArrayXb & b_vec,
                          const int dummy ) const
        {
              // recalculate fintess
              return  calc_fitness_from_wv_h(wv, h_vec, b_vec );
        }

  };
}
