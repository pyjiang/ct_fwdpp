#ifndef _CT_SINGLEPOP_
#define _CT_SINGLEPOP_


#include <utility>
#include <vector>
// #include <unordered_set>
#include <fwdpp/fwd_functional.hpp>

#include "ct_gamete.hpp"
#include "ct_sugar_singlepop.hpp"
#include <fwdpp/tags/gamete_tags.hpp>
#include <unordered_map>
#include <fwdpp/sugar/singlepop.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//#include "new_haploids.hpp"
#include "ct_eigen.hpp"


// define cell type single pop (haploid case!)
// change lookup table to unordered_map! (to have index info)
namespace CT
{

  using ct_gamete_t = ct_gamete<void>;
  using namespace Eigen;
  // using ct_gamete_scale_wij_t = ct_gamete_scale_wij<void,VectorXd>; // created on 1.23.17
  // using ct_gamete_robust_geno_t= ct_gamete_robust_geno<void,VectorXd>;
  using ct_gamete_wv_t = ct_gamete_wv<void,VectorXd>;




  // single_pop_w_h with b vector, and ct  [in order to include everything into pop!!]
  template<typename mutation_type,
           typename mcont,
           typename gcont,
           typename dipvector,
           typename mvector,
           typename ftvector,
           typename lookup_table_type >
  struct ct_singlepop_w_b : public sugar::singlepop<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>
  {
    ArrayXb b_vec;
    Cell_type ct; // the cell type struct that holds the current reference v, w, etc.

    using mcount_t = std::vector<uint_t>;

    // added on 1.12.17
    // for basic type, for gamete type ct_gamete_wv
    template<typename Vec>
    ct_singlepop_w_b(const uint_t & popsize,
                            double fitness,
                            const Vec & wv,
                            const ArrayXb & _b_vec,
                            const Cell_type & _ct) :  sugar::singlepop<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>(popsize,fitness,wv),
                                                      b_vec (_b_vec), ct(_ct)
    {

    }



    // created on 1.30.17
    // work with scale_vec (but do not change robsut bits ) : don't evolve
    // initial pop with different degrees of robustness
    template<typename Cell_type_count_vec,
             typename Scale_Vec>
    ct_singlepop_w_b(const uint_t & popsize,
                            double fitness,
                            const uint_t & robust_geno,
                            const Scale_Vec & scale_vec,
                            const Cell_type_count_vec & ct_count_vec,
                            const ArrayXb & _b_vec,
                            const Cell_type & _ct) :  sugar::singlepop<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>(popsize,fitness,robust_geno,scale_vec, ct_count_vec, _ct.wv),
                                                      b_vec (_b_vec), ct(_ct)
    {

    }

    // add on 2.2.17
    // do not pass Scale_Vec, rather, pass robust_geno_vec
    template<typename Cell_type_count_vec,
             typename Robust_Geno_Vec>
    ct_singlepop_w_b(const uint_t & popsize,
                            double fitness,
                            const uint & Lv,
                            const uint & robust_degrees,
                            const Robust_Geno_Vec & robust_geno_vec,
                            const Cell_type_count_vec & ct_count_vec,
                            int scale_function_flag,// dummy in this case (only needed when evolve robustness)
                            const ArrayXb & _b_vec,
                            const Cell_type & _ct) :  sugar::singlepop<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>(popsize,fitness,Lv, robust_degrees,robust_geno_vec, ct_count_vec,scale_function_flag,  _ct.wv),
                                                      b_vec (_b_vec), ct(_ct)
    {

    }




    // // initialize from initial parameters
    // // ct_gamete_with_h_t
    // template<typename Vec>
    // ct_singlepop_with_h_w_b(const uint_t & popsize,
    //                         const Vec & h_vec,
    //                         double fitness,
    //                         const int & robust_geno,
    //                         const int & robust_degrees,
    //                         const Vec & wv,
    //                         const ArrayXb & _b_vec,
    //                         const Cell_type & _ct) :  sugar::singlepop<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>(popsize,h_vec,fitness,robust_geno,robust_degrees,wv,0),
    //                                                   b_vec (_b_vec), ct(_ct)
    // {
    //
    // }

    // initialize from saved struct
    ct_singlepop_w_b( const gcont & gametes,
                             const mcont & mutations,
                             const mcount_t & mcounts,
                             const dipvector & haploids,
                             const ArrayXb & _b_vec,
                             const Cell_type & _ct) :  sugar::singlepop<mutation_type,mcont,gcont,dipvector,mvector,ftvector,lookup_table_type>(gametes,mutations,mcounts,haploids), b_vec(_b_vec), ct(_ct)
    {

    }

  };

  // template<typename mtype,
  //          typename haploid_t = std::size_t>
  // using singlepop_h_w_b_t = ct_singlepop_with_h_w_b <mtype,
  //                                           std::vector<mtype>,
  //                                           std::vector<ct_gamete_with_h_t>,
  //                                           std::vector<haploid_t>,
  //                                           std::vector<mtype>,
  //                                           std::vector<uint_t>,
  //                                           std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>
  //                                         >;

  // // created on 1.23.17
  // // for new type of evolving robustness, to change wij times v
  // template<typename mtype,
  //          typename haploid_t = std::size_t>
  // using singlepop_wij_b_t = ct_singlepop_w_b <mtype,
  //                                           std::vector<mtype>,
  //                                           std::vector<ct_gamete_scale_wij_t>,
  //                                           std::vector<haploid_t>,
  //                                           std::vector<mtype>,
  //                                           std::vector<uint_t>,
  //                                           std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>
  //                                         >;

  // for interface
  template<typename mtype,
           typename haploid_t = std::size_t >
  using popbase_basic_t = sugar::popbase<mtype,
                                  std::vector<mtype>,
                                  std::vector<ct_gamete_wv_t>,
                                  std::vector<haploid_t>,
                                  std::vector<mtype>,
                                  std::vector<uint_t>,
                                  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;

  // for interface
  // need to create instances of class that have the same variable structure as singlepop_basic_t, but is a base class it inherits from
  template<typename mtype,
           typename haploid_t = std::size_t>
  using ct_singlepop_basic_base = sugar::singlepop<mtype,
                                              std::vector<mtype>,
                                              std::vector<ct_gamete_wv_t>,
                                              std::vector<haploid_t>,
                                              std::vector<mtype>,
                                              std::vector<uint_t>,
                                              std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>
                                              >;


  // for basic gamete type
  template<typename mtype,
           typename haploid_t = std::size_t>
  using singlepop_basic_t = ct_singlepop_w_b <mtype,
                                            std::vector<mtype>,
                                            std::vector<ct_gamete_wv_t>,
                                            std::vector<haploid_t>,
                                            std::vector<mtype>,
                                            std::vector<uint_t>,
                                            std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>
                                            >;


}

#endif
