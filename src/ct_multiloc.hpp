#ifndef __CT_SUGAR_MULTILOC_HPP__
#define __CT_SUGAR_MULTILOC_HPP__

// for multi-locus

#include <vector>
#include <unordered_set>
#include <fwdpp/sugar/multiloc/multiloc.hpp>
#include <fwdpp/fwd_functional.hpp>

#include "ct_gamete.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unordered_map>
#include "ct_sugar_multiloc.hpp"

namespace CT
{
  /*!
    \brief Single population, multilocus simulation.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  // multi locus, gamete type is ct_gamete_with_h_t
  // using ct_gamete_robust_bits_t= ct_gamete_robust_geno_only<void,VectorXd>;
  using ct_gamete_wv_t= ct_gamete_wv<void,VectorXd>;
  using ct_gamete_geno_no_h_t = ct_gamete_robust_geno_only_wij<void>; // added on 1.30.17


  // for interface
  using ct_gamete_t = ct_gamete<void>;
  // using gamete_base_t = gamete_base<void>;




  // below 2 only added for interface
  // this is the base class, should shared by both gamete with h and without h
  template<typename mtype,
           typename haploid_t = std::pair<std::size_t,std::size_t> >
  using popbase_two_t = sugar::popbase<mtype,
                                  std::vector<mtype>,
                                  std::vector<ct_gamete_wv_t>,
                                  std::vector<haploid_t>,
                                  std::vector<mtype>,
                                  std::vector<uint_t>,
                                  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;

  // template<typename mtype,
  //          typename haploid_t = std::pair<std::size_t,std::size_t> >
  // using popbase_2loc_t = sugar::popbase_2loc<mtype,
  //                                         std::vector<mtype>,
  //                                         std::vector<ct_gamete_robust_bits_t>,
  //                                         std::vector<ct_gamete_wv_t>,
  //                                         std::vector<haploid_t>,
  //                                         std::vector<mtype>,
  //                                         std::vector<uint_t>,
  //                                         std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;



    // // using haploid instead of diploid (but have two locus random segregation)
    // template<typename mtype,
    //        typename haploid_t = std::pair<std::size_t,std::size_t> >
    // using multiloc_t = sugar::multiloc<mtype,
    // 		                                  std::vector<mtype>,
    // 		                                  std::vector<ct_gamete_robust_bits_t>,
    //                                       std::vector<ct_gamete_wv_t>,
    // 		                                  std::vector<haploid_t>,
    // 		                                  std::vector<mtype>,
    // 		                                  std::vector<uint_t>,
    // 		                                  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;

  // no h
    template<typename mtype,
             typename haploid_t = std::pair<std::size_t,std::size_t> >
    using popbase_2loc_no_h_t = sugar::popbase_2loc<mtype,
                                            std::vector<mtype>,
                                            std::vector<ct_gamete_geno_no_h_t>,
                                            std::vector<ct_gamete_wv_t>,
                                            std::vector<haploid_t>,
                                            std::vector<mtype>,
                                            std::vector<uint_t>,
                                            std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;

  // no h
    template<typename mtype,
           typename haploid_t = std::pair<std::size_t,std::size_t> >
    using multiloc_no_h_t = sugar::multiloc<mtype,
                                          std::vector<mtype>,
                                          std::vector<ct_gamete_geno_no_h_t>,
                                          std::vector<ct_gamete_wv_t>,
                                          std::vector<haploid_t>,
                                          std::vector<mtype>,
                                          std::vector<uint_t>,
                                          std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;



  // create a new struct that has two pop and also b vector. [for input output]
  template<typename mutation_type,
           typename mcont,
           typename gcont1,
           typename gcont2,
           typename dipvector,
           typename mvector,
           typename ftvector,
           typename lookup_table_type>
  struct twolocus_w_b: public sugar::multiloc <mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type>
  {
    ArrayXb b_vec;
    Cell_type ct; // the cell type struct that holds the current reference v, w, etc.

    using mcount_t = std::vector<uint_t>;
    //
    // // twolocus with h
    // // initialize from initial parameters
    // template <typename Vec > // gamete type
    // twolocus_w_b(const uint_t & popsize,
    //              const Vec & h_vec,
    //              double fitness, // will be a dummy variable (although calculated, but not pass to initialize pop)
    //              const int & robust_geno,
    //              const int & robust_degrees,
    //              const Vec & wv,
    //              const ArrayXb & _b_vec,
    //              const Cell_type & _ct):
    //              sugar::multiloc <mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type> (popsize, h_vec, robust_geno, robust_degrees, wv), b_vec(_b_vec), ct(_ct)
    // {
    // }

    // modified on 5.16.17
    // add scale
    // modifed on 2.1.17
    // remove parameter scale
    // modifed on 1.31.17
    // for two locus without h
    // added a parameter about scale_function
    // population with fixed scale value
    template <typename Vec > // gamete type
    twolocus_w_b(const uint_t & popsize,
                 double fitness, // will be a dummy variable (although calculated, but not pass to initialize pop)
                //  const uint & robust_geno,
                //  const int & robust_degrees,
                 double scale, // need to have this parameter
                //  int scale_func_flag,
                 const Vec & wv,
                 const ArrayXb & _b_vec,
                 const Cell_type & _ct):
                 sugar::multiloc <mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type> (popsize, fitness, scale, wv), b_vec(_b_vec), ct(_ct)
    {
    }

    // added on 2.1.17
    // for two locus with different scales
    // added scale_vec
    // for two locus without h
    // added a parameter about scale_function
    template <typename Cell_type_count_vec,
              typename Robust_Geno_Vec> // gamete type
    twolocus_w_b(const uint_t & popsize,
                //  double fitness, // will be a dummy variable (although calculated, but not pass to initialize pop)
                 const uint & robust_bits,
                 const uint & robust_degrees,
                 const Robust_Geno_Vec & robust_geno_vec,
                 const Cell_type_count_vec & ct_count_vec,
                 int scale_func_flag,
                 int robust_flag,
                 const ArrayXb & _b_vec,
                 const Cell_type & _ct):
                 sugar::multiloc <mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type> (popsize, robust_bits, robust_degrees, robust_geno_vec, ct_count_vec, scale_func_flag, robust_flag, _ct.wv), b_vec(_b_vec), ct(_ct)
    {
    }


    // two locus with h
    // initialize from saved state
    twolocus_w_b(const gcont1 & gametes1,
                     const gcont2 & gametes2,
                     const mcont & mutations,
                     const mcount_t & mcounts,
                     const dipvector & haploids,
                     const ArrayXb & _b_vec,
                     const Cell_type & _ct):
                     sugar::multiloc <mutation_type,mcont,gcont1,gcont2,dipvector,mvector,ftvector,lookup_table_type> (gametes1,gametes2,mutations,mcounts,haploids), b_vec(_b_vec), ct(_ct)
     {
     }

  };

  // template<typename mtype,
  //          typename haploid_t = std::pair<std::size_t,std::size_t> >
  // using twolocus_w_b_t = twolocus_w_b<mtype,
  //                                     std::vector<mtype>,
  //                                     std::vector<ct_gamete_robust_bits_t>,
  //                                     std::vector<ct_gamete_wv_t>,
  //                                     std::vector<haploid_t>,
  //                                     std::vector<mtype>,
  //                                     std::vector<uint_t>,
  //                                     std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;

  // two locus type without h
  template<typename mtype,
           typename haploid_t = std::pair<std::size_t,std::size_t> >
  using twolocus_no_h_t = twolocus_w_b<mtype,
	                                  std::vector<mtype>,
	                                  std::vector<ct_gamete_geno_no_h_t>,
                                    std::vector<ct_gamete_wv_t>,
	                                  std::vector<haploid_t>,
	                                  std::vector<mtype>,
	                                  std::vector<uint_t>,
	                                  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>>;
}
#endif
