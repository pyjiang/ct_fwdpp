#ifndef _CALC_MEAN_ROBUST_
#define _CALC_MEAN_ROBUST_

// new file which calculates the mean robustnes of genotype at every generation for the population
#include "ct_singlepop.hpp"
//#include "new_haploids.hpp"

#include <fstream>
namespace CT {


 //  // concrete functions with defined types
 //  // ct_gamete_with_h_t, haploid (just int)
 // template <typename gcont_t>
 // double get_robust_alpha(const gcont_t & gametes, const int & hap)
 // {
 //   return gametes[hap].get_robust_alpha();
 // }

 // for test (because get in consistent result for alpha and robust geno )
 template <typename gcont_t>
 int get_robust_geno(const gcont_t & gametes, const int & hap)
 {
   return gametes[hap].get_robust_geno();
 }

 // // ct_gamete_robust_geno_t, haploid (with Haploid type )
 // template <typename gcont_t>
 // double get_robust_alpha(const gcont_t & gametes, const CT::Haploid<VectorXd> & hap)
 // {
 //   return gametes[hap.index].get_robust_alpha();
 // }


  // // this function applies to ct_gamete_robust_geno_t
  // template <typename gcont_t,
  //           typename haploids_geno_t>
  // double calc_mean_robust(const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr)
  // {
  //   double robust_sum=0;
  //   for (uint_t i = 0 ; i < 2*N_curr ; ++i)
  //   {
  //     robust_sum += get_robust_alpha(gametes,haploids[i]);
  //   }
  //   robust_sum/= (2*N_curr);
  //   return robust_sum;
  // }

  // // added on 1.23.17
  // // this new robustness is by the ratio of the current alpha to initial alpha value
  // template <typename gcont_t,
  //           typename haploids_geno_t>
  // double calc_mean_robust(const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr, double initial_alpha)
  // {
  //   double robust_sum=0;
  //   for (uint_t i = 0 ; i < 2*N_curr ; ++i)
  //   {
  //     robust_sum += get_robust_alpha(gametes,haploids[i]);
  //   }
  //   robust_sum/= (2*N_curr);
  //   robust_sum /= initial_alpha ; // calcaulte the ratio of the current alpha to the initial alpha
  //   return robust_sum;
  // }

  // // added on 1.25.17
  // // calculate mean robust geno ()
  // // use robust geno to calculate(to debug)
  // template <typename gcont_t,
  //           typename haploids_geno_t>
  // double calc_mean_robust_geno(const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr)
  // {
  //   double robust_sum=0;
  //   for (uint_t i = 0 ; i < 2*N_curr ; ++i)
  //   {
  //     robust_sum += get_robust_geno(gametes,haploids[i]);
  //   }
  //   robust_sum/= (2*N_curr);
  //   return robust_sum;
  // }

  // // added on 1.25.17
  // // will output of each gametes's alpha at each generation!!
  // // option1
  // template <typename gcont_t,
  //           typename haploids_geno_t>
  // void output_all_individual_alpha(std::ofstream & out, const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr)
  // {
  //   for (uint_t i = 0 ; i < 2*N_curr ; ++i)
  //   {
  //     out << gametes[haploids[i]].alpha << "\t";
  //   }
  //   out << std::endl;
  // }


  // struct approx_equal
  // {
  //   using result_type = bool;
  //   template<typename T>
  //  inline bool operator()(const T a , const  T b ) const
  //  {
  //    if(abs(a-b) < 1e-6)
  //    {
  //      return true;
  //    }
  //    return false;
  //  }
  // };




  // // option 2
  // // calculate count of each gamete and return the count vector
  // // for the case of non-evolved alpha
  // // for some reason, std::map does not work properly for this countm but std::unordered_map works
  // template <typename gcont_t,
  //           typename haploids_geno_t>
  // std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>  count_alpha_dist(const gcont_t & gametes, const haploids_geno_t & haploids, const uint_t & N_curr)
  // {
  //   // make a dictionary that will record count
  // std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>   alpha_count_dict ;
  //
  //   for (uint_t i = 0 ; i < 2*N_curr ; ++i)
  //   {
  //     alpha_count_dict[gametes[haploids[i]].alpha]++;
  //   }
  //
  //   return alpha_count_dict;
  // }

  // get initial scale count
  template <typename gcont_t,
            typename haploids_geno_t>
  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>  count_scale_dist(const gcont_t & gametes, const haploids_geno_t & haploids, const uint_t & N_curr)
  {
    // make a dictionary that will record count
  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>   scale_count_dict ;
  double scale;

    for (uint_t i = 0 ; i < 2*N_curr ; ++i)
    {
      scale = gametes[haploids[i]].get_robust_alpha_ratio();
      scale_count_dict[scale]++;
    }

    return scale_count_dict;
  }

  // modified on 2.1.17
  // scale equals initial scale  times later scale
  // two locus info stored in haploids
  template <typename gcont_t,
            typename haploids_geno_t>
  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>  count_scale_dist_two_loc(const gcont_t & gametes, const haploids_geno_t & haploids, const uint_t & N_curr)
  {
    // make a dictionary that will record count
  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>   scale_count_dict ;
  double scale;

    for (uint_t i = 0 ; i < 2*N_curr ; ++i)
    {
      scale = gametes[haploids[i].first].get_robust_alpha_ratio();
      scale_count_dict[scale]++;
    }

    return scale_count_dict;
  }



// modified on 2.22.17
// will also pass the total number of allele counts in the population, and check if the total counts in the dictionary is the same as total allele count
// for output alpha dist
// in the order of alpha list
// note: because [] in unordered_map is not a constant function, cannot pass alpha_count_dict as a const parameter
template <typename dict,
          typename alpha_cont_t>
void output_alpha_dist(std::ofstream & out, dict & alpha_count_dict, const alpha_cont_t alpha_list , const int twoN)
{
  int count =0;
  for (auto & alpha : alpha_list)
  {
    count += alpha_count_dict[alpha];
    out << alpha_count_dict[alpha] << "\t";
  }
  out << std::endl;

  assert(count == twoN);
  // loop to the second to last element
  // auto it =alpha_count_dict.begin();
  // out << it-> first << "\t" << it-> second ;
  // ++it;
  //
  // for(; it !=  alpha_count_dict.end(); ++it )
  // {
  //   out << "," << it-> first << "\t" << it-> second ;
  // }
  // out << std::endl;
  // then output the last one
  // out << it-> first << "\t" << it-> second << std::endl;
}



//
// // for calculating and output alpha dist
// template <typename gcont_t,
//           typename haploids_geno_t,
//           typename alpha_cont_t>
// void calc_and_output_alpha_dist(std::ofstream & out, const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr, const alpha_cont_t alpha_list)
// {
//   std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>  alpha_count_dict = count_alpha_dist(gametes, haploids, N_curr);
//   output_alpha_dist (out, alpha_count_dict, alpha_list, 2*N_curr);
// }



template <typename gcont_t,
          typename haploids_geno_t,
          typename alpha_cont_t>
void calc_and_output_scale_dist(std::ofstream & out, const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr, const alpha_cont_t scale_list)
{
  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>  scale_count_dict = count_scale_dist(gametes, haploids, N_curr);
  output_alpha_dist (out, scale_count_dict, scale_list, 2*N_curr);
}

// for two locus case
template <typename gcont_t,
          typename haploids_geno_t,
          typename alpha_cont_t>
void calc_and_output_scale_dist_two_loc(std::ofstream & out, const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr, const alpha_cont_t scale_list)
{
  std::unordered_map<double,int,std::hash<double>,KTfwd::equal_eps>  scale_count_dict = count_scale_dist_two_loc(gametes, haploids, N_curr);
  output_alpha_dist (out, scale_count_dict, scale_list, 2*N_curr);
}



 //  // created on 12.21.16
 // // dealing with two locus (add a dummy variable)
 //  template <typename gcont_t,
 //            typename haploids_geno_t>
 //  double calc_mean_robust(const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr, int dummy)
 //  {
 //    double robust_sum=0;
 //    for (uint_t i = 0 ; i < 2*N_curr ; ++i)
 //    {
 //      robust_sum += get_robust_alpha(gametes,haploids[i].first);
 //    }
 //    robust_sum/= (2*N_curr);
 //    return robust_sum;
 //  }

// create on 1.2.17
// another way to calculate mean robustness
// also give the reference v
template <typename gcont_t,
          typename haploids_geno_t,
          typename mcont_t>
void check_robust_bits(const gcont_t & gametes,const haploids_geno_t & haploids, const uint_t & N_curr, VectorXi v, const mcont_t & mutations, int robust_bits, int generation)
{
  // first get the number of 1s in ref v
  int count = v.sum();

  int index, changed_bits;
  // get smutations in each robust bits
  for (uint_t i = 0 ; i < 2*N_curr ; ++i)
  {
    index = haploids[i].first;
    changed_bits= muts_to_robust(gametes[index].smutations.cbegin(),gametes[index].smutations.cend(),mutations,v,robust_bits,1);
    if (count + changed_bits != gametes[index].robust_geno)
      std::cout << generation << " " << count + changed_bits << " " << gametes[index].robust_geno << std::endl;
  }

}


  // //
  // // added on 12.8.16
  // // added calculation of standard deviation of h in every individual over time (to track the change in variance)
  // // the third parameter var is actually the variance of the h term without a parameter
  // template <typename haploid_t>
  // double calc_mean_h_sd_relative(const haploid_t & haploids, uint_t twoN, double initial_sd)
  // {
  //   double h_mean_sd= 0;
  //   for (uint_t i=0; i< twoN;i++ )
  //   {
  //     h_mean_sd+= calc_h_sd(haploids[i].h_vec);
  //   }
  //   h_mean_sd /= twoN;
  //   return h_mean_sd/initial_sd;
  // }
  //



}
#endif
