#ifndef _CT_GAMETE_
#define _CT_GAMETE_
#include <fwdpp/forward_types.hpp>
#include <cmath>

#include "ct_eigen.hpp" // change to hpp


namespace CT
{
  using namespace KTfwd;

  // basic gamete type
  template<typename TAG>
  struct ct_gamete: public gamete_base<TAG>
  {
    int flag; // whether it needs to recalculate fitness or not (if flag==0, do not need to recalculate)
    double fitness; // gamete will store fitness (if it does not need to re-calculate)


    ct_gamete(const uint_t & _icount, const int & _flag=1) noexcept //initizialize with gamete count and  flag
      : gamete_base<TAG>(_icount),flag(_flag)
    {
    }

    // template <typename fit_func>
    // ct_gamete(const uint_t & _icount, const fit_func & ff) noexcept //initizialize without flag (only use this Constructor to initialize the first time, need to initialize fitness!)
    //   : gamete_base<TAG>(_icount)
    // {
    //   fitness =  ff();
    //   flag=0;
    // }

    // if specify flag, will set as flag
    ct_gamete(const uint_t & _icount, const int & _flag, const double & _fitness) noexcept //initizialize with flag and fitness (from saved state)
      : gamete_base<TAG>(_icount),flag(_flag),fitness(_fitness)
    {
      // flag=1; // allow the saved gametes to recalculate fitness
    }

    // maybe obselete
    // // if not specify flag, set flag to 1
    // ct_gamete(const uint_t & _icount, const double & _fitness) noexcept //initizialize with flag and fitness (from saved state)
    //   : gamete_base<TAG>(_icount),fitness(_fitness)
    // {
    //   flag=1; // allow the saved gametes to recalculate fitness
    // }

    template <typename T>
    ct_gamete(const uint_t & _icount, T && n, T && s) noexcept // initialize with global vector neutral and selected
      : gamete_base<TAG> (_icount,n,s)
    {
      flag=1;
    }

    ct_gamete(const uint_t & _icount, const ct_gamete & g ) noexcept // initialize a new gamete from a previous gamete(will pass the smutation list )
      : gamete_base<TAG> (_icount,g.mutations,g.smutations), fitness(g.fitness)
    {

    }


    //write some get and set function for fitness and flag
    int get_flag() const
    {
      return flag;
    }

    void set_flag(int new_flag)
    {
      flag= new_flag;
    }

    double get_fitness() const
    {
      return fitness;
    }

    void set_fitness(double new_fit)
    {
      fitness= new_fit;
    }

    // 11.2.16
    // create copy functions
    // only copy mutations
    void copy(const ct_gamete & g)
    {
      this->mutations=g.mutations;
      this->smutations=g.smutations;
      fitness = g.fitness; // added on 11.22.16 (need to initialize fitness to record old fitness value of gamete)
    }

  };

  // created on 1.23.17
  // create new robust genotype, but will scale the index in wij instead of changing h
  template<typename TAG>
  struct ct_gamete_robust_geno_only_wij: public ct_gamete<TAG>
  {
    int robust_geno; // stores the robust genotypes in the current gamete
    double scale ; // this scale measures the degree of robustness change by mutations in robust geno

    // do not use static
    int scale_func_flag; // indicates which scale function to use


    // // function pointer of calc_scale
    // static  void (ct_gamete_robust_geno_only_wij<TAG>::*update_scale_func)(int);
    // define function pointer (do not use static)
    void (ct_gamete_robust_geno_only_wij<TAG>:: *update_scale_func)(int);
    // static FPTR update_scale_func;


    ct_gamete_robust_geno_only_wij(const uint_t & _icount, const uint_t & _robust_geno) noexcept
      : ct_gamete<TAG>(_icount),robust_geno(_robust_geno)
    {
      scale =1;
      //this -> flag =1; // default flag is 1, when calling ct_gamete<TAG>(_icount)
    }

    // modified on 5.23.17
    // change the scale to the first parameter, just to distinguish from the other constructor which initilize with count and fitness
    // modified on 5.16.17
    // do not pass robust_geno, just scale
    // used to initialize gamete with give scale (when not need to evolve robustness)
    ct_gamete_robust_geno_only_wij(const double & _scale, const uint_t & _icount) noexcept
      : ct_gamete<TAG>(_icount)
    {

      scale = _scale; // modified on 4.24.17
      // this -> flag =1;
      robust_geno = 0; // just put in as a dummy variable
    }



    // modifed on 8.8.17
    // add parameter of robust_bits
    void initialize_calc_function_ptr_and_rob(const uint_t &  _robust_geno, const int robust_bits, const int robust_degrees)
    {
      // flag=8 will not use bits additively but one genotype represent one alpha value for robustness
      // so the way to initialize mutations are different
      if (scale_func_flag ==8 || scale_func_flag ==16 ) //added on 8.8.17
     {
        if (scale_func_flag ==8 )
           update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_8_degree;
        else
           update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_16_degree;
         initialize_robustness_v2(_robust_geno, robust_bits, robust_degrees);
     }
     else{

       if (scale_func_flag ==0 ) // modified on 2.26.17, if flag=0, use 2 based as well!!
       {
         update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_additive_2base;
       }
       else if(scale_func_flag==7) // it is the additive case
       {
          update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_8_degree;
       }
       else if (scale_func_flag ==15)
       {
          update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_16_degree;
       }
       initialize_robustness(_robust_geno, robust_degrees);
     }
   }


     void initialize_calc_function_ptr()
     {
       // flag=8 will not use bits additively but one genotype represent one alpha value for robustness
       // so the way to initialize mutations are different
       if (scale_func_flag ==8  ) //added on 8.8.17
      {
          update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_8_degree;
      }
      else if (scale_func_flag ==16)
      {
         update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_16_degree;
      }
      else
      {
        if (scale_func_flag ==0 ) // modified on 2.26.17, if flag=0, use 2 based as well!!
        {
          update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_additive_2base;
        }
        else if(scale_func_flag==7) // it is the additive case
        {
           update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_8_degree;
        }
        else if(scale_func_flag ==15)
        {
          update_scale_func = &ct_gamete_robust_geno_only_wij<TAG>::update_scale_16_degree;
        }
      }

    }


    // // created on 1.30.17
    // // we can also initialize gamete initial_alpha to the given gamma value
    // ct_gamete_robust_geno_only_wij(const uint_t & _icount, const int & _robust_geno, const double gamma, double scale=1) noexcept
    //   : ct_gamete<TAG>(_icount),robust_geno(_robust_geno)
    // {
    //   // when initialize, the origianl alpha is set to prev_alpha
    //   initialize_alpha(gamma, scale);
    //   // initial_alpha = alpha ;
    //   this -> flag =1;
    // }



    // initizialize with a previous gamete
    ct_gamete_robust_geno_only_wij(const uint_t & _icount, const ct_gamete_robust_geno_only_wij & g ) noexcept // initialize a new gamete from a previous gamete(will pass the smutation list )
      : ct_gamete<TAG> (_icount,g.mutations,g.smutations), robust_geno(g.robust_geno),  scale (g.scale), scale_func_flag(g.scale_func_flag)
    {
      initialize_calc_function_ptr();
      // this -> flag=1; // already specified in the base class
    }

    // initialize with only count, flag, fitness (called when read in from file)
    ct_gamete_robust_geno_only_wij(const uint_t & _icount, const int flag, const double fitness ) noexcept
      : ct_gamete<TAG> (_icount,flag,fitness)
    {
    }

    // initialize with count and fitness (just to be compatible with read in function)
    ct_gamete_robust_geno_only_wij(const uint_t & _icount, double fitness ) noexcept
      : ct_gamete<TAG> (_icount, fitness)
    {
    }

    // new
    // initialize with count and fitness (just to be compatible with read in function)
    ct_gamete_robust_geno_only_wij(const uint_t & _icount ) noexcept
      : ct_gamete<TAG> (_icount)
    {
    }

    // this is needed when reading in gametes and initialize with count and flag
    ct_gamete_robust_geno_only_wij(const uint_t & _icount, const int flag) noexcept
      : ct_gamete<TAG> (_icount,flag)
    {
    }

    // note: because the mutation bits have been changed already
    // copy function (when a mutation occurs on a pre-existing gamete)
    template <typename ct_gamete_t>
    void copy(const ct_gamete_t & g)
    {
      ct_gamete<TAG>::copy(g);
      robust_geno= g.robust_geno;
      scale = g.scale ;
      scale_func_flag = g.scale_func_flag;
      initialize_calc_function_ptr();
      this -> flag=1;
    }


    void initialize_scale(double _scale)
    {
      scale = _scale ;
    }


    // added on 2.21.17
    // still additive, but use 2-based, instead of exp based
    void update_scale_additive_2base(int robust_degrees)
    {
        scale = calc_scale_2base (robust_degrees, robust_geno);

    }

    // added on 8.8.17
    void update_scale_8_degree(int robust_degrees)
    {
        scale = calc_scale_8_degree(robust_degrees, robust_geno);
    }

    void update_scale_16_degree(int robust_degrees)
    {
        scale = calc_scale_16_degree(robust_degrees, robust_geno);
    }


    // also used to initialize scale and robust_geno from initial robust_geno
    void update_robust_geno(const int & new_robust_geno, int robust_degrees)
    {
      robust_geno = new_robust_geno;

      // use function pointer
      (this->*update_scale_func)(robust_degrees);
    }

    // added on 2.6.17
    // for the gametes that haven't initialized the scale functions
    void initialize_robustness_all(const int _scale_func_flag, const int init_robust_geno, const int robust_bits, const int robust_degrees)
    {
      scale_func_flag = _scale_func_flag;
      initialize_calc_function_ptr_and_rob(init_robust_geno, robust_bits, robust_degrees);
    }


    void initialize_robustness(const int init_robust_geno, int robust_degrees)
    {
      update_robust_geno(init_robust_geno, robust_degrees);
      initialize_mutations();
    }

    // for flag=8
    void initialize_robustness_v2(const int init_robust_geno, const int robust_bits, const int robust_degrees)
    {
      update_robust_geno(init_robust_geno, robust_degrees);
      initialize_mutations_v2(robust_bits);
    }


    // note: this initialize mutations only works for additive bits case!!
    // when initilize with robust_geno, should also change smutations in gamete as well
    void initialize_mutations()
    {
      for(uint i=0; i< robust_geno; i++)
      {
        this->smutations.push_back(i);
      }
    }

    // added on 8.8.17
    // initialize mutation where each genotype value represent one scale, so each mutation will change one bit in genotype
    // basically converts integer to an array of bits
    void initialize_mutations_v2(int robust_bits)
    {
      int mask = 1;
      int val;
      for (uint i=0; i< robust_bits; i++)
      {
        val =  robust_geno & mask ;
        if (val)
        {
          this->smutations.push_back(i);
        }
        mask <<= 1;
      }
    }



    void update_robust_geno_flag(const int & changed_bits, const int & robust_degrees, int robust_flag)
    {
      if(changed_bits)
      {
        if (robust_flag==0)
          update_robust_geno(robust_geno ^ changed_bits, robust_degrees);
        else if(robust_flag==1)
          update_robust_geno(robust_geno + changed_bits, robust_degrees);
      }
    }

    void update_robust_geno_flag_and_set(const int & changed_bits, const int & robust_degrees, int flag)
    {
      update_robust_geno_flag(changed_bits,robust_degrees,flag);
      this-> flag=1; // should set flag=1 (still needs to recalculate fitness, just everything else is updated )
    }


    double get_robust_geno() const
    {
      return robust_geno;
    }

    double get_robust_alpha_ratio() const
    {
      // return alpha/initial_alpha ;
      // return initial_scale* scale ;
      return scale;
    }


  };



  // gametes for basic type simulations
  template<typename TAG,
           typename Vec>
  struct ct_gamete_wv:  public ct_gamete<TAG>
  {
    Vec wv; // added wv to this class

    // constructor for basic gamete type
    ct_gamete_wv(const uint_t & _icount, double fitness, const Vec & _wv ) noexcept
      : ct_gamete<TAG>(_icount,fitness), wv(_wv)
    {
    }

    // constructor ?when did it called
    ct_gamete_wv(const uint_t & _icount, const Vec & _wv) noexcept //initizialize with flag
      : ct_gamete<TAG>(_icount), wv(_wv)
    {
    }

    ct_gamete_wv(const uint_t & _icount, const ct_gamete_wv & g) noexcept
      : ct_gamete<TAG> (_icount,g.mutations,g.smutations), wv (g.wv)
    {
    }

    ct_gamete_wv(const uint_t & _icount, double fitness ) noexcept //initizialize with flag
      : ct_gamete<TAG>(_icount,fitness)
    {
    }



    // constructor for reading from file
    ct_gamete_wv(const uint_t & _icount, const int flag, const double fitness  ) noexcept //initizialize with flag
      : ct_gamete<TAG>(_icount,flag, fitness)
    {
    }

    ct_gamete_wv(const uint_t & _icount, const int flag ) noexcept //initizialize with flag
      : ct_gamete<TAG>(_icount,flag)
    {
    }

    ct_gamete_wv(const uint_t & _icount ) noexcept //initizialize with flag
      : ct_gamete<TAG>(_icount)
    {
    }

    void update_wv(const Vec & diff_vec)
    {
      wv += diff_vec;
    }

    // added on 1.24.17
    void update_wv(double scale)
    {
      wv = update_wv_scale(scale, wv);
    }

    // first add the diff_vec (for the current mutations )
    // and then rescale it if robustness changes
    void update_wv(const Vec & diff_vec, double scale)
    {
      update_wv(diff_vec);
      if (scale!=1)
        update_wv(scale);
    }


    void update_all(const Vec & diff_vec)
    {
      update_wv(diff_vec);
    }

    void copy(const ct_gamete_wv & g)
    {
      ct_gamete<TAG>::copy(g);
      wv=g.wv;
    }


  };




}
#endif
