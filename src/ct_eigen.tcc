#ifndef _CT_EIGEN_TCC_
#define _CT_EIGEN_TCC_

#include "ct_eigen.hpp"
// #include "eigen_IO.hpp"
// #include "evol_alpha_range.hpp"

//extern double LOW_SCALE_BOUND;
//extern double HIGH_SCALE_BOUND;

// calculate ref_val
VectorXd calc_cell_type_z_val(const MatrixXd & w, const VectorXi & v, const VectorXd & h)
{
  VectorXd z_valprev=calc_cell_type_wv(w,v);
  VectorXd z_val = calc_cell_type_z_val_from_wv(z_valprev,h);
  return z_val ;
  // ArrayXb res = (z_val > 0 ) ;
  // return res ;
}

// calculating intermediate values of w times v
VectorXd calc_cell_type_wv(const MatrixXd & w, const VectorXi & v)
{
  VectorXd wv =  w* v.cast <double> ();
  return wv;
}

VectorXd calc_cell_type_z_val_from_wv(const VectorXd & wv, const VectorXd & h)
{
  VectorXd z_val = wv -h;
  return z_val;
}

// created on 11.16.16
double calc_fitness_from_wv_h(const VectorXd & wv, const VectorXd & h, const ArrayXb & b_vec)
{
  VectorXd z_val= calc_cell_type_z_val_from_wv(wv,h);
  ArrayXb z = (z_val.array() > 0 ) ;
  return eval_fitness(z,b_vec);
}


ArrayXb calc_cell_type(const MatrixXd & w, const VectorXi & v, const VectorXd & h)
{
  VectorXd z_val= calc_cell_type_z_val(w,v,h);
  ArrayXb res = (z_val.array() > 0 ) ;
  return res;
}




// // added on 1.12.17
// // do not pass gsl_rng, but initiate it inside the function
// // this function is to simulate random cell types repeatedly until a desired z is generated
// Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, double c, double pv1, double ga)
// {
//   // first initialize r
//   // ref: https://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generator-Examples.html#Random-Number-Generator-Examples
//   gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
//
//   // then call the previous function
//   Cell_type ct= calc_cell_type_from_para_indeph(nz,Lv,c,pv1,ga,r);
//   gsl_rng_free (r);
//   return ct;
// }


// created on 1.12.17
// pass GSLrng type in the argument
Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, double c, double pv1, double ga, GSLrng & r)
{
  return calc_cell_type_from_para_indeph(nz,Lv,c,pv1,ga,r.get());
}




// modified on 12.7.16
// no longer use rvz as parameter (used to keep it to be compatible with previous version)
// calculating cell types, return a vector of z
Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, double c, double pv1, double ga, const gsl_rng * r)
{
  // uint32_t Lv=int(nz*rvz); // length of v vector
  VectorXi v(Lv);
  // first generate random numbers
  double v_prob_vec[]={1-pv1,pv1};
  gsl_ran_discrete_t * rd_v_pt= gsl_ran_discrete_preproc(2,v_prob_vec);
  // generate random vector for v
  for (uint32_t i=0; i<Lv; i++ )
  {
    v(i)=gsl_ran_discrete(r,rd_v_pt);
  }
  // generate random number for w matrix
  MatrixXd w(nz,Lv);
  double  range_w[]={ga, -ga, 0};
  double w_prob_vec[]={c/2, c/2, 1-c};
  gsl_ran_discrete_t * rd_w_pt= gsl_ran_discrete_preproc(3,w_prob_vec);
  for (uint32_t i=0; i<nz; i++)
  {
    for (uint32_t j=0; j<Lv; j++)
    {
      w(i,j)= range_w[gsl_ran_discrete(r,rd_w_pt)];
    }
  }
  // initialize h vector from a normal distribution
  VectorXd h(nz);
  double sigma=sqrt(Lv*pv1*c);
  for (uint32_t i=0; i<nz; i++)
  {
    h(i)= gsl_ran_gaussian(r,sigma);
  }

  // then calculate z
  // ArrayXb z= calc_cell_type(w,v,h);
  // ArrayXd z_val= calc_cell_type_z_val(w,v,h);

  // will get wv and z_val
  VectorXd wv = calc_cell_type_wv(w,v);
  VectorXd z_val = calc_cell_type_z_val_from_wv(wv,h);

  ArrayXb z = (z_val.array() > 0 ) ;

  // modifed on 1.23.17
  // add gamma
  Cell_type ct(w,v,h,z_val,wv,z, ga); // change wv to array
  return ct;
  // return z for the moment
  //return z;
}




// modified on 2.1.17
// use gamma to initialize w, and the standard way to initialize h (as default)
// use robust_bits only to initilize in v (to leave space for robust bits)
// remove robust_degrees in parameter
// modifed on 12.27.16
// modified on 12.6.16
// added a flag, indicating which way of robustness to code (flag=0, original method, 5 bits, code as 32 degrees; flag=1, using long bits: eg. 1000, code additively)
// modified on 11.28.16
// created on 11.21.16
// specify robust value for robust bits
Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, uint32_t robust_bits, int flag,  double c, double pv1, double ga, const gsl_rng * r, uint32_t robust_val)
{
  VectorXi v=VectorXi::Zero(Lv); // initialize to be 0, fixed a bug when exporting to python
  // first generate random numbers
  double v_prob_vec[]={1-pv1,pv1};
  gsl_ran_discrete_t * rd_v_pt= gsl_ran_discrete_preproc(2,v_prob_vec);
  // generate random vector for v
  // convert robust_val to the first few bits in v and then generate the rest based on random numbers
  if (flag==0)
    set_robust_geno(v,robust_val,robust_bits);
  else
  {
    assert(robust_val<robust_bits);
    set_robust_geno_v2(v,robust_val);
  }

  // to generate random numbers of robust_bits digit, but disgard them (to be the same as the other version for test )
  for (uint32_t i=0; i< robust_bits; i++)
  {
    gsl_ran_discrete(r,rd_v_pt);
  }


  for (uint32_t i=robust_bits; i<Lv; i++ )
  {
    v(i)=gsl_ran_discrete(r,rd_v_pt);
  }

  // generate random number for w matrix
  uint32_t func_geno_len= Lv- robust_bits;

  MatrixXd w(nz,func_geno_len);
  double  range_w[]={ga, -ga, 0};
  double w_prob_vec[]={c/2, c/2, 1-c};
  gsl_ran_discrete_t * rd_w_pt= gsl_ran_discrete_preproc(3,w_prob_vec);
  for (uint32_t i=0; i<nz; i++)
  {
    for (uint32_t j=0; j<func_geno_len; j++)
    {
      w(i,j)= range_w[gsl_ran_discrete(r,rd_w_pt)];
    }
  }

  // initialize h vector from a normal distribution
  VectorXd h(nz);
  double var= func_geno_len*pv1*c;
  resample_h_from_alpha(h,1,var,r);
  // resample_h(h,robust_val,robust_degrees, var,r); //original code
  // generate_rand_h(h, sqrt(var),r);// for test

  // then calculate z
  // ArrayXb z= calc_cell_type(w,v,h);
  // ArrayXd z_val= calc_cell_type_z_val(w,v,h);

  // will get wv and z_val
  VectorXd wv = calc_cell_type_wv(w,v,func_geno_len);
  VectorXd z_val = calc_cell_type_z_val_from_wv(wv,h);

  ArrayXb z = (z_val.array() > 0 ) ;

  // pass gamma
  Cell_type ct(w,v,h,z_val,wv,z,robust_val,ga); // change wv to array
  return ct;
  // return z for the moment
}


// created on 1.23.17
// cell type of scaled gamma
Cell_type scaled_ct(const Cell_type & ct, double scale_factor)
{
  // scale w matrix
  MatrixXd new_w = ct.w*scale_factor;
  // calculate new h
  VectorXd new_h  = calc_cell_type_wv( new_w - ct.w, ct.v ) + ct.h;
  // calculate new wv
  VectorXd new_wv =  calc_cell_type_wv(new_w, ct.v) ;
  Cell_type new_ct (new_w, ct.v, new_h, ct.z_val, new_wv, ct.z , ct.gamma * scale_factor );
  return new_ct;
}





// calculating intermediate values of w times v
// func_geno_len is when robustness is encoded in one big locus as the functional gneotype. the func_geno_len is the genotype that directly regulate phenotype
VectorXd calc_cell_type_wv(const MatrixXd & w, const VectorXi & v, uint32_t func_geno_len)
{
  VectorXd wv =  w* v.tail(func_geno_len).cast <double> ();
  return wv;
}

// fixed on 11.27.16
// where v(0)-v(5) represent robust_geno, v(0) is the highest digit
// get robust_geno from cell type
int get_robust_geno(const VectorXi & v, uint32_t robust_bits)
{
  int robust_geno=0;
  for (uint32_t i=0; i< robust_bits; i++)
  {
    robust_geno += v(i) << (robust_bits - i-1);
  }
  return robust_geno;
}

// created on 12.6.16
// code it additively
int get_robust_geno_v2(const VectorXi & v, uint32_t robust_bits)
{
  int robust_geno=0;
  for (uint32_t i=0; i< robust_bits; i++)
  {
    robust_geno += v(i);
  }
  return robust_geno;
}

// convert robust_val to robust geno in v
void set_robust_geno(VectorXi & v, uint32_t robust_val, uint32_t robust_bits)
{
  uint32_t mask = 1 << (robust_bits -1) ;
  for(uint32_t i=0; i< robust_bits; i++)
  {
    v(i)= (robust_val & mask) >> (robust_bits-i-1) ;
    mask >>= 1;
  }
}

// created on 12.6.16
// additive version (will just set the first robust_val bits to 1)
void set_robust_geno_v2(VectorXi & v, uint32_t robust_val)
{
  for(uint32_t i=0; i< robust_val; i++)
  {
    v(i)=1;
  }
}



// created on 1.23.17
// add parameters of left_bound and right_bound
double calc_alpha(uint32_t robust_degrees, uint32_t robust_geno, double left_bound, double right_bound)
{
  double range= right_bound- left_bound;
  double alpha= exp(left_bound + range/(robust_degrees-1)* robust_geno); // alpha ranges from 0.05 to 20, and the robust_geno separate it linearly on log scale
  return alpha;
}



// calculate scale factor
// make sure this scale is symmetric along 1
double calc_scale(uint32_t robust_degrees, uint32_t robust_geno )
{
  double scale = exp( LOW_SCALE_BOUND + ( HIGH_SCALE_BOUND - LOW_SCALE_BOUND) / (robust_degrees-1) * robust_geno);
  return scale ;
}

// added on 2.21.17
// new way to scale, using two-based instead of exp-based, also high bound and low bound should be integer
double calc_scale_2base(uint32_t robust_degrees, uint32_t robust_geno )
{
  double digit = LOW_SCALE_BOUND + ( HIGH_SCALE_BOUND - LOW_SCALE_BOUND ) / (robust_degrees-1) * robust_geno ;
  double scale = pow(2,digit);
  return scale ;
}


// added on 8.8.17
double calc_scale_8_degree(uint32_t robust_degrees, uint32_t robust_geno)
{
    double scale;
    assert(robust_degrees ==8);
    assert(robust_geno <=7 );
    switch(robust_geno)
    {
        case 0: scale = 0.5;
                break;
        case 1: scale = 1;
                break;
        case 2: scale = 2;
                break;
        case 3: scale = 2.5;
                break;
        case 4: scale = 3;
                break;
        case 5: scale = 3.5;
                break;
        case 6: scale = 4;
                break;
        case 7: scale = 8;
                break;
    }
    return scale;
}

// added on 8.8.17
double calc_scale_16_degree(uint32_t robust_degrees, uint32_t robust_geno)
{
    double scale;
    assert(robust_degrees == 16);
    assert( robust_geno <=15 );
    switch(robust_geno)
    {
        case 0: scale = 0.25;
                break;
        case 1: scale = 0.5;
                break;
        case 2: scale = 0.75;
                break;
        case 3: scale = 1;
                break;
        case 4: scale = 2;
                break;
        case 5: scale = 2.5;
                break;
        case 6: scale = 3;
                break;
        case 7: scale = 3.5;
                break;
        case 8: scale = 4;
                break;
        case 9: scale = 6;
                break;
        case 10: scale = 8;
                break;
        case 11: scale = 16;
                break;
        case 12: scale = 32;
                break;
        case 13: scale = 64;
                break;
        case 14: scale = 128;
                break;
        case 15: scale = 256;
                break;
    }
    return scale;
}





// given alpha as a parameter
void resample_h_from_alpha(VectorXd & h,double alpha,double var, const gsl_rng * r)
{
  double sigma= sqrt(var)*alpha;
  generate_rand_h(h, sigma,r);
}

// version 2
void resample_h_from_alpha(VectorXd & h,double alpha,double var,const double * norm_rand_tb, const gsl_rng * r)
{
  double sigma= sqrt(var)*alpha;
  generate_rand_h(h, sigma,norm_rand_tb,r);
}


// version 2, generate random numbers from table
// sigma is the standard deviation
void generate_rand_h(VectorXd & h, double sigma,const double * norm_rand_tb,const gsl_rng * r)
{
  uint32_t nz= h.size();
  for (uint32_t i=0; i<nz; i++)
  {
    h(i)= norm_rand_tb[gsl_rng_uniform_int(r,MAXSIZE)]*sigma ;
  }
}

void generate_rand_h(VectorXd & h, double sigma,const gsl_rng * r)
{
  uint32_t nz= h.size();
  for (uint32_t i=0; i<nz; i++)
  {
    h(i)= gsl_ran_gaussian(r,sigma);
  }
}

double eval_fitness(const ArrayXb & z, const ArrayXb & b_vec)
{
  double fitness = (z==b_vec).count()/ double(z.size());
  return fitness;
}


// added on 12.27.16
// calculate fitness based on w,v,and b (with length of robust_bits in v )
double calc_fitness(const MatrixXd & w, const VectorXi & v, uint32_t func_geno_len, const VectorXd & h, const ArrayXb & b_vec)
{
  VectorXd wv = calc_cell_type_wv(w,v,func_geno_len);
  return calc_fitness_from_wv_h(wv,h,b_vec);
}


// modified on 1.30.17
// add ratio
// modified on 12.20.16
// add relative_index, indicating the beginning of the range of the index [for multi locus]
// define the detailed calculation
// need to loop through all the columns that in the range and then times the value
// return the result of matrix multiplication of  submatrix  and subvector
// indicies are the positions where a mutation occur
// use iterator to loop through the mutations (given iterator begin and end )
template<typename iterator_t,
         typename mcont_t>
VectorXd recalc_mat(const MatrixXd & w, const VectorXi & v, iterator_t begin, iterator_t end, const mcont_t & mutations, double ratio, int relative_index=0)
{
  VectorXd diff_vec = VectorXd::Zero(w.rows());
  for(iterator_t i=begin;i!=end;i++)
  {
    int index= std::round( (mutations[*i].pos-relative_index) * v.size()) ; // fixed a bug on 10.17.16, use std::round instead of int!!
    // std::cout << index <<std::endl;
   //  VectorXd wcol= w.col(index);
    diff_vec.noalias() += ratio * w.col(index)*(1-2*v(index)); // to avoid temporary value
  }
  return diff_vec;
}


// moved recalc_mat from ct_fitness_policy.cc

// modified on 12.6.16
// added flag for robust bits calculation
// recalc_mat with robust_bits
// in this function, change the type of robust_bits to int (for calculation), but elsewhere, robust_bits are saved as uint32
// note: v is the v of g (previous to this round of mutation), not common v
template<typename iterator_t,
         typename mcont_t>
VectorXd recalc_mat(const MatrixXd & w, const VectorXi & v, const int & robust_bits, const iterator_t begin, const iterator_t end, const mcont_t  & mutations,  int & changed_bits, int flag)
{
  VectorXd diff_vec = VectorXd::Zero(w.rows());
  uint32_t func_geno_len= v.size() - robust_bits;
 //  int changed_bits =0; // store the value

  int index;
  for(iterator_t i=begin;i!=end;i++)
  {
    int ori_index= std::round( mutations[*i].pos * v.size()) ;// fixed a bug on 10.17.16, use std::round instead of int!!
    if((ori_index-robust_bits)>=0 )
    {
      index = ori_index-robust_bits;
      // std::cout << index <<std::endl;
     //  VectorXd wcol= w.col(index);
      diff_vec.noalias() += w.col(index)*(1-2*v.tail(func_geno_len)(index));
    }
    else // the change is in the robust_bits
    {
      index = ori_index;
      if(flag==0)
        //changed_bits += 1<< (robust_bits-index-1);
        changed_bits += 1 << index;
      else if(flag==1)
        changed_bits +=  pow(-1,v(index)); // if v==0, add 1; if v==1, minus 1
    }
  }

  return diff_vec;
}

// created on 1.23.17
// new function of recalc_mat, which scales the mutation with changed robustness in coef
template<typename iterator_t,
         typename mcont_t>
VectorXd recalc_mat(const MatrixXd & w, const VectorXi & v, const int & robust_bits, const iterator_t begin, const iterator_t end, const mcont_t  & mutations,  double ratio, int & changed_bits, int flag)
{
  VectorXd diff_vec = VectorXd::Zero(w.rows());
  uint32_t func_geno_len= v.size() - robust_bits;
 //  int changed_bits =0; // store the value

  int index;
  for(iterator_t i=begin;i!=end;i++)
  {
    int ori_index= std::round( mutations[*i].pos * v.size()) ;// fixed a bug on 10.17.16, use std::round instead of int!!
    if((ori_index-robust_bits)>=0 )
    {
      index = ori_index-robust_bits;
      // std::cout << index <<std::endl;
     //  VectorXd wcol= w.col(index);
      diff_vec.noalias() += ratio * w.col(index)*(1-2*v.tail(func_geno_len)(index));
    }
    else // the change is in the robust_bits
    {
      index = ori_index;
      if(flag==0)
        //changed_bits += 1<< (robust_bits-index-1);//
        changed_bits +=  1<< index;
      else if(flag==1)
        changed_bits +=  pow(-1,v(index)); // if v==0, add 1; if v==1, minus 1
    }
  }

  return diff_vec;
}




// added on 12.17.16
// convert the list of mutations in robust bits  to changed bits
template<typename iterator_t,
         typename mcont_t>
int muts_to_robust(const iterator_t begin, const iterator_t end, const mcont_t  & mutations, const VectorXi & v,  int robust_bits, int flag)
{
  int changed_bits=0;
  if(flag==0)
  {
    for(iterator_t i=begin;i!=end;i++)
    {
      int ori_index= std::round( mutations[*i].pos * robust_bits) ;
      //changed_bits += 1<< (robust_bits-ori_index-1); //old, change on 8.8.17
      changed_bits += 1<< ori_index;
    }
  }
  else if(flag==1)
  {
    for(iterator_t i=begin;i!=end;i++)
    {
      int ori_index= std::round( mutations[*i].pos * robust_bits) ;
      changed_bits +=  pow(-1,v(ori_index)); // if v==0, add 1; if v==1, minus 1
    }
  }
  return changed_bits;
}


// // calculate standard deviation of h vector
// double calc_h_sd(const VectorXd & h)
// {
//   double h_mean = h.mean();
//   double var=0;
//   uint32_t len = h.size();
//   for(uint32_t i=0; i< len ; i++)
//   {
//     var+= pow(h(i)-h_mean,2);
//   }
//   return sqrt(var/(len-1));
// }


// modifed on 2.15.17
// add an extra parameter diffk (when type=2), indicates how many different bits are from initial b to initial z
// created on 12.19.16
// generating initial b vector (now write in function to be called easily)
void generate_b_vec(ArrayXb & b_vec, const ArrayXb & z, int type, const gsl_rng * r, int diffk ) //diffk set to 0 default 
{
  if (type==0) //type=0, stabilizing selection
  {
    b_vec = z;
  }
  else if (type ==1) // random init
  {
    // random b vector
    double b_prob_vec[]={0.5,0.5};
    gsl_ran_discrete_t * rd_b_pt= gsl_ran_discrete_preproc(2,b_prob_vec); // equal probability of 0 or 1
    for (uint i=0; i< b_vec.size(); i ++)
    {
      b_vec(i)= gsl_ran_discrete(r,rd_b_pt);
    }
  }
  else if (type==2)
  {
    b_vec = z;
    generate_b_dis(b_vec, z, diffk, r);
  }
}

// modifed on 4.18.17
// add another method to initialize b vector, by randomly generating a genotype, and given the interaction of the current cell type (will return the gneotype)
VectorXi generate_b_vec(ArrayXb & b_vec, double pv1, const MatrixXd & w, const VectorXd & h,  const gsl_rng * r )
{
  // generate a random v vector
  uint32_t Lv = w.cols();
  VectorXi rand_v (Lv);
  double v_prob_vec[]={1-pv1,pv1};
  gsl_ran_discrete_t * rd_v_pt= gsl_ran_discrete_preproc(2,v_prob_vec);
  // generate random vector for v
  for (uint32_t i=0; i<Lv; i++ )
  {
    rand_v(i)=gsl_ran_discrete(r,rd_v_pt);
  }

  // then calcualte z from rand_v, w, h, and assign to b
  b_vec = calc_cell_type(w,rand_v,h);
  return rand_v;
}


// created on 1.5.17
// wrote function that generate b with cerain distance of z
void generate_b_dis(ArrayXb & b_vec, const ArrayXb & z, int dis, const gsl_rng * r)
{
  uint nz = z.size();
  int * z_idx = new int [nz];
  int * rand_idx = new int [dis];

  for (uint i=0; i<nz; i++)
    z_idx[i]= i;

  gsl_ran_choose(r, rand_idx, dis, z_idx, nz, sizeof(int)); // random choose index in vector to change

  for (int i=0; i< dis; i++)
  {
    b_vec(rand_idx[i])= 1 - b_vec(rand_idx[i]);
  }
}

// first v is Destination, copy from v1 and v2
void copy_back(VectorXi & v, const VectorXi & v1, const VectorXi & v2)
{
  uint len1=v1.size();
  uint len2=v2.size();
  // for (uint i=0; i <len1; i++)
  //   v(i)= v1(i);
  // for (uint i=len1; i< (len2+len1); i++)
  //   v(i)= v2(i-len1);
  v.head(len1)=v1; // this function works as left parameter as well!! (works the same as looping assignment)
  v.tail(len2)=v2;

}




VectorXd update_wv_scale(double scale, const VectorXd & wv)
{
  return scale*wv ;
}



#endif
