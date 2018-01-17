#ifndef _CT_EIGEN_
#define _CT_EIGEN_

// using Eigen library to calculate cell type and see speed
#include <iostream>
#include <Eigen/Dense>
#include <gsl/gsl_randist.h>
#include <stdint.h>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>
#include <functional>

#include "gsl_funcs.hpp" // needs definition of GSLrng


double LOW_SCALE_BOUND = -1;
double HIGH_SCALE_BOUND = 1;

#define MAXSIZE 1000000

// using Eigen::MatrixXd;
// using Eigen::VectorXi;
// using Eigen::VectorXd;
// using Eigen::ArrayXi;
// using Eigen::ArrayXd;

using namespace Eigen;

typedef Array<bool,Dynamic,1> ArrayXb;
// typedef Matrix<bool,Dynamic,1> VectorXb;

// modified on 1.23.17
// add gamma to the parameter!
// modified on 10.31.16
// will add intermediate value of w*v into cell type object as well!!
//define cell type struct  (for easy returns)
// contains w, v, h, z
struct Cell_type
{
  MatrixXd w;
  VectorXi v;
  VectorXd h;
  VectorXd z_val;
  VectorXd wv; // w times v
  ArrayXb z;
  double gamma;
  int robust_geno; // add robust geno to cell type
  // double var; // add variance (defaul variance for h) // obselete (remove var)
  Cell_type(const MatrixXd & _w,const VectorXi & _v,const VectorXd & _h,const VectorXd & _z_val,const VectorXd & _wv, const ArrayXb & _z, double _gamma, int _robust_geno): w(_w),v(_v),h(_h),z_val(_z_val),wv(_wv),z(_z),gamma(_gamma),robust_geno(_robust_geno){}
  Cell_type(const MatrixXd & _w,const VectorXi & _v,const VectorXd & _h,const VectorXd & _z_val,const VectorXd & _wv, const ArrayXb & _z): w(_w),v(_v),h(_h),z_val(_z_val),wv(_wv),z(_z){}
  Cell_type(const MatrixXd & _w,const VectorXi & _v,const VectorXd & _h,const VectorXd & _z_val,const VectorXd & _wv, const ArrayXb & _z, double _gamma): w(_w),v(_v),h(_h),z_val(_z_val),wv(_wv),z(_z),gamma(_gamma){}
  Cell_type(){} // empty constructor
  Cell_type(const Cell_type & ct) : w(ct.w), v(ct.v), h(ct.h), z_val(ct.z_val), wv(ct.wv), z(ct.z), gamma(ct.gamma),robust_geno(ct.robust_geno){}
};

VectorXd calc_cell_type_wv(const MatrixXd & w, const VectorXi & v);
VectorXd calc_cell_type_wv(const MatrixXd & w, const VectorXi & v, uint32_t func_geno_len);
VectorXd calc_cell_type_z_val_from_wv(const VectorXd & wv, const VectorXd & h);
VectorXd calc_cell_type_z_val(const MatrixXd & w, const VectorXi & v, const VectorXd & h);

ArrayXb calc_cell_type(const MatrixXd & w, const VectorXi & v, const VectorXd & h);
void generate_rand_h(VectorXd & h, double sigma,const gsl_rng * r);
void generate_rand_h(VectorXd & h, double sigma,const double * norm_rand_tb,const gsl_rng * r);
// void generate_rand_h(VectorXd & h, double sigma,const gsl_rng * r);


// void resample_h(VectorXd & h,int robust_geno,double var, const double * norm_rand_tb, const gsl_rng * r);
void resample_h_from_alpha(VectorXd & h,double alpha,double var, const gsl_rng * r);
void resample_h_from_alpha(VectorXd & h,double alpha,double var,const double * norm_rand_tb, const gsl_rng * r);
template<typename iterator_t, typename mcont_t> VectorXd recalc_mat(const MatrixXd & w, const VectorXi & v, const int & robust_bits, const iterator_t begin, const iterator_t end, const mcont_t  & mutations,  int & changed_bits, int flag) ;
template<typename iterator_t, typename mcont_t> VectorXd recalc_mat(const MatrixXd & w, const VectorXi & v, const int & robust_bits, const iterator_t begin, const iterator_t end, const mcont_t  & mutations,  double ratio, int & changed_bits, int flag) ;
double eval_fitness(const ArrayXb & z, const ArrayXb & b_vec) ;
double calc_fitness_from_wv_h(const VectorXd & wv, const VectorXd & h, const ArrayXb & b_vec);

int get_robust_geno(const VectorXi & v, uint32_t robust_bits);
double calc_alpha(uint32_t robust_degrees, uint32_t robust_geno);
int get_robust_geno_v2(const VectorXi & v, uint32_t robust_bits);
void set_robust_geno(VectorXi & v, uint32_t robust_val, uint32_t robust_bits);
void set_robust_geno_v2(VectorXi & v, uint32_t robust_val);
Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, uint32_t robust_bits, int flag,  double c, double pv1, double ga, const gsl_rng * r, uint32_t robust_val);
// Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, uint32_t robust_bits, int flag, uint32_t robust_degrees, double c, double pv1,  const gsl_rng * r, uint32_t robust_val);
Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, double c, double pv1, double ga, const gsl_rng * r);
// double calc_h_sd(const VectorXd & h);
template<typename iterator_t, typename mcont_t> int muts_to_robust(const iterator_t begin, const iterator_t end, const mcont_t  & mutations,const VectorXi & v,  int robust_bits, int flag);
void generate_b_vec(ArrayXb & b_vec, const ArrayXb & z, int type, const gsl_rng * r, int diffk =0);
void copy_back(VectorXi & v, const VectorXi & v1, const VectorXi & v2);
// void copy_back(VectorXi & v, const VectorXi & v2);
void generate_b_dis(ArrayXb & b_vec, const ArrayXb & z, int dis, const gsl_rng * r);
VectorXi generate_b_vec(ArrayXb & b_vec, double pv1, const MatrixXd & w, const VectorXd & h,  const gsl_rng * r );


VectorXd update_h_vec(const int & diff, int robust_degrees, const VectorXd & h_vec);
double calc_fitness(const MatrixXd & w, const VectorXi & v, uint32_t func_geno_len, const VectorXd & h, const ArrayXb & b_vec);

Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, double c, double pv1, double ga);
Cell_type  calc_cell_type_from_para_indeph(uint32_t nz, uint32_t Lv, double c, double pv1, double ga, GSLrng & r);
// VectorXd update_h_vec(double scale, const VectorXd  & wv, int robust_bits, const VectorXd & h_vec);
VectorXd update_wv_scale(double scale, const VectorXd & wv);
// #include"ct_eigen.tcc"

double calc_scale(uint32_t robust_degrees, uint32_t robust_geno);
template<typename iterator_t,
         typename mcont_t>
VectorXd recalc_mat(const MatrixXd & w, const VectorXi & v, iterator_t begin, iterator_t end, const mcont_t & mutations, double ratio, int relative_index=0);
// double calc_scale_threshold(uint32_t robust_degrees, uint32_t robust_geno );


double calc_scale_2base(uint32_t robust_degrees, uint32_t robust_geno );
double calc_scale_8_degree(uint32_t robust_degrees, uint32_t robust_geno);
double calc_scale_16_degree(uint32_t robust_degrees, uint32_t robust_geno);

#endif
