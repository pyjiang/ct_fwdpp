//  -*- C++ -*-
#ifndef __CT_IO_TCC__
#define __CT_IO_TCC__

#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/IOhelp.hpp>
#include "ct_internal_IOhelp.hpp"
#include <fwdpp/sugar/GSLrng_t.hpp>

using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

// modified from IO.tcc, will call functions under CT_internal namespace rather than KTfwd
// also needs to write global varibales, like v, w, h, z etc

namespace CT
{

  struct diploidIOplaceholder
    {
      using result_type = void;
      /*!
  Does nothing
      */
      template< typename iterator, typename streamtype >
      inline result_type operator()( iterator , streamtype & ) const
      {
  //Does nothing!
      }
    };



  // modifed on 1.5.17
  // rather than the original functions, now write functors for write_binary_pop_hap,
  // so that when write different types of populations, less redundant code will be written for higher level functions


  struct write_binary_pop_hap
  {

    //Binary I/O for individual-based simulation
    // basic haploid case
    template< typename gcont_t,
        typename mcont_t,
        typename dipvector_t,
        typename mutation_writer_type,
        typename ostreamtype,
        typename diploid_writer_t = diploidIOplaceholder,
        typename cell_type,
        typename Vec>
    void operator() ( const gcont_t & gametes,
  			  const mcont_t & mutations,
  			  const dipvector_t & haploids,
  			  const mutation_writer_type & mw,
          const std::vector<uint_t> & mcounts, // counts of mutations
          const cell_type & ct, // global cell type object
          const std::size_t & gen, // current stopped generation number
          const double & mu, // mutation rate
          const double & recomb, // recombination rate
          const std::string & out_common_name, // common name of output files
          const Vec & b_vec, // b vector
  			  ostreamtype & buffer,
  			  const diploid_writer_t & dw = diploid_writer_t() ) const
    {
        CT_internal::write_mutations()( mutations,mw,buffer ); // modified on 2.22.17
        CT_internal::write_haplotypes()( gametes, buffer );
        CT_internal::write_mcounts() (mcounts, buffer); // write mcounts
        CT_internal::write_cell_types()(ct,buffer); // write global cell type object
        buffer.write(reinterpret_cast<const char *>(&gen),sizeof(std::size_t)); // write current generation
        // gsl_rng_fwrite(buffer,r); // write gsl state to buffer

        buffer.write(reinterpret_cast<const char *>(&mu),sizeof(double)); // write mutation rate
        buffer.write(reinterpret_cast<const char *>(&recomb),sizeof(double)); // write recombination rate

        std::size_t str_size= out_common_name.size();
        buffer.write(reinterpret_cast<const char *>(&str_size),sizeof(std::size_t));// write output length of common name
        buffer.write(out_common_name.c_str(),str_size); // write output common name
        CT_internal::write_eigen_matrix(b_vec,buffer); // write b vector

        std::size_t NDIPS = haploids.size();
        buffer.write( reinterpret_cast<char *>(&NDIPS), sizeof(std::size_t) );

        for(const auto & hap : haploids )
        {
        	buffer.write(reinterpret_cast<const char *>(&hap),sizeof(std::size_t));
        	dw(hap,buffer);
        }
    }


    // modified on 12.6.16
    // added flag of robustness, need to write to state file as well
    // this function will be compatible with  robust bits in v
    // without re-writing the previous function
    template< typename gcont_t,
        typename mcont_t,
        typename dipvector_t,
        typename mutation_writer_type,
        typename ostreamtype,
        typename diploid_writer_t = diploidIOplaceholder,
        typename cell_type,
        typename Vec>
    void operator() ( const gcont_t & gametes,
          const mcont_t & mutations,
          const dipvector_t & haploids,
          const mutation_writer_type & mw,
          const std::vector<uint_t> & mcounts, // counts of mutations
          const cell_type & ct, // global cell type object
          const std::size_t & gen, // current stopped generation number
          const double & mu, // mutation rate
          const double & recomb, // recombination rate
          const std::string & out_common_name, // common name of output files
          const Vec & b_vec, // b vector
          const uint32_t robust_bits,
          const int robust_flag, // flag of robustness
          ostreamtype & buffer,
          const diploid_writer_t & dw = diploid_writer_t() ) const
      {
          this->operator()(gametes,mutations,haploids,mw,mcounts,ct,gen,mu,recomb,out_common_name,b_vec,buffer,dw);
          buffer.write(reinterpret_cast<const char *>(&robust_bits),sizeof(uint32_t));
          buffer.write(reinterpret_cast<const char *>(&robust_flag),sizeof(int)); // write robustness flag
      }



      // for all two locus cases
      // created on 12.21.16
      template< typename gcont1_t,
          typename gcont2_t,
          typename mcont_t,
          typename dipvector_t,
          typename mutation_writer_type,
          typename ostreamtype,
          typename diploid_writer_t = diploidIOplaceholder ,
          typename cell_type,
          typename Vec>
      void operator() ( const gcont1_t & gametes1,
            const gcont2_t & gametes2,
            const mcont_t & mutations,
            const dipvector_t & haploids,
            const mutation_writer_type & mw,
            const std::vector<uint_t> & mcounts, // counts of mutations
            const cell_type & ct, // global cell type object
            const std::size_t & gen, // current stopped generation number
            const double & mu1, // mutation rate of gamete1
            const double & mu2, // mutation rate of gamete2
            const double & recomb, // recombination rate
            const std::string & out_common_name, // common name of output files
            const Vec & b_vec, // b vector
            const uint32_t robust_bits,
            const int robust_flag, // flag of robustness
            ostreamtype & buffer,
            const diploid_writer_t & dw = diploid_writer_t() ) const
        {
            // actully pass the function to diploid writer, because here haploids has a pair of value
            write_binary_pop(gametes1,mutations,haploids,mw,mcounts,ct,gen,mu1,recomb,out_common_name,b_vec,buffer,dw);

            buffer.write(reinterpret_cast<const char *>(&robust_bits),sizeof(uint32_t));
            buffer.write(reinterpret_cast<const char *>(&robust_flag),sizeof(int)); // write robustness flag

            CT_internal::write_haplotypes()( gametes2, buffer );
            buffer.write(reinterpret_cast<const char *>(&mu2),sizeof(double));
        }



  };



// modify read functions as well! (as functor)
struct read_binary_pop_hap
{
  // modified on 12.5.16
  // add extra parameter indicating whether to keep gamete flag or not
  // basic type
    template< typename gcont_t,
        typename mcont_t,
        typename dipvector_t,
        typename mutation_reader_type,
        typename istreamtype,
        typename diploid_reader_t = diploidIOplaceholder,
        typename cell_type,
        typename Vec>
    void operator () (  gcont_t & gametes,
                        mcont_t & mutations,
                        dipvector_t & haploids,
                        const mutation_reader_type & mr,
                        std::vector<uint_t> & mcounts, // counts of mutations
                        cell_type & ct, // global cell type object
                        std::size_t * gen_pt, // current stopped generation number
                        double * mu_pt, // mutation rate
                        double * recomb_pt, // recombination rate
                        std::string  out_common_name, // common name of output files
                        Vec & b_vec, // b vector
                        istreamtype & in,
                        int keep_flag, // if keep_flag==1, will keep gamete flag ; if ==0, will not keep
                        const diploid_reader_t & dr = diploid_reader_t()  ) const
    {
      gametes.clear();
      mutations.clear();
      haploids.clear();
      CT_internal::read_mutations()(mutations,mr,in);

      if(keep_flag==1)
        CT_internal::read_haplotypes()(gametes,in,keep_flag);
      else
        CT_internal::read_haplotypes()(gametes,in);

      CT_internal::read_mcounts()(mcounts,in); // read in mcounts
      CT_internal::read_cell_types()(ct,in); // read in cell types
      std::size_t NDIPS,c;

      fwdpp_internal::scalar_reader<std::size_t>()(in,gen_pt); // read in previous generation
      // gsl_rng_fread(in,r); // read in random number state
      fwdpp_internal::scalar_reader<double>()(in,mu_pt); // read in mutation rate
      fwdpp_internal::scalar_reader<double>()(in,recomb_pt); // read in recombination rate
      std::size_t str_size;
      fwdpp_internal::scalar_reader<decltype(str_size)>()(in,&str_size); // read in the size
      out_common_name.resize(str_size); // resize of string
      in.read(&out_common_name[0],str_size); // read in common name

      CT_internal::read_eigen_matrix(b_vec,in); // read in b vector

      fwdpp_internal::scalar_reader<decltype(NDIPS)>()(in,&NDIPS);
      haploids.resize(NDIPS);
      for( auto & hap : haploids )
      {
        fwdpp_internal::scalar_reader<decltype(c)>()(in,&c);
        hap = c;
        dr(hap,in);
      }
    }


    // robust bits
    template< typename gcont_t,
        typename mcont_t,
        typename dipvector_t,
        typename mutation_reader_type,
        typename istreamtype,
        typename diploid_reader_t = diploidIOplaceholder ,
        typename cell_type,
        typename Vec>
    void operator () (  gcont_t & gametes,
          mcont_t & mutations,
          dipvector_t & haploids,
          const mutation_reader_type & mr,
          std::vector<uint_t> & mcounts, // counts of mutations
          cell_type & ct, // global cell type object
          std::size_t * gen_pt, // current stopped generation number
          double * mu_pt, // mutation rate
          double * recomb_pt, // recombination rate
          std::string  out_common_name, // common name of output files
          Vec & b_vec, // b vector
          uint32_t * p_robust_bits,
          int * p_robust_flag,
          istreamtype & in,
          int keep_flag, // if keep_flag==1, will keep gamete flag ; if ==0, will not keep
          const diploid_reader_t & dr  = diploid_reader_t()) const
      {
          this-> operator() (gametes,mutations,haploids,mr,mcounts,ct,gen_pt,mu_pt,recomb_pt,out_common_name,b_vec,in,keep_flag,dr);
          fwdpp_internal::scalar_reader<uint32_t>()(in,p_robust_bits);
          fwdpp_internal::scalar_reader<int>()(in,p_robust_flag);
      }



      // for all two locus cases
      // added on 12.22.16
      template< typename gcont1_t,
          typename gcont2_t,
    	    typename mcont_t,
    	    typename dipvector_t,
    	    typename mutation_reader_type,
    	    typename istreamtype,
    	    typename diploid_reader_t  = diploidIOplaceholder,
          typename cell_type,
          typename Vec>
      void operator () (  gcont1_t & gametes1,
            gcont2_t & gametes2,
    			  mcont_t & mutations,
    			  dipvector_t & haploids,
    			  const mutation_reader_type & mr,
            std::vector<uint_t> & mcounts, // counts of mutations
            cell_type & ct, // global cell type object
            std::size_t * gen_pt, // current stopped generation number
            double * mu1_pt, // mutation rate of locus 1
            double * mu2_pt, // mutation rate of locus 2
            double * recomb_pt, // recombination rate
            std::string  out_common_name, // common name of output files
            Vec & b_vec, // b vector
            uint32_t * p_robust_bits,
            int * p_robust_flag,
    			  istreamtype & in,
    			  const diploid_reader_t & dr = diploid_reader_t()) const
        {
            // use diploid reader
            read_binary_pop(gametes1,mutations,haploids,mr,mcounts,ct,gen_pt,mu1_pt,recomb_pt,out_common_name,b_vec,in,dr);
            fwdpp_internal::scalar_reader<uint32_t>()(in,p_robust_bits);
            fwdpp_internal::scalar_reader<int>()(in,p_robust_flag);

            CT_internal::read_haplotypes()(gametes2,in);
            fwdpp_internal::scalar_reader<double>()(in,mu2_pt);
        }

};





  // diploid case
  //Single-locus sims, single pop
  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t,
      typename cell_type,
      typename Vec>
  inline void write_binary_pop ( const gcont_t & gametes,
			  const mcont_t & mutations,
			  const dipvector_t & diploids,
			  const mutation_writer_type & mw,
        const std::vector<uint_t> & mcounts, // counts of mutations
        const cell_type & ct, // global cell type object
        const std::size_t & gen, // current stopped generation number
        const double & mu, // mutation rate
        const double & recomb, // recombination rate
        const std::string & out_common_name, // common name of output files
        const Vec & b_vec, // b vector
			  ostreamtype & buffer,
			  const diploid_writer_t & dw )
  {
    KTfwd::fwdpp_internal::write_mutations()( mutations,mw,buffer );
    CT_internal::write_haplotypes()( gametes, buffer );
    CT_internal::write_mcounts() (mcounts, buffer); // write mcounts
    CT_internal::write_cell_types()(ct,buffer); // write global cell type object
    buffer.write(reinterpret_cast<const char *>(&gen),sizeof(std::size_t)); // write current generation
    // gsl_rng_fwrite(buffer,r); // write gsl state to buffer

    buffer.write(reinterpret_cast<const char *>(&mu),sizeof(double)); // write mutation rate
    buffer.write(reinterpret_cast<const char *>(&recomb),sizeof(double)); // write recombination rate

    std::size_t str_size= out_common_name.size();
    buffer.write(reinterpret_cast<const char *>(&str_size),sizeof(std::size_t));// write output length of common name
    buffer.write(out_common_name.c_str(),str_size); // write output common name
    CT_internal::write_eigen_matrix(b_vec,buffer); // write b vector

    std::size_t NDIPS = diploids.size();
    buffer.write( reinterpret_cast<char *>(&NDIPS), sizeof(std::size_t) );

    for(const auto & dip : diploids )
    {
    	buffer.write(reinterpret_cast<const char *>(&dip.first),sizeof(std::size_t));
    	buffer.write(reinterpret_cast<const char *>(&dip.second),sizeof(std::size_t));
    	dw(dip,buffer);
    }
  }

  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t,
      typename cell_type,
      typename Vec>
  inline void read_binary_pop (  gcont_t & gametes,
			  mcont_t & mutations,
			  dipvector_t & diploids,
			  const mutation_reader_type & mr,
        std::vector<uint_t> & mcounts, // counts of mutations
        cell_type & ct, // global cell type object
        std::size_t * gen_pt, // current stopped generation number
        double * mu_pt, // mutation rate
        double * recomb_pt, // recombination rate
        std::string  out_common_name, // common name of output files
        Vec & b_vec, // b vector
			  istreamtype & in,
			  const diploid_reader_t & dr)
  {
    gametes.clear();
    mutations.clear();
    diploids.clear();
    fwdpp_internal::read_mutations()(mutations,mr,in);
    CT_internal::read_haplotypes()(gametes,in);
    CT_internal::read_mcounts()(mcounts,in); // read in mcounts
    CT_internal::read_cell_types()(ct,in); // read in cell types
    std::size_t NDIPS,c;

    fwdpp_internal::scalar_reader<std::size_t>()(in,gen_pt); // read in previous generation
    // gsl_rng_fread(in,r); // read in random number state
    fwdpp_internal::scalar_reader<double>()(in,mu_pt); // read in mutation rate
    fwdpp_internal::scalar_reader<double>()(in,recomb_pt); // read in recombination rate
    std::size_t str_size;
    fwdpp_internal::scalar_reader<decltype(str_size)>()(in,&str_size); // read in the size
    out_common_name.resize(str_size); // resize of string
    in.read(&out_common_name[0],str_size); // read in common name

    CT_internal::read_eigen_matrix(b_vec,in); // read in b vector

    fwdpp_internal::scalar_reader<decltype(NDIPS)>()(in,&NDIPS);
    diploids.resize(NDIPS);
    for( auto & dip : diploids )
    {
    	fwdpp_internal::scalar_reader<decltype(c)>()(in,&c);
    	dip.first = c;
    	fwdpp_internal::scalar_reader<decltype(c)>()(in,&c);
    	dip.second = c;
    	dr(dip,in);
    }
    //  std::cout << "test";
  }




}//ns KTfwd

#endif
