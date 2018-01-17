#ifndef __CT_INTERNAL_IOHELP_HPP__
#define __CT_INTERNAL_IOHELP_HPP__

/*
  Mechanics of data serialization
  The various write_binary_pop and read_binary pop
  functions rely on these implementations
*/
#include <zlib.h>

#include "ct_gamete.hpp"
// #include "noise.hpp"

namespace CT {
  namespace CT_internal {

    using ct_gamete_t=ct_gamete<void>;
    // using ct_gamete_robust_geno_t=ct_gamete_robust_geno<void,VectorXd>;
    // using ct_gamete_robust_geno_only_t=ct_gamete_robust_geno_only<void,VectorXd>;
    using ct_gamete_wv_t=ct_gamete_wv<void,VectorXd>;
    // using ct_gamete_scale_wij_t = ct_gamete_scale_wij<void,VectorXd>;
    using ct_gamete_geno_no_h_t = ct_gamete_robust_geno_only_wij<void>; // added on 1.30.17



    template<typename T>
    struct scalar_reader
    {
      template<typename streamtype>
      inline void operator()( streamtype & i, T * __t ) const
      {
	       i.read( reinterpret_cast<char*>(__t), sizeof(T) );
      }
      inline void operator()( gzFile & gzin, T * __t ) const
      {
	       gzread(gzin,__t,sizeof(T));
      }
    };

    struct write_mutations
    {
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_writer_type,
		typename ostreamtype>
      void operator()( const list_type< mutation_type, list_type_allocator > & mutations,
		       const mutation_writer_type & mw,
		       ostreamtype & buffer) const
      {
	std::size_t MUTNO = mutations.size();
	buffer.write( reinterpret_cast<char *>(&MUTNO), sizeof(std::size_t) );
	//write the mutation data to the buffer
	for(const auto & m : mutations ) mw(m,buffer);
      }
    };

    // reference: http://stackoverflow.com/questions/25389480/how-to-write-read-an-eigen-matrix-from-binary-file

    template <typename Matrix,
              typename ostreamtype>
    void write_eigen_matrix (const Matrix & matrix, ostreamtype & out)
    {
      typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
      out.write((char*) (&rows), sizeof(typename Matrix::Index));
      out.write((char*) (&cols), sizeof(typename Matrix::Index));
      out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    };

    // 11.18.16 modified
    // note: the written file is different from previous version (order is different, but the total info is the same!)
    struct write_haplotypes
    {
      template<typename gamete_type,
	             typename... gamete_cont_t_details,
	             template<typename,typename... > class gamete_cont_t,
	             typename ostreamtype>
      void operator()( const gamete_cont_t< gamete_type, gamete_cont_t_details... > & gametes,
		       ostreamtype & buffer) const
      {
	       std::size_t N = gametes.size();
	       buffer.write( reinterpret_cast< char * >(&N), sizeof(std::size_t) );
         for( const auto & g : gametes )
         {
        	    buffer.write(reinterpret_cast<const char *>(&g.n),sizeof(decltype(g.n)));
              buffer.write(reinterpret_cast<const char *>(&g.flag),sizeof(decltype(g.flag))); // write flag

              std::size_t nm = g.mutations.size();
        	    buffer.write(reinterpret_cast<const char *>(&nm),sizeof(std::size_t));
        	    if(nm)
        	      {
        		        buffer.write(reinterpret_cast<const char *>(&g.mutations[0]),
        			      nm*sizeof(typename gamete_type::mutation_container::value_type));
        	      }
        	    nm = g.smutations.size();
        	    buffer.write(reinterpret_cast<const char *>(&nm),sizeof(std::size_t));
        	    if(nm)
        	      {
        		        buffer.write(reinterpret_cast<const char *>(&g.smutations[0]),
        			      nm*sizeof(typename gamete_type::mutation_container::value_type));
        	      }
              // if gamete type is ct_gamete_with_h_t
              // will write h_vec and robust_geno
              write_extra(g,buffer);
        	  }
      }


      // template<typename gamete_type,
      //         typename ostreamtype>
      // typename std::enable_if<std::is_same<gamete_type,ct_gamete_robust_geno_only_t>::value>::type
      // write_extra(const gamete_type & g, ostreamtype & buffer ) const
      // {
      //   write_eigen_matrix(g.h_vec, buffer);
      //   buffer.write(reinterpret_cast<const char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
      //   buffer.write(reinterpret_cast<const char *>(&g.alpha),sizeof(decltype(g.alpha)));
      // }

      // added on 1.30.17
      template<typename gamete_type,
              typename ostreamtype>
      typename std::enable_if<std::is_same<gamete_type,ct_gamete_geno_no_h_t>::value>::type
      write_extra(const gamete_type & g, ostreamtype & buffer ) const
      {
        // write_eigen_matrix(g.h_vec, buffer);
        buffer.write(reinterpret_cast<const char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
        // buffer.write(reinterpret_cast<const char *>(&g.initial_scale),sizeof(decltype(g.initial_scale)));
        buffer.write(reinterpret_cast<const char *>(&g.scale),sizeof(decltype(g.scale)));
        // buffer.write(reinterpret_cast<const char *>(&g.alpha),sizeof(decltype(g.alpha)));
        buffer.write(reinterpret_cast<const char *>(&g.scale_func_flag),sizeof(decltype(g.scale_func_flag))); // store scale_func_flag
      }

      // // added on 1.23.17
      // // for ct_gamete_scale_wij_t
      // template<typename gamete_type,
      //         typename ostreamtype>
      // typename std::enable_if<std::is_same<gamete_type,ct_gamete_scale_wij_t>::value>::type
      // write_extra(const gamete_type & g, ostreamtype & buffer ) const
      // {
      //   buffer.write(reinterpret_cast<const char *>(&g.fitness),sizeof(decltype(g.fitness))); // write fitness
      //   buffer.write(reinterpret_cast<const char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
      //   write_eigen_matrix(g.wv, buffer);
      //   // buffer.write(reinterpret_cast<const char *>(&g.initial_alpha),sizeof(decltype(g.initial_alpha)));
      //   // buffer.write(reinterpret_cast<const char *>(&g.initial_scale),sizeof(decltype(g.initial_scale)));
      //   // buffer.write(reinterpret_cast<const char *>(&g.alpha),sizeof(decltype(g.alpha)));
      //   buffer.write(reinterpret_cast<const char *>(&g.scale),sizeof(decltype(g.scale)));
      // }


      // modified on 1.14.17
      // must write out fitness (for compatible with the default type)
      template<typename gamete_type,
              typename ostreamtype>
      typename std::enable_if<std::is_same<gamete_type,ct_gamete_wv_t>::value>::type
      write_extra(const gamete_type & g, ostreamtype & buffer ) const
      {
        buffer.write(reinterpret_cast<const char *>(&g.fitness),sizeof(decltype(g.fitness))); // write fitness
        write_eigen_matrix(g.wv, buffer);
      }


      //
      template<typename gamete_type,
              typename ostreamtype>
      typename std::enable_if<std::is_same<gamete_type,ct_gamete_t>::value>::type
      write_extra(const gamete_type & g, ostreamtype & buffer ) const
      {
        buffer.write(reinterpret_cast<const char *>(&g.fitness),sizeof(decltype(g.fitness))); // write fitness
      }

      // template<typename gamete_type,
      //         typename ostreamtype>
      // typename std::enable_if<std::is_same<gamete_type,ct_gamete_robust_geno_t>::value>::type
      // write_extra(const gamete_type & g, ostreamtype & buffer ) const
      // {
      //   buffer.write(reinterpret_cast<const char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
      //   buffer.write(reinterpret_cast<const char *>(&g.alpha),sizeof(decltype(g.alpha)));
      //   write_eigen_matrix(g.wv, buffer);
      // }

    };

    struct write_mcounts
    {
      template <typename mcont_t,
                typename ostreamtype>
      void operator () (const mcont_t & mcounts, ostreamtype & buffer)
      {
        uint_t len= mcounts.size();
        buffer.write( reinterpret_cast< char * >(&len), sizeof(uint_t) );
        for (const auto & count : mcounts )
        {
          buffer.write(reinterpret_cast<const char *>(&count),sizeof(decltype(count)));
        }
      }
    };






    // modified this function on 11.7.16
    // because cell type is modified that it stores extra values
    // so this function will write the extra values
    struct write_cell_types
    {
      template <typename cell_type,
                typename ostreamtype>
      void operator() (const cell_type & ct, ostreamtype & buffer) const
      {
        write_eigen_matrix(ct.v, buffer); //write v
        write_eigen_matrix(ct.w, buffer); //write w
        write_eigen_matrix(ct.h, buffer); //write h
        // write robust_geno and var
        buffer.write(reinterpret_cast<const char *>(&ct.robust_geno),sizeof(decltype(ct.robust_geno)));
        // buffer.write(reinterpret_cast<const char *>(&ct.var),sizeof(decltype(ct.var)));
        buffer.write(reinterpret_cast<const char *>(&ct.gamma),sizeof(decltype(ct.gamma)));
      }
    };


    template <typename Matrix,
              typename istreamtype>
    void read_eigen_matrix (Matrix & matrix, istreamtype & in)
    {
      typename Matrix::Index rows=0, cols=0;
      in.read((char*) (&rows),sizeof(typename Matrix::Index));
      in.read((char*) (&cols),sizeof(typename Matrix::Index));
      matrix.resize(rows, cols);
      in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
    };

    // modified this function on 11.7.16
    // because cell type is modified that it stores extra values
    // so this function will write the extra values
    struct read_cell_types
    {
      template <typename cell_type,
                typename istreamtype>
      void operator() (cell_type & ct, istreamtype & in) const
      {
        read_eigen_matrix(ct.v, in);
        read_eigen_matrix(ct.w, in);
        read_eigen_matrix(ct.h, in);
        in.read( reinterpret_cast<char *>(&ct.robust_geno), sizeof(decltype(ct.robust_geno)) );
        // in.read( reinterpret_cast<char *>(&ct.var), sizeof(decltype(ct.var)) );
        in.read( reinterpret_cast<char *>(&ct.gamma), sizeof(decltype(ct.gamma)) );
      }
    };

    struct read_mcounts
    {
      template <typename mcont_t,
                typename istreamtype>
      void operator () (mcont_t & mcounts, istreamtype & in)
      {
        uint_t len;
        in.read( reinterpret_cast<char *>(&len), sizeof(decltype(len)) );

        uint_t count;

        for (uint_t i=0; i <len; i++)
        {
          in.read( reinterpret_cast<char *>(&count), sizeof(decltype(count)) );
          mcounts.emplace_back(count);
        }
      }
    };



    struct read_mutations
    {
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_reader,
		typename istreamtype >
      void operator()(list_type< mutation_type, list_type_allocator > & mutations,
		      const mutation_reader & mr,
		      istreamtype & in) const
      {
	std::size_t NMUTS;
	in.read( reinterpret_cast<char *>(&NMUTS), sizeof(decltype(NMUTS)) );
	for(uint_t i = 0 ; i < NMUTS ; ++i)
	  {
	    mutations.emplace_back(mr(in));
	  }
      }

      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_reader>
      void
      operator()(list_type< mutation_type, list_type_allocator > & mutations,
		 const mutation_reader & mr,
		 gzFile gzin) const
      {
	std::size_t NMUTS;
	gzread( gzin, &NMUTS, sizeof(decltype(NMUTS)) );
	for(uint_t i = 0 ; i < NMUTS ; ++i)
	  {
	    mutations.emplace_back( mr(gzin) );
	  }
      }
    };

    // modified on 12.5.16
    // given an extra parameter keep_flag, specifying whether will set the gamete flag to 1 or to keep its original in the state file
    struct read_haplotypes
    {
      // keep gamete flag
      template< typename gamete_type,
		            typename list_type_allocator,
		            template<typename,typename> class list_type,
		            typename istreamtype >
      void  operator()(list_type< gamete_type, list_type_allocator > & gametes,
		                  istreamtype & in,
                      int keep_flag) const
      {
      	std::size_t NHAPS;
      	in.read( reinterpret_cast<char *>(&NHAPS),sizeof(decltype(NHAPS)) );
      	uint_t N;
        int flag;
        // double fitness;
      	std::size_t nm;
      	for( uint_t i = 0 ; i < NHAPS ; ++i )
      	{
      	    in.read(reinterpret_cast<char *>(&N),sizeof(decltype(N)));
            in.read(reinterpret_cast<char *>(&flag),sizeof(decltype(flag))); // read flag

      	    // gamete_type g(N,flag,fitness); // initialize gamete

            // modified on 1.14.17 (no fitness initialized here!)
            gamete_type g(N,flag); // initialize gamete
      	    in.read(reinterpret_cast<char *>(&nm),sizeof(decltype(nm)));
      	    if(nm)
      	      {
      		        g.mutations.resize(nm);
      		        in.read(reinterpret_cast<char *>(&g.mutations[0]),nm*sizeof(typename gamete_type::mutation_container::value_type));
      	      }
      	    in.read(reinterpret_cast<char *>(&nm),sizeof(decltype(nm)));
      	    if(nm)
      	      {
      		        g.smutations.resize(nm);
      		        in.read(reinterpret_cast<char *>(&g.smutations[0]),nm*sizeof(typename gamete_type::mutation_container::value_type));
      	      }

            // if gamete type is ct_gamete_with_h
            read_extra(g,in);

      	    gametes.emplace_back(std::move(g));
      	  }
        }

        // default will not use flag info in the gamete (set to 1)
        template< typename gamete_type,
  		            typename list_type_allocator,
  		            template<typename,typename> class list_type,
  		            typename istreamtype >
        void  operator()(list_type< gamete_type, list_type_allocator > & gametes,
  		                  istreamtype & in) const
        {
        	std::size_t NHAPS;
        	in.read( reinterpret_cast<char *>(&NHAPS),sizeof(decltype(NHAPS)) );
        	uint_t N;
          int flag;
          // double fitness;
        	std::size_t nm;
        	for( uint_t i = 0 ; i < NHAPS ; ++i )
        	{
        	    in.read(reinterpret_cast<char *>(&N),sizeof(decltype(N)));
              in.read(reinterpret_cast<char *>(&flag),sizeof(decltype(flag))); // read flag

        	    // gamete_type g(N,fitness); // initialize gamete (need to change flag)
              gamete_type g(N);
        	    in.read(reinterpret_cast<char *>(&nm),sizeof(decltype(nm)));
        	    if(nm)
        	      {
        		        g.mutations.resize(nm);
        		        in.read(reinterpret_cast<char *>(&g.mutations[0]),nm*sizeof(typename gamete_type::mutation_container::value_type));
        	      }
        	    in.read(reinterpret_cast<char *>(&nm),sizeof(decltype(nm)));
        	    if(nm)
        	      {
        		        g.smutations.resize(nm);
        		        in.read(reinterpret_cast<char *>(&g.smutations[0]),nm*sizeof(typename gamete_type::mutation_container::value_type));
        	      }

              // if gamete type is ct_gamete_with_h
              read_extra(g,in);

        	    gametes.emplace_back(std::move(g));
        	  }
          }



          // template<typename gamete_type,
          //         typename istreamtype>
          // typename std::enable_if<std::is_same<gamete_type,ct_gamete_robust_geno_only_t>::value>::type
          // read_extra(gamete_type & g, istreamtype & in ) const
          // {
          //   // in.read(reinterpret_cast<char *>(&g.fitness),sizeof(decltype(g.fitness))); // read fitness
          //   read_eigen_matrix(g.h_vec, in);
          //   in.read(reinterpret_cast<char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
          //   in.read(reinterpret_cast<char *>(&g.alpha),sizeof(decltype(g.alpha)));
          // }



          // modified on 2.6.17
          // initialize function pointer
          // added on 1.30.17
          template<typename gamete_type,
                  typename istreamtype>
          typename std::enable_if<std::is_same<gamete_type,ct_gamete_geno_no_h_t>::value>::type
          read_extra(gamete_type & g, istreamtype & in ) const
          {
            // in.read(reinterpret_cast<char *>(&g.fitness),sizeof(decltype(g.fitness))); // read fitness
            // read_eigen_matrix(g.h_vec, in);
            in.read(reinterpret_cast<char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
            // in.read(reinterpret_cast<char *>(&g.initial_scale),sizeof(decltype(g.initial_scale)));
            in.read(reinterpret_cast<char *>(&g.scale),sizeof(decltype(g.scale)));
            // in.read(reinterpret_cast<char *>(&g.alpha),sizeof(decltype(g.alpha)));
            in.read(reinterpret_cast<char *>(&g.scale_func_flag),sizeof(decltype(g.scale_func_flag)));
            g.initialize_calc_function_ptr();
          }


          // // added on 1.23.17
          // // for gamete type
          // template<typename gamete_type,
          //         typename istreamtype>
          // typename std::enable_if<std::is_same<gamete_type,ct_gamete_scale_wij_t>::value>::type
          // read_extra(gamete_type & g, istreamtype & in ) const
          // {
          //   in.read(reinterpret_cast<char *>(&g.fitness),sizeof(decltype(g.fitness))); // read fitness
          //   in.read(reinterpret_cast<char *>(&g.robust_geno),sizeof(decltype(g.robust_geno)));
          //   read_eigen_matrix(g.wv, in);
          //   // in.read(reinterpret_cast<char *>(&g.initial_alpha),sizeof(decltype(g.initial_alpha)));
          //   // in.read(reinterpret_cast<char *>(&g.initial_scale),sizeof(decltype(g.initial_scale)));
          //   in.read(reinterpret_cast<char *>(&g.scale),sizeof(decltype(g.scale)));
          // }

          // modified on 1.14.17
          // added read fitness to be compatible with basic type
          template<typename gamete_type,
                  typename istreamtype>
          typename std::enable_if<std::is_same<gamete_type,ct_gamete_wv_t>::value>::type
          read_extra(gamete_type & g, istreamtype & in ) const
          {
            in.read(reinterpret_cast<char *>(&g.fitness),sizeof(decltype(g.fitness))); // read fitness
            read_eigen_matrix(g.wv, in);
          }


          template<typename gamete_type,
                  typename istreamtype>
          typename std::enable_if<std::is_same<gamete_type,ct_gamete_t>::value>::type
          read_extra(gamete_type & g, istreamtype & in ) const
          {
            in.read(reinterpret_cast<char *>(&g.fitness),sizeof(decltype(g.fitness))); // read fitness
          }


      template< typename gamete_type,
		typename list_type_allocator,
		template<typename,typename> class list_type>
      void operator()(list_type< gamete_type, list_type_allocator > & gametes,
		      gzFile gzin) const
      {
	std::size_t NHAPS;
	gzread( gzin,&NHAPS,sizeof(decltype(NHAPS)) );
	uint_t N;
	std::size_t nm;
	for( uint_t i = 0 ; i < NHAPS ; ++i )
	  {
	    gzread( gzin,&N,sizeof(decltype(N)) );
	    gamete_type g(N);
	    gzread( gzin,&nm,sizeof(decltype(nm)) );
	    if(nm)
	      {
		g.mutations.resize(nm);
		gzread( gzin,&g.mutations[0],nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    gzread( gzin,&nm,sizeof(decltype(nm)) );
	    if(nm)
	      {
		g.smutations.resize(nm);
		gzread( gzin,&g.smutations[0],nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    gametes.emplace_back(std::move(g));
	  }
      }
    };
  }
}

#endif
