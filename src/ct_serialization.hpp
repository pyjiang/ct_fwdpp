// write my own mutation serialization because it is from mutation_base
#ifndef _CT_SERIAL_
#define _CT_SERIAL_

#include  <type_traits>
// #include  <fwdpp/forward_types.hpp>
#include  <iostream>
#include "ct_mutation.hpp"


namespace CT
{
  struct mutation_writer
  {
    using result_type = void;

    // activate only for mutation_base type
    template<typename mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,KTfwd::mutation_base>::value,result_type>::type
    operator()( const mutation_t &m,
		std::ostream & buffer) const
    {
      buffer.write( reinterpret_cast<const char *>(&m.pos),sizeof(double));
    }

    template<typename mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,CT::ct_mutation>::value,result_type>::type
    operator()( const mutation_t &m,
		std::ostream & buffer) const
    {
      buffer.write( reinterpret_cast<const char *>(&m.pos),sizeof(double));
      buffer.write( reinterpret_cast<const char *>(&m.scale),sizeof(double));
    }
  };


  template<typename mutation_t>
  struct mutation_reader
  {
    using result_type = mutation_t;

    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<U,KTfwd::mutation_base>::value,result_type>::type
    operator()( std::istream & in ) const
    {
      double pos;
      in.read( reinterpret_cast<char *>(&pos),sizeof(double));
      return result_type(pos,false); // do not set as neutral mutations
    }

    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<U,CT::ct_mutation>::value,result_type>::type
    operator()( std::istream & in ) const
    {
      double pos, scale;
      in.read( reinterpret_cast<char *>(&pos),sizeof(double));
      in.read( reinterpret_cast<char *>(&scale),sizeof(double));
      return result_type(pos,scale,false); // do not set as neutral mutations
    }
  };
}
#endif
