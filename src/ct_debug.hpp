#ifndef __CT_DEBUG_HPP__
#define __CT_DEBUG_HPP__

// added assertion functions to make some checks

// added on 3.9.17
// modifed version to check gamete flag  (only check those present in haploids)
template<typename gcont_t,
         typename haploids_t>
bool check_gamete_flags(const gcont_t & gametes, const haploids_t & haploids)
{
  for(const auto & hap : haploids)
  {
    if (gametes[hap].flag!=0)
    {
      return false;
    }
  }
  return true ;
}



// make sure all gamete flags are 0 (their fitness have been calculated)
template<typename gcont_t>
bool check_gamete_flags(const gcont_t & gametes)
{
  for(const auto & g :gametes )
  {
    if (g.flag!=0)
    {
      return false;
    }
  }
  return true;
}


// check if the total number of gamete equals 2N
template<typename gcont_t>
bool check_gamete_count(const gcont_t & gametes, int total )
{
  int count=0;
  for(const auto & g :gametes )
  {
    count += g.n;
  }
  if(count != total ) return false;
  return true;
}

#endif
