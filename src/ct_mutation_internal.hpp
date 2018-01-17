#ifndef __CT_INTERNAL_MUTATION_HPP__
#define __CT_INTERNAL_MUTATION_HPP__

#include <algorithm>
#include <cassert>
#include "ct_mut_model.hpp"
#include "ct_eigen.hpp"

#include "ct_singlepop.hpp" // defines different types of pop

namespace CT {
  namespace CT_internal
  {
    /*!
      Mechanics of adding a new mutation.

      Ensures that it is always added into a vector
      in less-than-sorted order according to position.

      The insertion position is found via a binary search,
      and the call to emplace will typically result in a
      call to memmove or memcopy, depending on the system/compiler.
    */
    // will search if the mutation exists in the list or not!
    // changed the return type of the mutation to mut_pos_indx_t (of an index and a flag, not just index!)
    template< typename mcont_t,
	            typename gamete_type>
    void add_new_mutation( const mut_pos_indx_t idx_obj,
			                     const mcont_t & mutations,
			                     gamete_type & new_gamete )
    {
      uint_t idx= idx_obj.index;

      assert(idx<mutations.size());
      //Establish a pointer to the container to be modified
      auto mc = (mutations[idx].neutral) ? &new_gamete.mutations : &new_gamete.smutations;

      if (idx_obj.flag ==0) // if does not exist, just add the index to the gamete!
      {
        //Insert new mutation at position retaining sort-by-position order in the gamete!!
        mc->emplace(std::upper_bound(mc->begin(),
                                    mc->end(),mutations[idx].pos,
                                    [&mutations](const double & __value,const std::size_t __mut) noexcept {
                                    assert(__mut<mutations.size());
                                    return __value < mutations[__mut].pos;}),
                                    idx );
      }
      else// if it exists, will find out whether the index exists in current gamete or not! if exists, will remove the link, otherwise will add the link
      {
        // assert( std::find(mc->cbegin(),mc->cend(),idx) == mc->end() );
        auto iter_pos=std::find(mc->begin(),mc->end(),idx); //the original code is using cbegin and cend, so iter_pos is a const iterator!! therefore cannot used by erase to modify the Vector!!!
                                                            // when changed to begin() and end(), it solves the problem!

        if (iter_pos==mc->end()) // does not exists in the gamete mutation list, add the new mutation
        {
          //Insert new mutation at position retaining sort-by-position order
          mc->emplace(std::upper_bound(mc->begin(),
    				                          mc->end(),mutations[idx].pos,
    				                          [&mutations](const double & __value,const std::size_t __mut) noexcept {
    				                          assert(__mut<mutations.size());
    				                          return __value < mutations[__mut].pos;}),
    		                              idx );
          // do not modify mutation number now (this is useless, because # of gamete will change and will influence the count at the end of gen)
          // then mutation count +1
          //mutations[idx].xtra+=1;

          //Check post-condition in debug mode...
          assert(std::is_sorted(mc->cbegin(),mc->cend(),
    			                      [&mutations](const std::size_t i, const std::size_t j)noexcept
    			                      {
    			                        return mutations[i].pos<mutations[j].pos;
    			                      }));

        }
        else //exists in the list, will remove the mutation
        {
          mc->erase(iter_pos);
          // then mutation count -1
          //mutations[idx].xtra-=1;
        }
      }

    }

    // template<typename mmodel,
    // 	     typename gamete_type,
    // 	     typename mlist_type,
    // 	     typename queue_t>
    // inline typename std::result_of<mmodel(queue_t &,mlist_type &)>::type
    // mmodel_dispatcher( const mmodel & m, gamete_type & , mlist_type & mutations, queue_t & recycling_bin)
    // /*!
    //   Run-time dispatcher for mutation model
    // */
    // {
    //   return m(recycling_bin,mutations);
    // }

    template<typename mmodel,
	     typename gamete_type,
    	     typename mlist_type,
    	     typename queue_t>
    inline typename std::result_of<mmodel(queue_t &,mlist_type &)>::type
    mmodel_dispatcher( const mmodel & m, gamete_type & , mlist_type & mutations, queue_t & recycling_bin)
    /*!
      Run-time dispatcher for mutation model
    */
    {
      return m(recycling_bin,mutations);
    }

    /*!
      Apply mutation model N times to a new gamete.
      Updates mutation list
    */

    // struct idx_Comp
    // {
    //   mut_pos_indx_t idx;
    //   bool operator()(const mut_pos_indx_t & idx_obj2)
    //   {
    //     return (idx.index ==idx_obj2.index );
    //   }
    // };



    // modified on 12.20.16
    // add relative_index indicating where v starts
    // put common functionality in this one
    template<typename queue_type,
             typename mutation_model,
             typename mlist_type,
             typename gamete_type,
             typename Vec>
    std::vector<uint_t> generate_N_mutations_recycle( queue_type & recycling_bin,
                                  const mutation_model & mmodel,
                                  const unsigned & n,
                                  mlist_type & mutations,
                                  gamete_type & g,
                                  const Vec & ref_v,
                                  Vec & v_of_g,
                                  int relative_index=0 )
    {
      assert(gamete_is_sorted_n(g,mutations));
      assert(gamete_is_sorted_s(g,mutations));
      // create a vector of int that stores all the new mutations created for the current gamete
      std::vector<mut_pos_indx_t> new_muts_obj;
      std::vector<uint_t> new_muts;


      // Vec v_of_g = ref_v; // first copy from ref_v

      auto old_g_muts(g.smutations); // store old mutations already present in g

      // note: omiting below is not a bug!!
      // because v_of_g is not exactly the genotype of current v, but mostly ref_v, and certain bits are flipped if and only if mutations occur at the current generation
      // flipps the previous mutation back which has been considered at the end!

      // // add on 1.24.17
      // // PJ found a bug in the code
      // // v_of_g should be the v of the current genotype before new mutations
      // int ori_index;
      // for(auto iter= old_g_muts.begin(); iter!= old_g_muts.end(); iter++)
      // {
      //   ori_index= std::round( (mutations[*iter].pos-relative_index) * ref_v.size()) ;
      //   v_of_g[ori_index]= 1- v_of_g[ori_index];
      // }


      for( unsigned i = 0 ; i < n ; ++i )
    	{
         auto  idx_obj= mmodel_dispatcher(mmodel,g,mutations,recycling_bin);

         // only add new_muts to the list if there's no double mutation at the same position
         std::vector<CT::mut_pos_indx_t>::iterator p=std::find_if(new_muts_obj.begin(),new_muts_obj.end(),[idx_obj](const mut_pos_indx_t & idx) { return idx_obj.index == idx.index ;});
         if( p == new_muts_obj.cend())
         {
           new_muts_obj.push_back(idx_obj);
         }
         else // already present in the list, remove it
         {
           new_muts_obj.erase(p);
         }
    	}

      for (auto iter= new_muts_obj.cbegin(); iter!= new_muts_obj.cend(); iter++)
      {
        uint_t index = iter->index;
        new_muts.push_back(index);
        add_new_mutation(*iter,mutations,g);
        // fix a bug: add extra step: find positions that has been filpped twice, change
        auto find_iter = std::find(old_g_muts.cbegin(),old_g_muts.cend(),index);
        if (find_iter!=old_g_muts.cend() ) // if there's a mutation arising at this generation that flips back
        {
          int ori_index= std::round( (mutations[*find_iter].pos-relative_index) * ref_v.size()) ;
          v_of_g[ori_index]= 1- v_of_g[ori_index];
        }
      }
      return new_muts;
    }

    // created on 1.24.17
    // for test for two gamete type
    template<typename queue_type,
             typename mutation_model,
             typename mlist_type,
             typename gamete_type1,
             typename gamete_type2,
             typename Vec>
    std::vector<uint_t> generate_N_mutations_recycle( queue_type & recycling_bin,
                                                      const mutation_model & mmodel,
                                                      const unsigned & n,
                                                      mlist_type & mutations,
                                                      gamete_type1 & g1,
                                                      gamete_type2 & g2,
                                                      const Vec & ref_v,
                                                      Vec & v_of_g,
                                                      int relative_index=0 )
    {
      assert(gamete_is_sorted_n(g1,mutations));
      assert(gamete_is_sorted_s(g1,mutations));

      assert(gamete_is_sorted_n(g2,mutations));
      assert(gamete_is_sorted_s(g2,mutations));

      // create a vector of int that stores all the new mutations created for the current gamete
      std::vector<mut_pos_indx_t> new_muts_obj;
      std::vector<uint_t> new_muts;


      // Vec v_of_g = ref_v; // first copy from ref_v

      auto old_g_muts(g1.smutations); // store old mutations already present in old_g_muts

      // note: omiting below is not a bug!!
      // because v_of_g is not exactly the genotype of current v, but mostly ref_v, and certain bits are flipped if and only if mutations occur at the current generation
      // flipps the previous mutation back which has been considered at the end!

      // // add on 1.24.17
      // // PJ found a bug in the code
      // // v_of_g should be the v of the current genotype before new mutations
      // int ori_index;
      // for(auto iter= old_g_muts.begin(); iter!= old_g_muts.end(); iter++)
      // {
      //   ori_index= std::round( (mutations[*iter].pos-relative_index) * ref_v.size()) ;
      //   v_of_g[ori_index]= 1- v_of_g[ori_index];
      // }


      for( unsigned i = 0 ; i < n ; ++i )
      {
         auto  idx_obj= mmodel_dispatcher(mmodel,g1,mutations,recycling_bin);

         // only add new_muts to the list if there's no double mutation at the same position
         std::vector<CT::mut_pos_indx_t>::iterator p=std::find_if(new_muts_obj.begin(),new_muts_obj.end(),[idx_obj](const mut_pos_indx_t & idx) { return idx_obj.index == idx.index ;});
         if( p == new_muts_obj.cend())
         {
           new_muts_obj.push_back(idx_obj);
         }
         else // already present in the list, remove it
         {
           new_muts_obj.erase(p);
         }
      }

      for (auto iter= new_muts_obj.cbegin(); iter!= new_muts_obj.cend(); iter++)
      {
        uint_t index = iter->index;
        new_muts.push_back(index);

        add_new_mutation(*iter,mutations,g1);
        add_new_mutation(*iter,mutations,g2);

        // fix a bug: add extra step: find positions that has been filpped twice, change
        auto find_iter = std::find(old_g_muts.cbegin(),old_g_muts.cend(),index);
        if (find_iter!=old_g_muts.cend() ) // if there's a mutation arising at this generation that flips back
        {
          int ori_index= std::round( (mutations[*find_iter].pos-relative_index) * ref_v.size()) ;
          v_of_g[ori_index]= 1- v_of_g[ori_index];
        }
      }
      return new_muts;
    }






    // modified on 12.17.16
    // will write add_N_mutations_recycle not as single function, but as
    // operator(), so that it can be called by higher level function seemlessly
    struct add_N_mutations_recycle
    {
      // comment out previous version
      // modified on 1.23.17
      // specify gamete type : ct_gamete_with_h
      // modifed on 12.6.16
      // for gamete type : ct_gamete_with_h
      // added robustness flag and robust_degrees
      // modified on 11.16.16, add return value (vector of indicies that new mutations arised on the current gamete)
      // this version will update the wv vector in the gamete when new mutations occur (Note: this is only compatible without recombination!)
      // template<typename queue_type,
  	  //          typename mutation_model,
      // 	       typename mlist_type,
      //          typename Mat,
      //          typename Vec>




      // // created on 1.23.17
      // // for the new gamete type: ct_gamete_scale_wij
      // template<typename queue_type,
  	  //          typename mutation_model,
      // 	       typename mlist_type,
      //          typename Mat,
      //          typename Vec>
      // void operator ()( queue_type & recycling_bin,
      //             const mutation_model & mmodel,
      //             const unsigned & n,
      //             mlist_type & mutations,
      //             ct_gamete_scale_wij_t & g,
      //             const Mat & ref_w,
      //             const Vec & ref_v,
      //             const int & robust_bits,
      //             const int & flag, //robustness flag
      //             const int & robust_degrees) const
      // {
      //   Vec v_of_g = ref_v; // first copy from ref_v
      //
      //   std::vector<uint_t> new_muts= generate_N_mutations_recycle(recycling_bin,mmodel,n,mutations,g, ref_v,v_of_g);
      //   int changed_bits=0;
      //
      //   // note: the generation when robustness changes, it will not reflect changed immediately, but will be reflected in the following generations
      //   // where mutation occur
      //   double ratio = g.get_robust_alpha_ratio();
      //
      //
      //   VectorXd diff_vec= recalc_mat(ref_w, v_of_g , robust_bits, new_muts.cbegin(), new_muts.cend(),mutations,ratio, changed_bits,flag );
      //   g.update_all(diff_vec,changed_bits, robust_degrees, flag);
      // }


      // modified on 2.21.17
      // add scale factor, to store into new mutations
      // modified on 1.30.17
      // add ratio as a parameter
      // modified on 12.20.16
      // add relative_index indicating starting range of the mutations
      // this should be work with gamemet save wv but do not save robust bits
      template<typename queue_type,
  	           typename mutation_model,
      	       typename mlist_type,
      	       typename gamete_type,
               typename Mat,
               typename Vec>
      void operator() ( queue_type & recycling_bin,
  				                          const mutation_model & mmodel,
  				                          const unsigned & n,
  				                          mlist_type & mutations,
  				                          gamete_type & g,
                                    double ratio, // the scaling factor
                                    const Mat & ref_w,
                                    const Vec & ref_v,
                                    int relative_index) const
      {
        Vec v_of_g = ref_v; // first copy from ref_v
        std::vector<uint_t> new_muts= generate_N_mutations_recycle(recycling_bin,mmodel,n,mutations,g, ref_v,v_of_g,relative_index);

        // // for test
        // // print out new muts pos
        // for (auto i=new_muts.cbegin();i!=new_muts.cend();i++)
        //   std::cout<< mutations[*i].pos << std::endl;

        // to set the scales in the new mutations with ratio
        for (auto idx : new_muts)
        {
          mutations[idx].scale = ratio ;
        }


        VectorXd diff_vec= recalc_mat(ref_w, v_of_g, new_muts.cbegin(), new_muts.cend(),mutations,ratio, relative_index );
        g.update_all(diff_vec);
      }


      // modified on 1.3.17
      // fixed a bug for robust bits (if it was a mutation that flips a previous mutation back)
      // should use v_of_g in muts_to_robust instead of ref_g!!
      // added on 12.17.16
      // for mutations just in robust bits
      template<typename queue_type,
  	           typename mutation_model,
      	       typename mlist_type,
      	       typename gamete_type,
               typename Vec>
      void operator ()( queue_type & recycling_bin,
                  const mutation_model & mmodel,
                  const unsigned & n,
                  mlist_type & mutations,
                  gamete_type & g,
                  const Vec & ref_v,
                  const int & robust_bits,
                  const int & flag, //robustness flag
                  const int & robust_degrees) const // important! this should be const! (if passed as a const functor when called!!), but the compiler would not say it explictly,
                                                    // just say error passing ...
      {
        Vec v_of_g = ref_v; // first copy from ref_v
        std::vector<uint_t> new_muts= generate_N_mutations_recycle(recycling_bin,mmodel,n,mutations,g, ref_v,v_of_g);

        // // for test
        // // print out new muts pos
        // for (auto i=new_muts.cbegin();i!=new_muts.cend();i++)
        //   std::cout<< mutations[*i].pos << std::endl;


        int changed_bits= muts_to_robust(new_muts.cbegin(),new_muts.cend(),mutations,v_of_g,robust_bits,flag);
        g.update_robust_geno_flag_and_set(changed_bits, robust_degrees, flag);

      }


    };


  }
}

#endif
