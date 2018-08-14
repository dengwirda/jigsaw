
/* 
------------------------------------------------------------
 * AABB-TYPE-K: data-types for AABB-trees in R^k. 
------------------------------------------------------------
 *
 * This program may be freely redistributed under the 
 * condition that the copyright notices (including this 
 * entire header) are not removed, and no compensation 
 * is received through use of the software.  Private, 
 * research, and institutional use is free.  You may 
 * distribute modified versions of this code UNDER THE 
 * CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
 * TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
 * ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
 * MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
 * NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
 * of this code as part of a commercial system is 
 * permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
 * AUTHOR.  (If you are not directly supplying this 
 * code to a customer, and you are instead telling them 
 * how they can obtain it for free, then you are not 
 * required to make any arrangement with me.) 
 *
 * Disclaimer:  Neither I nor: Columbia University, The
 * Massachusetts Institute of Technology, The 
 * University of Sydney, nor The National Aeronautics
 * and Space Administration warrant this code in any 
 * way whatsoever.  This code is provided "as-is" to be 
 * used at your own risk.
 *
------------------------------------------------------------
 *
 * Last updated: 10 September, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __AABB_TYPE_K__
#   define __AABB_TYPE_K__

    namespace geom_tree {

    
    class aabb_node_base_k
    {
/*------------ base-kind for user-defined aabb-node types */
    public  :
    } ;


    template <
    typename R,
    typename I,
    size_t   K
             >
    class aabb_item_rect_k
    { 
/*------------ k-dim. "rectangle" item type for AABB-tree */
    public  :
    
    typedef R                   real_type;
    typedef I                   iptr_type;
    
    iptr_type static constexpr _dims = K ;
    
    private :
    
    iptr_type                  _ipos ;
    
    real_type                  _pmin [ K];
    real_type                  _pmax [ K];

    public  :
/*---------------------------------------- "write" access */
    __inline_call real_type& pmin (
        iptr_type _ipos
        )
    {   return this->_pmin[ _ipos] ;
    }
    __inline_call real_type& pmax (
        iptr_type _ipos
        )
    {   return this->_pmax[ _ipos] ;
    }
    __inline_call iptr_type& ipos (
        )
    {   return this->_ipos;
    }
/*---------------------------------------- "const" access */
    __inline_call real_type const& pmin (
        iptr_type _ipos
        ) const
    {   return this->_pmin[ _ipos] ;
    }
    __inline_call real_type const& pmax (
        iptr_type _ipos
        ) const
    {   return this->_pmax[ _ipos] ;
    }
    __inline_call iptr_type const& ipos (
        ) const
    {   return this->_ipos;
    }
/*--------------------------- calc. centroid, length, etc */
    __inline_call real_type  pmid (
        iptr_type _ipos
        ) const
    {   return(this->_pmin[ _ipos] + 
               this->_pmax[ _ipos])*
              (real_type) +.5 ;
    }
    __inline_call real_type  plen (
        iptr_type _ipos
        ) const
    {   return this->_pmax[ _ipos] -
               this->_pmin[ _ipos] ;
    }
    
    } ;
    
    
    template <
    typename R,
    typename I,
    size_t   K
             >
    class aabb_item_node_k
    { 
/*---------------- k-dim. "point" item type for AABB-tree */
    public  :
    
    typedef R                   real_type;
    typedef I                   iptr_type;
    
    iptr_type static constexpr _dims = K ;
    
    private :
    
    iptr_type                  _ipos ;

    real_type                  _pval [ K];

    public  :
/*---------------------------------------- "write" access */
    __inline_call real_type& pval (
        iptr_type _ipos
        )
    {   return this->_pval[ _ipos] ;
    }
    __inline_call real_type& pmin (
        iptr_type _ipos
        )
    {   return this->_pval[ _ipos] ;
    }
    __inline_call real_type& pmax (
        iptr_type _ipos
        )
    {   return this->_pval[ _ipos] ;
    }
    __inline_call iptr_type& ipos (
        )
    {   return this->_ipos  ;
    }
/*---------------------------------------- "const" access */
    __inline_call real_type const& pval (
        iptr_type _ipos
        ) const
    {   return this->_pval[ _ipos] ;
    }
    __inline_call real_type const& pmin (
        iptr_type _ipos
        ) const
    {   return this->_pval[ _ipos] ;
    }
    __inline_call real_type const& pmax (
        iptr_type _ipos
        ) const
    {   return this->_pval[ _ipos] ;
    }
    __inline_call iptr_type const& ipos (
        ) const
    {   return this->_ipos;
    }
/*--------------------------- calc. centroid, length, etc */
    __inline_call real_type const& pmid (
        iptr_type _ipos
        ) const
    {   return this->_pval[ _ipos] ;
    }
    __inline_call real_type  plen (
        iptr_type
        ) const
    {   return   (real_type) +0.00 ;
    }
    
    } ;
    
    
    }
    
#   endif   //__AABB_TYPE_K__
    
    
    
