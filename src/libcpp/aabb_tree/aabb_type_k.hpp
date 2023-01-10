
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
 * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant
 * this code in any way whatsoever.  This code is
 * provided "as-is" to be used at your own risk.
 *
 * THE CONTRIBUTORS include:
 * (a) The University of Sydney
 * (b) The Massachusetts Institute of Technology
 * (c) Columbia University
 * (d) The National Aeronautics & Space Administration
 * (e) Los Alamos National Laboratory
 *
------------------------------------------------------------
 *
 * Last updated: 10 Sept., 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * d.engwirda@gmail.com
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
    typename R,  // reals typedef
    typename I,  // integer typedef
    size_t   K   // number of geometric dimensions
             >
    class aabb_item_rect_k
    {
/*------------ k-dim. "rectangle" item type for AABB-tree */
    public  :

    typedef R                   real_type;
    typedef I                   iptr_type;

    iptr_type static constexpr _dims = K ;
    bool_type static constexpr _is_a_rect = true ;

    private :

    iptr_type                  _ipos ;

    real_type                  _pmin [ K];
    real_type                  _pmax [ K];

    public  :
/*---------------------------------------- "write" access */
    __inline_call real_type& pmin (
        iptr_type _ppos
        )
    {   return this->_pmin[ _ppos] ;
    }
    __inline_call real_type& pmax (
        iptr_type _ppos
        )
    {   return this->_pmax[ _ppos] ;
    }
    __inline_call iptr_type& ipos (
        )
    {   return this->_ipos;
    }
/*---------------------------------------- "const" access */
    __inline_call real_type const& pmin (
        iptr_type _ppos
        ) const
    {   return this->_pmin[ _ppos] ;
    }
    __inline_call real_type const& pmax (
        iptr_type _ppos
        ) const
    {   return this->_pmax[ _ppos] ;
    }
    __inline_call iptr_type const& ipos (
        ) const
    {   return this->_ipos;
    }
/*--------------------------- calc. centroid, length, etc */
    __inline_call real_type  pmid (
        iptr_type _ppos
        ) const
    {   return(this->_pmin[ _ppos] +
               this->_pmax[ _ppos])*
              (real_type) +.5 ;
    }
    __inline_call real_type  plen (
        iptr_type _ppos
        ) const
    {   return this->_pmax[ _ppos] -
               this->_pmin[ _ppos] ;
    }

    } ;


    template <
    typename R,  // reals typedef
    typename I,  // integer typedef
    size_t   K   // number of geometric dimensions
             >
    class aabb_item_node_k
    {
/*---------------- k-dim. "point" item type for AABB-tree */
    public  :

    typedef R                   real_type;
    typedef I                   iptr_type;

    iptr_type static constexpr _dims = K ;
    bool_type static constexpr _is_a_rect = false;

    private :

    iptr_type                  _ipos ;

    real_type                  _pval [ K];

    public  :
/*---------------------------------------- "write" access */
    __inline_call real_type& pval (
        iptr_type _ppos
        )
    {   return this->_pval[ _ppos] ;
    }
    __inline_call real_type& pmin (
        iptr_type _ppos
        )
    {   return this->_pval[ _ppos] ;
    }
    __inline_call real_type& pmax (
        iptr_type _ppos
        )
    {   return this->_pval[ _ppos] ;
    }
    __inline_call iptr_type& ipos (
        )
    {   return this->_ipos  ;
    }
/*---------------------------------------- "const" access */
    __inline_call real_type const& pval (
        iptr_type _ppos
        ) const
    {   return this->_pval[ _ppos] ;
    }
    __inline_call real_type const& pmin (
        iptr_type _ppos
        ) const
    {   return this->_pval[ _ppos] ;
    }
    __inline_call real_type const& pmax (
        iptr_type _ppos
        ) const
    {   return this->_pval[ _ppos] ;
    }
    __inline_call iptr_type const& ipos (
        ) const
    {   return this->_ipos;
    }
/*--------------------------- calc. centroid, length, etc */
    __inline_call real_type const& pmid (
        iptr_type _ppos
        ) const
    {   return this->_pval[ _ppos] ;
    }
    __inline_call real_type  plen (
        iptr_type
        ) const
    {   return   (real_type)+0.E+0 ;
    }

    } ;


    }

#   endif   //__AABB_TYPE_K__



