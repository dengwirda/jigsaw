
    /*
    --------------------------------------------------------
     * MATRIX-UTIL: (small) matrix utilities. 
    --------------------------------------------------------
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
    --------------------------------------------------------
     *
     * Last updated: 07 January, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MATRIX_UTIL__
#   define __MATRIX_UTIL__

    namespace geometry
    {

    /*
    --------------------------------------------------------
     * small matrix utilities. 
    --------------------------------------------------------
     */
      
#   define __ij(__ir , __ic , __nr) \
              ((__ir)+(__ic)*(__nr))

    template <
    typename      real_type, 
    typename      size_type
             >
    __inline_call real_type det_2x2 (
        size_type  _la,
    __const_ptr  (real_type) _aa
        )
    {   return 
        _aa[__ij(0,0,_la)] *
        _aa[__ij(1,1,_la)] -
        _aa[__ij(0,1,_la)] *
        _aa[__ij(1,0,_la)] ;
    }

    template <
    typename      real_type, 
    typename      size_type
             >
    __inline_call real_type det_3x3 (
        size_type  _la,
    __const_ptr  (real_type) _aa
        )
    {   return 
        _aa[__ij(0,0,_la)] * (
        _aa[__ij(1,1,_la)] * 
        _aa[__ij(2,2,_la)] - 
        _aa[__ij(1,2,_la)] * 
        _aa[__ij(2,1,_la)] ) -

        _aa[__ij(0,1,_la)] * (
        _aa[__ij(1,0,_la)] * 
        _aa[__ij(2,2,_la)] - 
        _aa[__ij(1,2,_la)] * 
        _aa[__ij(2,0,_la)] ) +

        _aa[__ij(0,2,_la)] * (
        _aa[__ij(1,0,_la)] * 
        _aa[__ij(2,1,_la)] - 
        _aa[__ij(1,1,_la)] * 
        _aa[__ij(2,0,_la)] ) ;
    }

    template <
    typename      real_type,
    typename      size_type 
             >
    __inline_call void_type inv_2x2 (
        size_type  _la, 
    __const_ptr  (real_type) _aa,
        size_type  _lx,
    __write_ptr  (real_type) _xx,
        real_type &_da
        )
    {
        _da = det_2x2(_la, _aa) ;
    
        _xx[__ij(0,0,_lx)] = 
        _aa[__ij(1,1,_la)] ;
        _xx[__ij(1,1,_lx)] = 
        _aa[__ij(0,0,_la)] ;
        _xx[__ij(0,1,_lx)] = 
       -_aa[__ij(0,1,_la)] ;
        _xx[__ij(1,0,_lx)] = 
       -_aa[__ij(1,0,_la)] ;
    }

    template <
    typename      real_type,
    typename      size_type 
             >
    __inline_call void_type inv_3x3 (
        size_type  _la, 
    __const_ptr  (real_type) _aa,
        size_type  _lx,
    __write_ptr  (real_type) _xx,
        real_type &_da
        )
    {
        _da = det_3x3(_la, _aa) ;
    
        _xx[__ij(0,0,_lx)] =
        _aa[__ij(2,2,_la)] * 
        _aa[__ij(1,1,_la)] - 
        _aa[__ij(2,1,_la)] * 
        _aa[__ij(1,2,_la)] ;
    
        _xx[__ij(0,1,_lx)] =
        _aa[__ij(2,1,_la)] * 
        _aa[__ij(0,2,_la)] - 
        _aa[__ij(2,2,_la)] * 
        _aa[__ij(0,1,_la)] ;

        _xx[__ij(0,2,_lx)] =
        _aa[__ij(1,2,_la)] * 
        _aa[__ij(0,1,_la)] - 
        _aa[__ij(1,1,_la)] * 
        _aa[__ij(0,2,_la)] ;

        _xx[__ij(1,0,_lx)] =
        _aa[__ij(2,0,_la)] * 
        _aa[__ij(1,2,_la)] - 
        _aa[__ij(2,2,_la)] * 
        _aa[__ij(1,0,_la)] ;
    
        _xx[__ij(1,1,_lx)] =
        _aa[__ij(2,2,_la)] * 
        _aa[__ij(0,0,_la)] - 
        _aa[__ij(2,0,_la)] * 
        _aa[__ij(0,2,_la)] ;
    
        _xx[__ij(1,2,_lx)] =
        _aa[__ij(1,0,_la)] * 
        _aa[__ij(0,2,_la)] - 
        _aa[__ij(1,2,_la)] * 
        _aa[__ij(0,0,_la)] ;

        _xx[__ij(2,0,_lx)] =
        _aa[__ij(2,1,_la)] * 
        _aa[__ij(1,0,_la)] - 
        _aa[__ij(2,0,_la)] * 
        _aa[__ij(1,1,_la)] ;

        _xx[__ij(2,1,_lx)] =
        _aa[__ij(2,0,_la)] * 
        _aa[__ij(0,1,_la)] - 
        _aa[__ij(2,1,_la)] * 
        _aa[__ij(0,0,_la)] ;
    
        _xx[__ij(2,2,_lx)] =
        _aa[__ij(1,1,_la)] * 
        _aa[__ij(0,0,_la)] - 
        _aa[__ij(1,0,_la)] * 
        _aa[__ij(0,1,_la)] ;
    }


    }

#   endif//__MATRIX_UTIL__


