
    /*
    --------------------------------------------------------
     * VECT-BASE-K: basic vector operations in R^2/R^3.
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
     * Last updated: 06 June, 2016
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
#   pragma once

#   ifndef __VECT_BASE_K__
#   define __VECT_BASE_K__

    namespace geometry {
 
/*-------------------------------- form vector from coord */
    template <
    typename      data_type
             > 
    __inline_call void_type vector_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __write_ptr  (data_type) _vv
        )
    {
        _vv[0] = _p2[0] - _p1[0];
        _vv[1] = _p2[1] - _p1[1];
    }

    template <
    typename      data_type
             > 
    __inline_call void_type vector_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __write_ptr  (data_type) _vv
        )
    {
        _vv[0] = _p2[0] - _p1[0];
        _vv[1] = _p2[1] - _p1[1];
        _vv[2] = _p2[2] - _p1[2];
    }

/*-------------------------------- calc. (squared) length */
    template <
    typename      data_type
             > 
    __inline_call data_type lensqr_2d (
    __const_ptr  (data_type) _vv
        )
    {   return _vv[0] * _vv[0] + 
               _vv[1] * _vv[1] ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type lensqr_3d (
    __const_ptr  (data_type) _vv
        )
    {   return _vv[0] * _vv[0] + 
               _vv[1] * _vv[1] +
               _vv[2] * _vv[2] ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type lensqr_2d (
    __const_ptr  (data_type) _p1 ,
    __const_ptr  (data_type) _p2
        )
    {
        data_type _vv[ +2];
        vector_2d(_p1, _p2, _vv) ;

        return _vv[0] * _vv[0] + 
               _vv[1] * _vv[1] ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type lensqr_3d (
    __const_ptr  (data_type) _p1 ,
    __const_ptr  (data_type) _p2
        )
    {
        data_type _vv[ +3];
        vector_3d(_p1, _p2, _vv) ;

        return _vv[0] * _vv[0] + 
               _vv[1] * _vv[1] +
               _vv[2] * _vv[2] ;
    }

/*-------------------------------- calc. euclidean length */
    template <
    typename      data_type
             > 
    __inline_call data_type length_2d (
    __const_ptr  (data_type) _vv
        )
    {   return std::sqrt(lensqr_2d(_vv)) ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type length_2d (
    __const_ptr  (data_type) _p1 ,
    __const_ptr  (data_type) _p2
        )
    {   return 
        std::sqrt(lensqr_2d(_p1, _p2)) ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type length_3d (
    __const_ptr  (data_type) _vv
        )
    {   return std::sqrt(lensqr_3d(_vv)) ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type length_3d (
    __const_ptr  (data_type) _p1 ,
    __const_ptr  (data_type) _p2
        )
    {   return 
        std::sqrt(lensqr_3d(_p1, _p2)) ;
    }

/*-------------------------------- vector "cross" product */
    template <
    typename      data_type
             > 
    __inline_call void_type cross_2d (
    __const_ptr  (data_type) _v1 ,
    __const_ptr  (data_type) _v2 ,
        data_type &_cp
        )
    {   _cp    = _v1[0] * _v2[1] - 
                 _v1[1] * _v2[0] ;
    }

    template <
    typename      data_type
             > 
    __inline_call void_type cross_3d (
    __const_ptr  (data_type) _v1 ,
    __const_ptr  (data_type) _v2 ,
    __write_ptr  (data_type) _cp
        )
    {
        _cp[0] = _v1[1] * _v2[2] - 
                 _v1[2] * _v2[1] ;
        _cp[1] = _v1[2] * _v2[0] - 
                 _v1[0] * _v2[2] ;
        _cp[2] = _v1[0] * _v2[1] - 
                 _v1[1] * _v2[0] ;
    }

/*-------------------------------- vector "inner" product */
    template <
    typename      data_type
             > 
    __inline_call data_type dot_2d (
    __const_ptr  (data_type) _v1 ,
    __const_ptr  (data_type) _v2
        )
    {   return   _v1[0] * _v2[0] + 
                 _v1[1] * _v2[1] ;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type dot_3d (
    __const_ptr  (data_type) _v1 ,
    __const_ptr  (data_type) _v2
        )
    {   return   _v1[0] * _v2[0] + 
                 _v1[1] * _v2[1] + 
                 _v1[2] * _v2[2] ;
    }

/*-------------------------------- cosine between vectors */
    template <
    typename      data_type
             > 
    __inline_call data_type cosine_2d (
    __const_ptr  (data_type) _v1,
    __const_ptr  (data_type) _v2
        )
    {
        data_type _ll = length_2d(_v1) * 
                        length_2d(_v2) ;

        return dot_2d(_v1, _v2) / _ll;
    }

    template <
    typename      data_type
             > 
    __inline_call data_type cosine_3d (
    __const_ptr  (data_type) _v1,
    __const_ptr  (data_type) _v2
        )
    {
        data_type _ll = length_3d(_v1) * 
                        length_3d(_v2) ;

        return dot_3d(_v1, _v2) / _ll;
    }


    }

#   endif   //__VECT_BASE_K__



