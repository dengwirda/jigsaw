
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
     * Last updated: 30 April, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
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
    typename      real_type
             >
    __inline_call void_type vector_2d (
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __write_ptr  (real_type) _vv
        )
    {
        _vv[0] = _p2[0] - _p1[0];
        _vv[1] = _p2[1] - _p1[1];
    }

    template <
    typename      real_type
             >
    __inline_call void_type vector_3d (
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __write_ptr  (real_type) _vv
        )
    {
        _vv[0] = _p2[0] - _p1[0];
        _vv[1] = _p2[1] - _p1[1];
        _vv[2] = _p2[2] - _p1[2];
    }

    template <
    typename      real_type
             >
    __inline_call void_type vector_4d (
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __write_ptr  (real_type) _vv
        )
    {
        _vv[0] = _p2[0] - _p1[0];
        _vv[1] = _p2[1] - _p1[1];
        _vv[2] = _p2[2] - _p1[2];
        _vv[3] = _p2[3] - _p1[3];
    }

/*-------------------------------- calc. (squared) length */
    template <
    typename      real_type
             >
    __inline_call real_type lensqr_2d (
    __const_ptr  (real_type) _vv
        )
    {   return _vv[0] * _vv[0] +
               _vv[1] * _vv[1] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type lensqr_3d (
    __const_ptr  (real_type) _vv
        )
    {   return _vv[0] * _vv[0] +
               _vv[1] * _vv[1] +
               _vv[2] * _vv[2] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type lensqr_4d (
    __const_ptr  (real_type) _vv
        )
    {   return _vv[0] * _vv[0] +
               _vv[1] * _vv[1] +
               _vv[2] * _vv[2] +
               _vv[3] * _vv[3] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type lensqr_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2
        )
    {
        real_type _vv[ +2];
        vector_2d(_p1, _p2, _vv) ;

        return _vv[0] * _vv[0] +
               _vv[1] * _vv[1] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type lensqr_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2
        )
    {
        real_type _vv[ +3];
        vector_3d(_p1, _p2, _vv) ;

        return _vv[0] * _vv[0] +
               _vv[1] * _vv[1] +
               _vv[2] * _vv[2] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type lensqr_4d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2
        )
    {
        real_type _vv[ +4];
        vector_4d(_p1, _p2, _vv) ;

        return _vv[0] * _vv[0] +
               _vv[1] * _vv[1] +
               _vv[2] * _vv[2] +
               _vv[3] * _vv[3] ;
    }

/*-------------------------------- calc. euclidean length */
    template <
    typename      real_type
             >
    __inline_call real_type length_2d (
    __const_ptr  (real_type) _vv
        )
    {   return std::sqrt(lensqr_2d(_vv)) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type length_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2
        )
    {   return
        std::sqrt(lensqr_2d(_p1, _p2)) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type length_3d (
    __const_ptr  (real_type) _vv
        )
    {   return std::sqrt(lensqr_3d(_vv)) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type length_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2
        )
    {   return
        std::sqrt(lensqr_3d(_p1, _p2)) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type length_4d (
    __const_ptr  (real_type) _vv
        )
    {   return std::sqrt(lensqr_4d(_vv)) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type length_4d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2
        )
    {   return
        std::sqrt(lensqr_4d(_p1, _p2)) ;
    }

/*-------------------------------- scale onto unit vector */
    template <
    typename      real_type
             >
    __inline_call real_type normalise_2d (
    __write_ptr  (real_type) _vv
        )
    {
        real_type _ll = length_2d(_vv) ;
        _vv[0] /= _ll ;
        _vv[1] /= _ll ;

        return (  _ll ) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type normalise_3d (
    __write_ptr  (real_type) _vv
        )
    {
        real_type _ll = length_3d(_vv) ;
        _vv[0] /= _ll ;
        _vv[1] /= _ll ;
        _vv[2] /= _ll ;

        return (  _ll ) ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type normalise_4d (
    __write_ptr  (real_type) _vv
        )
    {
        real_type _ll = length_4d(_vv) ;
        _vv[0] /= _ll ;
        _vv[1] /= _ll ;
        _vv[2] /= _ll ;
        _vv[3] /= _ll ;

        return (  _ll ) ;
    }

/*-------------------------------- vector "cross" product */
    template <
    typename      real_type
             >
    __inline_call void_type cross_2d (
    __const_ptr  (real_type) _v1 ,
    __const_ptr  (real_type) _v2 ,
        real_type &_cp
        )
    {   _cp    = _v1[0] * _v2[1] -
                 _v1[1] * _v2[0] ;
    }

    template <
    typename      real_type
             >
    __inline_call void_type cross_3d (
    __const_ptr  (real_type) _v1 ,
    __const_ptr  (real_type) _v2 ,
    __write_ptr  (real_type) _cp
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
    typename      real_type
             >
    __inline_call real_type dot_2d (
    __const_ptr  (real_type) _v1 ,
    __const_ptr  (real_type) _v2
        )
    {   return   _v1[0] * _v2[0] +
                 _v1[1] * _v2[1] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type dot_3d (
    __const_ptr  (real_type) _v1 ,
    __const_ptr  (real_type) _v2
        )
    {   return   _v1[0] * _v2[0] +
                 _v1[1] * _v2[1] +
                 _v1[2] * _v2[2] ;
    }

    template <
    typename      real_type
             >
    __inline_call real_type dot_4d (
    __const_ptr  (real_type) _v1 ,
    __const_ptr  (real_type) _v2
        )
    {   return   _v1[0] * _v2[0] +
                 _v1[1] * _v2[1] +
                 _v1[2] * _v2[2] +
                 _v1[3] * _v2[3] ;
    }

/*-------------------------------- cosine between vectors */
    template <
    typename      real_type
             >
    __inline_call real_type cosine_2d (
    __const_ptr  (real_type) _v1,
    __const_ptr  (real_type) _v2
        )
    {
        real_type _ll = length_2d(_v1) *
                        length_2d(_v2) ;

        return dot_2d(_v1, _v2) / _ll;
    }

    template <
    typename      real_type
             >
    __inline_call real_type cosine_3d (
    __const_ptr  (real_type) _v1,
    __const_ptr  (real_type) _v2
        )
    {
        real_type _ll = length_3d(_v1) *
                        length_3d(_v2) ;

        return dot_3d(_v1, _v2) / _ll;
    }


    }

#   endif   //__VECT_BASE_K__



