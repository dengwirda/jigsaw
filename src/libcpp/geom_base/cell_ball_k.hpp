
    /*
    --------------------------------------------------------
     * CELL-BALL-K: various circumscribing ball calc.'s.
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
     * Last updated: 26 July, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __CELL_BALL_K__
#   define __CELL_BALL_K__

    namespace geometry {

    /*
    --------------------------------------------------------
     * MINI-BALL: min. enclosing balls.
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_2d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_3d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_2d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_3d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_3d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  bool_type  _in = false
        ) ;

    #define __circface12(__cc, __pi, __pj, __pq) \
        do {                                \
    /*---------- calc. face circumball */   \
        circ_ball_2d(__bf, __pi, __pj, true);    \
    /*---------- keep if min-enclosing */   \
        if (__bf[2] >= lensqr_2d(__bf, __pq))    \
        {                                   \
            if (__bf[2] < __cc[2])          \
            {                               \
                __cc[0] = __bf[0];          \
                __cc[1] = __bf[1];          \
                __cc[2] = __bf[2];          \
            }                               \
        }                                   \
        } while (false) ;

    #define __circface13(__cc, __pi, __pj, __pq) \
        do {                                \
    /*---------- calc. face circumball */   \
        circ_ball_3d(__bf, __pi, __pj, true);    \
    /*---------- keep if min-enclosing */   \
        if (__bf[3] >= lensqr_3d(__bf, __pq))    \
        {                                   \
            if (__bf[3] < __cc[3])          \
            {                               \
                __cc[0] = __bf[0];          \
                __cc[1] = __bf[1];          \
                __cc[2] = __bf[2];          \
                __cc[3] = __bf[3];          \
            }                               \
        }                                   \
        } while (false) ;

    #define __circface23(__cc, __pi, __pj, __pk, __pq) \
        do {                                \
    /*---------- calc. face circumball */   \
        circ_ball_3d(__bf, __pi, __pj, __pk, true); \
    /*---------- keep if min-enclosing */   \
        if (__bf[3] >= lensqr_3d(__bf, __pq))       \
        {                                   \
            if (__bf[3] < __cc[3])          \
            {                               \
                __cc[0] = __bf[0];          \
                __cc[1] = __bf[1];          \
                __cc[2] = __bf[2];          \
                __cc[3] = __bf[3];          \
            }                               \
        }                                   \
        } while (false) ;

    /*
    --------------------------------------------------------
     * CIRC-BALL: circumscribing balls.
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __inline_call void_type tria_norm_3d (  // forward dec's
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3,
    __write_ptr  (real_type) _nv
         ) ;

    template <
    typename      real_type
             >
    __inline_call void_type circ_ball_2d (
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in
        )
    {
        _cc[0] =
        (real_type).5*(_p1[0]+_p2[0]);
        _cc[1] =
        (real_type).5*(_p1[1]+_p2[1]);

        __unreferenced(_in);

        _cc[2] = (real_type)+0. ;
        _cc[2] = std::max (
        _cc[2] , lensqr_2d(_cc, _p1));
        _cc[2] = std::max (
        _cc[2] , lensqr_2d(_cc, _p2));
    }

    template <
    typename      real_type
             >
    __inline_call void_type circ_ball_3d (
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in
        )
    {
        _cc[0] =
        (real_type).5*(_p1[0]+_p2[0]);
        _cc[1] =
        (real_type).5*(_p1[1]+_p2[1]);
        _cc[2] =
        (real_type).5*(_p1[2]+_p2[2]);

        __unreferenced(_in);

        _cc[3] = (real_type)+0. ;
        _cc[3] = std::max (
        _cc[3] , lensqr_3d(_cc, _p1));
        _cc[3] = std::max (
        _cc[3] , lensqr_3d(_cc, _p2));
    }

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_2d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in
        )
    {
        real_type _xm[2*2] ;
        real_type _xi[2*2] ;
        _xm[__ij(0,0,2)] = _p2[0]-_p1[0] ;
        _xm[__ij(0,1,2)] = _p2[1]-_p1[1] ;

        _xm[__ij(1,0,2)] = _p3[0]-_p1[0] ;
        _xm[__ij(1,1,2)] = _p3[1]-_p1[1] ;

        real_type _xr[2*1] ;
        _xr[0] = (real_type)+.5 * (
        std::pow(_xm[__ij(0,0,2)], +2) +
        std::pow(_xm[__ij(0,1,2)], +2)
            ) ;

        _xr[1] = (real_type)+.5 * (
        std::pow(_xm[__ij(1,0,2)], +2) +
        std::pow(_xm[__ij(1,1,2)], +2)
            ) ;

        real_type _dd ;
        math::inv_2x2 (
            +2, _xm, +2, _xi, _dd )  ;
        _bb[0] = (
        _xi[__ij(0,0,2)] * _xr[0] +
        _xi[__ij(0,1,2)] * _xr[1] )  ;

        _bb[1] = (
        _xi[__ij(1,0,2)] * _xr[0] +
        _xi[__ij(1,1,2)] * _xr[1] )  ;

        double   _P1[2];
        copy_node_2d (_P1, _p1) ;

        double   _P2[2];
        copy_node_2d (_P2, _p2) ;

        double   _P3[2];
        copy_node_2d (_P3, _p3) ;

        _dd = +geompred::orient2d (
            _P1, _P2, _P3) ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;

        real_type _ee[2*1];
        real_type _db[2*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,2)] * _bb[0] +
        _xm[__ij(0,1,2)] * _bb[1] )  ;

        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,2)] * _bb[0] +
        _xm[__ij(1,1,2)] * _bb[1] )  ;

        _db[0] = (
        _xi[__ij(0,0,2)] * _ee[0] +
        _xi[__ij(0,1,2)] * _ee[1] )  ;

        _db[1] = (
        _xi[__ij(1,0,2)] * _ee[0] +
        _xi[__ij(1,1,2)] * _ee[1] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        }

        _bb[ 2]  = lensqr_2d(_bb) ;

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;

        if (!_in)  return  ;

        real_type __bf[3]  ;
        __circface12( _bb, _p1, _p2, _p3);
        __circface12( _bb, _p2, _p3, _p1);
        __circface12( _bb, _p3, _p1, _p2);
    }

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_3d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in
        )
    {
        real_type _nv[4*1] ;
        tria_norm_3d(_p1, _p2, _p3, _nv) ;

        real_type _xm[3*3] ;
        real_type _xi[3*3] ;
        _xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
        _xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;

        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
        _xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;

        _xm[__ij(2,0,3)] = _nv[0] ;
        _xm[__ij(2,1,3)] = _nv[1] ;
        _xm[__ij(2,2,3)] = _nv[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        std::pow(_xm[__ij(0,0,3)], +2) +
        std::pow(_xm[__ij(0,1,3)], +2) +
        std::pow(_xm[__ij(0,2,3)], +2)
            ) ;

        _xr[1] = (real_type)+.5 * (
        std::pow(_xm[__ij(1,0,3)], +2) +
        std::pow(_xm[__ij(1,1,3)], +2) +
        std::pow(_xm[__ij(1,2,3)], +2)
            ) ;

        _xr[2] = (real_type)+.0 ;

        real_type _dd ;
        math::inv_3x3 (
            +3, _xm, +3, _xi, _dd )  ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
        _xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;

        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
        _xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
        _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;

        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;

        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;

        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
        _xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;

        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
        _xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
        _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;
        }

        _bb[ 3]  = lensqr_3d(_bb) ;

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;

        if (!_in)  return  ;

        real_type __bf[4]  ;
        __circface13( _bb, _p1, _p2, _p3);
        __circface13( _bb, _p2, _p3, _p1);
        __circface13( _bb, _p3, _p1, _p2);
    }

    template <
    typename      real_type
             >
    __inline_call void_type quad_ball_2d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  bool_type  _in
        )
    {
        __unreferenced(_in);

        real_type _b1[3];
        circ_ball_2d(
            _b1, _p1, _p2, _p3, false) ;
        real_type _b2[3];
        circ_ball_2d(
            _b2, _p1, _p3, _p4, false) ;
        real_type _b3[3];
        circ_ball_2d(
            _b3, _p1, _p2, _p4, false) ;
        real_type _b4[3];
        circ_ball_2d(
            _b4, _p2, _p3, _p4, false) ;

        _bb[0] =
       (_b1[0]+_b2[0]+_b3[0]+_b4[0]) / +4. ;
        _bb[1] =
       (_b1[1]+_b2[1]+_b3[1]+_b4[1]) / +4. ;

        _bb[2] = std::pow(
        _b1[2]*_b2[2]*_b3[2]*_b4[2], 1./4.) ;
    }

    template <
    typename      real_type
             >
    __inline_call void_type quad_ball_3d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  bool_type  _in
        )
    {
        __unreferenced(_in);

        real_type _b1[4];
        circ_ball_3d(
            _b1, _p1, _p2, _p3, false) ;
        real_type _b2[4];
        circ_ball_3d(
            _b2, _p1, _p3, _p4, false) ;
        real_type _b3[4];
        circ_ball_3d(
            _b3, _p1, _p2, _p4, false) ;
        real_type _b4[4];
        circ_ball_3d(
            _b4, _p2, _p3, _p4, false) ;

        _bb[0] =
       (_b1[0]+_b2[0]+_b3[0]+_b4[0]) / +4. ;
        _bb[1] =
       (_b1[1]+_b2[1]+_b3[1]+_b4[1]) / +4. ;
        _bb[2] =
       (_b1[2]+_b2[2]+_b3[2]+_b4[2]) / +4. ;

        _bb[3] = std::pow(
        _b1[3]*_b2[3]*_b3[3]*_b4[3], 1./4.) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type circ_ball_3d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  bool_type  _in
        )
    {
        real_type _xm[3*3] ;
        real_type _xi[3*3] ;
        _xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
        _xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;

        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
        _xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;

        _xm[__ij(2,0,3)] = _p4[0]-_p1[0] ;
        _xm[__ij(2,1,3)] = _p4[1]-_p1[1] ;
        _xm[__ij(2,2,3)] = _p4[2]-_p1[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        std::pow(_xm[__ij(0,0,3)], +2) +
        std::pow(_xm[__ij(0,1,3)], +2) +
        std::pow(_xm[__ij(0,2,3)], +2)
            ) ;

        _xr[1] = (real_type)+.5 * (
        std::pow(_xm[__ij(1,0,3)], +2) +
        std::pow(_xm[__ij(1,1,3)], +2) +
        std::pow(_xm[__ij(1,2,3)], +2)
            ) ;

        _xr[2] = (real_type)+.5 * (
        std::pow(_xm[__ij(2,0,3)], +2) +
        std::pow(_xm[__ij(2,1,3)], +2) +
        std::pow(_xm[__ij(2,2,3)], +2)
            ) ;

        real_type _dd ;
        math::inv_3x3 (
            +3, _xm, +3, _xi, _dd )  ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
        _xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;

        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
        _xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
        _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        double   _P1[3];
        copy_node_3d (_P1, _p1) ;

        double   _P2[3];
        copy_node_3d (_P2, _p2) ;

        double   _P3[3];
        copy_node_3d (_P3, _p3) ;

        double   _P4[3];
        copy_node_3d (_P4, _p4) ;

        _dd = -geompred::orient3d (
            _P1, _P2, _P3, _P4) ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;

        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;

        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;

        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
        _xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;

        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
        _xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
        _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;
        }

        _bb[ 3]  = lensqr_3d(_bb) ;

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;

        if (!_in) return   ;

        real_type __bf[4]  ;
        __circface23(
            _bb, _p1, _p2, _p3, _p4) ;
        __circface23(
            _bb, _p1, _p2, _p4, _p3) ;
        __circface23(
            _bb, _p2, _p3, _p4, _p1) ;
        __circface23(
            _bb, _p3, _p1, _p4, _p2) ;
    }

    #undef  __circface12
    #undef  __circface13
    #undef  __circface23

    /*
    --------------------------------------------------------
     * MINI-BALL: min. enclosing balls.
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_2d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_3d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_2d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_3d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in = false
        ) ;

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_3d (  // forward dec's
    __write_ptr  (real_type) _cc ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  bool_type  _in = false
        ) ;

    #define __perpface12(__cc, __pi, __pj, __pq)    \
        do {                                \
    /*---------- calc. face circumball */   \
        perp_ball_2d(__bf, __pi, __pj, true);       \
    /*---------- keep if min-enclosing */   \
        if (__bf[2] >=                      \
                lensqr_2d(__bf, __pq) - __pq[2])    \
        {                                   \
            if (__bf[2] < __cc[2])          \
            {                               \
                __cc[0] = __bf[0];          \
                __cc[1] = __bf[1];          \
                __cc[2] = __bf[2];          \
            }                               \
        }                                   \
        } while (false) ;

    #define __perpface13(__cc, __pi, __pj, __pq)    \
        do {                                \
    /*---------- calc. face circumball */   \
        perp_ball_3d(__bf, __pi, __pj, true);       \
    /*---------- keep if min-enclosing */   \
        if (__bf[3] >=                      \
                lensqr_3d(__bf, __pq) - __pq[3])    \
        {                                   \
            if (__bf[3] < __cc[3])          \
            {                               \
                __cc[0] = __bf[0];          \
                __cc[1] = __bf[1];          \
                __cc[2] = __bf[2];          \
                __cc[3] = __bf[3];          \
            }                               \
        }                                   \
        } while (false) ;

    #define __perpface23(__cc, __pi, __pj, __pk, __pq)  \
        do {                                \
    /*---------- calc. face circumball */   \
        perp_ball_3d(__bf, __pi, __pj, __pk, true); \
    /*---------- keep if min-enclosing */   \
        if (__bf[3] >=                      \
                lensqr_3d(__bf, __pq) - __pq[3])    \
        {                                   \
            if (__bf[3] < __cc[3])          \
            {                               \
                __cc[0] = __bf[0];          \
                __cc[1] = __bf[1];          \
                __cc[2] = __bf[2];          \
                __cc[3] = __bf[3];          \
            }                               \
        }                                   \
        } while (false) ;

    /*
    --------------------------------------------------------
     * PERP-BALL: perpendicular (orthogonal) balls.
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_2d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in
        )
    {
        real_type _dd[3*1] ;
        _dd[0] = _p1[0] - _p2[0] ;
        _dd[1] = _p1[1] - _p2[1] ;
        _dd[2] = _p1[2] - _p2[2] ;

        real_type _dp =
            geometry::lensqr_2d(_dd) ;

        real_type _tt =
            (real_type).5 *
                (_dd[2] + _dp) / _dp ;

        if (_in == true)
        {
            _tt  =
        std::min((real_type)1., _tt) ;
            _tt  =
        std::max((real_type)0., _tt) ;
        }

        _bb[0] =
        _p1[0] -  _tt * _dd[0] ;
        _bb[1] =
        _p1[1] -  _tt * _dd[1] ;

        real_type _r1 =
        geometry::lensqr_2d(_p1, _bb);
        real_type _r2 =
        geometry::lensqr_2d(_p2, _bb);

        _r1 -= _p1[2] ;
        _r2 -= _p2[2] ;

        _bb[2] = std::max (_r1, _r2) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_3d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  bool_type  _in
        )
    {
        real_type _dd[4*1] ;
        _dd[0] = _p1[0] - _p2[0] ;
        _dd[1] = _p1[1] - _p2[1] ;
        _dd[2] = _p1[2] - _p2[2] ;
        _dd[3] = _p1[3] - _p2[3] ;

        real_type _dp =
            geometry::lensqr_3d(_dd) ;

        real_type _tt =
            (real_type).5 *
                (_dd[3] + _dp) / _dp ;

        if (_in == true)
        {
            _tt  =
        std::min((real_type)1., _tt) ;
            _tt  =
        std::max((real_type)0., _tt) ;
        }

        _bb[0] =
        _p1[0] -  _tt * _dd[0] ;
        _bb[1] =
        _p1[1] -  _tt * _dd[1] ;
        _bb[2] =
        _p1[2] -  _tt * _dd[2] ;

        real_type _r1 =
        geometry::lensqr_3d(_p1, _bb);
        real_type _r2 =
        geometry::lensqr_3d(_p2, _bb);

        _r1 -= _p1[3] ;
        _r2 -= _p2[3] ;

        _bb[3] = std::max (_r1, _r2) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_2d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in
        )
    {
        real_type _xm[2*2] ;
        real_type _xi[2*2] ;
        _xm[__ij(0,0,2)] = _p2[0]-_p1[0] ;
        _xm[__ij(0,1,2)] = _p2[1]-_p1[1] ;

        _xm[__ij(1,0,2)] = _p3[0]-_p1[0] ;
        _xm[__ij(1,1,2)] = _p3[1]-_p1[1] ;

        real_type _xr[2*1] ;
        _xr[0] = (real_type)+.5 * (
        std::pow(_xm[__ij(0,0,2)], +2) +
        std::pow(_xm[__ij(0,1,2)], +2)
            ) ;

        _xr[1] = (real_type)+.5 * (
        std::pow(_xm[__ij(1,0,2)], +2) +
        std::pow(_xm[__ij(1,1,2)], +2)
            ) ;

        real_type _w21 = _p2[2]-_p1[2] ;
        real_type _w31 = _p3[2]-_p1[2] ;

        _xr[0]-= (real_type)+.5 * _w21 ;
        _xr[1]-= (real_type)+.5 * _w31 ;

        real_type _dd ;
        math::inv_2x2 (
            +2, _xm, +2, _xi, _dd )  ;
        _bb[0] = (
        _xi[__ij(0,0,2)] * _xr[0] +
        _xi[__ij(0,1,2)] * _xr[1] )  ;

        _bb[1] = (
        _xi[__ij(1,0,2)] * _xr[0] +
        _xi[__ij(1,1,2)] * _xr[1] )  ;

        double   _P1[2];
        copy_node_2d (_P1, _p1) ;

        double   _P2[2];
        copy_node_2d (_P2, _p2) ;

        double   _P3[2];
        copy_node_2d (_P3, _p3) ;

        _dd = +geompred::orient2d (
            _P1, _P2, _P3) ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;

        real_type _ee[2*1];
        real_type _db[2*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,2)] * _bb[0] +
        _xm[__ij(0,1,2)] * _bb[1] )  ;

        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,2)] * _bb[0] +
        _xm[__ij(1,1,2)] * _bb[1] )  ;

        _db[0] = (
        _xi[__ij(0,0,2)] * _ee[0] +
        _xi[__ij(0,1,2)] * _ee[1] )  ;

        _db[1] = (
        _xi[__ij(1,0,2)] * _ee[0] +
        _xi[__ij(1,1,2)] * _ee[1] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        }

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;

        real_type _r1 =
        geometry::lensqr_2d(_p1, _bb) ;
        real_type _r2 =
        geometry::lensqr_2d(_p2, _bb) ;
        real_type _r3 =
        geometry::lensqr_2d(_p3, _bb) ;

        _r1     -= _p1[2]  ;
        _r2     -= _p2[2]  ;
        _r3     -= _p3[2]  ;

        _bb[2] =
        (_r1+_r2+_r3) / (real_type)3. ;

        if (!_in)  return  ;

        real_type __bf[3]  ;
        __perpface12(_bb, _p1, _p2, _p3) ;
        __perpface12(_bb, _p2, _p3, _p1) ;
        __perpface12(_bb, _p3, _p1, _p2) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_3d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  bool_type  _in
        )
    {
        real_type _nv[4*1] ;
        tria_norm_3d(_p1, _p2, _p3, _nv) ;

        real_type _xm[3*3] ;
        real_type _xi[3*3] ;
        _xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
        _xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;

        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
        _xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;

        _xm[__ij(2,0,3)] = _nv[0] ;
        _xm[__ij(2,1,3)] = _nv[1] ;
        _xm[__ij(2,2,3)] = _nv[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        std::pow(_xm[__ij(0,0,3)], +2) +
        std::pow(_xm[__ij(0,1,3)], +2) +
        std::pow(_xm[__ij(0,2,3)], +2)
            ) ;

        _xr[1] = (real_type)+.5 * (
        std::pow(_xm[__ij(1,0,3)], +2) +
        std::pow(_xm[__ij(1,1,3)], +2) +
        std::pow(_xm[__ij(1,2,3)], +2)
            ) ;

        real_type _w21 = _p2[3]-_p1[3] ;
        real_type _w31 = _p3[3]-_p1[3] ;

        _xr[0]-= (real_type)+.5 * _w21 ;
        _xr[1]-= (real_type)+.5 * _w31 ;

        _xr[2] = (real_type)+.0 ;

        real_type _dd ;
        math::inv_3x3 (
            +3, _xm, +3, _xi, _dd )  ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
        _xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;

        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
        _xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
        _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;

        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;

        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;

        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
        _xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;

        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
        _xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
        _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;
        }

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;

        real_type _r1 =
        geometry::lensqr_3d(_p1, _bb) ;
        real_type _r2 =
        geometry::lensqr_3d(_p2, _bb) ;
        real_type _r3 =
        geometry::lensqr_3d(_p3, _bb) ;

        _r1     -= _p1[3]  ;
        _r2     -= _p2[3]  ;
        _r3     -= _p3[3]  ;

        _bb[3] =
        (_r1+_r2+_r3) / (real_type)3. ;

        if (!_in)  return  ;

        real_type __bf[4]  ;
        __perpface13(_bb, _p1, _p2, _p3) ;
        __perpface13(_bb, _p2, _p3, _p1) ;
        __perpface13(_bb, _p3, _p1, _p2) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type perp_ball_3d (
    __write_ptr  (real_type) _bb ,
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  bool_type  _in
        )
    {
        real_type _xm[3*3] ;
        real_type _xi[3*3] ;
        _xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
        _xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;

        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
        _xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;

        _xm[__ij(2,0,3)] = _p4[0]-_p1[0] ;
        _xm[__ij(2,1,3)] = _p4[1]-_p1[1] ;
        _xm[__ij(2,2,3)] = _p4[2]-_p1[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        std::pow(_xm[__ij(0,0,3)], +2) +
        std::pow(_xm[__ij(0,1,3)], +2) +
        std::pow(_xm[__ij(0,2,3)], +2)
            ) ;

        _xr[1] = (real_type)+.5 * (
        std::pow(_xm[__ij(1,0,3)], +2) +
        std::pow(_xm[__ij(1,1,3)], +2) +
        std::pow(_xm[__ij(1,2,3)], +2)
            ) ;

        _xr[2] = (real_type)+.5 * (
        std::pow(_xm[__ij(2,0,3)], +2) +
        std::pow(_xm[__ij(2,1,3)], +2) +
        std::pow(_xm[__ij(2,2,3)], +2)
            ) ;

        real_type _w21 = _p2[3]-_p1[3] ;
        real_type _w31 = _p3[3]-_p1[3] ;
        real_type _w41 = _p4[3]-_p1[3] ;

        _xr[0]-= (real_type)+.5 * _w21 ;
        _xr[1]-= (real_type)+.5 * _w31 ;
        _xr[2]-= (real_type)+.5 * _w41 ;

        real_type _dd ;
        math::inv_3x3 (
            +3, _xm, +3, _xi, _dd )  ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
        _xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;

        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
        _xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
        _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        double   _P1[3];
        copy_node_3d (_P1, _p1) ;

        double   _P2[3];
        copy_node_3d (_P2, _p2) ;

        double   _P3[3];
        copy_node_3d (_P3, _p3) ;

        double   _P4[3];
        copy_node_3d (_P4, _p4) ;

        _dd = -geompred::orient3d (
            _P1, _P2, _P3, _P4) ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;

        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;

        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;

        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
        _xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;

        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
        _xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
        _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;
        }

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;

        real_type _r1 =
        geometry::lensqr_3d(_p1, _bb) ;
        real_type _r2 =
        geometry::lensqr_3d(_p2, _bb) ;
        real_type _r3 =
        geometry::lensqr_3d(_p3, _bb) ;
        real_type _r4 =
        geometry::lensqr_3d(_p4, _bb) ;

        _r1     -= _p1[3]  ;
        _r2     -= _p2[3]  ;
        _r3     -= _p3[3]  ;
        _r4     -= _p4[3]  ;

        _bb[3] =
        (_r1+_r2+_r3+_r4)/(real_type)4. ;

        if (!_in) return   ;

        real_type __bf[4]  ;
        __perpface23 (
            _bb, _p1 , _p2, _p3, _p4) ;
        __perpface23 (
            _bb, _p1 , _p2, _p4, _p3) ;
        __perpface23 (
            _bb, _p2 , _p3, _p4, _p1) ;
        __perpface23 (
            _bb, _p3 , _p1, _p4, _p2) ;
    }

    #undef  __perpface12
    #undef  __perpface13
    #undef  __perpface23

    /*
    --------------------------------------------------------
     * MASS-BALL: centre-of-mass balls.
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call void_type mass_ball_2d (
    __write_ptr  (real_type) _bb,
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]/= (real_type) +2. ;

        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]/= (real_type) +2. ;

        real_type _r1 =
        geometry::lensqr_2d(_bb, _p1);
        real_type _r2 =
        geometry::lensqr_2d(_bb, _p2);

        _bb[2] =
        (_r1+_r2)/(real_type)+2. ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type mass_ball_3d (
    __write_ptr  (real_type) _bb,
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]/= (real_type) +2. ;

        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]/= (real_type) +2. ;

        _bb[2] = _p1[2] ;
        _bb[2]+= _p2[2] ;
        _bb[2]/= (real_type) +2. ;

        real_type _r1 =
        geometry::lensqr_3d(_bb, _p1);
        real_type _r2 =
        geometry::lensqr_3d(_bb, _p2);

        _bb[3] =
        (_r1+_r2)/(real_type)+2. ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type mass_ball_2d (
    __write_ptr  (real_type) _bb,
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]+= _p3[0] ;
        _bb[0]/= (real_type) +3. ;

        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]+= _p3[1] ;
        _bb[1]/= (real_type) +3. ;

        real_type _r1 =
        geometry::lensqr_2d(_bb, _p1);
        real_type _r2 =
        geometry::lensqr_2d(_bb, _p2);
        real_type _r3 =
        geometry::lensqr_2d(_bb, _p3);

        _bb[2] = std::max ( _r3,
                 std::max ( _r1, _r2)
                          ) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type mass_ball_3d (
    __write_ptr  (real_type) _bb,
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]+= _p3[0] ;
        _bb[0]/= (real_type) +3. ;

        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]+= _p3[1] ;
        _bb[1]/= (real_type) +3. ;

        _bb[2] = _p1[2] ;
        _bb[2]+= _p2[2] ;
        _bb[2]+= _p3[2] ;
        _bb[2]/= (real_type) +3. ;

        real_type _r1 =
        geometry::lensqr_3d(_bb, _p1);
        real_type _r2 =
        geometry::lensqr_3d(_bb, _p2);
        real_type _r3 =
        geometry::lensqr_3d(_bb, _p3);

        _bb[3] = std::max ( _r3,
                 std::max ( _r1, _r2)
                          ) ;
    }

    template <
    typename      real_type
             >
    __normal_call void_type mass_ball_3d (
    __write_ptr  (real_type) _bb,
    __const_ptr  (real_type) _p1,
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3,
    __const_ptr  (real_type) _p4
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]+= _p3[0] ;
        _bb[0]+= _p4[0] ;
        _bb[0]/= (real_type) +4. ;

        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]+= _p3[1] ;
        _bb[1]+= _p4[1] ;
        _bb[1]/= (real_type) +4. ;

        _bb[2] = _p1[2] ;
        _bb[2]+= _p2[2] ;
        _bb[2]+= _p3[2] ;
        _bb[2]+= _p4[2] ;
        _bb[2]/= (real_type) +4. ;

        real_type _r1 =
        geometry::lensqr_3d(_bb, _p1);
        real_type _r2 =
        geometry::lensqr_3d(_bb, _p2);
        real_type _r3 =
        geometry::lensqr_3d(_bb, _p3);
        real_type _r4 =
        geometry::lensqr_3d(_bb, _p4);

        _bb[3] = std::max ( _r4,
                 std::max ( _r3,
                 std::max ( _r1, _r2)
                        ) ) ;
    }

    }

#   endif //__CELL_BALL_K__



