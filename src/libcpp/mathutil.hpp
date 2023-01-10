
    /*
    --------------------------------------------------------
     * MATH-UTIL: general-purpose math utilities.
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
    --------------------------------------------------------
     *
     * Last updated: 30 Mar., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MATH_UTIL__
#   define __MATH_UTIL__

    namespace math
    {

    /*
    --------------------------------------------------------
     * calc. polynomial roots.
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __normal_call bool_type polyroots (
        real_type _aa,      // aa*xx^2+bb*xx+cc=0
        real_type _bb,
        real_type _cc,
    __write_ptr  (real_type) _xx
        )
    {
        real_type _rt =
            std::numeric_limits
                <real_type>::epsilon() ;

        bool_type _real = false ;

        real_type _sq = _bb * _bb -
       (real_type)+4. * _aa * _cc ;

        if (_sq >= (real_type)+0.)  // real roots
        {
            _sq  = std::sqrt(_sq) ;

            _xx[0] = (-_bb + _sq) ;
            _xx[1] = (-_bb - _sq) ;

            real_type _xm = std::max (
                std::abs(_xx[0]),
                std::abs(_xx[1])) ;

            if (std::abs(_aa) >
                std::abs(_xm) * _rt)
            {
                _real =  true ;

                _aa *=(real_type)+2. ;

                _xx[0] /= _aa ;
                _xx[1] /= _aa ;
            }
            else
            if (std::abs(_bb) >
                std::abs(_cc) * _rt)
            {
                _real =  true ;

                _xx[0] = -_cc / _bb  ;
                _xx[1] = -_cc / _bb  ;
            }
        }

        return _real ;
    }

    /*
    --------------------------------------------------------
     * axis-aligned ellipsoid.
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __inline_call void_type toR3 (
      __const_ptr(real_type) _rrad ,
      __const_ptr(real_type) _apos ,
      __write_ptr(real_type) _ppos
        )
    {
    /*------------ helper: convert from S^2 to R^3 coord. */
        _ppos[0] = _rrad[0] *
            std::cos( _apos[0] ) *
            std::cos( _apos[1] ) ;

        _ppos[1] = _rrad[1] *
            std::sin( _apos[0] ) *
            std::cos( _apos[1] ) ;

        _ppos[2] = _rrad[2] *
            std::sin( _apos[1] ) ;
    }

    template <
        typename  real_type
             >
    __inline_call void_type toS2 (
      __const_ptr(real_type) _rrad ,
      __const_ptr(real_type) _ppos ,
      __write_ptr(real_type) _apos
        )
    {
    /*------------ helper: convert from R^3 to S^2 coord. */
        real_type _xmul =
            _ppos[0] * _rrad [1] ;
        real_type _ymul =
            _ppos[1] * _rrad [0] ;
        real_type _zrat =
            _ppos[2] / _rrad [2] ;

        _zrat = std::min(+1.,_zrat);
        _zrat = std::max(-1.,_zrat);

        _apos[0]= std::atan2(_ymul,
                             _xmul);
        _apos[1]= std::asin (_zrat);
    }

    /*
    --------------------------------------------------------
     * small matrix utilities.
    --------------------------------------------------------
     */

#   define __ij(__ir , __ic , __nr) \
              ((__ir)+(__ic)*(__nr))

    template <
        typename  real_type,
        typename  size_type
             >
    __inline_call real_type det_2x2 (
        size_type  _la,
    __const_ptr  (real_type) _aa
        )
    {   return
        _aa[__ij(0, 0, _la)] *
        _aa[__ij(1, 1, _la)] -
        _aa[__ij(0, 1, _la)] *
        _aa[__ij(1, 0, _la)] ;
    }

    template <
        typename  real_type,
        typename  size_type
             >
    __inline_call real_type det_3x3 (
        size_type  _la,
    __const_ptr  (real_type) _aa
        )
    {   return
        _aa[__ij(0, 0, _la)] * (
        _aa[__ij(1, 1, _la)] *
        _aa[__ij(2, 2, _la)] -
        _aa[__ij(1, 2, _la)] *
        _aa[__ij(2, 1, _la)] ) -

        _aa[__ij(0, 1, _la)] * (
        _aa[__ij(1, 0, _la)] *
        _aa[__ij(2, 2, _la)] -
        _aa[__ij(1, 2, _la)] *
        _aa[__ij(2, 0, _la)] ) +

        _aa[__ij(0, 2, _la)] * (
        _aa[__ij(1, 0, _la)] *
        _aa[__ij(2, 1, _la)] -
        _aa[__ij(1, 1, _la)] *
        _aa[__ij(2, 0, _la)] ) ;
    }

    template <
        typename  real_type,
        typename  size_type
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

        _xx[__ij(0, 0, _lx)] =
        _aa[__ij(1, 1, _la)] ;
        _xx[__ij(1, 1, _lx)] =
        _aa[__ij(0, 0, _la)] ;
        _xx[__ij(0, 1, _lx)] =
       -_aa[__ij(0, 1, _la)] ;
        _xx[__ij(1, 0, _lx)] =
       -_aa[__ij(1, 0, _la)] ;
    }

    template <
        typename  real_type,
        typename  size_type
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

        _xx[__ij(0, 0, _lx)] =
        _aa[__ij(2, 2, _la)] *
        _aa[__ij(1, 1, _la)] -
        _aa[__ij(2, 1, _la)] *
        _aa[__ij(1, 2, _la)] ;

        _xx[__ij(0, 1, _lx)] =
        _aa[__ij(2, 1, _la)] *
        _aa[__ij(0, 2, _la)] -
        _aa[__ij(2, 2, _la)] *
        _aa[__ij(0, 1, _la)] ;

        _xx[__ij(0, 2, _lx)] =
        _aa[__ij(1, 2, _la)] *
        _aa[__ij(0, 1, _la)] -
        _aa[__ij(1, 1, _la)] *
        _aa[__ij(0, 2, _la)] ;

        _xx[__ij(1, 0, _lx)] =
        _aa[__ij(2, 0, _la)] *
        _aa[__ij(1, 2, _la)] -
        _aa[__ij(2, 2, _la)] *
        _aa[__ij(1, 0, _la)] ;

        _xx[__ij(1, 1, _lx)] =
        _aa[__ij(2, 2, _la)] *
        _aa[__ij(0, 0, _la)] -
        _aa[__ij(2, 0, _la)] *
        _aa[__ij(0, 2, _la)] ;

        _xx[__ij(1, 2, _lx)] =
        _aa[__ij(1, 0, _la)] *
        _aa[__ij(0, 2, _la)] -
        _aa[__ij(1, 2, _la)] *
        _aa[__ij(0, 0, _la)] ;

        _xx[__ij(2, 0, _lx)] =
        _aa[__ij(2, 1, _la)] *
        _aa[__ij(1, 0, _la)] -
        _aa[__ij(2, 0, _la)] *
        _aa[__ij(1, 1, _la)] ;

        _xx[__ij(2, 1, _lx)] =
        _aa[__ij(2, 0, _la)] *
        _aa[__ij(0, 1, _la)] -
        _aa[__ij(2, 1, _la)] *
        _aa[__ij(0, 0, _la)] ;

        _xx[__ij(2, 2, _lx)] =
        _aa[__ij(1, 1, _la)] *
        _aa[__ij(0, 0, _la)] -
        _aa[__ij(1, 0, _la)] *
        _aa[__ij(0, 1, _la)] ;
    }


    }

#   endif//__MATH_UTIL__



