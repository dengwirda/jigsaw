
    /*
    --------------------------------------------------------
     * HFUN-CLIP-kD: kernels for eikonal solvers.
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
     * Last updated: 24 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_CLIP_K__
#   define __HFUN_CLIP_K__

    namespace mesh {

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type EIKONAL_edge_2d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
                  vals_type  _g1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type& _h2 ,
                  vals_type  _g2
        )
    {
        real_type _gg  = 
      ((real_type)_g1 + _g2) / (real_type)2.;

    /*---------------------- form "limited" extrap. to P2 */
        real_type _hn = (real_type)_h1 +
            _gg * std::sqrt(
        geometry::lensqr_2d(_p1, _p2)) ;

        vals_type _hm = (vals_type)_hn ;
        if (_hm < _h2)
        {
            _h2 = _hm; return  true ;
        }
        else/*no clip*/return false ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type EIKONAL_edge_3d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
                  vals_type  _g1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type& _h2 ,
                  vals_type  _g2
        )
    {
        real_type _gg  = 
      ((real_type)_g1 + _g2) / (real_type)2.;

    /*---------------------- form "limited" extrap. to P2 */
        real_type _hn = (real_type)_h1 +
            _gg * std::sqrt(
        geometry::lensqr_3d(_p1, _p2)) ;

        vals_type _hm = (vals_type)_hn ;
        if (_hm < _h2)
        {
            _h2 = _hm; return  true ;
        }
        else/*no clip*/return false ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call real_type EIKONAL_proj_2d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type  _h2 ,
    __const_ptr  (real_type) _pk ,
                  real_type  _gg ,
                  real_type  _g1 ,
                  real_type  _g2
        )
    {
    /*---------------------- form "limited" extrap. to PK */
        real_type const _zero = (real_type)0.;

        real_type _edge[2];
        real_type _src1[2], _src2[2];
        geometry::vector_2d(_p1, _p2, _edge) ;
        geometry::vector_2d(_p1, _pk, _src1) ;
        geometry::vector_2d(_p2, _pk, _src2) ;

        real_type _esqr = 
            geometry::lensqr_2d(_edge) ;
        real_type _1sqr = 
            geometry::lensqr_2d(_src1) ;
        real_type _2sqr = 
            geometry::lensqr_2d(_src2) ;

        real_type _1len = std::sqrt(_1sqr) ;
        real_type _2len = std::sqrt(_2sqr) ;

        real_type _hn = std::min( 
       (real_type)_h1 + std::min(_gg, _g1) *_1len,
       (real_type)_h2 + std::min(_gg, _g2) *_2len) ;

        real_type _dh = 
           (real_type)_h2 - (real_type)_h1 ;

        real_type _mm = 
            _gg * _gg * _esqr -  _dh * _dh ;
        if (_mm > (real_type)+0.)
        {
            real_type _proj = 
                geometry::dot_2d(_edge, _src1) ;

            real_type _area = std::abs(
                _edge[0] * _src1[1] - 
                _edge[1] * _src1[0] ) ;

            real_type _rm = std::sqrt(_mm) ;
            real_type _te = 
                _proj - _area * _dh / _rm  ;
            if (_te > _zero && _te < _esqr)
            {
                real_type _dd = 
                   (_dh * _proj + _area * _rm) ;
                _hn = std::min( _hn, 
                    (real_type) _h1 + _dd / _esqr) ;
            }
        }
        return _hn ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call real_type EIKONAL_proj_3d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type  _h2 ,
    __const_ptr  (real_type) _pk ,
                  real_type  _gg ,
                  real_type  _g1 ,
                  real_type  _g2
        )
    {
    /*---------------------- form "limited" extrap. to P3 */
        real_type const _zero = (real_type)0.;

        real_type _edge[3];
        real_type _src1[3], _src2[3];
        geometry::vector_3d(_p1, _p2, _edge) ;
        geometry::vector_3d(_p1, _pk, _src1) ;
        geometry::vector_3d(_p2, _pk, _src2) ;

        real_type _esqr = 
            geometry::lensqr_3d(_edge) ;
        real_type _1sqr = 
            geometry::lensqr_3d(_src1) ;
        real_type _2sqr = 
            geometry::lensqr_3d(_src2) ;

        real_type _1len = std::sqrt(_1sqr) ;
        real_type _2len = std::sqrt(_2sqr) ;

        real_type _hn = std::min( 
       (real_type)_h1 + std::min(_gg, _g1) *_1len,
       (real_type)_h2 + std::min(_gg, _g2) *_2len) ;

        real_type _dh = 
           (real_type)_h2 - (real_type)_h1 ;

        real_type _mm = 
            _gg * _gg * _esqr -  _dh * _dh ;
        if (_mm > (real_type)+0.)
        {
            real_type _proj = 
                geometry::dot_3d(_edge, _src1) ;

            real_type _avec[3];
            geometry::cross_3d(_edge, _src1, _avec);
            real_type _area = std::sqrt(
                geometry::lensqr_3d(_avec));

            real_type _rm = std::sqrt(_mm) ;
            real_type _te = 
                _proj - _area * _dh / _rm  ;
            if (_te > _zero && _te < _esqr)
            {
                real_type _dd = 
                   (_dh * _proj + _area * _rm) ;
                _hn = std::min( _hn, 
                    (real_type) _h1 + _dd / _esqr) ;
            }
        }
        return _hn ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type EIKONAL_tria_2d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
                  vals_type  _g1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type  _h2 ,
                  vals_type  _g2 ,
    __const_ptr  (real_type) _p3 ,
                  vals_type& _h3 ,
                  vals_type  _g3
        )
    {
        vals_type _h0 = std::min(_h1, _h2) ;
        if(_h3 <= _h0) return false ;

        real_type _gg  = 
      ((real_type)_g1 + _g2 + _g3) / (real_type)3.;
        
        if (_gg <= (real_type)+0.)
        {
            _h3 = _h0; return  true ;
        }

        real_type _G1 = 
      ((real_type)_g1 + _g3 ) / (real_type)2.;
        real_type _G2 = 
      ((real_type)_g2 + _g3 ) / (real_type)2.;

        vals_type _hm = 
           (vals_type)EIKONAL_proj_2d (
            _p1, _h1, _p2, _h2, _p3, 
            _gg, _G1, _G2) ;

        if (_hm < _h3)
        {
            _h3 = _hm; return  true ;
        }
        else/*no clip*/return false ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type EIKONAL_tria_3d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
                  vals_type  _g1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type  _h2 ,
                  vals_type  _g2 ,
    __const_ptr  (real_type) _p3 ,
                  vals_type& _h3 ,
                  vals_type  _g3
        )
    {
        vals_type _h0 = std::min(_h1, _h2) ;
        if(_h3 <= _h0) return false ;

        real_type _gg  = 
      ((real_type)_g1 + _g2 + _g3) / (real_type)3.;
        
        if (_gg <= (real_type)+0.)
        {
            _h3 = _h0; return  true ;
        }

        real_type _G1 = 
      ((real_type)_g1 + _g3 ) / (real_type)2.;
        real_type _G2 = 
      ((real_type)_g2 + _g3 ) / (real_type)2.;

        vals_type _hm = 
           (vals_type)EIKONAL_proj_3d (
            _p1, _h1, _p2, _h2, _p3, 
            _gg, _G1, _G2) ;

        if (_hm < _h3)
        {
            _h3 = _hm; return  true ;
        }
        else/*no clip*/return false ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type EIKONAL_tria_3d (
    __const_ptr  (real_type) _p1 ,
                  vals_type  _h1 ,
                  vals_type  _g1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type  _h2 ,
                  vals_type  _g2 ,
    __const_ptr  (real_type) _p3 ,
                  vals_type  _h3 ,
                  vals_type  _g3 ,
    __const_ptr  (real_type) _p4 ,
                  vals_type& _h4 ,
                  vals_type  _g4
        )
    {
    /*---------------------- form "limited" extrap. to P4 */
        vals_type _h0 = std::min({_h1, _h2, _h3}) ;
        if(_h4 <= _h0) return false ;

        real_type _gg  = 
      ((real_type)_g1 + _g2 + _g3 + _g4 ) / (real_type)4. ;
        
        if (_gg <= (real_type)+0.)
        {
            _h4 = _h0; return  true ;
        }

        __const_ptr(real_type) 
                  _pj[3] ={_p1, _p2, _p3} ;
        vals_type _hj[3] ={_h1, _h2, _h3} ;
        vals_type _gj[3] ={_g1, _g2, _g3} ;
        
        real_type _ge[3] ;
        for (int _ii = 0; _ii < 3; ++_ii)
            _ge [_ii] = 
          ((real_type)_gj[_ii]  + _g4) / (real_type)2.;

        real_type _hn = _h4;
        for (int _ii = 0; _ii < 3; ++_ii)
        {
            auto _i0 = (_ii + 0) % 3 ;
            auto _G0 = _gj[_i0] ;
            auto _i1 = (_ii + 1) % 3 ;
            auto _G1 = _gj[_i1] ;

            real_type _gk =
          ((real_type)_G0 + _G1 + _g4) / (real_type)3.;

            _gk = std::min(_gg, _gk) ;
            _hn = std::min(_hn, EIKONAL_proj_3d(
                _pj[_i0], _hj[_i0], 
                _pj[_i1], _hj[_i1],
                _p4, _gk, _ge[_i0], _ge[_i1])) ;
        }

        real_type _e2[3], _e3[3], _e4[3], _norm[3] ;
        geometry::vector_3d(_p1, _p2, _e2);
        geometry::vector_3d(_p1, _p3, _e3);
        geometry::vector_3d(_p1, _p4, _e4);
        geometry:: cross_3d(_e2, _e3, _norm);

        real_type _nsqr = geometry::lensqr_3d(_norm) ;
        if (_nsqr > (real_type)0.)
        {
            real_type _opp1[3], _opp2[3], _grad[3] ;
            geometry::cross_3d(_e3, _norm, _opp1);
            geometry::cross_3d(_norm, _e2, _opp2);

            real_type _inv = (real_type)+1. / _nsqr  ;
            real_type _dh1 = (real_type)_h2 - 
                             (real_type)_h1 ;
            real_type _dh2 = (real_type)_h3 - 
                             (real_type)_h1 ;
            for (int _ii = 0; _ii < 3; ++_ii)
            {
                _opp1[_ii]*= _inv ;
                _opp2[_ii]*= _inv ;
                _grad[_ii] = _dh1 * _opp1[_ii] +
                             _dh2 * _opp2[_ii] ;
            }

            real_type _mm = 
                _gg * _gg - geometry::lensqr_3d(_grad) ;
            if (_mm > (real_type)0.)
            {
                real_type _dist = std::abs(
                    geometry::dot_3d (_e4, _norm)) 
                        / std::sqrt(_nsqr) ;

                real_type _push = 
                    _dist / std::sqrt(_mm) ;

                real_type _es[3];
                for (int _ii = 0; _ii < 3; ++_ii)
                    _es [_ii] = 
                        _e4[_ii] - _push * _grad [_ii] ;

                real_type _uu = 
                    geometry::dot_3d (_es, _opp1) ;
                real_type _vv = 
                    geometry::dot_3d (_es, _opp2) ;
                real_type _ww = 
                    (real_type)1. - _uu - _vv ;

                if (_uu > (real_type)0. &&
                    _vv > (real_type)0. &&
                    _ww > (real_type)0. )
                {
                    real_type _offs [3] ;
                    for (int _ii = 0; _ii < 3; ++_ii)
                        _offs[_ii] = _uu * _e2[_ii] +
                                     _vv * _e3[_ii] - 
                                           _e4[_ii] ;

                    _hn = std::min( _hn,
                        _ww * _h1 + _uu * _h2 + _vv * _h3 + 
                            _gg * geometry::length_3d(_offs)) ;
                }
            }
        }

        vals_type _hm = (vals_type) _hn ;
        if (_hm < _h4)
        {
            _h4 = _hm; return  true;
        }
        else/*no clip*/return false;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_edge_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type  _g1 ,
                  vals_type  _g2
        )
    {
    /*---------------------- limit h-values within EDGE-2 */
        bool_type _clip = false ;

        if (_h1 >=(vals_type)+0. &&
            _h2 >=(vals_type)+0. )
        {            
            if (_si != +1)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_edge_2d  (
                    _p1, _o1, _g1 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;

            if (_si != +0)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_edge_2d  (
                    _p2, _o2, _g2 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_edge_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type  _g1 ,
                  vals_type  _g2
        )
    {
    /*---------------------- limit h-values within EDGE-2 */
        bool_type _clip = false ;

        if (_h1 >=(vals_type)+0. &&
            _h2 >=(vals_type)+0. )
        {            
            if (_si != +1)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_edge_3d  (
                    _p1, _o1, _g1 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;

            if (_si != +0)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_edge_3d  (
                    _p2, _o2, _g2 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_tria_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3
        )
    {
    /*---------------------- limit h-values within TRIA-3 */
        bool_type _clip = false ;

        if (_o1 >=(vals_type)+0. &&
            _o2 >=(vals_type)+0. &&
            _o3 >=(vals_type)+0. )
        {
            if (_si != +2)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_tria_2d  (
                    _p1, _o1, _g1 ,
                    _p2, _o2, _g2 ,
                    _p3, _h3, _g3 )  )
                    _clip =  true ;

            if (_si != +0)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_tria_2d  (
                    _p2, _o2, _g2 ,
                    _p3, _o3, _g3 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;

            if (_si != +1)
    /*--------------------------------- 3rd node ordering */
                if (EIKONAL_tria_2d  (
                    _p3, _o3, _g3 ,
                    _p1, _o1, _g1 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_tria_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3
        )
    {
    /*---------------------- limit h-values within TRIA-3 */
        bool_type _clip = false ;

        if (_o1 >=(vals_type)+0. &&
            _o2 >=(vals_type)+0. &&
            _o3 >=(vals_type)+0. )
        {
            if (_si != +2)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_tria_3d  (
                    _p1, _o1, _g1 ,
                    _p2, _o2, _g2 ,
                    _p3, _h3, _g3 )  )
                    _clip =  true ;

            if (_si != +0)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_tria_3d  (
                    _p2, _o2, _g2 ,
                    _p3, _o3, _g3 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;

            if (_si != +1)
    /*--------------------------------- 3rd node ordering */
                if (EIKONAL_tria_3d  (
                    _p3, _o3, _g3 ,
                    _p1, _o1, _g1 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;
        }

        return ( _clip ) ;
    }
    
    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_tria_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  int        _si ,  // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4
        )
    {
    /*---------------------- limit h-values within TRIA-4 */
        bool_type _clip = false ;

        if (_o1 >=(vals_type)+0. &&
            _o2 >=(vals_type)+0. &&
            _o3 >=(vals_type)+0. &&
            _o4 >=(vals_type)+0. )
        {
            if (_si != +3)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_tria_3d  (
                    _p1, _o1, _g1 ,
                    _p2, _o2, _g2 ,
                    _p3, _o3, _g3 ,
                    _p4, _h4, _g4 )  )
                    _clip =  true ;

            if (_si != +0)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_tria_3d  (
                    _p2, _o2, _g2 ,
                    _p3, _o3, _g3 ,
                    _p4, _o4, _g4 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;

            if (_si != +1)
    /*--------------------------------- 3rd node ordering */
                if (EIKONAL_tria_3d  (
                    _p3, _o3, _g3 ,
                    _p1, _o1, _g1 ,
                    _p4, _o4, _g4 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;
                    
            if (_si != +2)
    /*--------------------------------- 4th node ordering */
                if (EIKONAL_tria_3d  (
                    _p1, _o1, _g1 ,
                    _p2, _o2, _g2 ,
                    _p4, _o4, _g4 ,
                    _p3, _h3, _g3 )  )
                    _clip =  true ;
        }
        
        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_quad_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4
        )
    {
    /*---------------------- limit h-values within QUAD-4 */
        bool_type _clip = false ;
        bool_type _okay = false ;

        real_type _a1 =
        geometry::tria_area_2d(_p1, _p2, _p3) ;
        real_type _a2 =
        geometry::tria_area_2d(_p1, _p3, _p4) ;

        if (_a1*_a2 > (real_type)+0. )
        {
            _okay =  true ;
    /*--------------------------------- 1st tria ordering */
            if (eikonal_tria_2d (
                _p1, _p2, _p3 , _si - 0,
                _o1, _o2, _o3 ,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p1, _p3, _p4 , 
                _si == 1 ? -1 : _si - (_si > 1),
                _o1, _o3, _o4 , 
                _h1, _h3, _h4 ,
                _g1, _g3, _g4 ) )
                _clip =  true ;
        }

        real_type _a3 =
        geometry::tria_area_2d(_p1, _p2, _p4) ;
        real_type _a4 =
        geometry::tria_area_2d(_p2, _p3, _p4) ;

        if (_a3*_a4 > (real_type)+0. )
        {
            _okay =  true ;
    /*--------------------------------- 2nd tria ordering */
            if (eikonal_tria_2d (
                _p1, _p2, _p4 , 
                _si == 2 ? -1 : _si - (_si > 2),
                _o1, _o2, _o4 ,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p2, _p3, _p4 , _si - 1,
                _o2, _o3, _o4 ,
                _h2, _h3, _h4 ,
                _g2, _g3, _g4 ) )
                _clip =  true ;
        }

        if (!_okay)
        {
    /*--------------------------------- a degenerate quad */
            if (_h1 >=(vals_type)+0. &&
                _h2 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p1, _p2, 
                _si == 0 ? 0 : _si == 1 ? 1 : -1,
                _o1, _o2,
                _h1, _h2, _g1 , _g2) )
                _clip =  true ;

            if (_h2 >=(vals_type)+0. &&
                _h3 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p2, _p3, 
                _si == 1 ? 0 : _si == 2 ? 1 : -1,
                _o2, _o3,
                _h2, _h3, _g2 , _g3) )
                _clip =  true ;

            if (_h3 >=(vals_type)+0. &&
                _h4 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p3, _p4, 
                _si == 2 ? 0 : _si == 3 ? 1 : -1,
                _o3, _o4,
                _h3, _h4, _g3 , _g4) )
                _clip =  true ;

            if (_h4 >=(vals_type)+0. &&
                _h1 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p4, _p1, 
                _si == 3 ? 0 : _si == 0 ? 1 : -1,
                _o4, _o1,
                _h4, _h1, _g4 , _g1) )
                _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_grid_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4
        )
    {
    /*---------------------- limit h-values within GRID-4 */
        bool_type _clip = false ;
        {
    /*--------------------------------- 1st tria ordering */
            if (eikonal_tria_2d (
                _p1, _p2, _p3 , _si - 0,
                _o1, _o2, _o3 ,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p1, _p3, _p4 , 
                _si == 1 ? -1 : _si - (_si > 1),
                _o1, _o3, _o4 , 
                _h1, _h3, _h4 ,
                _g1, _g3, _g4 ) )
                _clip =  true ;
        
    /*--------------------------------- 2nd tria ordering */
            if (eikonal_tria_2d (
                _p1, _p2, _p4 , 
                _si == 2 ? -1 : _si - (_si > 2),
                _o1, _o2, _o4 ,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p2, _p3, _p4 , _si - 1,
                _o2, _o3, _o4 ,
                _h2, _h3, _h4 ,
                _g2, _g3, _g4 ) )
                _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_quad_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4
        )
    {
    /*---------------------- limit h-values within QUAD-4 */
        bool_type _clip = false ;
        bool_type _okay = false ;

        real_type _n1[3];
        geometry::tria_norm_3d(_p1, _p2, _p3,
                               _n1) ;
        real_type _n2[3];
        geometry::tria_norm_3d(_p1, _p3, _p4,
                               _n2) ;

        if (geometry::dot_3d(
                _n1, _n2) > (real_type)+0.)
        {
            _okay =  true ;
    /*--------------------------------- 1st tria ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p3 , _si - 0,
                _o1, _o2, _o3 ,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p1, _p3, _p4 , 
                _si == 1 ? -1 : _si - (_si > 1),
                _o1, _o3, _o4 , 
                _h1, _h3, _h4 ,
                _g1, _g3, _g4 ) )
                _clip =  true ;
        }

        real_type _n3[3];
        geometry::tria_norm_3d(_p1, _p2, _p4,
                               _n3) ;
        real_type _n4[3];
        geometry::tria_norm_3d(_p2, _p3, _p4,
                               _n4) ;

        if (geometry::dot_3d(
                _n3, _n4) > (real_type)+0.)
        {
            _okay =  true ;
    /*--------------------------------- 2nd tria ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p4 , 
                _si == 2 ? -1 : _si - (_si > 2),
                _o1, _o2, _o4 ,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p2, _p3, _p4 , _si - 1,
                _o2, _o3, _o4 ,
                _h2, _h3, _h4 ,
                _g2, _g3, _g4 ) )
                _clip =  true ;
        }

        if (!_okay)
        {
    /*--------------------------------- a degenerate quad */
            if (_h1 >=(vals_type)+0. &&
                _h2 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p1, _p2, 
                _si == 0 ? 0 : _si == 1 ? 1 : -1,
                _o1, _o2,
                _h1, _h2, _g1 , _g2) )
                _clip =  true ;

            if (_h2 >=(vals_type)+0. &&
                _h3 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p2, _p3, 
                _si == 1 ? 0 : _si == 2 ? 1 : -1,
                _o2, _o3,
                _h2, _h3, _g2 , _g3) )
                _clip =  true ;

            if (_h3 >=(vals_type)+0. &&
                _h4 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p3, _p4, 
                _si == 2 ? 0 : _si == 3 ? 1 : -1,
                _o3, _o4,
                _h3, _h4, _g3 , _g4) )
                _clip =  true ;

            if (_h4 >=(vals_type)+0. &&
                _h1 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p4, _p1, 
                _si == 3 ? 0 : _si == 0 ? 1 : -1,
                _o4, _o1,
                _h4, _h1, _g4 , _g1) )
                _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_grid_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4
        )
    {
    /*---------------------- limit h-values within GRID-4 */
        bool_type _clip = false ;
        {
    /*--------------------------------- 1st tria ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p3 , _si - 0,
                _o1, _o2, _o3 ,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p1, _p3, _p4 , 
                _si == 1 ? -1 : _si - (_si > 1),
                _o1, _o3, _o4 , 
                _h1, _h3, _h4 ,
                _g1, _g3, _g4 ) )
                _clip =  true ;
        
    /*--------------------------------- 2nd tria ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p4 , 
                _si == 2 ? -1 : _si - (_si > 2),
                _o1, _o2, _o4 ,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p2, _p3, _p4 , _si - 1,
                _o2, _o3, _o4 ,
                _h2, _h3, _h4 ,
                _g2, _g3, _g4 ) )
                _clip =  true ;
        }

        return ( _clip ) ;
    }

    __inline_call int eikonal_subindex (
        int _si, 
        int _i1, int _i2, int _i3, int _i4, int _i5 = -1
        )
    {
        return _si < 0 ?-1 :
            _si == _i1 ? 0 : _si == _i2 ? 1 :
            _si == _i3 ? 2 : _si == _i4 ? 3 :
            _si == _i5 ? 4 : -1 ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __inline_call bool_type eikonal_pyra_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
    __const_ptr  (real_type) _p5 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type  _o5 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type& _h5 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4 ,
                  vals_type  _g5
        )
    {
    /*---------------------- limit h-values within PYRA-5 */
        bool_type _clip = false ;
        {
    /*--------------------------------- 1st base ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p3, _p5,
                eikonal_subindex(_si, 0, 1, 2, 4),
                _o1, _o2, _o3, _o5,
                _h1, _h2, _h3, _h5,
                _g1, _g2, _g3, _g5) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p1, _p3, _p4, _p5,
                eikonal_subindex(_si, 0, 2, 3, 4),
                _o1, _o3, _o4, _o5,
                _h1, _h3, _h4, _h5,
                _g1, _g3, _g4, _g5) )
                _clip =  true ;
        
    /*--------------------------------- 2nd base ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p4, _p5,
                eikonal_subindex(_si, 0, 1, 3, 4),
                _o1, _o2, _o4, _o5,
                _h1, _h2, _h4, _h5,
                _g1, _g2, _g4, _g5) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p2, _p3, _p4, _p5,
                eikonal_subindex(_si, 1, 2, 3, 4),
                _o2, _o3, _o4, _o5,
                _h2, _h3, _h4, _h5,
                _g2, _g3, _g4, _g5) )
                _clip =  true ;
        }

        return ( _clip ) ;
    }
    
    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_wedg_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
    __const_ptr  (real_type) _p5 ,
    __const_ptr  (real_type) _p6 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type  _o5 ,
                  vals_type  _o6 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type& _h5 ,
                  vals_type& _h6 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4 ,
                  vals_type  _g5 ,
                  vals_type  _g6
        )
    {
    /*---------------------- limit h-values within WEDG-6 */
        bool_type _clip = false ;
        
    //  Similiarly to the hexa case, iterate on a regular 
    //  division about the cell centroid...
    
        real_type constexpr _6r = (real_type) +6. ;
        real_type _pc[3] = {
        (_p1[0] + _p2[0] + _p3[0] + _p4[0] + 
         _p5[0] + _p6[0] ) / _6r,
        (_p1[1] + _p2[1] + _p3[1] + _p4[1] + 
         _p5[1] + _p6[1] ) / _6r,
        (_p1[2] + _p2[2] + _p3[2] + _p4[2] + 
         _p5[2] + _p6[2] ) / _6r,
            } ;
            
        vals_type _hc = (vals_type)((
       (real_type)_o1 + _o2 + _o3 + _o4 + _o5 + _o6) / _6r);
        vals_type _gc = (vals_type)((
       (real_type)_g1 + _g2 + _g3 + _g4 + _g5 + _g6) / _6r);
         
         for (auto _iter = 2; _iter-- != 0; )
        {
            bool_type _diff = false;
            vals_type _oc = _hc;
    /*----------------------------- subcells about centre */
            if (eikonal_tria_3d (
                _p1, _p2, _p3, _pc,
                eikonal_subindex(_si, 0, 1, 2, 6),
                _o1, _o2, _o3, _oc,
                _h1, _h2, _h3, _hc,
                _g1, _g2, _g3, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
    
            if (eikonal_pyra_3d (
                _p1, _p2, _p5, _p4, _pc,
                eikonal_subindex(_si, 0, 1, 4, 3, 6),
                _o1, _o2, _o5, _o4, _oc,
                _h1, _h2, _h5, _h4, _hc,
                _g1, _g2, _g5, _g4, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p2, _p3, _p6, _p5, _pc,
                eikonal_subindex(_si, 1, 2, 5, 4, 6),
                _o2, _o3, _o6, _o5, _oc,
                _h2, _h3, _h6, _h5, _hc,
                _g2, _g3, _g6, _g5, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p3, _p1, _p4, _p6, _pc,
                eikonal_subindex(_si, 2, 0, 3, 5, 6),
                _o3, _o1, _o4, _o6, _oc,
                _h3, _h1, _h4, _h6, _hc,
                _g3, _g1, _g4, _g6, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_tria_3d (
                _p6, _p5, _p4, _pc,
                eikonal_subindex(_si, 5, 4, 3, 6),
                _o6, _o5, _o4, _oc,
                _h6, _h5, _h4, _hc,
                _g6, _g5, _g4, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (!_diff ) break;
        }

        return ( _clip ) ;
    }
    
    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_hexa_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
    __const_ptr  (real_type) _p5 ,
    __const_ptr  (real_type) _p6 ,
    __const_ptr  (real_type) _p7 ,
    __const_ptr  (real_type) _p8 ,
                  int        _si ,      // current source
                  vals_type  _o1 ,
                  vals_type  _o2 ,
                  vals_type  _o3 ,
                  vals_type  _o4 ,
                  vals_type  _o5 ,
                  vals_type  _o6 ,
                  vals_type  _o7 ,
                  vals_type  _o8 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type& _h3 ,
                  vals_type& _h4 ,
                  vals_type& _h5 ,
                  vals_type& _h6 ,
                  vals_type& _h7 ,
                  vals_type& _h8 ,
                  vals_type  _g1 ,
                  vals_type  _g2 ,
                  vals_type  _g3 ,
                  vals_type  _g4 ,
                  vals_type  _g5 ,
                  vals_type  _g6 ,
                  vals_type  _g7 ,
                  vals_type  _g8
        )
    {
    /*---------------------- limit h-values within HEXA-8 */
        bool_type _clip = false ;
    
    //  There are a very large number of hexa => tet cases:
    //  J. Pellerin, K. Verhetsel, J.F. Remacle (2018): 
    //  There are 174 Subdivisions of the Hexahedron 
    //  into Tetrahedra. 
    //  ACM Transactions on Graphics (TOG), 37(6), pp. 1-9.
    
    //  Instead, do a few iterations on a regular division
    //  about the cell centroid...
    
        real_type constexpr _8r = (real_type)+8.;
        real_type _pc[3] = {
        (_p1[0] + _p2[0] + _p3[0] + _p4[0] + 
         _p5[0] + _p6[0] + _p7[0] + _p8[0] ) / _8r,
        (_p1[1] + _p2[1] + _p3[1] + _p4[1] + 
         _p5[1] + _p6[1] + _p7[1] + _p8[1] ) / _8r,
        (_p1[2] + _p2[2] + _p3[2] + _p4[2] + 
         _p5[2] + _p6[2] + _p7[2] + _p8[2] ) / _8r,
            } ;
            
        vals_type _hc = (vals_type)((
       (real_type)_o1 + _o2 + _o3 + _o4 + 
                  _o5 + _o6 + _o7 + _o8 ) / _8r);
        vals_type _gc = (vals_type)((
       (real_type)_g1 + _g2 + _g3 + _g4 + 
                  _g5 + _g6 + _g7 + _g8 ) / _8r);
        
        for (auto _iter = 2; _iter-- != 0; )
        {
            bool_type _diff = false;
            vals_type _oc = _hc;
    /*----------------------------- subcells about centre */
            if (eikonal_pyra_3d (
                _p1, _p2, _p3, _p4, _pc,
                eikonal_subindex(_si, 0, 1, 2, 3, 8),
                _o1, _o2, _o3, _o4, _oc,
                _h1, _h2, _h3, _h4, _hc,
                _g1, _g2, _g3, _g4, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
                
            if (eikonal_pyra_3d (
                _p1, _p2, _p6, _p5, _pc,
                eikonal_subindex(_si, 0, 1, 5, 4, 8),
                _o1, _o2, _o6, _o5, _oc,
                _h1, _h2, _h6, _h5, _hc,
                _g1, _g2, _g6, _g5, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
                
            if (eikonal_pyra_3d (
                _p2, _p3, _p7, _p6, _pc,
                eikonal_subindex(_si, 1, 2, 6, 5, 8),
                _o2, _o3, _o7, _o6, _oc,
                _h2, _h3, _h7, _h6, _hc,
                _g2, _g3, _g7, _g6, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
                
            if (eikonal_pyra_3d (
                _p3, _p4, _p8, _p7, _pc,
                eikonal_subindex(_si, 2, 3, 7, 6, 8),
                _o3, _o4, _o8, _o7, _oc,
                _h3, _h4, _h8, _h7, _hc,
                _g3, _g4, _g8, _g7, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p4, _p1, _p5, _p8, _pc,
                eikonal_subindex(_si, 3, 0, 4, 7, 8),
                _o4, _o1, _o5, _o8, _oc,
                _h4, _h1, _h5, _h8, _hc,
                _g4, _g1, _g5, _g8, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p8, _p7, _p6, _p5, _pc,
                eikonal_subindex(_si, 7, 6, 5, 4, 8),
                _o8, _o7, _o6, _o5, _oc,
                _h8, _h7, _h6, _h5, _hc,
                _g8, _g7, _g6, _g5, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }

            if (!_diff ) break;
        }

        return ( _clip ) ;
    }


    }

#   endif//__HFUN_CLIP_K__



