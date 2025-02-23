
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
     * Last updated: 20 Oct., 2024
     *
     * Copyright 2013-2024
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
        vals_type _gg =
            (_g1 + _g2) / (vals_type)2. ;

    /*---------------------- form "limited" extrap. to P2 */
        vals_type _hn = (vals_type) (
            _h1 + _gg *  std::sqrt(
        geometry::lensqr_2d(_p1, _p2))) ;

        if (_hn < _h2)
        {
            _h2 = _hn; return  true ;
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
        vals_type _gg =
            (_g1 + _g2) / (vals_type)2. ;

    /*---------------------- form "limited" extrap. to P2 */
        vals_type _hn = (vals_type) (
            _h1 + _gg *  std::sqrt(
        geometry::lensqr_3d(_p1, _p2))) ;

        if (_hn < _h2)
        {
            _h2 = _hn; return  true ;
        }
        else/*no clip*/return false ;
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
        vals_type _gg =
            (_g1+_g2+_g3)/(vals_type)3. ;

    /*---------------------- form "limited" extrap. to P3 */
        real_type _d1[2];
        real_type _d3[2];
        geometry::
            vector_2d (_p1, _p2, _d1) ;
        geometry::
            vector_2d (_p3, _p1, _d3) ;

        vals_type _dh =_h2- _h1;

        real_type _aa =
        geometry::   dot_2d(_d1, _d3) ;
        real_type _bb =
        geometry::lensqr_2d(_d1) ;
        real_type _cc =
        geometry::lensqr_2d(_d3) ;

        real_type _rr =
            std::pow(_dh/_gg, +2);
        real_type _at =
            _bb * _bb - _rr * _bb;
        real_type _bt =
            +2. * _aa * (_bb-_rr);
        real_type _ct =
            _aa * _aa - _rr * _cc;

        vals_type _hn =
        std::numeric_limits<vals_type>::infinity() ;

        real_type _tt[2] ;
        if (!math::
        polyroots(_at, _bt, _ct, _tt) )
        {
    /*-------------------------- test flow along boundary */
        return EIKONAL_edge_2d (
                _p1, _h1, _g1,
                _p3, _h3, _g3)

             | EIKONAL_edge_2d (
                _p2, _h2, _g2,
                _p3, _h3, _g3) ;
        }
        else
        {
    /*-------------------------- test flow about interior */
        _tt[0] = std::min(
            (real_type)+1.,_tt[0]) ;
        _tt[0] = std::max(
            (real_type)+0.,_tt[0]) ;

        _tt[1] = std::min(
            (real_type)+1.,_tt[1]) ;
        _tt[1] = std::max(
            (real_type)+0.,_tt[1]) ;

        real_type _vp[2] = {
        _d3[0] +  _tt[0] * _d1[0],
        _d3[1] +  _tt[0] * _d1[1],
            } ;
        real_type _vm[2] = {
        _d3[0] +  _tt[1] * _d1[0],
        _d3[1] +  _tt[1] * _d1[1],
            } ;

        real_type _lp = std::sqrt(
            geometry::lensqr_2d(_vp)) ;

        vals_type _hp =(vals_type)(
        _h1 + _tt[0]*_dh + _gg*_lp) ;

        real_type _lm = std::sqrt(
            geometry::lensqr_2d(_vm)) ;

        vals_type _hm =(vals_type)(
        _h1 + _tt[1]*_dh + _gg*_lm) ;

        _hn = std::min(_hp, _hm) ;

        if (_hn < _h3)
        {
            _h3 = _hn; return  true ;
        }
        else/*no clip*/return false ;

        }
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
        vals_type _gg =
            (_g1+_g2+_g3)/(vals_type)3. ;

    /*---------------------- form "limited" extrap. to P3 */
        real_type _d1[3];
        real_type _d3[3];
        geometry::
            vector_3d (_p1, _p2, _d1) ;
        geometry::
            vector_3d (_p3, _p1, _d3) ;

        vals_type _dh =_h2- _h1;

        real_type _aa =
        geometry::   dot_3d(_d1, _d3) ;
        real_type _bb =
        geometry::lensqr_3d(_d1) ;
        real_type _cc =
        geometry::lensqr_3d(_d3) ;

        real_type _rr =
            std::pow(_dh/_gg, +2);
        real_type _at =
            _bb * _bb - _rr * _bb;
        real_type _bt =
            +2. * _aa * (_bb-_rr);
        real_type _ct =
            _aa * _aa - _rr * _cc;

        vals_type _hn =
        std::numeric_limits<vals_type>::infinity() ;

        real_type _tt[2] ;
        if (!math::
        polyroots(_at, _bt, _ct, _tt) )
        {
    /*-------------------------- test flow along boundary */
        return EIKONAL_edge_3d (
                _p1, _h1, _g1,
                _p3, _h3, _g3)

             | EIKONAL_edge_3d (
                _p2, _h2, _g2,
                _p3, _h3, _g3) ;
        }
        else
        {
    /*-------------------------- test flow about interior */
        _tt[0] = std::min(
            (real_type)+1.,_tt[0]) ;
        _tt[0] = std::max(
            (real_type)+0.,_tt[0]) ;

        _tt[1] = std::min(
            (real_type)+1.,_tt[1]) ;
        _tt[1] = std::max(
            (real_type)+0.,_tt[1]) ;

        real_type _vp[3] = {
        _d3[0] +  _tt[0] * _d1[0],
        _d3[1] +  _tt[0] * _d1[1],
        _d3[2] +  _tt[0] * _d1[2],
            } ;
        real_type _vm[3] = {
        _d3[0] +  _tt[1] * _d1[0],
        _d3[1] +  _tt[1] * _d1[1],
        _d3[2] +  _tt[1] * _d1[2],
            } ;

        real_type _lp = std::sqrt(
            geometry::lensqr_3d(_vp)) ;

        vals_type _hp =(vals_type)(
        _h1 + _tt[0]*_dh + _gg*_lp) ;

        real_type _lm = std::sqrt(
            geometry::lensqr_3d(_vm)) ;

        vals_type _hm =(vals_type)(
        _h1 + _tt[1]*_dh + _gg*_lm) ;

        _hn = std::min(_hp, _hm) ;

        if (_hn < _h3)
        {
            _h3 = _hn; return  true ;
        }
        else/*no clip*/return false ;

        }
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
        vals_type _gg =
            (_g1+_g2+_g3+_g4)/(vals_type)4. ;

    /*---------------------- form "limited" extrap. to P4 */
        real_type _pa[3]; 
        real_type _pb[3]; 
        real_type _pc[3];
        geometry::
            vector_3d (_p1, _p3, _pa) ;
        geometry::
            vector_3d (_p2, _p3, _pb) ;
        geometry::
            vector_3d (_p3, _p4, _pc) ;
            
        vals_type _ha =_h1- _h3;
        vals_type _hb =_h2- _h3;
            
        real_type _aa = 
        geometry::   dot_3d(_pa, _pa) ;
        real_type _ab = 
        geometry::   dot_3d(_pa, _pb) ;
        real_type _ac = 
        geometry::   dot_3d(_pa, _pc) ;
        real_type _bb = 
        geometry::   dot_3d(_pb, _pb) ;
        real_type _bc = 
        geometry::   dot_3d(_pb, _pc) ;
        real_type _cc = 
        geometry::   dot_3d(_pc, _pc) ;
        
    /*---------------------- local gradient descent soln. */
        vals_type _hm, _hn =
        std::numeric_limits<vals_type>::infinity() ;
        
        real_type _ss = (real_type) .5;
        real_type _uu = (real_type) 1. / 3. ;
        real_type _vv = (real_type) 1. / 3. ;
        real_type _ww = (real_type) 1. - _uu - _vv ;
 
    //  init. from best of base midpoint, and adj. faces
        real_type _pi[3] = {
        _p1[0]*_uu+_p2[0]*_vv+_p3[0]*_ww ,
        _p1[1]*_uu+_p2[1]*_vv+_p3[1]*_ww ,
        _p1[2]*_uu+_p2[2]*_vv+_p3[2]*_ww
            }  ;
        
        _hn =  (vals_type) (
            _h1*_uu + _h2*_vv + _h3*_ww + 
        _gg * geometry::length_3d(_pi, _p4) ) ;
        
        if (EIKONAL_tria_3d(_p1, _h1, _g1,
                            _p2, _h2, _g2,
                            _p4, _hn, _g4)  )
        {
            _uu = (real_type)0.;
            _vv = (real_type).5;
        }
                        
        if (EIKONAL_tria_3d(_p2, _h2, _g2,
                            _p3, _h3, _g3,
                            _p4, _hn, _g4)  )
        {
            _uu = (real_type).5;
            _vv = (real_type)0.;
        }
                  
        if (EIKONAL_tria_3d(_p3, _h3, _g3,
                            _p1, _h1, _g1,
                            _p4, _hn, _g4)  )
        {
            _uu = (real_type).5;
            _vv = (real_type).5;   
        }
        
        for (auto _iter = 8; _iter-- != +0; )
        {
    /*---------------------- analytical gradients + soln. */
        real_type _dd = _aa * _uu * _uu + 
                        _bb * _vv * _vv +
                   2. * _ab * _uu * _vv +
                   2. * _ac * _uu +
                   2. * _bc * _vv + _cc ;
                   
        real_type _du =
                   2. * _aa * _uu + 
                   2. * _ab * _vv + 2. *_ac ;           
        real_type _dv =
                   2. * _bb * _vv + 
                   2. * _ab * _uu + 2. *_bc ;
        
        _dd = std::sqrt(_dd);           
        _du = _ha + .5* _gg * _du / _dd ;
        _dv = _hb + .5* _gg * _dv / _dd ;
        
    /*---------------------- local line-search relaxation */
        bool _up = false;
        for (auto _line = 4; _line-- != +0; )
        {
        real_type _tt = +1. ;
        real_type _un = _uu - _ss * _du ;
        real_type _vn = _vv - _ss * _dv ;
        
    /*---------------------- keep [u,v] in base of tetra. */
        if (1. - _un - _vn < 0.)
            _tt = std::min(
        _tt, (_uu + _vv - 1.) / (_du + _dv) / _ss) ;
        
        if (_un < 0.)
            _tt = std::min(
        _tt, (0. + _uu) / (_uu - _un));
            
        if (_un > 1.)
            _tt = std::min(
        _tt, (1. - _uu) / (_un - _uu));
                
        if (_vn < 0.)
            _tt = std::min(
        _tt, (0. + _vv) / (_vv - _vn));
            
        if (_vn > 1.)
            _tt = std::min(
        _tt, (1. - _vv) / (_vn - _vv));
                
        _ss = _ss * std::max(0., _tt) ;
        
    /*---------------------- eval. extrap. from new [u,v] */
        real_type _um = _uu - _ss * _du ;
        real_type _vm = _vv - _ss * _dv ;
        
        real_type _dm = _aa * _um * _um + 
                        _bb * _vm * _vm +
                   2. * _ab * _um * _vm +
                   2. * _ac * _um +
                   2. * _bc * _vm + _cc ;
                
        _hm = (vals_type)(
            _ha * _um + _hb * _vm + _h3
                + _gg * std::sqrt(_dm)) ;
            
        if (_hm < _hn || _ss <= 0.)
        {
    /*---------------------- move on to new better optima */
            _uu = _um; _vv = _vm; _hn = _hm; 
            _up = true; break ;
        }
        else { _ss *= (real_type).5; }
        }       
    /*---------------------- done if stuck or on boundary */
        if (!_up || _ss <= 0.) break ;
    /*---------------------- update length for next outer */
        _ss *= (real_type)+2.;
        }
    /*---------------------- update if extrap. does limit */
        if (_hn < _h4)
        {
            _h4 = _hn; return  true;
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
                  vals_type  _hb ,      // current anchor
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
            if (_h2 > _hb)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_edge_2d  (
                    _p1, _h1, _g1 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;

            if (_h1 > _hb)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_edge_2d  (
                    _p2, _h2, _g2 ,
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
                  vals_type  _hb ,      // current anchor
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
            if (_h2 > _hb)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_edge_3d  (
                    _p1, _h1, _g1 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;

            if (_h1 > _hb)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_edge_3d  (
                    _p2, _h2, _g2 ,
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
                  vals_type  _hb ,      // current anchor
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

        if (_h1 >=(vals_type)+0. &&
            _h2 >=(vals_type)+0. &&
            _h3 >=(vals_type)+0. )
        {
            if (_h3 > _hb)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_tria_2d  (
                    _p1, _h1, _g1 ,
                    _p2, _h2, _g2 ,
                    _p3, _h3, _g3 )  )
                    _clip =  true ;

            if (_h1 > _hb)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_tria_2d  (
                    _p2, _h2, _g2 ,
                    _p3, _h3, _g3 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;

            if (_h2 > _hb)
    /*--------------------------------- 3rd node ordering */
                if (EIKONAL_tria_2d  (
                    _p3, _h3, _g3 ,
                    _p1, _h1, _g1 ,
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
                  vals_type  _hb ,      // current anchor
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

        if (_h1 >=(vals_type)+0. &&
            _h2 >=(vals_type)+0. &&
            _h3 >=(vals_type)+0. )
        {
            if (_h3 > _hb)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_tria_3d  (
                    _p1, _h1, _g1 ,
                    _p2, _h2, _g2 ,
                    _p3, _h3, _g3 )  )
                    _clip =  true ;

            if (_h1 > _hb)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_tria_3d  (
                    _p2, _h2, _g2 ,
                    _p3, _h3, _g3 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;

            if (_h2 > _hb)
    /*--------------------------------- 3rd node ordering */
                if (EIKONAL_tria_3d  (
                    _p3, _h3, _g3 ,
                    _p1, _h1, _g1 ,
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
                  vals_type  _hb ,  // current anchor
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

        if (_h1 >=(vals_type)+0. &&
            _h2 >=(vals_type)+0. &&
            _h3 >=(vals_type)+0. &&
            _h4 >=(vals_type)+0. )
        {
            if (_h4 > _hb)
    /*--------------------------------- 1st node ordering */
                if (EIKONAL_tria_3d  (
                    _p1, _h1, _g1 ,
                    _p2, _h2, _g2 ,
                    _p3, _h3, _g3 ,
                    _p4, _h4, _g4 )  )
                    _clip =  true ;

            if (_h1 > _hb)
    /*--------------------------------- 2nd node ordering */
                if (EIKONAL_tria_3d  (
                    _p2, _h2, _g2 ,
                    _p3, _h3, _g3 ,
                    _p4, _h4, _g4 ,
                    _p1, _h1, _g1 )  )
                    _clip =  true ;

            if (_h2 > _hb)
    /*--------------------------------- 3rd node ordering */
                if (EIKONAL_tria_3d  (
                    _p3, _h3, _g3 ,
                    _p1, _h1, _g1 ,
                    _p4, _h4, _g4 ,
                    _p2, _h2, _g2 )  )
                    _clip =  true ;
                    
            if (_h3 > _hb)
    /*--------------------------------- 4th node ordering */
                if (EIKONAL_tria_3d  (
                    _p1, _h1, _g1 ,
                    _p2, _h2, _g2 ,
                    _p4, _h4, _g4 ,
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
                  vals_type  _hb ,      // current anchor
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
                _p1, _p2, _p3 , _hb,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p1, _p3, _p4 , _hb,
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
                _p1, _p2, _p4 , _hb,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p2, _p3, _p4 , _hb,
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
                _p1, _p2, _hb ,
                _h1, _h2, _g1 , _g2) )
                _clip =  true ;

            if (_h2 >=(vals_type)+0. &&
                _h3 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p2, _p3, _hb ,
                _h2, _h3, _g2 , _g3) )
                _clip =  true ;

            if (_h3 >=(vals_type)+0. &&
                _h4 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p3, _p4, _hb ,
                _h3, _h4, _g3 , _g4) )
                _clip =  true ;

            if (_h4 >=(vals_type)+0. &&
                _h1 >=(vals_type)+0. )
            if (eikonal_edge_2d (
                _p4, _p1, _hb ,
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
                  vals_type  _hb ,      // current anchor
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
                _p1, _p2, _p3 , _hb,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p1, _p3, _p4 , _hb,
                _h1, _h3, _h4 ,
                _g1, _g3, _g4 ) )
                _clip =  true ;
        
    /*--------------------------------- 2nd tria ordering */
            if (eikonal_tria_2d (
                _p1, _p2, _p4 , _hb,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_2d (
                _p2, _p3, _p4 , _hb,
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
                  vals_type  _hb ,      // current anchor
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
                _p1, _p2, _p3 , _hb,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p1, _p3, _p4 , _hb,
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
                _p1, _p2, _p4 , _hb,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p2, _p3, _p4 , _hb,
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
                _p1, _p2, _hb ,
                _h1, _h2, _g1 , _g2) )
                _clip =  true ;

            if (_h2 >=(vals_type)+0. &&
                _h3 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p2, _p3, _hb ,
                _h2, _h3, _g2 , _g3) )
                _clip =  true ;

            if (_h3 >=(vals_type)+0. &&
                _h4 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p3, _p4, _hb ,
                _h3, _h4, _g3 , _g4) )
                _clip =  true ;

            if (_h4 >=(vals_type)+0. &&
                _h1 >=(vals_type)+0. )
            if (eikonal_edge_3d (
                _p4, _p1, _hb ,
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
                  vals_type  _hb ,      // current anchor
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
                _p1, _p2, _p3 , _hb,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p1, _p3, _p4 , _hb,
                _h1, _h3, _h4 ,
                _g1, _g3, _g4 ) )
                _clip =  true ;
        
    /*--------------------------------- 2nd tria ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p4 , _hb,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p2, _p3, _p4 , _hb,
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
    __inline_call bool_type eikonal_pyra_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
    __const_ptr  (real_type) _p4 ,
    __const_ptr  (real_type) _p5 ,
                  vals_type  _hb ,      // current anchor
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
                _p1, _p2, _p3, _p5, _hb,
                _h1, _h2, _h3, _h5,
                _g1, _g2, _g3, _g5) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p1, _p3, _p4, _p5, _hb,
                _h1, _h3, _h4, _h5,
                _g1, _g3, _g4, _g5) )
                _clip =  true ;
        
    /*--------------------------------- 2nd base ordering */
            if (eikonal_tria_3d (
                _p1, _p2, _p4, _p5, _hb,
                _h1, _h2, _h4, _h5,
                _g1, _g2, _g4, _g5) )
                _clip =  true ;

            if (eikonal_tria_3d (
                _p2, _p3, _p4, _p5, _hb,
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
                  vals_type  _hb ,      // current anchor
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
            
        vals_type constexpr _6v = (vals_type) +6. ;
        vals_type _hc = (
         _h1 + _h2 + _h3 + _h4 + _h5 + _h6) / _6v ;
        vals_type _gc = (
         _g1 + _g2 + _g3 + _g4 + _g5 + _g6) / _6v ;
         
         for (auto _iter = 4; _iter-- != 0; )
        {
            bool_type _diff = false;
    /*----------------------------- subcells about centre */
            if (eikonal_tria_3d (
                _p1, _p2, _p3, _pc, _hb,
                _h1, _h2, _h3, _hc,
                _g1, _g2, _g3, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
    
            if (eikonal_pyra_3d (
                _p1, _p2, _p5, _p4, _pc, _hb,
                _h1, _h2, _h5, _h4, _hc,
                _g1, _g2, _g5, _g4, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p2, _p3, _p6, _p5, _pc, _hb,
                _h2, _h3, _h6, _h5, _hc,
                _g2, _g3, _g6, _g5, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p3, _p1, _p4, _p6, _pc, _hb,
                _h3, _h1, _h4, _h6, _hc,
                _g3, _g1, _g4, _g6, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_tria_3d (
                _p6, _p5, _p4, _pc, _hb,
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
                  vals_type  _hb ,      // current anchor
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
            
        vals_type constexpr _8v = (vals_type)+8.;
        vals_type _hc = (
         _h1 + _h2 + _h3 + _h4 + _h5 + _h6 + _h7 + _h8
            ) / _8v ;
        vals_type _gc = (
         _g1 + _g2 + _g3 + _g4 + _g5 + _g6 + _g7 + _g8
            ) / _8v ;
        
        for (auto _iter = 4; _iter-- != 0; )
        {
            bool_type _diff = false;
    /*----------------------------- subcells about centre */
            if (eikonal_pyra_3d (
                _p1, _p2, _p3, _p4, _pc, _hb,
                _h1, _h2, _h3, _h4, _hc,
                _g1, _g2, _g3, _g4, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
                
            if (eikonal_pyra_3d (
                _p1, _p2, _p6, _p5, _pc, _hb,
                _h1, _h2, _h6, _h5, _hc,
                _g1, _g2, _g6, _g5, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
                
            if (eikonal_pyra_3d (
                _p2, _p3, _p7, _p6, _pc, _hb,
                _h2, _h3, _h7, _h6, _hc,
                _g2, _g3, _g7, _g6, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
                
            if (eikonal_pyra_3d (
                _p3, _p4, _p8, _p7, _pc, _hb,
                _h3, _h4, _h8, _h7, _hc,
                _g3, _g4, _g8, _g7, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p4, _p1, _p5, _p8, _pc, _hb,
                _h4, _h1, _h5, _h8, _hc,
                _g4, _g1, _g5, _g8, _gc) )
            {
                _diff =  true ; _clip =  true ;
            }
            
            if (eikonal_pyra_3d (
                _p8, _p7, _p6, _p5, _pc, _hb,
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



