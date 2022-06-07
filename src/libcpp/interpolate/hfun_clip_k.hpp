
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 30 Mar., 2022
     *
     * Copyright 2013-2022
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
    __normal_call bool_type EIKONAL_edge_2d (
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
    __normal_call bool_type EIKONAL_edge_3d (
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
            std::numeric_limits
                <vals_type>::infinity() ;

        real_type _tt[2] ;
        if (!math::
        polyroots(_at, _bt, _ct, _tt) )
        {
    /*-------------------------- test flow along boundary */
        return EIKONAL_edge_2d (
                _p1, _h1, _gg,
                _p3, _h3, _gg)

             | EIKONAL_edge_2d (
                _p2, _h2, _gg,
                _p3, _h3, _gg) ;
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
            std::numeric_limits
                <vals_type>::infinity() ;

        real_type _tt[2] ;
        if (!math::
        polyroots(_at, _bt, _ct, _tt) )
        {
    /*-------------------------- test flow along boundary */
        return EIKONAL_edge_3d (
                _p1, _h1, _gg,
                _p3, _h3, _gg)

             | EIKONAL_edge_3d (
                _p2, _h2, _gg,
                _p3, _h3, _gg) ;
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
    __normal_call bool_type eikonal_edge_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type  _g1 ,
                  vals_type  _g2
        )
    {
    /*---------------------- limit h-values within EDGE-2 */
        bool_type _clip = false ;

        vals_type _h0 = std::min(_h1,_h2) ;

        if (_h2 > _h0)
    /*--------------------------------- 1st node ordering */
            if (EIKONAL_edge_2d (
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ) )
                _clip =  true ;

        if (_h1 > _h0)
    /*--------------------------------- 2nd node ordering */
            if (EIKONAL_edge_2d (
                _p2, _h2, _g2 ,
                _p1, _h1, _g1 ) )
                _clip =  true ;

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_edge_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
                  vals_type& _h1 ,
                  vals_type& _h2 ,
                  vals_type  _g1 ,
                  vals_type  _g2
        )
    {
    /*---------------------- limit h-values within EDGE-2 */
        bool_type _clip = false ;

        vals_type _h0 = std::min(_h1,_h2) ;

        if (_h2 > _h0)
    /*--------------------------------- 1st node ordering */
            if (EIKONAL_edge_3d (
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ) )
                _clip =  true ;

        if (_h1 > _h0)
    /*--------------------------------- 2nd node ordering */
            if (EIKONAL_edge_3d (
                _p2, _h2, _g2 ,
                _p1, _h1, _g1 ) )
                _clip =  true ;

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_tria_2d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
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

        vals_type _h0 =
        std::min( _h3, std::min(_h1,_h2)) ;

        if (_h3 > _h0)
    /*--------------------------------- 1st node ordering */
            if (EIKONAL_tria_2d (
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ,
                _p3, _h3, _g3 ) )
                _clip =  true ;

        if (_h1 > _h0)
    /*--------------------------------- 2nd node ordering */
            if (EIKONAL_tria_2d (
                _p2, _h2, _g2 ,
                _p3, _h3, _g3 ,
                _p1, _h1, _g1 ) )
                _clip =  true ;

        if (_h2 > _h0)
    /*--------------------------------- 3rd node ordering */
            if (EIKONAL_tria_2d (
                _p3, _h3, _g3 ,
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ) )
                _clip =  true ;

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_tria_3d (
    __const_ptr  (real_type) _p1 ,
    __const_ptr  (real_type) _p2 ,
    __const_ptr  (real_type) _p3 ,
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

        vals_type _h0 =
        std::min( _h3, std::min(_h1,_h2)) ;

        if (_h3 > _h0)
    /*--------------------------------- 1st node ordering */
            if (EIKONAL_tria_3d (
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ,
                _p3, _h3, _g3 ) )
                _clip =  true ;

        if (_h1 > _h0)
    /*--------------------------------- 2nd node ordering */
            if (EIKONAL_tria_3d (
                _p2, _h2, _g2 ,
                _p3, _h3, _g3 ,
                _p1, _h1, _g1 ) )
                _clip =  true ;

        if (_h2 > _h0)
    /*--------------------------------- 3rd node ordering */
            if (EIKONAL_tria_3d (
                _p3, _h3, _g3 ,
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ) )
                _clip =  true ;

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
            if (_h1 >=(real_type)+0. &&
                _h2 >=(real_type)+0. &&
                _h3 >=(real_type)+0. )
            if (eikonal_tria_2d (
                _p1, _p2, _p3 ,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (_h1 >=(real_type)+0. &&
                _h3 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (eikonal_tria_2d (
                _p1, _p3, _p4 ,
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
            if (_h1 >=(real_type)+0. &&
                _h2 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (eikonal_tria_2d (
                _p1, _p2, _p4 ,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (_h2 >=(real_type)+0. &&
                _h3 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (eikonal_tria_2d (
                _p2, _p3, _p4 ,
                _h2, _h3, _h4 ,
                _g2, _g3, _g4 ) )
                _clip =  true ;
        }

        if (!_okay)
        {
    /*--------------------------------- a degenerate quad */
            if (_h1 >=(real_type)+0. &&
                _h2 >=(real_type)+0. )
            if (EIKONAL_edge_2d (
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ) )
                _clip =  true ;

            if (_h2 >=(real_type)+0. &&
                _h3 >=(real_type)+0. )
            if (EIKONAL_edge_2d (
                _p2, _h2, _g2 ,
                _p3, _h3, _g3 ) )
                _clip =  true ;

            if (_h3 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (EIKONAL_edge_2d (
                _p3, _h3, _g3 ,
                _p4, _h4, _g4 ) )
                _clip =  true ;

            if (_h4 >=(real_type)+0. &&
                _h1 >=(real_type)+0. )
            if (EIKONAL_edge_2d (
                _p4, _h4, _g4 ,
                _p1, _h1, _g1 ) )
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
            if (_h1 >=(real_type)+0. &&
                _h2 >=(real_type)+0. &&
                _h3 >=(real_type)+0. )
            if (eikonal_tria_3d (
                _p1, _p2, _p3 ,
                _h1, _h2, _h3 ,
                _g1, _g2, _g3 ) )
                _clip =  true ;

            if (_h1 >=(real_type)+0. &&
                _h3 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (eikonal_tria_3d (
                _p1, _p3, _p4 ,
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
            if (_h1 >=(real_type)+0. &&
                _h2 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (eikonal_tria_3d (
                _p1, _p2, _p4 ,
                _h1, _h2, _h4 ,
                _g1, _g2, _g4 ) )
                _clip =  true ;

            if (_h2 >=(real_type)+0. &&
                _h3 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (eikonal_tria_3d (
                _p2, _p3, _p4 ,
                _h2, _h3, _h4 ,
                _g2, _g3, _g4 ) )
                _clip =  true ;
        }

        if (!_okay)
        {
    /*--------------------------------- a degenerate quad */
            if (_h1 >=(real_type)+0. &&
                _h2 >=(real_type)+0. )
            if (EIKONAL_edge_3d (
                _p1, _h1, _g1 ,
                _p2, _h2, _g2 ) )
                _clip =  true ;

            if (_h2 >=(real_type)+0. &&
                _h3 >=(real_type)+0. )
            if (EIKONAL_edge_3d (
                _p2, _h2, _g2 ,
                _p3, _h3, _g3 ) )
                _clip =  true ;

            if (_h3 >=(real_type)+0. &&
                _h4 >=(real_type)+0. )
            if (EIKONAL_edge_3d (
                _p3, _h3, _g3 ,
                _p4, _h4, _g4 ) )
                _clip =  true ;

            if (_h4 >=(real_type)+0. &&
                _h1 >=(real_type)+0. )
            if (EIKONAL_edge_3d (
                _p4, _h4, _g4 ,
                _p1, _h1, _g1 ) )
                _clip =  true ;
        }

        return ( _clip ) ;
    }

    template <
        typename  real_type ,
        typename  vals_type
             >
    __normal_call bool_type eikonal_tria_3d (
    __const_ptr  (real_type) /*_p1*/ ,
    __const_ptr  (real_type) /*_p2*/ ,
    __const_ptr  (real_type) /*_p3*/ ,
    __const_ptr  (real_type) /*_p4*/ ,
                  vals_type& /*_h1*/ ,
                  vals_type& /*_h2*/ ,
                  vals_type& /*_h3*/ ,
                  vals_type& /*_h4*/ ,
                  vals_type  /*_g1*/ ,
                  vals_type  /*_g2*/ ,
                  vals_type  /*_g3*/ ,
                  vals_type  /*_g4*/
        )
    {
        return false ;
    }


    }

#   endif//__HFUN_CLIP_K__



