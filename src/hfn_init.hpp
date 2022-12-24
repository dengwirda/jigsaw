
    /*
    --------------------------------------------------------
     * HFN-INIT: init HFUN data via GEOM/CFG., etc.
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
     * Last updated: 06 Jun., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFN_INIT__
#   define __HFN_INIT__

    /*
    --------------------------------------------------------
     * INIT-HFUN: default HFUN data via BBOX.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type init_hfun (
        jcfg_data &_jcfg,
        jlog_data &_jlog,
        geom_data &_geom,
        hfun_data &_hfun
        )
    {
        iptr_type _errv = __no_error ;

        __unreferenced(_jlog) ;

    /*--------------------------------- find GEOM scaling */
        real_type _scal = (real_type) +1. ;
        real_type _blen = (real_type) +0. ;

        if (_geom._ndim == +2 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _blen  = (real_type) +0. ;
            _blen +=
        _geom._euclidean_mesh_2d._bmax[0] -
        _geom._euclidean_mesh_2d._bmin[0] ;
            _blen +=
        _geom._euclidean_mesh_2d._bmax[1] -
        _geom._euclidean_mesh_2d._bmin[1] ;

            _blen /= (real_type) +2. ;
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _blen  = (real_type) +0. ;
            _blen +=
        _geom._euclidean_mesh_3d._bmax[0] -
        _geom._euclidean_mesh_3d._bmin[0] ;
            _blen +=
        _geom._euclidean_mesh_3d._bmax[1] -
        _geom._euclidean_mesh_3d._bmin[1] ;
            _blen +=
        _geom._euclidean_mesh_3d._bmax[2] -
        _geom._euclidean_mesh_3d._bmin[2] ;

            _blen /= (real_type) +3. ;
        }
        else
        if (_geom._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            _blen  = (real_type) +0. ;
            _blen +=
        _geom._ellipsoid_mesh_3d._bmax[0] -
        _geom._ellipsoid_mesh_3d._bmin[0] ;
            _blen +=
        _geom._ellipsoid_mesh_3d._bmax[1] -
        _geom._ellipsoid_mesh_3d._bmin[1] ;
            _blen +=
        _geom._ellipsoid_mesh_3d._bmax[2] -
        _geom._ellipsoid_mesh_3d._bmin[2] ;

            _blen /= (real_type) +3. ;
        }

    /*--------------------------------- push GEOM scaling */
        real_type _hmin , _hmax;
        if (_jcfg._hfun_scal ==
        jcfg_data::hfun_scal::relative)
        {
            _scal = _blen ;          // swap to abs scale
            _hmax =
             std::min( +1.0 , _jcfg._hfun_hmax);
            _hmin =
             std::max( +0.0 , _jcfg._hfun_hmin);
        }
        else
        {
            _hmax =
             std::min(_blen , _jcfg._hfun_hmax);
            _hmin =
             std::max( +0.0 , _jcfg._hfun_hmin);
        }

        if (_hfun._ndim == +0)
        {
    /*--------------------------------- constant-value-kd */
            _hfun.
        _constant_value_kd._hval = _scal*_hmax ;

        }
        else
        if (_hfun._ndim == +2)
        {
        if (_hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
        for (auto _iter  = _hfun.
            _euclidean_mesh_2d._hval.head();
                  _iter != _hfun.
            _euclidean_mesh_2d._hval.tend();
                ++_iter  )
        {
           *_iter*= (fp32_t) _scal ;

           *_iter = std::min(
           *_iter, (fp32_t) (_scal * _hmax)) ;
           *_iter = std::max(
           *_iter, (fp32_t) (_scal * _hmin)) ;
        }

        }
        if (_hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-2d */
        for (auto _iter  = _hfun.
            _euclidean_grid_2d._hmat.head();
                  _iter != _hfun.
            _euclidean_grid_2d._hmat.tend();
                ++_iter  )
        {
           *_iter*= (fp32_t) _scal ;

           *_iter = std::min(
           *_iter, (fp32_t) (_scal * _hmax)) ;
           *_iter = std::max(
           *_iter, (fp32_t) (_scal * _hmin)) ;
        }

        }
        }
        else
        if (_hfun._ndim == +3)
        {
        if (_hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
        for (auto _iter  = _hfun.
            _euclidean_mesh_3d._hval.head();
                  _iter != _hfun.
            _euclidean_mesh_3d._hval.tend();
                ++_iter  )
        {
           *_iter*= (fp32_t) _scal ;

           *_iter = std::min(
           *_iter, (fp32_t) (_scal * _hmax)) ;
           *_iter = std::max(
           *_iter, (fp32_t) (_scal * _hmin)) ;
        }

        }
        if (_hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-3d */
        for (auto _iter  = _hfun.
            _euclidean_grid_3d._hmat.head();
                  _iter != _hfun.
            _euclidean_grid_3d._hmat.tend();
                ++_iter  )
        {
           *_iter*= (fp32_t) _scal ;

           *_iter = std::min(
           *_iter, (fp32_t) (_scal * _hmax)) ;
           *_iter = std::max(
           *_iter, (fp32_t) (_scal * _hmin)) ;
        }

        }
        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
        for (auto _iter  = _hfun.
            _ellipsoid_mesh_3d._hval.head();
                  _iter != _hfun.
            _ellipsoid_mesh_3d._hval.tend();
                ++_iter  )
        {
           *_iter*= (fp32_t) _scal ;

           *_iter = std::min(
           *_iter, (fp32_t) (_scal * _hmax)) ;
           *_iter = std::max(
           *_iter, (fp32_t) (_scal * _hmin)) ;
        }

        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_grid)
        {
    /*--------------------------------- ellipsoid-grid-3d */
        for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._hmat.head();
                  _iter != _hfun.
            _ellipsoid_grid_3d._hmat.tend();
                ++_iter  )
        {
           *_iter*= (fp32_t) _scal ;

           *_iter = std::min(
           *_iter, (fp32_t) (_scal * _hmax)) ;
           *_iter = std::max(
           *_iter, (fp32_t) (_scal * _hmin)) ;
        }

        }


        return (  _errv ) ;
    }


#   endif   //__HFN_INIT__



