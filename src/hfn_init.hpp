
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The 
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any 
     * way whatsoever.  This code is provided "as-is" to be 
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 14 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
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
        
        if (_jcfg._hfun_scal == 
        jcfg_data::hfun_scal::relative)
        {
        
        if (_geom._ndim == +2 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _scal  = (real_type) +0. ;
            _scal += 
        _geom._euclidean_mesh_2d._bmax[0] -
        _geom._euclidean_mesh_2d._bmin[0] ;
            _scal += 
        _geom._euclidean_mesh_2d._bmax[1] -
        _geom._euclidean_mesh_2d._bmin[1] ;
            
            _scal /= (real_type) +2. ;
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _scal  = (real_type) +0. ;
            _scal += 
        _geom._euclidean_mesh_3d._bmax[0] -
        _geom._euclidean_mesh_3d._bmin[0] ;
            _scal += 
        _geom._euclidean_mesh_3d._bmax[1] -
        _geom._euclidean_mesh_3d._bmin[1] ;
            _scal += 
        _geom._euclidean_mesh_3d._bmax[2] -
        _geom._euclidean_mesh_3d._bmin[2] ;
            
            _scal /= (real_type) +3. ;
        }
        else
        if (_geom._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _scal  = (real_type) +0. ;
            _scal += 
        _geom._ellipsoid_mesh_3d._bmax[0] -
        _geom._ellipsoid_mesh_3d._bmin[0] ;
            _scal += 
        _geom._ellipsoid_mesh_3d._bmax[1] -
        _geom._ellipsoid_mesh_3d._bmin[1] ;
            _scal += 
        _geom._ellipsoid_mesh_3d._bmax[2] -
        _geom._ellipsoid_mesh_3d._bmin[2] ;
            
            _scal /= (real_type) +3. ;
        }
        
        }
         
    /*--------------------------------- push GEOM scaling */   
        if (_hfun._ndim == +0)
        {
    /*--------------------------------- constant-value-kd */
            _hfun._constant_value_kd.
                _hval = _scal*_jcfg._hfun_hmax ;
        }
        else 
        if (_hfun._ndim == +2)
        {
        if (_hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
        for (auto _node  = _hfun.
        _euclidean_mesh_2d._mesh._set1.head();
                  _node != _hfun.
        _euclidean_mesh_2d._mesh._set1.tend();
                ++_node  )
        {
            if (_node->mark() < 0) continue;
            
            _node->hval() = 
                _scal  * _node->hval();
            
            _node->hval() = 
                std::min(_node->hval(), 
                    _scal *_jcfg._hfun_hmax) ;
                
            _node->hval() = 
                std::max(_node->hval(), 
                    _scal *_jcfg._hfun_hmin) ;
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
           *_iter = _scal *  *_iter;
            
           *_iter = std::min(*_iter, 
                    _scal *_jcfg._hfun_hmax) ;
                
           *_iter = std::max(*_iter, 
                    _scal *_jcfg._hfun_hmin) ;
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
        for (auto _node  = _hfun.
        _euclidean_mesh_3d._mesh._set1.head();
                  _node != _hfun.
        _euclidean_mesh_3d._mesh._set1.tend();
                ++_node  )
        {
            if (_node->mark() < 0) continue;

            _node->hval() = 
                _scal  * _node->hval();
            
            _node->hval() = 
                std::min(_node->hval(), 
                    _scal *_jcfg._hfun_hmax) ;
                
            _node->hval() = 
                std::max(_node->hval(), 
                    _scal *_jcfg._hfun_hmin) ;
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
           *_iter = _scal *  *_iter;
            
           *_iter = std::min(*_iter, 
                    _scal *_jcfg._hfun_hmax) ;
                
           *_iter = std::max(*_iter, 
                    _scal *_jcfg._hfun_hmin) ;
        }
    
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
           *_iter = _scal *  *_iter;
            
           *_iter = std::min(*_iter, 
                    _scal *_jcfg._hfun_hmax) ;
                
           *_iter = std::max(*_iter, 
                    _scal *_jcfg._hfun_hmin) ;
        }
    
        }
    
     
        return (  _errv ) ;
    }


#   endif   //__HFN_INIT__



