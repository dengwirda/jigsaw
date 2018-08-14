
    /*
    --------------------------------------------------------
     * GEOM-MESH-ELLIPSOID-3: ellipsoidal geom. in R^3.
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
     * Last updated: 14 March, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __GEOM_MESH_ELLIPSOID_3__
#   define __GEOM_MESH_ELLIPSOID_3__

    namespace mesh {

    template <
    typename R ,
    typename I 
             >
    class geom_mesh_ellipsoid_3d
    {
    public  :
    
    /*-------------- analytic ellipsoidal geometry in R^3 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    
    typedef geom_mesh_ellipsoid_3d  <
                real_type ,
                iptr_type    >      geom_type ;
                
    public  :
    
    containers::
        fixed_array<real_type,3>    _bmin ;
    containers::
        fixed_array<real_type,3>    _bmax ;
  
    /*--------------- (x/a)**2 + (y/b)**2 + (z/c)**2 = 1. */
        
    real_type                       _radA ;
    real_type                       _radB ;
    real_type                       _radC ;

    public  :
    
    /*
    --------------------------------------------------------
     * INIT-GEOM: init. geometry data structures.
    --------------------------------------------------------
     */
     
    template <
        typename  geom_opts
             >
    __normal_call void_type init_geom (
        geom_opts &_opts
        )
    {
        __unreferenced(_opts) ;

        this->_bmin[0] = -this->_radA ;
        this->_bmin[1] = -this->_radB ;
        this->_bmin[2] = -this->_radC ;
    
        this->_bmax[0] = +this->_radA ;
        this->_bmax[1] = +this->_radB ;
        this->_bmax[2] = +this->_radC ;
        
    //!! blah something about coastlines...
        
    }
    
    /*
    --------------------------------------------------------
     * SEED-FEAT: init. "seed" vertex set on geom.
    --------------------------------------------------------
     */
     
    template <
        typename  mesh_type ,
        typename  geom_opts
             >
    __normal_call void_type seed_feat (
        mesh_type &_mesh ,
        geom_opts &_opts
        )
    {   
        __unreferenced(_opts) ;

        real_type  _ppos[3] ;
        iptr_type  _inod;
        _ppos[0] = (real_type) +0.0E+0;
        _ppos[1] = (real_type) +0.0E+0;
        _ppos[2] = (real_type) +0.0E+0;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +4;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +0;     
    }       
        
    /*
    --------------------------------------------------------
     * SEED-MESH: init. "seed" vertex set on geom.
    --------------------------------------------------------
     */
    
    template <
        typename  mesh_type ,
        typename  geom_opts
             >
    __normal_call void_type seed_mesh (
        mesh_type &_mesh ,
        geom_opts &_opts
        )
    {   
        __unreferenced(_opts) ;

        real_type  _pi = 
       (real_type)std::atan(1.0) * 4. ;
        
        if (_mesh._tria.
                _nset.count() <= +8 )
        {
        real_type  _ppos[3] ;
        iptr_type  _inod;
        _ppos[0] = this->_radA * 
        std::sin(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)+0.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.0) *
        std::sin(_pi*(real_type)+0.0) ;
        _ppos[2] = this->_radC * 
        std::cos(_pi*(real_type)+0.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +2;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::sin(_pi*(real_type)+1.0) *
        std::cos(_pi*(real_type)+0.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+1.0) *
        std::sin(_pi*(real_type)+0.0) ;
        _ppos[2] = this->_radC * 
        std::cos(_pi*(real_type)+1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +2;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::sin(_pi*(real_type)+0.5) *
        std::cos(_pi*(real_type)+0.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.5) *
        std::sin(_pi*(real_type)+0.0) ;
        _ppos[2] = this->_radC * 
        std::cos(_pi*(real_type)+0.5) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +2;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2; 
        
        _ppos[0] = this->_radA * 
        std::sin(_pi*(real_type)+0.5) *
        std::cos(_pi*(real_type)+0.5) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.5) *
        std::sin(_pi*(real_type)+0.5) ;
        _ppos[2] = this->_radC * 
        std::cos(_pi*(real_type)+0.5) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +2;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2; 
        
        _ppos[0] = this->_radA * 
        std::sin(_pi*(real_type)+0.5) *
        std::cos(_pi*(real_type)+1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.5) *
        std::sin(_pi*(real_type)+1.0) ;
        _ppos[2] = this->_radC * 
        std::cos(_pi*(real_type)+0.5) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +2;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;      
    
        _ppos[0] = this->_radA * 
        std::sin(_pi*(real_type)+0.5) *
        std::cos(_pi*(real_type)+1.5) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.5) *
        std::sin(_pi*(real_type)+1.5) ;
        _ppos[2] = this->_radC * 
        std::cos(_pi*(real_type)+0.5) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +2;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2; 
        }
        
    }
    
    /*
    --------------------------------------------------------
     * LINE-SURF: helper!
    --------------------------------------------------------
     */
    
    __normal_call bool_type line_surf (
        real_type *_ipos ,
        real_type *_jpos ,
        real_type &_ttii ,
        real_type &_ttjj
        )
    {

    /*- calc. intersection of a line & spheroidal surface */
    
        bool_type _find = false ;
    
        real_type _mm[3] = {
        _jpos[0] * (real_type)+.5 +
        _ipos[0] * (real_type)+.5 ,
        _jpos[1] * (real_type)+.5 +
        _ipos[1] * (real_type)+.5 ,
        _jpos[2] * (real_type)+.5 +
        _ipos[2] * (real_type)+.5
            } ;
        real_type _dd[3] = {
        _jpos[0] * (real_type)+.5 -
        _ipos[0] * (real_type)+.5 ,
        _jpos[1] * (real_type)+.5 -
        _ipos[1] * (real_type)+.5 ,
        _jpos[2] * (real_type)+.5 -
        _ipos[2] * (real_type)+.5
            } ;    
            
        real_type _ax = 
            _dd [0] / this->_radA ;
        real_type _ay = 
            _dd [1] / this->_radB ;
        real_type _az = 
            _dd [2] / this->_radC ;
            
        real_type _aa = _ax * _ax
                      + _ay * _ay
                      + _az * _az ;
        
        real_type _bx = 
            _mm [0] * _dd[0] / 
                (this->_radA * this->_radA) ;
        real_type _by = 
            _mm [1] * _dd[1] / 
                (this->_radB * this->_radB) ;
        real_type _bz = 
            _mm [2] * _dd[2] / 
                (this->_radC * this->_radC) ;
                      
        real_type _bb = 
            _bx * (real_type) +2. +
            _by * (real_type) +2. +
            _bz * (real_type) +2. ; 
                      
        real_type _cx = 
            _mm [0] / this->_radA ;
        real_type _cy = 
            _mm [1] / this->_radB ;
        real_type _cz = 
            _mm [2] / this->_radC ;
            
        real_type _cc = _cx * _cx
                      + _cy * _cy
                      + _cz * _cz 
                   - (real_type) 1.  ;
                      
        real_type _tsqr = _bb*_bb - 
            (real_type)4.*_aa*_cc ;
            
        if (_tsqr >= (real_type)+0.)
        {
            _tsqr = std::sqrt(_tsqr) ;

            _find = true ;

            _ttii = (-_bb + _tsqr) / 
            ((real_type)+2. * _aa) ;
            
            _ttjj = (-_bb - _tsqr) / 
            ((real_type)+2. * _aa) ;
        }
    
        return ( _find ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: compute poly-edge intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  poly_list ,
        typename  hits_func
             >
    __normal_call bool_type intersect (
        real_type *_ppos ,
        real_type *_nvec ,
        poly_list &_poly ,
        hits_func &_hfun
        )
    {
        __unreferenced(_ppos) ;
        __unreferenced(_nvec) ;
        __unreferenced(_poly) ;
        __unreferenced(_hfun) ;
    
        return ( false ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: compute ball-edge intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        real_type *_cmid ,
        real_type  _rsiz ,
        hits_func &_hfun
        )
    {
    
        __unreferenced(_cmid) ;
        __unreferenced(_rsiz) ;
        __unreferenced(_hfun) ;
    
        return ( false ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: compute line-face intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        real_type *_ipos ,
        real_type *_jpos ,
        hits_func &_hfun
        )
    {
        bool_type _find = false ;

        real_type _ttaa, _ttbb;
        if (line_surf(_ipos, _jpos, _ttaa, _ttbb))
        {
            real_type _pmid[3] = {
            _jpos[0] * (real_type)+.5 +
            _ipos[0] * (real_type)+.5 ,
            _jpos[1] * (real_type)+.5 +
            _ipos[1] * (real_type)+.5 ,
            _jpos[2] * (real_type)+.5 +
            _ipos[2] * (real_type)+.5
                } ;
            real_type _pdel[3] = {
            _jpos[0] * (real_type)+.5 -
            _ipos[0] * (real_type)+.5 ,
            _jpos[1] * (real_type)+.5 -
            _ipos[1] * (real_type)+.5 ,
            _jpos[2] * (real_type)+.5 -
            _ipos[2] * (real_type)+.5
                } ;
        
            real_type  _apos[3] = {
            _pmid[0] + _ttaa * _pdel[0] ,
            _pmid[1] + _ttaa * _pdel[1] ,
            _pmid[2] + _ttaa * _pdel[2]
                } ;
            real_type  _bpos[3] = {
            _pmid[0] + _ttbb * _pdel[0] ,
            _pmid[1] + _ttbb * _pdel[1] ,
            _pmid[2] + _ttbb * _pdel[2]
                } ;

            char_type _hits = 
                geometry::null_hits ;
            char_type _feat = +2;
            char_type _topo = +2;
            iptr_type _itag = +0;

            if (_ttaa >= (real_type)-1. &&
                _ttaa <= (real_type)+1. )
            {
                _find  =  true  ; 
                _hfun  ( _apos, _hits , 
                         _feat, _topo , 
                         _itag) ;
            }
            
            if (_ttbb >= (real_type)-1. &&
                _ttbb <= (real_type)+1. )
            {
                _find  =  true  ; 
                _hfun  ( _bpos, _hits , 
                         _feat, _topo , 
                         _itag) ;
            }
        }

        return ( _find ) ;
    }
      
    /*
    --------------------------------------------------------
     * INTERSECT: compute disk-face intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        real_type *_cmid ,
        real_type *_nvec ,
        real_type  _rsiz ,
        real_type *_ball ,
        hits_func &_hfun
        )
    {
        __unreferenced(_nvec) ;
        
        bool_type _find = false ;

        real_type  _bvec[3] = {
        _ball[0] - _cmid[0] ,
        _ball[1] - _cmid[1] ,
        _ball[2] - _cmid[2]
            } ;
        real_type  _blen = 
        geometry::length_3d(_bvec) ;
        _bvec[0]/= _blen ;
        _bvec[1]/= _blen ;
        _bvec[2]/= _blen ;

        _bvec[0]*= _rsiz ;
        _bvec[1]*= _rsiz ;
        _bvec[2]*= _rsiz ;

        real_type  _tpos[3] = {
        _cmid[0] + _bvec[0] ,
        _cmid[1] + _bvec[1] ,
        _cmid[2] + _bvec[2] ,
            } ;
            
        real_type  _circ[3] = {
            (real_type) +0. ,
            (real_type) +0. ,
            (real_type) +0.
            } ;

        real_type _ttaa;
        real_type _ttbb;
        if (line_surf(_circ, _tpos, _ttaa, _ttbb))
        {
            real_type _pmid[3] = {
            _tpos[0] * (real_type)+.5 +
            _circ[0] * (real_type)+.5 ,
            _tpos[1] * (real_type)+.5 +
            _circ[1] * (real_type)+.5 ,
            _tpos[2] * (real_type)+.5 +
            _circ[2] * (real_type)+.5
                } ;
            real_type _pdel[3] = {
            _tpos[0] * (real_type)+.5 -
            _circ[0] * (real_type)+.5 ,
            _tpos[1] * (real_type)+.5 -
            _circ[1] * (real_type)+.5 ,
            _tpos[2] * (real_type)+.5 -
            _circ[2] * (real_type)+.5
                } ;
        
            real_type  _apos[3] = {
            _pmid[0] + _ttaa * _pdel[0] ,
            _pmid[1] + _ttaa * _pdel[1] ,
            _pmid[2] + _ttaa * _pdel[2]
                } ;
            real_type  _bpos[3] = {
            _pmid[0] + _ttbb * _pdel[0] ,
            _pmid[1] + _ttbb * _pdel[1] ,
            _pmid[2] + _ttbb * _pdel[2]
                } ;

            char_type _hits = 
                geometry::null_hits ;
            char_type _feat = +2;
            char_type _topo = +2;
            iptr_type _itag = +0;

            if (_ttaa > (real_type)-1.)
            {
                _find  =  true  ; 
                _hfun  ( _apos, _hits , 
                         _feat, _topo , 
                         _itag) ;
            }
            
            if (_ttbb > (real_type)-1.)
            {
                _find  =  true  ; 
                _hfun  ( _bpos, _hits , 
                         _feat, _topo , 
                         _itag) ;
            }
        }

        return ( _find ) ;
    }
    
    /*
    --------------------------------------------------------
     * IS-INSIDE: TRUE if point is "inside" geometry.
    --------------------------------------------------------
     */
    
    __normal_call bool_type is_inside (
        real_type *_ppos
        )
    {
        real_type _xx = 
            _ppos[0] / this->_radA ;
        real_type _yy = 
            _ppos[1] / this->_radB ;
        real_type _zz = 
            _ppos[2] / this->_radC ;
            
        return ( _xx * _xx 
               + _yy * _yy 
               + _zz * _zz 
              <= (real_type) +1. ) ;
    }
                
    } ;
    
    
    }

#   endif   //__GEOM_MESH_ELLIPSOID_3__
    
    
    
