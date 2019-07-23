
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
     * Last updated: 10 July, 2019
     *
     * Copyright 2013-2019
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
    typename I ,
    typename A = allocators::basic_alloc
             >
    class geom_mesh_ellipsoid_3d : 
            public  geom_base_3d<R, I>
    {
    public  :
    
    /*-------------- analytic ellipsoidal geometry in R^3 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;
    
    typedef geom_mesh_ellipsoid_3d  <
                real_type ,
                iptr_type ,
                allocator    >      geom_type ;

    typedef geom_base_3d     <
                real_type ,
                iptr_type    >      base_type ;
    
    typedef typename 
            base_type::line_type    line_type ;
    typedef typename 
            base_type::flat_type    flat_type ;
    typedef typename 
            base_type::disc_type    disc_type ;
    typedef typename 
            base_type::ball_type    ball_type ;
            
            
    class node_type: public tria_complex_node_3<I,R>
        {
    /*------------------------------------ loc. node type */
        public  :       
        iptr_type                     _itag ;

        char_type                     _fdim ;
        char_type                     _feat ;
        char_type                     _topo ;
        
        public  :
    /*------------------------------------ "write" access */
        __inline_call iptr_type&       itag (
            )
        {   return  this->_itag ;
        }
        __inline_call char_type&       feat (
            )
        {   return  this->_feat ;
        }
        __inline_call char_type&       fdim (
            )
        {   return  this->_fdim ;
        }
        __inline_call char_type&       topo (
            )
        {   return  this->_topo ;
        }
    /*------------------------------------ "const" access */
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        __inline_call char_type const& feat (
            ) const
        {   return  this->_feat ;
        }
        __inline_call char_type const& fdim (
            ) const
        {   return  this->_fdim ;
        }
        __inline_call char_type const& topo (
            ) const
        {   return  this->_topo ;
        }
        
        } ;
        
    class edge_type: public tria_complex_edge_2<I>
        {
    /*------------------------------------ loc. edge type */
        public  :
        iptr_type                     _itag ;
    
        char_type                     _feat ;
        char_type                     _topo ;
        
        public  :
    /*------------------------------------ "write" access */
        __inline_call iptr_type&       itag (
            )
        {   return  this->_itag ;
        }
        __inline_call char_type&       topo (
            )
        {   return  this->_topo ;
        }
        __inline_call char_type&       feat (
            )
        {   return  this->_feat ;
        }
    /*------------------------------------ "const" access */
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        __inline_call char_type const& topo (
            ) const
        {   return  this->_topo ;
        }
        __inline_call char_type const& feat (
            ) const
        {   return  this->_feat ;
        }
        
        } ;
        
    typedef mesh::tria_complex_1<
                node_type,
                edge_type,
                allocator       >   mesh_type ; 
            
    typedef geom_tree::aabb_node_base_k
                                    tree_node ;

    typedef geom_tree::aabb_item_rect_k <
                    float,
                iptr_type, +3   >   tree_item ;
                
    typedef geom_tree::aabb_tree<
                tree_item, +3   ,
                tree_node,
                allocator       >   tree_type ;

    iptr_type static constexpr _nbox = +4 ;
                
    public  :
    
    containers::
        fixed_array<real_type,3>   _bmin ;
    containers::
        fixed_array<real_type,3>   _bmax ;
        
    mesh_type                      _mesh ;

    tree_type                      _ebox ;
  
    /*--------------- (x/a)**2 + (y/b)**2 + (z/c)**2 = 1. */
        
    real_type                      _radA = 1. ;
    real_type                      _radB = 1. ;
    real_type                      _radC = 1. ;
    
    real_type                      _rEPS ;

    public  :
    
    __inline_call void_type toR3 (
      __const_ptr(real_type) _apos ,
      __write_ptr(real_type) _ppos        
        ) const
    {
    /*------------ helper: convert from S^2 to R^3 coord. */
        _ppos[0] = this->_radA *
            std::cos( _apos[0] ) * 
            std::cos( _apos[1] ) ;

        _ppos[1] = this->_radB *
            std::sin( _apos[0] ) * 
            std::cos( _apos[1] ) ;
        
        _ppos[2] = this->_radC *
            std::sin( _apos[1] ) ;
    }

    __inline_call void_type toS2 (
      __const_ptr(real_type) _ppos ,
      __write_ptr(real_type) _apos        
        ) const
    {
    /*------------ helper: convert from R^3 to S^2 coord. */
        real_type _xmul = 
            _ppos[0] * this->_radB ;
        real_type _ymul = 
            _ppos[1] * this->_radA ;
        real_type _zrat = 
            _ppos[2] / this->_radC ;

        _zrat = std::min(+1.,_zrat);
        _zrat = std::max(-1.,_zrat);

        _apos[0]= std::atan2(_ymul, 
                             _xmul);
        _apos[1]= std::asin (_zrat);
    }

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

    /*--------------------------- form AABB for full obj. */
        this->_bmin[0] = -this->_radA ;
        this->_bmin[1] = -this->_radB ;
        this->_bmin[2] = -this->_radC ;
    
        this->_bmax[0] = +this->_radA ;
        this->_bmax[1] = +this->_radB ;
        this->_bmax[2] = +this->_radC ;  
        
    /*--------------------------- convert to R^3 coord.'s */
        for (auto _iter  = 
             this->_mesh._set1.head() ;
                  _iter !=
             this->_mesh._set1.tend() ;
                ++_iter  )
        {
            if (_iter->mark() < 0)
                continue ;
            
            real_type _apos[2];
            real_type _ppos[3];

            _apos[0] = _iter->pval(0) ;
            _apos[1] = _iter->pval(1) ;

            toR3(_apos, _ppos);

            _iter->pval(0) = _ppos[0] ;
            _iter->pval(1) = _ppos[1] ;
            _iter->pval(2) = _ppos[2] ;
        }

    /*--------------------------- form rel.-tol. for prj. */
        this->_rEPS = std::sqrt  (
            std::numeric_limits
                <real_type>::epsilon()) ;
                
        real_type _rBAR ; 
        _rBAR  = (real_type) +0. ;
        _rBAR += this->_radA ;
        _rBAR += this->_radB ;
        _rBAR += this->_radC ;
        _rBAR /= (real_type) +3. ;
  
        this->_rEPS *= _rBAR ;
  
    /*--------------------------- init. AABB for arc-seg. */   
        containers::
        block_array<tree_item> _bbox;
       
        iptr_type  _inum = +0;
       
        for (auto _iter  = 
             this->_mesh._set2.head() ;
                  _iter !=
             this->_mesh._set2.tend() ;
                ++_iter, ++_inum )
        {
            if (_iter->mark() < 0)
                continue ;

            iptr_type _enod[ +2] ;
            _enod[0] = _iter->node(0) ;
            _enod[1] = _iter->node(1) ;
      
            tree_item _tdat ;      
            _tdat.ipos() = _inum ;
        
            make_aabb (
           &this->_mesh.
            _set1[_enod[ 0 ]].pval(0) ,
           &this->_mesh.
            _set1[_enod[ 1 ]].pval(0) ,
           &_tdat .pmin( 0 ), 
           &_tdat .pmax( 0 )) ;
            
            _bbox.push_tail(_tdat) ;
        }
        
        this->_ebox.load (
            _bbox.head(), 
            _bbox.tend(),this->_nbox) ;
    }
    
    /*
    --------------------------------------------------------
     * HAVE-FEAT: TRUE if has k-dim. FEAT.
    --------------------------------------------------------
     */
    
    __inline_call bool_type have_feat (
        iptr_type _fdim
        )
    {
        if (_fdim == +2)
        return ( true  ) ;
        else
        if (_fdim == +1)
        return ! this->_ebox.empty () ;
        else
        return ( false ) ;
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
        
        real_type  _la = 
       (real_type)std::atan(0.5) * 1. ;
        
        real_type  _lo = 2.*_pi / 10. ; 
        
        if (_mesh._tria.
                _nset.count() <= +8 )
        {
    /*--------------------------- init. reg.-icosahedron */
        real_type  _ppos[3] ;
        iptr_type  _inod;
        _ppos[0] = this->_radA * 
        std::cos(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)+0.5) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)+0.5) ;
        _ppos[2] = this->_radC * 
        std::sin(_pi*(real_type)+0.5) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::cos(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)-0.5) ;
        _ppos[1] = this->_radB * 
        std::sin(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)-0.5) ;
        _ppos[2] = this->_radC * 
        std::sin(_pi*(real_type)-0.5) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+0.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+0.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)+1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+1.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+1.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)-1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+2.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+2.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)+1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
     
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+3.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+3.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)-1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
     
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+4.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+4.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)+1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
     
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+5.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+5.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)-1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;   
        
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+6.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+6.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)+1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;     
        
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+7.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+7.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)-1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
        
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+8.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+8.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)+1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
   
        _ppos[0] = this->_radA * 
        std::cos(_lo*(real_type)+9.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB * 
        std::sin(_lo*(real_type)+9.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC * 
        std::sin(_la*(real_type)-1.0) ;
        _mesh.
        _tria.push_node(_ppos, _inod) ;
        _mesh.
        _tria.node(_inod)->fdim() = +0;
        _mesh.
        _tria.node(_inod)->feat() = +0;
        _mesh.
        _tria.node(_inod)->topo() = +2;
         
        }
        
    }
    
    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */
    
    __normal_call void_type make_aabb (
        real_type *_apos,
        real_type *_bpos,
        float     *_rmin,
        float     *_rmax
        )
    {
    /*- build an AABB that encloses a spheroidal arc-seg. */
        real_type _rEPS = this->_rEPS ;
    
        _rmin[0] = std::min(
            _apos[0], _bpos[0]) ;
        _rmin[1] = std::min(
            _apos[1], _bpos[1]) ;
        _rmin[2] = std::min(
            _apos[2], _bpos[2]) ;
            
        _rmax[0] = std::max(
            _apos[0], _bpos[0]) ;
        _rmax[1] = std::max(
            _apos[1], _bpos[1]) ;
        _rmax[2] = std::max(
            _apos[2], _bpos[2]) ;
            
        float     _rmid[3] = {
       (float)     +.5 * _rmin[0] +
       (float)     +.5 * _rmax[0] ,
       (float)     +.5 * _rmin[1] +
       (float)     +.5 * _rmax[1] ,
       (float)     +.5 * _rmin[2] +
       (float)     +.5 * _rmax[2] ,
            } ;
            
        float     _rlen =   +0. ;
        _rlen = std::max (
        _rlen , _rmax[0]-_rmin[0]);
        _rlen = std::max (
        _rlen , _rmax[1]-_rmin[1]);
        _rlen = std::max (
        _rlen , _rmax[2]-_rmin[2]);
    
        _rlen*= (float)     +.5 ;
        _rlen+= _rEPS ;
    
        _rmin[0] = std::min(
        _rmin[0], _rmid[0]-_rlen) ;
        _rmin[1] = std::min(
        _rmin[1], _rmid[1]-_rlen) ;
        _rmin[2] = std::min(
        _rmin[2], _rmid[2]-_rlen) ;
        
        _rmax[0] = std::max(
        _rmax[0], _rmid[0]+_rlen) ;
        _rmax[1] = std::max(
        _rmax[1], _rmid[1]+_rlen) ;
        _rmax[2] = std::max(
        _rmax[2], _rmid[2]+_rlen) ;
    }
    
    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
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
     * HELPERS: predicates for projection to geom. 
    --------------------------------------------------------
     */
    
    __normal_call void_type proj_surf (
        real_type *_psrc ,
        real_type *_pprj
        )
    {
    /*--------------------------- great-ellipse projector */
        real_type _zero[3] ;
        _zero[0] = (real_type) +.0 ;
        _zero[1] = (real_type) +.0 ;
        _zero[2] = (real_type) +.0 ;
                
        real_type _ttaa, _ttbb ;
        if (line_surf(
                _zero, _psrc, _ttaa, _ttbb) )
        {
        real_type _pmid[3] = {
        _psrc[0] * (real_type) +.5 +
        _zero[0] * (real_type) +.5 ,
        _psrc[1] * (real_type) +.5 +
        _zero[1] * (real_type) +.5 ,
        _psrc[2] * (real_type) +.5 +
        _zero[2] * (real_type) +.5
            } ;
        real_type _pdel[3] = {
        _psrc[0] * (real_type) +.5 -
        _zero[0] * (real_type) +.5 ,
        _psrc[1] * (real_type) +.5 -
        _zero[1] * (real_type) +.5 ,
        _psrc[2] * (real_type) +.5 -
        _zero[2] * (real_type) +.5
            } ;
    
        if (_ttaa > (real_type)-1.)
        {
        _pprj[0] = 
        _pmid[0] + _ttaa*_pdel[0] ;
        _pprj[1] = 
        _pmid[1] + _ttaa*_pdel[1] ;
        _pprj[2] = 
        _pmid[2] + _ttaa*_pdel[2] ;
        }
        else
        if (_ttbb > (real_type)-1.)
        {
        _pprj[0] = 
        _pmid[0] + _ttbb*_pdel[0] ;
        _pprj[1] = 
        _pmid[1] + _ttbb*_pdel[1] ;
        _pprj[2] = 
        _pmid[2] + _ttbb*_pdel[2] ;
        }
        }
    }
    
    __normal_call void_type proj_curv (
        real_type *_psrc ,
        real_type *_pprj ,
        real_type *_apos ,
        real_type *_bpos
        )
    {
    /*--------------------------- great-ellipse projector */
        real_type _zero[3] ;
        _zero[0] = (real_type) +.0 ;
        _zero[1] = (real_type) +.0 ;
        _zero[2] = (real_type) +.0 ;
        
        real_type _vnrm[3] ;
        geometry::tria_norm_3d (
            _zero , _apos, 
            _bpos , _vnrm) ;
            
        real_type _ptmp[3] ;
        geometry::proj_flat_3d (
            _psrc , _zero,
            _vnrm , _ptmp) ;
            
        proj_surf(_ptmp, _pprj);
        
        real_type _anrm[3] ;
        real_type _bnrm[3] ;
        geometry::tria_norm_3d (
            _zero , _apos, 
            _pprj , _anrm) ;
            
        geometry::tria_norm_3d (
            _zero , _pprj, 
            _bpos , _bnrm) ;
            
        real_type _asgn = 
        geometry::dot_3d(_vnrm, _anrm);
        
        real_type _bsgn = 
        geometry::dot_3d(_vnrm, _bnrm);
 
        if (_asgn < (real_type)+0.||
            _bsgn < (real_type)+0.)
        {
        
            real_type _alen = 
                geometry::
            lensqr_3d(_apos, _pprj);
            
            real_type _blen = 
                geometry::
            lensqr_3d(_bpos, _pprj);
        
            if (_alen <= _blen)
            {
            _pprj[0] = _apos[0];
            _pprj[1] = _apos[1];
            _pprj[2] = _apos[2];
            }
            else
            {
            _pprj[0] = _bpos[0];
            _pprj[1] = _bpos[1];
            _pprj[2] = _bpos[2];
            }        
        }    
    }
    
    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type ball_test (
        ball_type &_ball ,
        edge_type &_edge ,
        real_type *_apos ,
        real_type *_bpos ,
        hits_func &_hfun ,
        iptr_type &_hnum
        )
    {
    /*--------------------------- bisect curve about ball */
        bool_type _okay = false ;
        
        real_type _pprj[3] ;
            proj_curv( _ball. _pmid , 
                _pprj, 
                _apos, _bpos) ;
        
        _okay = _okay | 
            ball_kern( _ball, _edge , 
                _apos, _pprj, 0 , 
                _hfun, _hnum) ;
           
        _okay = _okay | 
            ball_kern( _ball, _edge , 
                _pprj, _bpos, 0 , 
                _hfun, _hnum) ;
    
        return  _okay ;  
    }
    
    template <
        typename  hits_func
             >
    __normal_call bool_type ball_kern (
        ball_type &_ball ,
        edge_type &_edge ,
        real_type *_apos ,
        real_type *_bpos ,
        iptr_type  _call ,
        hits_func &_hfun ,
        iptr_type &_hnum
        )
    {
    /*- calc. intersection of a ball & spheroidal arc-seg */
        geometry::hits_type 
            _htmp = geometry::face_hits ;
        
        if (_call++ > +32) return false ;
        
    /*--------------------------- call linear intersector */
        real_type _ppos[3] ;
        real_type _qpos[3] ;
        iptr_type _inum = 
       (iptr_type)geometry::ball_line_3d (
            _ball._pmid, 
            _ball._rrad, 
            _apos,_bpos, 
            _ppos,_qpos )  ;
        
        if (_inum >= (iptr_type) +1)
        {           
        
        if (_inum == (iptr_type) +1)
        {
    /*--------------------------- call hit output functor */
        real_type _pprj[3] ;
        proj_surf(_ppos, _pprj) ;

        real_type _plen =
            geometry::
        lensqr_3d(_ppos, _pprj) ;
    
        if (_plen < this->_rEPS * 
                    this->_rEPS )
        {
            _hfun(&_pprj[0], _htmp ,
                _edge.feat() ,
                _edge.topo() ,
                _edge.itag() )  ;
        
            _hnum += +1 ;
            
            return true ;
        }
        }
        else
        if (_inum == (iptr_type) +2)
        {
    /*--------------------------- call hit output functor */
        real_type _pprj[3] ;
        real_type _qprj[3] ;
        proj_surf(_ppos, _pprj) ;
        proj_surf(_qpos, _qprj) ;

        real_type _plen =
            geometry::
        lensqr_3d(_ppos, _pprj) ;
        
        real_type _qlen =
            geometry::
        lensqr_3d(_qpos, _qprj) ;
    
        real_type _xlen = 
        std::max (_plen, _qlen) ;
    
        if (_xlen < this->_rEPS * 
                    this->_rEPS )
        {
            _hfun(&_pprj[0], _htmp ,
                _edge.feat() ,
                _edge.topo() ,
                _edge.itag() )  ;
                
            _hfun(&_qprj[0], _htmp ,
                _edge.feat() ,
                _edge.topo() ,
                _edge.itag() )  ;
        
            _hnum += +2 ;
            
            return true ;
        }
        }
       
    /*--------------------------- recursive arc bisection */
        bool_type _okay = false ;

        real_type _cpos[3] = {
        _apos[0] * (real_type) +.5 +
        _bpos[0] * (real_type) +.5 ,
        _apos[1] * (real_type) +.5 +
        _bpos[1] * (real_type) +.5 ,
        _apos[2] * (real_type) +.5 +
        _bpos[2] * (real_type) +.5 ,
            } ;

        real_type _cprj[3] ;
        proj_surf(_cpos, _cprj) ; 
        
        _okay = _okay | 
            ball_kern( _ball, _edge, 
                _apos, _cprj, _call, 
                _hfun, _hnum) ;
           
        _okay = _okay | 
            ball_kern( _ball, _edge, 
                _cprj, _bpos, _call, 
                _hfun, _hnum) ;
    
        return  _okay ;   
        
        }
    
        return  false ;
    }

    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type flat_kern (
        flat_type &_flat ,
        edge_type &_edge ,
        real_type *_apos ,
        real_type *_bpos ,
        iptr_type  _call ,
        hits_func &_hfun ,
        iptr_type &_hnum
        )
    {
    /*- calc. intersection of a flat & spheroidal arc-seg */
        geometry::hits_type 
            _htmp = geometry::face_hits ;
        
        if (_call++ > +32) return false ;
        
    /*--------------------------- call linear intersector */
        real_type _xpos[3] ;
        if (geometry::line_flat_3d (
            _flat._ppos, 
            _flat._nvec, 
            _apos,_bpos, 
            _xpos, true) )
        {           
        
    /*--------------------------- call hit output functor */
        real_type _xprj[3] ;
        proj_surf(_xpos, _xprj) ;

        real_type _xlen =
            geometry::
        lensqr_3d(_xpos, _xprj) ;
    
        if (_xlen < this->_rEPS * 
                    this->_rEPS )
        {
            _hfun(&_xprj[0], _htmp ,
                _edge.feat() ,
                _edge.topo() ,
                _edge.itag() )  ;
            
            _hnum +=   +1 ;
            
            return   true ;
        }
       
    /*--------------------------- recursive arc bisection */
        bool_type _okay = false ;

        real_type _cpos[3] = {
        _apos[0] * (real_type) +.5 +
        _bpos[0] * (real_type) +.5 ,
        _apos[1] * (real_type) +.5 +
        _bpos[1] * (real_type) +.5 ,
        _apos[2] * (real_type) +.5 +
        _bpos[2] * (real_type) +.5 ,
            } ;

        real_type _cprj[3] ;
        proj_surf(_cpos, _cprj) ; 
        
        _okay = _okay | 
            flat_kern( _flat, _edge, 
                _apos, _cprj, _call, 
                _hfun, _hnum) ;
           
        _okay = _okay | 
            flat_kern( _flat, _edge, 
                _cprj, _bpos, _call, 
                _hfun, _hnum) ;
    
        return  _okay ;   
        
        }
    
        return  false ;
    }

    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type disc_kern (
        disc_type &_disc ,
        real_type *_apos ,
        real_type *_bpos ,
        iptr_type  _call ,
        hits_func &_hfun
        )
    {
    /*- calc. intersection of a disc & spheroidal surface */  
        if (is_inside(_apos) !=
            is_inside(_bpos) )
        {
        if (_call++ > +3 )
        {
    /*--------------------------- call linear intersector */
        line_type _ldat;
        _ldat._ipos[0] = _apos[0] ;
        _ldat._ipos[1] = _apos[1] ;
        _ldat._ipos[2] = _apos[2] ;
        _ldat._jpos[0] = _bpos[0] ;
        _ldat._jpos[1] = _bpos[1] ;
        _ldat._jpos[2] = _bpos[2] ;
        
        return intersect(_ldat, _hfun);
        
        }
        else
        {   
    /*--------------------------- recursive arc bisection */   
        bool_type _okay = false ;

        real_type _cpos[4] = {
       (real_type) +.5 * _apos[0] +
       (real_type) +.5 * _bpos[0] ,
       (real_type) +.5 * _apos[1] +
       (real_type) +.5 * _bpos[1] ,
       (real_type) +.5 * _apos[2] +
       (real_type) +.5 * _bpos[2]
            } ;

        real_type _cdir[4] = {
        _cpos[0] - _disc._pmid[0] ,
        _cpos[1] - _disc._pmid[1] ,
        _cpos[2] - _disc._pmid[2] ,
            } ;

        _cdir[3] = 
        geometry::length_3d(_cdir);
        _cdir[0]/= _cdir[3] ;
        _cdir[1]/= _cdir[3] ;
        _cdir[2]/= _cdir[3] ;

        _cpos[0] = _disc._pmid[0] +
        _cdir[0] * _disc._rrad ;
        _cpos[1] = _disc._pmid[1] +
        _cdir[1] * _disc._rrad ;
        _cpos[2] = _disc._pmid[2] +
        _cdir[2] * _disc._rrad ;

        _okay = _okay | 
            disc_kern( _disc, 
        _apos , _cpos, _call, _hfun) ;
           
        _okay = _okay | 
            disc_kern( _disc, 
        _cpos , _bpos, _call, _hfun) ;
    
        return  _okay ;
            
        }
        }
        
        return  false ;    
    }

    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */

    template <
    typename   hits_func
             >
    class flat_intersect
        {
    /*------------------ flat-geom intersection predicate */
        public  :
        geom_type     &_geom ;
        flat_type     &_flat ;
        
        hits_func     &_hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        public  :
        flat_intersect operator = (
            flat_intersect & 
            )                   =   delete ;
        flat_intersect operator = (
            flat_intersect&& 
            )                   =   delete ;

        public  :
    /*----------------------- construct using _src. obj. */
        __normal_call flat_intersect (
            flat_type &_fsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc )  ,
                _flat( _fsrc )  , 
                _hfun( _hsrc )
        {
            this->_hnum    = +  0  ;
            this->_find    = false ;
        }
    /*----------------------- all intersection about node */
        __normal_call  void_type operator()  (
            typename  
            tree_type::item_data *_iptr
            )
        {
            for ( ; _iptr != nullptr; 
                        _iptr = _iptr->_next )
            {
        /*--------------- flat/edge intersection test */
            iptr_type  _epos = 
                _iptr->_data.ipos() ;
            
            iptr_type _enod[2];
            _enod[0] =_geom.
                _mesh._set2[_epos].node(0) ;
            _enod[1] =_geom.
                _mesh._set2[_epos].node(1) ;
        
        /*--------------- call output function on hit */
            _geom .flat_kern  ( 
                 this->_flat,
                 this->_geom.
                _mesh ._set2[_epos] ,
                &this->_geom.
            _mesh._set1[_enod[ 0]].pval(0) , 
                &this->_geom.
            _mesh._set1[_enod[ 1]].pval(0) ,
                    + 0 ,
                 this->_hfun,
                 this->_hnum) ;
            
            this->_find = 
            this->_find | (this->_hnum!=0) ;
            }
        }
        
        } ;
        
    template <
    typename   hits_func
             >
    class ball_intersect
        {
    /*------------------ ball-geom intersection predicate */
        public  :
        geom_type     &_geom ;
        ball_type     &_ball ;
 
        real_type      _rmin[3] ;
        real_type      _rmax[3] ;
        
        hits_func     &_hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        public  :
        ball_intersect operator = (
            ball_intersect & 
            )                   =   delete ;
        ball_intersect operator = (
            ball_intersect&& 
            )                   =   delete ;

        public  :
    /*----------------------- construct using _src. obj. */
        __normal_call ball_intersect (
            ball_type &_bsrc ,
            real_type *_psrc ,
            real_type *_qsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc )  ,
                _ball( _bsrc )  , 
                _hfun( _hsrc )
        {
            this->_rmin[0] = _psrc[0] ;
            this->_rmin[1] = _psrc[1] ;
            this->_rmin[2] = _psrc[2] ;
            
            this->_rmax[0] = _qsrc[0] ;
            this->_rmax[1] = _qsrc[1] ;
            this->_rmax[2] = _qsrc[2] ;
        
            this->_hnum    = +  0  ;
            this->_find    = false ;
        }
    /*----------------------- all intersection about node */
        __normal_call  void_type operator()  (
            typename  
            tree_type::item_data *_iptr
            )
        {
            for ( ; _iptr != nullptr; 
                        _iptr = _iptr->_next )
            {
        /*--------------- ball/edge intersection test */
            iptr_type  _epos = 
                _iptr->_data.ipos() ;
            
            iptr_type _enod[2];
            _enod[0] =_geom.
                _mesh._set2[_epos].node(0) ;
            _enod[1] =_geom.
                _mesh._set2[_epos].node(1) ;
            
        /*--------------- call output function on hit */
            _geom. ball_test  ( 
                 this->_ball,
                 this->_geom.
                _mesh ._set2[_epos] ,
                &this->_geom.
            _mesh._set1[_enod[ 0]].pval(0) , 
                &this->_geom.
            _mesh._set1[_enod[ 1]].pval(0) ,
                 this->_hfun,
                 this->_hnum) ;
            
            this->_find = 
            this->_find | (this->_hnum!=0) ;
            }
        }
        
        } ;
    
    /*
    --------------------------------------------------------
     * INTERSECT: find FLAT/1-GEOM. intersections.
    --------------------------------------------------------
     */
    
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        flat_type &_flat ,
        hits_func &_hfun
        )
    {
    /*------------------ tree-flat intersection predicate */
        typedef 
        geom_tree::aabb_pred_flat_3 <
                 float , 
             iptr_type >    tree_pred ; 

    /*------------------ dual-face intersection predicate */
        typedef 
        flat_intersect <
             hits_func >    hits_pred ;
             
        float           _PPOS[3] ;
        _PPOS[0] =      _flat. _ppos[0] ;
        _PPOS[1] =      _flat. _ppos[1] ;
        _PPOS[2] =      _flat. _ppos[2] ;

        float           _NVEC[3] ;
        _NVEC[0] =      _flat. _nvec[0] ;
        _NVEC[1] =      _flat. _nvec[1] ;
        _NVEC[2] =      _flat. _nvec[2] ;

        float           _RMIN[3] ;
        _RMIN[0] =      _flat. _rmin[0] ;
        _RMIN[1] =      _flat. _rmin[1] ;
        _RMIN[2] =      _flat. _rmin[2] ;

        float           _RMAX[3] ;
        _RMAX[0] =      _flat. _rmax[0] ;
        _RMAX[1] =      _flat. _rmax[1] ;
        _RMAX[2] =      _flat. _rmax[2] ;

    /*------------------ call actual intersection testing */
        tree_pred _pred(_PPOS, _NVEC,
                        _RMIN, _RMAX) ;     
        hits_pred _func(_flat,
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: find BALL/1-GEOM. intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        ball_type &_ball ,
        hits_func &_hfun
        )
    {
    /*------------------ tree-ball intersection predicate */
        typedef 
        geom_tree::aabb_pred_ball_3 <
                 float , 
             iptr_type >    tree_pred ; 
            
    /*------------------ ball-line intersection predicate */
        typedef 
        ball_intersect <
             hits_func >    hits_pred ;

        float           _PMID[3] ;
        _PMID[0] =      _ball. _pmid[0] ;
        _PMID[1] =      _ball. _pmid[1] ;
        _PMID[2] =      _ball. _pmid[2] ;

        float           _RRAD;
        _RRAD    =      _ball. _rrad;

    /*------------------ call actual intersection testing */      
        real_type _rmin[3] = {
        _ball._pmid[0] -_ball. _rrad,
        _ball._pmid[1] -_ball. _rrad,
        _ball._pmid[2] -_ball. _rrad
            } ;
        real_type _rmax[3] = {
        _ball._pmid[0] +_ball. _rrad,
        _ball._pmid[1] +_ball. _rrad,
        _ball._pmid[2] +_ball. _rrad
            } ;
      
        tree_pred _pred(_PMID, _RRAD) ;
        hits_pred _func(_ball,
                        _rmin, _rmax,
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: find LINE/2-GEOM. intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        line_type &_line ,
        hits_func &_hfun
        )
    {
        bool_type _find = false ;

        real_type  _ipos[3] ;
        _ipos[0] = _line._ipos[ 0];
        _ipos[1] = _line._ipos[ 1];
        _ipos[2] = _line._ipos[ 2];

        real_type  _jpos[3] ;
        _jpos[0] = _line._jpos[ 0];
        _jpos[1] = _line._jpos[ 1];
        _jpos[2] = _line._jpos[ 2];

        real_type _ttaa, _ttbb;
        if (line_surf(_ipos, _jpos, _ttaa, _ttbb))
        {
        
        real_type  _pmid[3] = {
        _jpos[0] * (real_type)+.5 +
        _ipos[0] * (real_type)+.5 ,
        _jpos[1] * (real_type)+.5 +
        _ipos[1] * (real_type)+.5 ,
        _jpos[2] * (real_type)+.5 +
        _ipos[2] * (real_type)+.5
            } ;
        real_type  _pdel[3] = {
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
        
        if (_ttaa >= (real_type)-1.)
        if (_ttaa <= (real_type)+1.)
        {
            _find  =  true  ; 
            _hfun  ( _apos, _hits , 
                     _feat, _topo , 
                     _itag) ;
        }
        
        if (_ttbb >= (real_type)-1.)
        if (_ttbb <= (real_type)+1.)
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
     * INTERSECT: find DISC/2-GEOM. intersections.
    --------------------------------------------------------
     */
     
    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        disc_type &_disc ,
        real_type *_sbal ,
        hits_func &_hfun
        )
    {
        bool_type _find = false ;

        __unreferenced(_sbal) ;

        real_type  _circ[3] = {
            (real_type) +0. ,
            (real_type) +0. ,
            (real_type) +0. } ;

        real_type  _pdir[4] = {
        _disc._pmid[0] - _circ[0] ,
        _disc._pmid[1] - _circ[1] ,
        _disc._pmid[2] - _circ[2] ,
            (real_type) +0. } ;
            
        real_type  _ddir[4] ;
        geometry::cross_3d(
        _pdir, _disc._nvec, _ddir) ;
        
        _pdir[3] = 
        geometry::length_3d(_pdir) ;
        _pdir[0]/= _pdir[3] ;
        _pdir[1]/= _pdir[3] ;
        _pdir[2]/= _pdir[3] ;
        
        _ddir[3] = 
        geometry::length_3d(_ddir) ;
        _ddir[0]/= _ddir[3] ;
        _ddir[1]/= _ddir[3] ;
        _ddir[2]/= _ddir[3] ;
        
        real_type  _apos[3] = {
        _disc._pmid[0] - 
            _disc._rrad *_pdir[0] ,
        _disc._pmid[1] - 
            _disc._rrad *_pdir[1] ,
        _disc._pmid[2] - 
            _disc._rrad *_pdir[2] 
            } ;
            
        real_type  _bpos[3] = {
        _disc._pmid[0] + 
            _disc._rrad *_ddir[0] ,
        _disc._pmid[1] + 
            _disc._rrad *_ddir[1] ,
        _disc._pmid[2] + 
            _disc._rrad *_ddir[2] 
            } ;
            
        real_type  _cpos[3] = {
        _disc._pmid[0] +
            _disc._rrad *_pdir[0] ,
        _disc._pmid[1] + 
            _disc._rrad *_pdir[1] ,
        _disc._pmid[2] + 
            _disc._rrad *_pdir[2] 
            } ;
            
        real_type  _dpos[3] = {
        _disc._pmid[0] - 
            _disc._rrad *_ddir[0] ,
        _disc._pmid[1] - 
            _disc._rrad *_ddir[1] ,
        _disc._pmid[2] - 
            _disc._rrad *_ddir[2] 
            } ;
        
        _find = _find | 
            disc_kern(_disc, 
        _apos , _bpos, +0, _hfun) ;
        
        _find = _find | 
            disc_kern(_disc, 
        _bpos , _cpos, +0, _hfun) ;
        
        _find = _find | 
            disc_kern(_disc, 
        _cpos , _dpos, +0, _hfun) ;
        
        _find = _find | 
            disc_kern(_disc, 
        _dpos , _apos, +0, _hfun) ;
     
        return  _find ;
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
    
    
    
