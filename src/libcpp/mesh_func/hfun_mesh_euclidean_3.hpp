
    /*
    --------------------------------------------------------
     * HFUN-MESH-EUCLIDEAN-kD: unstructured H(X) in R^k.
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
     * Last updated: 24 July, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_MESH_EUCLIDEAN_3__
#   define __HFUN_MESH_EUCLIDEAN_3__

    namespace mesh {

    template <
    typename R ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_mesh_euclidean_3d 
        : public hfun_base_kd <I, R> 
    {
    public  :
    
    /*------------------------- euclidean size-fun in R^3 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;
    
    typedef typename
            allocator::size_type    uint_type ;

    typedef hfun_mesh_euclidean_3d  <
            real_type ,
            iptr_type >             hfun_type ;

    typedef typename  hfun_base_kd  <
            iptr_type , 
            real_type >::hint_type  hint_type ;

    class node_type: public tria_complex_node_3<I,R>
        {
    /*------------------------------------ loc. node type */
        public  :       
        real_type                     _hval ;
        
        public  :
    /*------------------------------------ "write" access */
        __inline_call real_type&       hval (
            )
        {   return  this->_hval ;
        }
    /*------------------------------------ "const" access */
        __inline_call real_type const& hval (
            ) const
        {   return  this->_hval ;
        }
        
        } ;
        
    typedef tria_complex_edge_2<I>  edge_type ;
    typedef tria_complex_tria_3<I>  tri3_type ;
    typedef tria_complex_tria_4<I>  tri4_type ;
       
    typedef mesh::tria_complex_3<
                node_type,
                edge_type,
                tri3_type,
                tri4_type,
                allocator       >   mesh_type ;

    typedef geom_tree::aabb_node_base_k
                                    tree_node ;

    typedef geom_tree::aabb_item_rect_k <
                real_type,
                iptr_type,
                + 3             >   tree_item ;
                
    typedef geom_tree::aabb_tree<
                tree_item,
                + 3 ,
                tree_node,
                allocator       >   tree_type ; 
              
    public  :              
    
    containers::
        fixed_array<real_type,3>   _bmin ;
    containers::
        fixed_array<real_type,3>   _bmax ;
   
    real_type                      _vtol ;
   
    mesh_type                      _mesh ;
    tree_type                      _tree ;
  
    public  :

    /*
    --------------------------------------------------------
     * INIT-HFUN: make tree - check signs.
    --------------------------------------------------------
     */  
     
    __normal_call void_type init (
        )
    {
        class tria_pred
            {
    /*-------------------- TRUE if tree should index tria */
            public  :
            __inline_call 
                bool_type operator () (
                tri4_type const& _tdat
                ) const
            {   return _tdat.mark() >= 0 ;
            }
            } ;
    
    /*----------------------------- init. aabb at -+ inf. */
        for(auto _idim = 3; _idim-- != 0 ; )
        {
            this->_bmin[_idim] = 
                +std::numeric_limits<
                    real_type>::infinity() ;
                    
            this->_bmax[_idim] = 
                -std::numeric_limits<
                    real_type>::infinity() ;
        }
  
    /*----------------------------- calc. aabb for inputs */
        for (auto  _iter  = 
             this->_mesh._set1.head() ;
                   _iter !=
             this->_mesh._set1.tend() ;
                 ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
            this->_bmin[0] = std::min (
            this->_bmin[0] , 
            _iter->pval(0) ) ;
            this->_bmin[1] = std::min (
            this->_bmin[1] , 
            _iter->pval(1) ) ;
            this->_bmin[2] = std::min (
            this->_bmin[2] , 
            _iter->pval(2) ) ;
        
            this->_bmax[0] = std::max (
            this->_bmax[0] , 
            _iter->pval(0) ) ;
            this->_bmax[1] = std::max (
            this->_bmax[1] , 
            _iter->pval(1) ) ;
            this->_bmax[2] = std::max (
            this->_bmax[2] , 
            _iter->pval(2) ) ;
            }
        }
  
        real_type static const _RTOL = 
            std::pow (
        std::numeric_limits<real_type>
            ::epsilon(),(real_type).8) ;
            
        iptr_type static
        constexpr _NBOX=(iptr_type)+8  ;

        real_type _BTOL[3] ;
        _BTOL[0] = this->_bmax[0] - 
                   this->_bmin[0] ;
        _BTOL[1] = this->_bmax[1] - 
                   this->_bmin[1] ;
        _BTOL[2] = this->_bmax[2] - 
                   this->_bmin[2] ;
  
        real_type _vbox =_BTOL[0] * 
                         _BTOL[1] *
                         _BTOL[2] ;
  
        _BTOL[0]*=_RTOL ;
        _BTOL[1]*=_RTOL ;
        _BTOL[2]*=_RTOL ;

        this->_vtol = _vbox*_RTOL ;
  
    /*-------------------- make aabb-tree and init. bbox. */
        aabb_mesh(this->_mesh._set1 , 
                  this->_mesh._set4 , 
                  this->_tree,_BTOL ,
                 _NBOX , tria_pred()) ;
    
   /*--------------------- flip node order - orientations */    
        for (auto _tpos  = 
             this->_mesh._set4.head() ;
                  _tpos != 
             this->_mesh._set4.tend() ;
                ++_tpos  )
        {
            if (_tpos->mark() >= +0)
            {
            real_type _tvol = 
                geometry::tetra_vol_3d (
           &this->_mesh._set1[
                _tpos->node(0)].pval(0),
           &this->_mesh._set1[
                _tpos->node(1)].pval(0),
           &this->_mesh._set1[
                _tpos->node(2)].pval(0),
           &this->_mesh._set1[
                _tpos->node(3)].pval(0)
                ) ;
        
            if (_tvol < (real_type) +0.)
            std::swap(_tpos->node (0),
                      _tpos->node (1)) ;
            }
        }
    }
    
    /*
    --------------------------------------------------------
     * FIND-TRIA: scan for enclosing tria.
    --------------------------------------------------------
     */  
     
    class find_tria
        {
        public  :
        real_type            *_ppos ;
        
        mesh_type            *_mesh ;
        
        real_type             _vtol ;
        bool_type             _find ;
        iptr_type             _tpos ;
        
        public  :
    
    /*------------------------ make a tree-tria predicate */
        __inline_call find_tria (
            real_type*_psrc = nullptr ,
            mesh_type*_msrc = nullptr ,
            real_type _vsrc = 
                (real_type) + 0.0
            ) : _ppos(_psrc) ,
                _mesh(_msrc) ,
                _vtol(_vsrc) ,
                _find(false) ,
                _tpos(   -1)   {}

    /*------------------------ call pred. on tree matches */
        __inline_call 
            void_type operator () (
            typename  
        tree_type::item_data *_iptr
            )
        {
            if (this->_find)  return;
        
            for ( ; _iptr != nullptr; 
                    _iptr = _iptr->_next)
            {
                iptr_type  _tpos = 
                    _iptr->_data.ipos() ;
                
                if (find_pred ( _ppos ,
                        _mesh , _tpos ,
                        _vtol ) )
                {
                    this->_find =  true ;
                    this->_tpos = _tpos ;
                    
                    break ;
                }
            }
        }
        
        } ;
   
    typedef geom_tree::aabb_pred_node_k <
            real_type , 
            iptr_type , 3 > tree_pred ;
   
    /*
    --------------------------------------------------------
     * TRIA-PRED: TRUE if PPOS is in TPOS.
    --------------------------------------------------------
     */
     
    __static_call
    __normal_call bool_type find_pred (
        real_type*_ppos ,
        mesh_type*_mesh ,
        iptr_type _tpos ,
        real_type _vtol
        )
    {
        for (auto _fpos = +4; _fpos-- != +0; )
        {
            iptr_type _fnod[ 4];
            mesh_type::tri4_type::
                face_node(_fnod, _fpos, 3, 2);
            _fnod[0] = _mesh->
            _set4[_tpos].node(_fnod[0]);
            _fnod[1] = _mesh->
            _set4[_tpos].node(_fnod[1]);
            _fnod[2] = _mesh->
            _set4[_tpos].node(_fnod[2]);
            
            real_type _tvol = 
                geometry::tetra_vol_3d (
               &_mesh->
               _set1[ _fnod[0]].pval(0),
               &_mesh->
               _set1[ _fnod[1]].pval(0),
               &_mesh->
               _set1[ _fnod[2]].pval(0),
               _ppos) ;
            
            if (_tvol < -_vtol) return false ;
        }
        
        return  true ;
    }
   
    /*
    --------------------------------------------------------
     * HINT: check for valid index.
    --------------------------------------------------------
     */
   
    __inline_call bool_type hint_okay (
        hint_type _hint
        )
    {   
    /*------------------------- test whether hint is valid */
        return _hint >= (iptr_type)0
         && _hint < (iptr_type) 
            this->_mesh._set4.count()
         && this->_mesh.
           _set4 [_hint].mark() >= 0 ;
    }
   
    /*
    --------------------------------------------------------
     * EVAL: eval. size-fun. value.
    --------------------------------------------------------
     */
    
    __normal_call real_type eval (
        real_type *_ppos ,
        hint_type &_hint
        )
    /*------------------------ find tria + linear interp. */
    {
        real_type _hval = 
    +std::numeric_limits<real_type>::infinity() ;
    
        if (hint_okay(_hint))
        {
    /*------------------------ test whether hint is valid */
            if(!find_pred( _ppos, 
                   &_mesh, _hint,
                this->_vtol ) )
            {
            _hint =  this->null_hint();
            }
        }
        else
        {
    /*------------------------ hint is definitely invalid */
            _hint =  this->null_hint();
        }

        if (_hint == this->null_hint())
        {
    /*------------------------ scan to find bounding tria */
            tree_pred _pred (_ppos) ;
            find_tria _func (_ppos,
                            &_mesh,
                       this->_vtol) ;
        
            this->
           _tree.find(_pred, _func) ;
        
           _hint = _func._find ? 
                   _func._tpos : 
            hfun_type::null_hint () ;
        }
    
        if (_hint != this->null_hint())
        {
    /*------------------------ linear interp. within tria */
        real_type _hsum = (real_type)+.0 ;
        real_type _vsum = (real_type)+.0 ;

        for(auto _fpos = 4; _fpos-- != 0; )
        {
            iptr_type  _fnod [4] ;
            tri4_type::
            face_node(_fnod, _fpos, 3, 2) ;
            _fnod[0] = this->_mesh.
            _set4[_hint].node(_fnod[0]);
            _fnod[1] = this->_mesh.
            _set4[_hint].node(_fnod[1]);
            _fnod[2] = this->_mesh.
            _set4[_hint].node(_fnod[2]);
            _fnod[3] = this->_mesh.
            _set4[_hint].node(_fnod[3]);
            
            real_type _tvol = 
                geometry::tetra_vol_3d (
               &this->_mesh.
               _set1[ _fnod[0]].pval(0),
               &this->_mesh.
               _set1[ _fnod[1]].pval(0),
               &this->_mesh.
               _set1[ _fnod[2]].pval(0),
               _ppos) ;

            _hsum += _tvol * this-> 
            _mesh._set1[_fnod[3]].hval() ; 
            
            _vsum += _tvol ;
        }


        _hval = _hsum / _vsum ;
        
        }
        else
        {
    /*------------------------- inv-dist. weighted extrap */
        typename tree_type::item_data 
             *_nptr = nullptr ;    
    
        this->_tree.
             near(_ppos, _nptr) ;   // find nearest tria.
    
        if (_nptr  != nullptr )
        {      
    /*------------------------- extrap. using tria. nodes */
        iptr_type  _tpos = 
            _nptr->_data.ipos() ;
            
        iptr_type  _inod = 
          _mesh._set4[_tpos].node(0) ;
        iptr_type  _jnod = 
          _mesh._set4[_tpos].node(1) ;
        iptr_type  _knod = 
          _mesh._set4[_tpos].node(2) ;
        iptr_type  _lnod = 
          _mesh._set4[_tpos].node(3) ;
            
        real_type _isqr  = 
            geometry::lensqr_3d(_ppos, 
       &_mesh._set1[_inod].pval( +0));
        real_type _jsqr  = 
            geometry::lensqr_3d(_ppos, 
       &_mesh._set1[_jnod].pval( +0));
        real_type _ksqr  = 
            geometry::lensqr_3d(_ppos, 
       &_mesh._set1[_knod].pval( +0));
        real_type _lsqr  = 
            geometry::lensqr_3d(_ppos, 
       &_mesh._set1[_lnod].pval( +0));
    
    /*------------------------- extrap. as inv-dist. fun. */   
        real_type _wsum = (real_type)0. ;
        real_type _hsum = (real_type)0. ;
       
        real_type _iwsc = 
            (real_type)+1. / _isqr ;
        real_type _jwsc = 
            (real_type)+1. / _jsqr ;
        real_type _kwsc = 
            (real_type)+1. / _ksqr ;
        real_type _lwsc = 
            (real_type)+1. / _lsqr ;
       
        _wsum +=  _iwsc + _jwsc + 
                  _kwsc + _lwsc ;
 
        _hsum +=  _iwsc * 
            _mesh._set1[_inod].hval();
        _hsum +=  _jwsc * 
            _mesh._set1[_jnod].hval();
        _hsum +=  _kwsc * 
            _mesh._set1[_knod].hval();
        _hsum +=  _lwsc * 
            _mesh._set1[_lnod].hval();
      
      
        _hval = _hsum / _wsum ;
        
        }
        
        }  
          
    /*------------------------- size-fun interp. to ppos. */
        return  _hval ;  
    }
    
    } ;


    }

#   endif   //__HFUN_MESH_EUCLIDEAN_3__

    
    
