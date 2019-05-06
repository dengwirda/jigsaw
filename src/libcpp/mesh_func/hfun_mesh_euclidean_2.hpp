
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
     * Last updated: 30 April, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_MESH_EUCLIDEAN_2__
#   define __HFUN_MESH_EUCLIDEAN_2__

    namespace mesh {

    template <
    typename R ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_mesh_euclidean_2d
        : public hfun_base_kd <I, R> 
    {
    public  :
    
    /*------------------------- euclidean size-fun in R^2 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef hfun_mesh_euclidean_2d  <
            real_type ,
            iptr_type >             hfun_type ;

    typedef typename  hfun_base_kd  <
            iptr_type , 
            real_type >::hint_type  hint_type ;
            

    class node_type: public tria_complex_node_2<I,R>
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
       
    typedef mesh::tria_complex_2<
                node_type,
                edge_type,
                tri3_type,
                allocator       >   mesh_type ;

    typedef geom_tree::aabb_node_base_k
                                    tree_node ;

    typedef geom_tree::aabb_item_rect_k <
                real_type,
                iptr_type,
                + 2             >   tree_item ;
                
    typedef geom_tree::aabb_tree<
                tree_item,
                + 2 ,
                tree_node,
                allocator       >   tree_type ; 
                
    typedef geom_tree::aabb_pred_node_2 <
                real_type, 
                iptr_type       >   tree_pred ;
              
    public  :              
    
    containers::
        fixed_array<real_type,2>   _bmin ;
    containers::
        fixed_array<real_type,2>   _bmax ;
   
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
                tri3_type const& _tdat
                ) const
            {   return _tdat.mark() >= 0 ;
            }
            } ;
    
    /*----------------------------- init. aabb at -+ inf. */
        for(auto _idim = 2; _idim-- != 0 ; )
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
        
            this->_bmax[0] = std::max (
            this->_bmax[0] , 
            _iter->pval(0) ) ;
            this->_bmax[1] = std::max (
            this->_bmax[1] , 
            _iter->pval(1) ) ;
            }
        }
  
        real_type static const _RTOL = 
            std::pow (
        std::numeric_limits<real_type>
            ::epsilon(),(real_type).8) ;
            
        iptr_type static
        constexpr _NBOX=(iptr_type)+4  ;

        real_type _BTOL[2] ;
        _BTOL[0] = this->_bmax [ 0 ] - 
                   this->_bmin [ 0 ] ;
        _BTOL[1] = this->_bmax [ 1 ] - 
                   this->_bmin [ 1 ] ;
  
        _BTOL[0]*=_RTOL ;
        _BTOL[1]*=_RTOL ;

    /*-------------------- make aabb-tree and init. bbox. */
        aabb_mesh( this->_mesh._set1 , 
                   this->_mesh._set3 , 
                   this->_tree,_BTOL ,
                  _NBOX , tria_pred()) ;
    }
    
    /*
    --------------------------------------------------------
     * FIND-TRIA: scan for nearest tria.
    --------------------------------------------------------
     */  
    
    class find_tria
        {
        public  :
        real_type              *_ppos ;
        
        mesh_type              *_mesh ;
        
        bool_type               _find ;
        iptr_type               _tpos ;
        
        public  :
    
    /*------------------------ make a tree-tria predicate */
        __inline_call find_tria (
            real_type*_psrc = nullptr ,
            mesh_type*_msrc = nullptr
            ) : _ppos(_psrc) ,
                _mesh(_msrc) ,
                _find(false) ,
                _tpos(   -1)   {}

    /*------------------------ call pred. on tree matches */
        __inline_call 
            void_type operator () (
                typename  
            tree_type::item_data  *_iptr
            )
        {
            if (this->_find) return ;
        
            for ( ; _iptr != nullptr; 
                    _iptr = _iptr->_next)
            {
                real_type  _qtmp[+2];
                iptr_type  _tpos = 
                    _iptr->_data.ipos() ;
                
                if (near_pred ( _ppos ,
                        _qtmp ,*_mesh , 
                        _tpos ) )
                {
    /*------------------------ is fully inside: finished! */
                this->_find =  true ;
                this->_tpos = _tpos ;
                
                break ;

                }
            }
        }
        
        } ;
 
    class near_tria
        {
        public  :
        real_type              *_ppos ;
        real_type              *_qpos ;

        real_type               _dsqr ;
        
        mesh_type              *_mesh ;
        
        bool_type               _find ;
        iptr_type               _tpos ;
        
        public  :
    
    /*------------------------ make a tree-tria predicate */
        __inline_call near_tria (
            real_type*_psrc = nullptr ,
            real_type*_qsrc = nullptr ,
            mesh_type*_msrc = nullptr
            ) : _ppos(_psrc) ,
                _qpos(_qsrc) ,
           _dsqr(+std::numeric_limits
           <real_type>::infinity()) ,
                _mesh(_msrc) ,
                _find(false) ,
                _tpos(   -1)   {}

    /*------------------------ call pred. on tree matches */
        __inline_call 
            real_type operator () (
                typename  
            tree_type::item_data  *_iptr
            )
        {
            if (this->_find) return +0. ;
        
            for ( ; _iptr != nullptr; 
                    _iptr = _iptr->_next)
            {
                real_type  _qtmp[+2];
                iptr_type  _tpos = 
                    _iptr->_data.ipos() ;
                
                if (near_pred ( _ppos ,
                        _qtmp ,*_mesh , 
                        _tpos ) )
                {
    /*------------------------ is fully inside: finished! */
                this->_dsqr = 
                    (real_type) +0. ;

                this->_find =  true ;
                this->_tpos = _tpos ;
                
                break ;

                }
                else
                {
    /*------------------------ projected match: keep best */
                real_type _dtmp = 
            geometry::lensqr_2d(_ppos, _qtmp);

                if (_dtmp<_dsqr )
                {
                _qpos[0] = _qtmp[0] ;
                _qpos[1] = _qtmp[1] ;

                this->_dsqr = _dtmp ;
                this->_tpos = _tpos ;
                }

                }
            }

            return ( this->_dsqr )  ;
        }
        
        } ;
   
    /*
    --------------------------------------------------------
     * NEAR-PRED: TRUE if PPOS is in TPOS.
    --------------------------------------------------------
     */
     
    __static_call
    __normal_call bool_type near_pred (
        real_type*_ppos ,
        real_type*_qpos ,
        mesh_type&_mesh ,
        iptr_type _tpos
        )
    {
        geometry::hits_type _hits ;
        if (geometry::proj_tria_2d(_ppos, 
           &_mesh._set1[
            _mesh._set3[
            _tpos].node(0)].pval(0) ,
           &_mesh._set1[
            _mesh._set3[
            _tpos].node(1)].pval(0) ,
           &_mesh._set1[
            _mesh._set3[
            _tpos].node(2)].pval(0) ,
            _qpos,_hits) )
        {
            return (_hits ==
                 geometry::tria_hits) ;
        }
        else
        {
            return ( false ) ;
        } 
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
    /*------------------------ test whether hint is valid */
        return _hint >= (iptr_type)0
         && _hint < (iptr_type) 
            this->_mesh._set3.count()
         && this->_mesh.
           _set3 [_hint].mark() >= 0 ;
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
        real_type _QPOS[ 2] ;
        _QPOS[0] = _ppos[0] ;
        _QPOS[1] = _ppos[1] ;        

        real_type _hval = 
    +std::numeric_limits<real_type>::infinity() ;
    
        if (hint_okay(_hint))
        {
    /*------------------------ test whether hint is valid */
            if(!near_pred( _ppos, 
                    _QPOS, _mesh, 
                    _hint)  )
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
                            &_mesh) ;
        
            this->
           _tree.find(_pred, _func) ;
        
           _hint  = _func._find ? 
                    _func._tpos : 
            hfun_type::null_hint () ;
        }

        if (_hint == this->null_hint())
        {
    /*------------------------ outside: find nearest tria */
            near_tria _func (_ppos,
                      _QPOS,&_mesh) ;
        
            this->
           _tree.near(_ppos, _func) ;
        
           _hint =  _func._tpos ;
        }
    
        if (_hint != this->null_hint())
        {
    /*------------------------ linear interp. within tria */
        real_type _hsum = (real_type)+.0 ;
        real_type _vsum = (real_type)+.0 ;

        for(auto _fpos = 3; _fpos-- != 0; )
        {
            iptr_type  _fnod [3] ;
            tri3_type::
            face_node(_fnod, _fpos, 2, 1) ;

            _fnod[0] = this->_mesh.
            _set3[_hint].node(_fnod[0]);
            _fnod[1] = this->_mesh.
            _set3[_hint].node(_fnod[1]);
            _fnod[2] = this->_mesh.
            _set3[_hint].node(_fnod[2]);
            
            real_type _tvol = 
                geometry::tria_area_2d (
               &this->_mesh.
               _set1[_fnod[0]].pval(0) ,
               &this->_mesh.
               _set1[_fnod[1]].pval(0) ,
               _QPOS) ;

            _hsum += _tvol * this-> 
            _mesh._set1[_fnod[2]].hval() ;
            
            _vsum += _tvol ;
        }


        _hval = _hsum / _vsum ;

        }

    /*------------------------- size-fun interp. to ppos. */
        return  _hval ;
    }
    
    } ;


    }

#   endif   //__HFUN_MESH_EUCLIDEAN_2__

    
    
