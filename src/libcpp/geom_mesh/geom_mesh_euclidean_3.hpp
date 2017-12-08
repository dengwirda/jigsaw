
    /*
    --------------------------------------------------------
     * GEOM-MESH-EUCLIDEAN-3: euclidean geom. in R^3.
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
     * Last updated: 16 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __GEOM_MESH_EUCLIDEAN_3__
#   define __GEOM_MESH_EUCLIDEAN_3__

    namespace mesh {
 
    template <
    typename R ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class geom_mesh_euclidean_3d
    {
    public  :
    /*---------- discrete geometry in R^3 via nodes/faces */
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef geom_mesh_euclidean_3d  <
                real_type ,
                iptr_type ,
                allocator    >      geom_type ;

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
        
    class tri3_type: public tria_complex_tria_3<I>
        {
    /*------------------------------------ loc. face type */
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
	            node_type::
	                _dims       >   tree_item ;
	typedef geom_tree::aabb_tree<
		        tree_item,
		        node_type::
		            _dims,
                tree_node,
                allocator       >   tree_type ;

    iptr_type static constexpr _nbox = + 4;

    public  :
    
    containers::
        fixed_array<real_type,3>    _bmin ;
    containers::
        fixed_array<real_type,3>    _bmax ;

    mesh_type                       _tria ;

    tree_type                       _ebox ;
    tree_type                       _tbox ;

    public  :
    
    /*
    --------------------------------------------------------
     * construct geometry from alloc. etc.
    --------------------------------------------------------
     */
     
    __normal_call geom_mesh_euclidean_3d (
        allocator const& _asrc = allocator ()
        ) : _tria(_asrc) ,
            _ebox(_asrc) ,
            _tbox(_asrc) 
        {
        }
            
    /*
    --------------------------------------------------------
     * NODE-FEAT: calc. node feature type.
    --------------------------------------------------------
     */
     
    template <
        typename  list_type
             >
    __normal_call void_type feat_list (
        list_type &_lsrc ,
        list_type &_ldst
        )
    {
        for (auto _iter  = _lsrc.head() ;
                  _iter != _lsrc.tend() ;
                ++_iter  )
        {
        if  (this->_tria .
            _set2[*_iter].
                feat() != null_feat)
        {
            _ldst.push_tail (*_iter)  ;
        }
        }
    }
     
    template <
        typename  list_type ,
        typename  geom_opts
             >
    __normal_call char_type node_feat (
        iptr_type *_node ,
        list_type &_aset ,
        geom_opts &_opts
        )
    {
        char_type _feat = 
           (_aset.count() == +0) ||
           (_aset.count() == +2)
                ? null_feat : soft_feat ;
    
        real_type _DtoR = 
       (real_type) +3.1415926536 / 180.0;
    
        real_type _phi1 = 
       (real_type)+180. - _opts.phi1();
        real_type _eta1 = 
       (real_type)+  0. + _opts.eta1();
    
        real_type _hard = 
            std::cos( _phi1 * _DtoR) ;
        real_type _soft = 
            std::cos( _eta1 * _DtoR) ;
  
        __unreferenced(_node) ;
  
        for (auto _ipos  = _aset.head() ;
                  _ipos != _aset.tend() ;
                ++_ipos  )
        for (auto _jpos  = _ipos+1 ;
                  _jpos != _aset.tend() ;
                ++_jpos  )
        {
            iptr_type _iedg = * _ipos ;
            iptr_type _jedg = * _jpos ;
          
            iptr_type _inod[2] ;
            _inod[0] = this->_tria.
                _set2[_iedg]. node(0) ;
            _inod[1] = this->_tria.
                _set2[_iedg]. node(1) ;
            
            iptr_type _jnod[2] ;
            _jnod[0] = this->_tria.
                _set2[_jedg]. node(0) ;
            _jnod[1] = this->_tria.
                _set2[_jedg]. node(1) ;

            real_type _ivec[3] ;
            geometry::vector_3d(
               &this->_tria.
               _set1[ _inod[0]].pval(0) ,
               &this->_tria.
               _set1[ _inod[1]].pval(0) ,
               _ivec) ;
            
            real_type _jvec[3] ;
            geometry::vector_3d(
               &this->_tria.
               _set1[ _jnod[0]].pval(0) ,
               &this->_tria.
               _set1[ _jnod[1]].pval(0) ,
               _jvec) ;
                  
            real_type _acos = 
                geometry::cosine_3d(
                    _ivec , _jvec) ;
                
            if (_inod[0] == _jnod[1] ||
                _inod[1] == _jnod[0] )
                _acos *= (real_type)+1.;
            else
                _acos *= (real_type)-1.;
                
            if (_acos <= _hard)
            {
                _feat  = 
            std::max(_feat, hard_feat) ;
            }
            else
            if (_acos <= _soft)
            {          
                _feat  = 
            std::max(_feat, soft_feat) ;
            } 
        }    
        
        return  _feat  ;   
    }
    
    template <
        typename  list_type ,
        typename  geom_opts
             >
    __normal_call char_type node_feat (
        iptr_type *_node ,
        list_type &_eadj ,
        list_type &_fadj ,
        geom_opts &_opts
        )
    {
        char_type _feat =   null_feat ;
    
        real_type _DtoR = 
       (real_type) +3.1415926536 / 180.0;
    
        real_type _phi1 = 
       (real_type)+180. - _opts.phi1();
        real_type _eta1 = 
       (real_type)+  0. + _opts.eta1();
    
        real_type _hard = 
            std::cos( _phi1 * _DtoR) ;
        real_type _soft = 
            std::cos( _eta1 * _DtoR) ;
  
        __unreferenced(_node) ;
  
        for (auto _epos  = _eadj.head() ;
                  _epos != _eadj.tend() ;
                ++_epos  )
        for (auto _fpos  = _fadj.head() ;
                  _fpos != _fadj.tend() ;
                ++_fpos  )
        {
            iptr_type _iedg = * _epos ;
            iptr_type _ifac = * _fpos ;
          
            iptr_type _enod[2] ;
            _enod[0] = this->_tria.
                _set2[_iedg]. node(0) ;
            _enod[1] = this->_tria.
                _set2[_iedg]. node(1) ;
            
            iptr_type _fnod[3] ;
            _fnod[0] = this->_tria.
                _fset[_ifac]. node(0) ;
            _fnod[1] = this->_tria.
                _fset[_ifac]. node(1) ;
            _fnod[2] = this->_tria.
                _fset[_ifac]. node(2) ;

            iptr_type _same = +0 ;
            if (_fnod[0] == _enod[0])
                _same += +1 ;
            if (_fnod[0] == _enod[1])
                _same += +1 ;
            if (_fnod[1] == _enod[0])
                _same += +1 ;
            if (_fnod[1] == _enod[1])
                _same += +1 ;
            if (_fnod[2] == _enod[0])
                _same += +1 ;
            if (_fnod[2] == _enod[1])
                _same += +1 ;
            if (_same >= +2) continue ;

            real_type _evec[3] ;
            geometry::vector_3d(
               &this->_tria.
               _set1[ _enod[0]].pval(0) ,
               &this->_tria.
               _set1[ _enod[1]].pval(0) ,
               _evec) ;
            
            real_type _fvec[3] ;
            
            
            
            //!! actually, need to project edge onto tria
            
            
                  
            real_type _acos = 
                geometry::cosine_3d(
                    _evec , _fvec) ;
                    
            _acos =std::abs(_acos) ;
                
            if (_acos <= _hard)
            {
                _feat  = 
            std::max(_feat, hard_feat) ;
            }
            else
            if (_acos <= _soft)
            {          
                _feat  = 
            std::max(_feat, soft_feat) ;
            } 
        }    
        
        return  _feat  ;   
    }
    
    /*
    ------------------------------------------------------------
     * EDGE-FEAT: calc. edge feature type.
    ------------------------------------------------------------
     */
     
    template <
        typename  list_type ,
        typename  geom_opts
             >
    __normal_call char_type edge_feat (
        iptr_type *_enod ,
        list_type &_aset ,
        geom_opts &_opts
        )
    {
        char_type _feat = 
           (_aset.count() == +2)
                ? null_feat : soft_feat ;
        
        real_type _DtoR = 
       (real_type) +3.1415926536 / 180.0;
    
        real_type _phi2 = 
       (real_type)+180. - _opts.phi2();
        real_type _eta2 = 
       (real_type)+  0. + _opts.eta2();
    
        real_type _hard = 
            std::cos( _phi2 * _DtoR) ;
        real_type _soft = 
            std::cos( _eta2 * _DtoR) ;
    
        for (auto _ipos  = _aset.head() ;
                  _ipos != _aset.tend() ;
                ++_ipos  )
        for (auto _jpos  = _ipos+1 ;
                  _jpos != _aset.tend() ;
                ++_jpos  )
        {
            iptr_type _itri = * _ipos ;
            iptr_type _jtri = * _jpos ;
            
            iptr_type _inod[3];
            iptr_type _jnod[3];
            
            iptr_type _iloc ;
            for (_iloc = 3; _iloc-- != 0; )
            {
                tri3_type::face_node (
                    _inod, _iloc, 2, 1) ;
                    
                _inod[0] = this->_tria.
                    _set3[_itri]. 
                        node(_inod[0]);
                _inod[1] = this->_tria.
                    _set3[_itri]. 
                        node(_inod[1]);
                _inod[2] = this->_tria.
                    _set3[_itri]. 
                        node(_inod[2]);
                        
                iptr_type _same =  +0 ;
                if (_inod[0]==_enod[0])
                    _same += +1 ;
                if (_inod[0]==_enod[1])
                    _same += +1 ;
                if (_inod[1]==_enod[0])
                    _same += +1 ;
                if (_inod[1]==_enod[1])
                    _same += +1 ;
                    
                if (_same == +2 ) break ;
            }
            
            iptr_type _jloc ;
            for (_jloc = 3; _jloc-- != 0; )
            {
                tri3_type::face_node (
                    _jnod, _jloc, 2, 1) ;
                    
                _jnod[0] = this->_tria.
                    _set3[_jtri]. 
                        node(_jnod[0]);
                _jnod[1] = this->_tria.
                    _set3[_jtri]. 
                        node(_jnod[1]);
                _jnod[2] = this->_tria.
                    _set3[_jtri]. 
                        node(_jnod[2]);
                        
                iptr_type _same =  +0 ;
                if (_jnod[0]==_enod[0])
                    _same += +1 ;
                if (_jnod[0]==_enod[1])
                    _same += +1 ;
                if (_jnod[1]==_enod[0])
                    _same += +1 ;
                if (_jnod[1]==_enod[1])
                    _same += +1 ;
                    
                if (_same == +2 ) break ;
            }
            
            real_type _ivec[3];
            geometry::tria_norm_3d (
            &this->_tria.
            _set1 [_inod[0]].pval(0) ,
            &this->_tria.
            _set1 [_inod[1]].pval(0) ,
            &this->_tria.
            _set1 [_inod[2]].pval(0) ,
            _ivec) ;
                   
            real_type _jvec[3];
            geometry::tria_norm_3d (
            &this->_tria.
            _set1 [_jnod[0]].pval(0) ,
            &this->_tria.
            _set1 [_jnod[1]].pval(0) ,
            &this->_tria.
            _set1 [_jnod[2]].pval(0) ,
            _jvec) ;
            
            real_type _acos = 
                geometry::cosine_3d(
                    _ivec , _jvec) ;
 
            if (_inod[0] == _jnod[1] &&
                _inod[1] == _jnod[0] )
                _acos *= (real_type)+1.;
            else
                _acos *= (real_type)-1.;
           
            if (_acos <= _hard)
            {
                _feat  = 
            std::max(_feat, hard_feat) ;
            }
            else
            if (_acos <= _soft)
            {          
                _feat  = 
            std::max(_feat, soft_feat) ;
            }
        }
        
        return  _feat  ;
    }

    /*
    --------------------------------------------------------
     * SCAN-NODE: scan nodes and find features.
    --------------------------------------------------------
     */
     
    template <
        typename  geom_opts
             >
    __normal_call void_type scan_node (
        geom_opts &_opts
        )
    {
        containers::
            array<iptr_type>  _eadj ;
        containers::
            array<iptr_type>  _fadj ;

        containers::
            array<iptr_type>  _set2 ;
        
        for (auto _npos  = 
             this->_tria._set1.head() ;
                  _npos != 
             this->_tria._set1.tend() ;
                ++_npos  )
        {
            if (_npos->mark() >= +0)
            {
            _eadj.set_count (0);
            _fadj.set_count (0);
            
            _set2.set_count (0);
            
             this->_tria.node_edge (
                &_npos->node(0), _eadj) ;
                
             this->_tria.node_tri3 (
                &_npos->node(0), _fadj) ;

             feat_list(_eadj, _set2);

            _npos->topo () =  
               (char_type)_set2.count() ;
            
            _npos->feat () = node_feat(
                &_npos->node(0), _set2  , 
                    _opts) ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * SCAN-EDGE: scan edges and find features.
    --------------------------------------------------------
     */
     
    template <
        typename  geom_opts
             >
    __normal_call void_type scan_edge (
        geom_opts &_opts
        )
    {
        containers::
            array<iptr_type>  _fadj ;
        
        for (auto _epos  = 
             this->_tria._set2.head() ;
                  _epos != 
             this->_tria._set2.tend() ;
                ++_epos  )
        {
            if (_epos->mark() >= +0)
            {
            _fadj.set_count (0);
            
             this->_tria.edge_tri3 (
                &_epos->node(0), _fadj) ;
                
            _epos->topo () =  
               (char_type)_fadj.count() ;

            _epos->feat () = edge_feat(
                &_epos->node(0), _fadj  , 
                    _opts) ;        
            }
        }
    }

    /*
    --------------------------------------------------------
     * FIND-FEAT: scan geometry and find features.
    --------------------------------------------------------
     */
     
    template <
        typename  geom_opts
             >
    __normal_call void_type find_feat (
        geom_opts &_opts
        )
    {
    /*---------------------------------- init. geom feat. */
        for (auto _iter  = 
             this->_tria._set1.head() ;
                  _iter != 
             this->_tria._set1.tend() ;
                ++_iter  )
        {
            _iter->fdim () = +0  ;
            _iter->feat () = null_feat ;
            _iter->topo () = +2  ;
        }
        for (auto _iter  = 
             this->_tria._set2.head() ;
                  _iter != 
             this->_tria._set2.tend() ;
                ++_iter  )
        {
            _iter->feat () = null_feat ;
            _iter->topo () = +2  ;
        }
        for (auto _iter  = 
             this->_tria._set3.head() ;
                  _iter != 
             this->_tria._set3.tend() ;
                ++_iter  )
        {
            _iter->feat () = null_feat ;
            _iter->topo () = +2  ;
        }
  
        if (_opts.feat ())
        {
    /*---------------------------------- find sharp feat. */
            scan_edge(_opts) ;
            scan_node(_opts) ;
        }
        
        for (auto _iter  = 
             this->_tria._set3.head() ;
                  _iter != 
             this->_tria._set3.tend() ;
                ++_iter  )
        {
    /*----------------------------- assign nodes to trias */
            this->_tria.
           _set1[_iter->node(0)].fdim() = 2;
            this->_tria.
           _set1[_iter->node(1)].fdim() = 2;
            this->_tria.
           _set1[_iter->node(2)].fdim() = 2;
        }
        
        for (auto _iter  = 
             this->_tria._set2.head() ;
                  _iter != 
             this->_tria._set2.tend() ;
                ++_iter  )
        {
    /*----------------------------- assign nodes to edges */
            if (_iter->feat() !=  null_feat)
            {
            this->_tria.
           _set1[_iter->node(0)].fdim() = 1;
            this->_tria.
           _set1[_iter->node(1)].fdim() = 1;
            }
        }
        
        for (auto _iter  = 
             this->_tria._set1.head() ;
                  _iter != 
             this->_tria._set1.tend() ;
                ++_iter  )
        {
            if (_iter->feat() !=  null_feat)
            {
    /*----------------------------- assign nodes to feat. */
                _iter->fdim()  = 0 ;
            }
        }
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
        class edge_pred
            {
        /*---------------- TRUE if tree should index edge */
            public  :
            __inline_call 
                bool_type operator () (
                edge_type const& _edat
                ) const
            {   return _edat.feat() != 0 
                    && _edat.mark() >= 0 ;
            }
            } ;
            
        class tri3_pred
            {
        /*---------------- TRUE if tree should index tria */
            public  :
            __inline_call 
                bool_type operator () (
                tri3_type const& _tdat
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
             this->_tria._set1.head() ;
                   _iter !=
             this->_tria._set1.tend() ;
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
  
        real_type _BTOL[3] ;
        _BTOL[0] = this->_bmax[ 0] - 
                   this->_bmin[ 0] ;
        _BTOL[1] = this->_bmax[ 1] - 
                   this->_bmin[ 1] ;
        _BTOL[2] = this->_bmax[ 2] - 
                   this->_bmin[ 2] ;
  
        _BTOL[0]*= _RTOL ;
        _BTOL[1]*= _RTOL ;
        _BTOL[2]*= _RTOL ;
  
    /*-------------------- find sharp "features" in geom. */    
        find_feat (_opts);
        
    /*-------------------- make aabb-tree and init. bbox. */
        aabb_mesh(this->_tria._set1, 
                  this->_tria._set2, 
                  this->_ebox,_BTOL,
                  this->_nbox,edge_pred() 
                 ) ;
                  
        aabb_mesh(this->_tria._set1, 
                  this->_tria._set3, 
                  this->_tbox,_BTOL,
                  this->_nbox,tri3_pred() 
                 ) ;
    }

    /*
    --------------------------------------------------------
     * SEED-MESH: setup initial node set.
    --------------------------------------------------------
     */
     
    template <
        typename  mesh_type ,
        typename  size_type ,
        typename  geom_opts
             >
    __normal_call void_type seed_mesh (
        mesh_type &_mesh ,
        size_type &_size ,
        geom_opts &_opts
        )
    {
        __unreferenced(_size) ;
    
    /*------------------------- push set of feature nodes */
        for (auto _iter  = 
                   this->_tria._set1.head() ;
                  _iter != 
                   this->_tria._set1.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0 )
            {
            if (_iter->feat() != null_feat ||
                _iter->itag() <= -1 )
            {
            iptr_type _node = -1 ;
            if (_mesh._tria.push_node (
               &_iter->pval(0), _node))
            {
                _mesh._tria.node
                    (_node)->fdim() 
                        = _iter->fdim() ;
                        
                _mesh._tria.node
                    (_node)->feat() 
                        = _iter->feat() ;
                        
                _mesh._tria.node
                    (_node)->topo() 
                        = _iter->topo() ;          
            }
            }
            }
        }
    
    /*------------------------- well-distributed sampling */
        while (_mesh._tria.
               _nset.count() < (std::size_t)_opts.seed())
        {
            typename geom_type::
                     mesh_type::
                     node_list::_write_it _best;
            
            real_type _dmax  = (real_type) +.0 ;
            char_type _fdim  ;
            
            for (_fdim = 1; _fdim != 4; ++_fdim)
            {
            for (auto _ipos  = 
                 this->_tria._set1.head() ;
                      _ipos != 
                 this->_tria._set1.tend() ;
                    ++_ipos  )
            {
                if (_ipos->fdim() == _fdim)
                {
                real_type _dmin  = 
            +std::numeric_limits<real_type>::infinity() ;

                for (auto _jpos  = 
                    _mesh._tria._nset.head() ;
                          _jpos != 
                    _mesh._tria._nset.tend() ;
                        ++_jpos  )
                {
                    real_type _dist = 
                        geometry::lensqr_3d(
                        &_ipos->pval(+0), 
                        &_jpos->pval(+0)) ;

                    _dmin = std::min(_dmin, _dist);
                }

                if (_dmax < _dmin)
                {
                    _dmax = _dmin;
                    _best = _ipos;
                }
                }
            }
            if (_dmax > (real_type)0.) break ;
            }
 
            if (_dmax > (real_type)0.) 
            {
    /*------------------------- add current furthest node */
                iptr_type _node = -1;
                if (_mesh._tria.push_node(
                   &_best->pval(0), _node) )
                {
                    _mesh._tria.node
                        (_node)->fdim() 
                            = _best->fdim() ;
                            
                    _mesh._tria.node
                        (_node)->feat() 
                            = _best->feat() ;
                            
                    _mesh._tria.node
                        (_node)->topo() 
                            = _best->topo() ;           
                }
            }
            else break  ;
        }

    }

    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */

    template <
    typename   hits_func
             >
    class line_tria_pred
        {
    /*------------------ line-tria intersection predicate */
        public  :
        real_type      _ipos[3] ;
        real_type      _jpos[3] ;
        
        geom_type&     _geom ;
        
        hits_func&     _hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        bool_type      _exct ;
        
        geometry::
        hits_type      _hits ;
        
        public  :
        
        line_tria_pred operator = (
            line_tria_pred & 
            )                   =     delete ;
        line_tria_pred operator = (
            line_tria_pred&& 
            )                   =     delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call line_tria_pred (
            real_type *_isrc ,
            real_type *_jsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits = geometry::null_hits ;
        
            this->_exct = false;
        
            this->_hnum =     0;
            this->_find = false;
        
            for (auto _idim = 3; _idim-- != 0; )
            {
                this->
               _ipos[_idim] = _isrc[_idim];
                this->
               _jpos[_idim] = _jsrc[_idim];
            }
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
                geometry::hits_type 
                    _hits = geometry::null_hits;
            
            /*--------------- line-tria intersection test */
                iptr_type  _tpos = 
                    _iptr->_data.ipos() ;
    
                iptr_type  _tnod[3];
                _tnod[0] =_geom.
                    _tria._set3[_tpos].node(0) ;
                _tnod[1] =_geom.
                    _tria._set3[_tpos].node(1) ;
                _tnod[2] =_geom.
                    _tria._set3[_tpos].node(2) ;
                
                real_type  _xpos[3];
                _hits = geometry::line_tria_3d (
                   & this->_ipos[0], 
                   & this->_jpos[0],
                   &_geom ._tria.
                    _set1 [_tnod[0]].pval(0),
                   &_geom ._tria.
                    _set1 [_tnod[1]].pval(0),
                   &_geom ._tria.
                    _set1 [_tnod[2]].pval(0),
                    _xpos, true, 2 );

                if(_hits != geometry::null_hits)
                {
            /*--------------- call output function on hit */
                    this->_hfun (_xpos, _hits ,
                        _geom._tria .
                        _set3[_tpos].feat() ,
                        _geom._tria .
                        _set3[_tpos].topo() ,
                        _geom._tria .
                        _set3[_tpos].itag() ) ;
                    
                    this->_hnum+= +1 ;
                    
                    this->_find = true ;
                    this->_hits =_hits ;
                    
                if (_hits != geometry::face_hits)
                    this->_exct = true ;
                }
            }
        }
        
        } ;
        
    template <
    typename   hits_func
             >
    class tria_line_pred
        {
    /*------------------ tria-line intersection predicate */
        public  :
        real_type      _ipos[3] ;
        real_type      _jpos[3] ;
        real_type      _kpos[3] ;
        
        geom_type&     _geom ;
        
        hits_func&     _hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        bool_type      _exct ;
        
        geometry::
        hits_type      _hits ;
        
        public  :
        tria_line_pred operator = (
            tria_line_pred & 
            )                   =   delete ;
        tria_line_pred operator = (
            tria_line_pred&& 
            )                   =   delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call tria_line_pred (
            real_type *_isrc ,
            real_type *_jsrc ,
            real_type *_ksrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits = geometry::null_hits ;
  
            this->_exct = false;
  
            this->_hnum =     0;
            this->_find = false;
            
            for (auto _idim = 3; _idim-- != 0; )
            {
                this->
               _ipos[_idim] = _isrc[_idim];
                this->
               _jpos[_idim] = _jsrc[_idim];
                this->
               _kpos[_idim] = _ksrc[_idim];
            }
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
                geometry::hits_type 
                    _hits = geometry::null_hits ;
            
            /*--------------- line-tria intersection test */
                iptr_type  _epos = 
                    _iptr->_data.ipos() ;
                
                iptr_type  _enod[2];
                _enod[0] =_geom.
                    _tria._set2[_epos].node(0) ;
                _enod[1] =_geom.
                    _tria._set2[_epos].node(1) ;
                
                real_type  _xpos[3];
                _hits = geometry::line_tria_3d (
                   &_geom ._tria.
                    _set1 [_enod[0]].pval(0),
                   &_geom ._tria.
                    _set1 [_enod[1]].pval(0),
                   & this->_ipos[0], 
                   & this->_jpos[0],
                   & this->_kpos[0],
                    _xpos, true, 1) ;
  
                if(_hits != geometry::null_hits)
                {
            /*--------------- call output function on hit */
                    this->_hfun (_xpos, _hits ,
                        _geom._tria .
                        _set2[_epos].feat() ,
                        _geom._tria .
                        _set2[_epos].topo() ,
                        _geom._tria .
                        _set2[_epos].itag() ) ;
                    
                    this->_hnum+= +1 ;
                    
                    this->_find = true ;
                    this->_hits =_hits ;
                    
                if(_hits != geometry::face_hits)
                    this->_exct = true ;
                }
            }
        }
        
        } ;

    template <
    typename   poly_list,
    typename   hits_func
             >
    class poly_line_pred
        {
    /*------------------ poly-line intersection predicate */
        public  :
        poly_list&     _poly ;

        geom_type&     _geom ;
        
        hits_func&     _hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        geometry::
        hits_type      _hits ;
        
        public  :
        poly_line_pred operator = (
            poly_line_pred & 
            )                   =   delete ;
        poly_line_pred operator = (
            poly_line_pred&& 
            )                   =   delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call poly_line_pred (
            poly_list &_psrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _poly( _psrc),
                _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits = geometry::null_hits ;

            this->_hnum =     0 ;
            this->_find = false ;
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
                geometry::hits_type 
                    _hits = geometry::null_hits ;
            
            /*--------------- line-tria intersection test */
                iptr_type  _epos = 
                    _iptr->_data.ipos() ;
                
                iptr_type  _enod[2];
                _enod[0] =_geom.
                    _tria._set2[_epos].node(0) ;
                _enod[1] =_geom.
                    _tria._set2[_epos].node(1) ;
            
                real_type  _xpos[3];
                _hits = geometry::line_poly_3d (
                   &_geom._tria.
                    _set1[_enod[0]].pval(0) ,
                   &_geom._tria.
                    _set1[_enod[1]].pval(0) ,
                    this->_poly, _xpos) ;

                if(_hits != geometry::null_hits)
                {
            /*--------------- call output function on hit */
                    this->_hfun (_xpos, _hits ,
                        _geom._tria .
                        _set2[_epos].feat() ,
                        _geom._tria .
                        _set2[_epos].topo() ,
                        _geom._tria .
                        _set2[_epos].itag() ) ;
                    
                    this->_hnum+= +1 ;
                    
                    this->_find = true ;
                    this->_hits =_hits ;
                }
            }
        }
        
        } ;


    template <
    typename   hits_func
             >
    class ball_line_pred
        {
    /*------------------ tria-line intersection predicate */
        public  :
        real_type      _ball[3] ;
        real_type      _rsiz ;
        
        geom_type&     _geom ;
        
        hits_func&     _hfun ;
        
        geometry::
        hits_type      _hits ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        public  :
        ball_line_pred operator = (
            ball_line_pred & 
            )                   =   delete ;
        ball_line_pred operator = (
            ball_line_pred&& 
            )                   =   delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call ball_line_pred  (
            real_type *_cmid ,
            real_type  _rsiz ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits    =geometry::null_hits ;
        
            this->_find    = false ;
            this->_hnum    = + 0   ;
        
            this->_ball[0] = _cmid[0] ;
            this->_ball[1] = _cmid[1] ;
            this->_ball[2] = _cmid[2] ;
            
            this->_rsiz    = _rsiz ;
        }
    /*----------------------- all intersection about node */
        __normal_call  void_type operator()  (
            typename  
            tree_type::item_data*_iptr
            )
        {
            for ( ; _iptr != nullptr; 
                        _iptr = _iptr->_next )
            {
                geometry::hits_type _hits = 
                    geometry::face_hits ;
                
            /*--------------- ball_line intersection test */
                iptr_type _epos = 
                    _iptr->_data.ipos() ;
                
                iptr_type _enod[2];
                _enod[0] =_geom.
                    _tria._set2[_epos].node(0) ;
                _enod[1] =_geom.
                    _tria._set2[_epos].node(1) ;
                
                real_type _ipos[3];
                real_type _jpos[3];
	            size_t    _nhit = 
	                geometry::ball_line_3d (
		             this->_ball, 
                     this->_rsiz, 
		           &_geom ._tria.
		            _set1 [_enod[0]].pval(0) ,
                   &_geom ._tria.
                    _set1 [_enod[1]].pval(0) , 
		            _ipos, _jpos ) ;

                switch (_nhit)
                {
            /*--------------- call output function on hit */
                case +2 :
                    {
                this->_hfun (_jpos, _hits ,
                    _geom._tria .
                    _set2[_epos].feat() ,
                    _geom._tria .
                    _set2[_epos].topo() ,
                    _geom._tria .
                    _set2[_epos].itag() ) ;
                   
                this->_hnum += +1;
                    }
            /*--------------- call output function on hit */
                case +1 :
                    {
                this->_hfun (_ipos, _hits ,
                    _geom._tria .
                    _set2[_epos].feat() ,
                    _geom._tria .
                    _set2[_epos].topo() ,
                    _geom._tria .
                    _set2[_epos].itag() ) ;
                   
                this->_hnum += +1;
                
                this->_find = true ;
                this->_hits =_hits ;
                    }
                }
            }
        }
        
        } ;
  
  
    template <
    typename   hits_func
             >
    class circ_tria_pred
        {
    /*------------------ circ-tria intersection predicate */
        public  :
        real_type      _cmid[3] ;
        real_type      _nvec[3] ;
        real_type      _rsiz ;
        
        geom_type&     _geom ;
        
        hits_func&     _hfun ;
        
        bool_type      _find ;

        public  :
        circ_tria_pred operator = (
            circ_tria_pred & 
            )                   =   delete ;
        circ_tria_pred operator = (
            circ_tria_pred&& 
            )                   =   delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call circ_tria_pred (
            real_type *_csrc ,
            real_type *_nsrc ,
            real_type  _rsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc), 
                _find( false)
        {
            for (auto _idim = 3; _idim-- != 0; )
            {
                this->
               _cmid[_idim] = _csrc[_idim];
                this->
               _nvec[_idim] = _nsrc[_idim];
            }
            this->_rsiz     = _rsrc ;
        }
    /*----------------------- all intersection about node */
        __normal_call  void_type operator()  (
            typename  
            tree_type::item_data*_iptr
            )
        {
            for ( ; _iptr != nullptr; 
                        _iptr = _iptr->_next )
            {
                geometry::hits_type _hits = 
                    geometry::face_hits ;
                
            /*--------------- tria-flat intersection test */
                iptr_type  _tpos = 
                    _iptr->_data.ipos() ;
                    
                iptr_type _tnod[3];
                _tnod[0] =_geom.
                    _tria._set3[_tpos].node(0) ;
                _tnod[1] =_geom.
                    _tria._set3[_tpos].node(1) ;
                _tnod[2] =_geom.
                    _tria._set3[_tpos].node(2) ;
                
                real_type  _apos[3] ;
	            real_type  _bpos[3] ;
                real_type  _ipos[3] ;
	            real_type  _jpos[3] ;
                size_t     _nh = 0;
	            if (geometry::tria_flat_3d (
		             this->_cmid , 
                     this->_nvec , 
		           &_geom ._tria.
		            _set1 [_tnod[0]].pval(0) ,
                   &_geom ._tria.
                    _set1 [_tnod[1]].pval(0) ,
                   &_geom ._tria.
                    _set1 [_tnod[2]].pval(0) , 
                    _apos, _bpos))
            /*--------------- circ-tria intersection test */
	            _nh=geometry::ball_line_3d (
		             this->_cmid , 
                     this->_rsiz , 
		            _apos, _bpos , 
		            _ipos, _jpos ) ;

                switch (_nh)
                {
            /*--------------- call output function on hit */
                case +2 :
                    {
                this->_hfun (_jpos, _hits ,
                    _geom._tria .
                    _set3[_tpos].feat() ,
                    _geom._tria .
                    _set3[_tpos].topo() ,
                    _geom._tria .
                    _set3[_tpos].itag() ) ;
                    }
            /*--------------- call output function on hit */ 
                case +1 :
                    {
                this->_hfun (_ipos, _hits ,
                    _geom._tria .
                    _set3[_tpos].feat() ,
                    _geom._tria .
                    _set3[_tpos].topo() ,
                    _geom._tria .
                    _set3[_tpos].itag() ) ;
            
                this->_find = true ;
                    }
                }
            }
        }
        
        } ;

    /*
    --------------------------------------------------------
     * INTERSECT: compute poly-edge intersections.
    --------------------------------------------------------
     */
     
    template <
    typename      poly_list,
    typename      hits_func
             >
    __normal_call bool_type intersect (
        poly_list &_poly,
        hits_func &_hfun
        )
    {
    /*------------------ tree-line intersection predicate */
        typedef 
        geom_tree::aabb_pred_rect_k <
             real_type, 
             iptr_type, 
             +3        >    tree_pred ; 

    /*------------------ tria-line intersection predicate */
        typedef 
        poly_line_pred <
             poly_list ,
             hits_func >    hits_pred ; 

    /*------------------ call actual intersection testing */
        if (_poly.empty()) return false ;

        real_type _rmin[3] ;
        real_type _rmax[3] ;
        _rmin [ +0] = _poly[ +0][ +0];
        _rmin [ +1] = _poly[ +0][ +1];
        _rmin [ +2] = _poly[ +0][ +2];
        
        _rmax [ +0] = _poly[ +0][ +0];
        _rmax [ +1] = _poly[ +0][ +1];
        _rmax [ +2] = _poly[ +0][ +2];

        for(auto _iter  = _poly.head() ;
                 _iter != _poly.tend() ;
               ++_iter  )
        {
        _rmin [ +0] =  std::min (
        _rmin [ +0] , (*_iter)[ +0]) ;
        _rmin [ +1] =  std::min (
        _rmin [ +1] , (*_iter)[ +1]) ;
        _rmin [ +2] =  std::min (
        _rmin [ +2] , (*_iter)[ +2]) ;

        _rmax [ +0] =  std::max (
        _rmax [ +0] , (*_iter)[ +0]) ;
        _rmax [ +1] =  std::max (
        _rmax [ +1] , (*_iter)[ +1]) ;
        _rmax [ +2] =  std::max (
        _rmax [ +2] , (*_iter)[ +2]) ;
        }

        tree_pred _pred(_rmin, _rmax) ;     
        hits_pred _func(_poly, 
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: compute tria-edge intersections.
    --------------------------------------------------------
     */
     
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        real_type *_ipos,
        real_type *_jpos,
        real_type *_kpos,
        hits_func &_hfun
        )
    {
    /*------------------ tree-line intersection predicate */
        typedef 
        geom_tree::aabb_pred_rect_k <
             real_type, 
             iptr_type, 
             +3        >    tree_pred ; 

    /*------------------ tria-line intersection predicate */
        typedef 
        tria_line_pred <
             hits_func >    hits_pred ; 

    /*------------------ call actual intersection testing */ 
        real_type _rmin[3] ;
        real_type _rmax[3] ;
        iptr_type _idim;         
        for (_idim = +3; _idim-- != +0; )
        {
        _rmin [_idim] = _ipos[_idim];
        _rmin [_idim] =  std::min (
        _rmin [_idim] , _jpos[_idim]) ;
        _rmin [_idim] =  std::min (
        _rmin [_idim] , _kpos[_idim]) ;
        
        _rmax [_idim] = _ipos[_idim];
        _rmax [_idim] =  std::max (
        _rmax [_idim] , _jpos[_idim]) ;
        _rmax [_idim] =  std::max (
        _rmax [_idim] , _kpos[_idim]) ;
        }

        tree_pred _pred(_rmin, _rmax) ;     
        hits_pred _func(_ipos, _jpos, 
                        _kpos,
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: compute ball-edge intersections.
    --------------------------------------------------------
     */
     
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        real_type *_cmid,
        real_type  _rsiz,
        hits_func &_hfun
        )
    {
    /*------------------ tree-rect intersection predicate */
        typedef 
        geom_tree::aabb_pred_rect_k <
             real_type, 
             iptr_type, 
             +3        >    tree_pred ; 
            
    /*------------------ ball-line intersection predicate */
        typedef 
        ball_line_pred <
             hits_func >    hits_pred ;

    /*------------------ call actual intersection testing */      
        real_type _bmin[ 3] = {
            _cmid[ 0] - _rsiz ,
            _cmid[ 1] - _rsiz ,
            _cmid[ 2] - _rsiz
            } ;
        real_type _bmax[ 3] = {
            _cmid[ 0] + _rsiz ,
            _cmid[ 1] + _rsiz ,
            _cmid[ 2] + _rsiz
            } ;
      
        tree_pred _pred(_bmin, _bmax) ;
        hits_pred _func(_cmid, _rsiz,
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: compute line-face intersections.
    --------------------------------------------------------
     */
    
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        real_type *_ipos,
        real_type *_jpos,
        hits_func &_hfun
        )
    {
    /*------------------ tree-line intersection predicate */
        typedef 
        geom_tree::aabb_pred_line_k <
             real_type, 
             iptr_type, 
             +3        >    tree_pred ; 

    /*------------------ tria-line intersection predicate */
        typedef 
        line_tria_pred <
             hits_func >    hits_pred ; 

    /*------------------ call actual intersection testing */ 
        tree_pred _pred(_ipos, _jpos) ;
        hits_pred _func(_ipos, _jpos, 
                        *this, _hfun) ;

        this->_tbox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: compute disk-face intersections.
    --------------------------------------------------------
     */
     
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        real_type *_cmid,
        real_type *_nvec,
        real_type  _rsiz,
        real_type *_ball,
        hits_func &_hfun
        )
    {
        __unreferenced(_ball);

    /*------------------ tree-rect intersection predicate */
        typedef 
        geom_tree::aabb_pred_rect_k <
             real_type, 
             iptr_type, 
             +3        >    tree_pred ; 

    /*------------------ circ_tria intersection predicate */
        typedef 
        circ_tria_pred <
             hits_func >    hits_pred ;

    /*------------------ call actual intersection testing */ 
        real_type _rmin[ 3] = {
            _cmid[ 0] - _rsiz ,
            _cmid[ 1] - _rsiz ,
            _cmid[ 2] - _rsiz
            } ;
        real_type _rmax[ 3] = {
            _cmid[ 0] + _rsiz ,
            _cmid[ 1] + _rsiz ,
            _cmid[ 2] + _rsiz
            } ;
            
        tree_pred _pred(_rmin, _rmax) ;
        hits_pred _func(_cmid, _nvec, 
                        _rsiz, 
                        *this, _hfun) ;

        this->_tbox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }

    /*
    --------------------------------------------------------
     * IS-INSIDE: TRUE if point is "inside" geometry.
    --------------------------------------------------------
     */
    
    class null_pred
    {
    public  :
    __inline_call void_type operator()(
        real_type*_ppos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag 
        ) 
    {
        __unreferenced( _ppos) ;
        __unreferenced( _hits) ;
        __unreferenced( _feat) ;
        __unreferenced( _topo) ;
        __unreferenced( _itag) ;
    }
    } ;
    
    __normal_call bool_type is_inside (
        real_type *_ppos
        )
    {
    /*--------------------------- check against tree root */
        if (_ppos[0] < this->_bmin[0] ||
            _ppos[1] < this->_bmin[1] ||
            _ppos[2] < this->_bmin[2] )
            return ( false ) ;
            
        if (_ppos[0] > this->_bmax[0] ||
            _ppos[1] > this->_bmax[1] ||
            _ppos[2] > this->_bmax[2] )
            return ( false ) ;

    /*--------------------------- calc. axis-aligned dir. */  
        iptr_type _vdim = (iptr_type)+0;
        iptr_type _sign = (iptr_type)+0;
        
        real_type _vlen = 
            -std::numeric_limits
                <real_type>::infinity();
        real_type _near = 
            +std::numeric_limits
                <real_type>::infinity();

        for(iptr_type _idim = +3; _idim-- != +0; )
        {
            real_type _blen  = 
                _bmax[_idim] - 
                _bmin[_idim] ;
            
            _vlen = 
            std::max (_vlen, _blen) ;
            
        /*------------------------------ distance to bbox */
            real_type _dmin  = 
                _ppos[_idim] - 
                _bmin[_idim] ;
            real_type _dmax  = 
                _ppos[_idim] - 
                _bmax[_idim] ;

            _dmin = std::abs(_dmin) ;
            _dmax = std::abs(_dmax) ;

        /*------------------------------ keep min. so far */
            if (_near > _dmin)
            {
                _near = _dmin;
                _sign = -1 ;
                _vdim = _idim;
            }
            if (_near > _dmax)
            {
                _near = _dmax;
                _sign = +1 ;
                _vdim = _idim;
            }
        }

        _vlen *= (real_type)+2.;

    /*--------------------------- calc. "is-inside" state */
        for(iptr_type _iter =+0; _iter++ != +8; )
        {
            real_type _rvec[ +4] ;
            _rvec[0] = (real_type)+.0 ;
            _rvec[1] = (real_type)+.0 ;
            _rvec[2] = (real_type)+.0 ;
        
        /*----------------------- linear search direction */
            _rvec[_vdim] =  _sign ;

            if (_iter > +1)
            {
        /*----------------------- random search direction */
            _rvec[0] = 
        ((real_type)std::rand()) / RAND_MAX ;
            _rvec[1] = 
        ((real_type)std::rand()) / RAND_MAX ;
            _rvec[2] = 
        ((real_type)std::rand()) / RAND_MAX ;
            
            _rvec[0]-= (real_type)+.5 ;
            _rvec[1]-= (real_type)+.5 ;
            _rvec[2]-= (real_type)+.5 ;
            }

            _rvec[3] = 
            geometry::length_3d(_rvec);
            _rvec[0]*= 
                _vlen /_rvec[ +3] ;
            _rvec[1]*= 
                _vlen /_rvec[ +3] ;
            _rvec[2]*= 
                _vlen /_rvec[ +3] ;

            real_type  _rpos[ +3] = {
            _ppos[0] + _rvec[ +0] ,
            _ppos[1] + _rvec[ +1] ,
            _ppos[2] + _rvec[ +2] } ;

        /*-------------- tree-line intersection predicate */
            typedef 
            geom_tree::aabb_pred_line_k <
                 real_type, 
                 iptr_type, 
                 +3        >    tree_pred ; 

        /*-------------- tria-line intersection predicate */
            typedef 
            line_tria_pred <
                 null_pred >    hits_pred ;

            null_pred _hfun;
            tree_pred _pred(_ppos, _rpos) ;
            hits_pred _func(_ppos, _rpos, 
                            *this, _hfun) ;

            this->_tbox.find(_pred,_func) ;

            if (_func._exct) continue ;

            return (_func._hnum%+2) != +0 ;
        }

        return ( false ) ;
    }

    } ;
    
    
    }
    
#   endif   //__GEOM_MESH_EUCLIDEAN_3__



