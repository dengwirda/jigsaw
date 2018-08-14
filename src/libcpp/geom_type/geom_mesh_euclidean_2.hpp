
    /*
    --------------------------------------------------------
     * GEOM-MESH-EUCLIDEAN-2: euclidean geom. in R^2.
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
     * Last updated: 31 July, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __GEOM_MESH_EUCLIDEAN_2__
#   define __GEOM_MESH_EUCLIDEAN_2__

    namespace mesh {

    template <
    typename R ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class geom_mesh_euclidean_2d
    {
    public  :
    
    /*---------- discrete geometry in R^2 via nodes/faces */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef geom_mesh_euclidean_2d  <
                real_type ,
                iptr_type ,
                allocator    >      geom_type ;

    class node_type: public tria_complex_node_2<I,R>
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

    iptr_type static constexpr _nbox = +4;

    public  :
    
    containers::
        fixed_array<real_type,2>   _bmin ;
    containers::
        fixed_array<real_type,2>   _bmax ;

    mesh_type                      _tria ;

    tree_type                      _ebox ;

    public  :
    
    /*
    --------------------------------------------------------
     * construct geometry from alloc. etc.
    --------------------------------------------------------
     */
     
    __normal_call geom_mesh_euclidean_2d (
        allocator const& _asrc = allocator ()
        ) : _tria(_asrc) ,
            _ebox(_asrc)
        {
        }
            
    /*
    --------------------------------------------------------
     * NODE-FEAT: calc. node feature type.
    --------------------------------------------------------
     */
     
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
           (_aset.count() == +2)
                ? null_feat : soft_feat ;
    
        real_type _DtoR = 
       (real_type) +3.1415926536 / 180. ;
    
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

            real_type _ivec[2] ;
            geometry::vector_2d(
               &this->_tria.
               _set1[ _inod[0]].pval(0) ,
               &this->_tria.
               _set1[ _inod[1]].pval(0) ,
               _ivec) ;
            
            real_type _jvec[2] ;
            geometry::vector_2d(
               &this->_tria.
               _set1[ _jnod[0]].pval(0) ,
               &this->_tria.
               _set1[ _jnod[1]].pval(0) ,
               _jvec) ;
                  
            real_type _acos = 
                geometry::cosine_2d(
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
        containers::array<iptr_type> _eadj ;
    
    /*---------------------------------- init. geom feat. */
        for (auto _iter  = 
             this->_tria._set1.head() ;
                  _iter != 
             this->_tria._set1.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
                _iter->fdim () = +0  ;
                _iter->feat () = null_feat ;
                _iter->topo () = +2  ;
            }
        }
        
        for (auto _iter  = 
             this->_tria._set2.head() ;
                  _iter != 
             this->_tria._set2.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
                _iter->feat () = null_feat ;
                _iter->topo () = +2  ;
            }
        }
    
        if (_opts.feat() )
        {
    /*---------------------------------- find sharp feat. */
        for (auto _iter  = 
             this->_tria._set1.head() ;
                  _iter != 
             this->_tria._set1.tend() ;
                ++_iter  )
        {
    /*---------------------------------- find sharp 0-dim */
            if (_iter->mark() >= +0)
            {
            _eadj.set_count (0);
            
             this->_tria.node_edge (
                &_iter->node(0), _eadj) ;
                
            _iter->topo () =  
               (char_type)_eadj.count() ;

            _iter->feat () = node_feat(
                &_iter->node(0), _eadj, 
                    _opts) ;
            }
        }
        
        }
        
        for (auto _iter  = 
             this->_tria._set2.head() ;
                  _iter != 
             this->_tria._set2.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
    /*----------------------------- assign nodes to edges */
            this->_tria._set1[
                _iter->node(0)].fdim() = 1 ;
            this->_tria._set1[
                _iter->node(1)].fdim() = 1 ;
           }
        }
        
        for (auto _iter  = 
             this->_tria._set1.head() ;
                  _iter != 
             this->_tria._set1.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
            if (_iter->feat() != null_feat )
            {
    /*----------------------------- assign nodes to feat. */
                _iter->fdim()  = +0;
            }
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
            {   return _edat.mark() >= 0 ;
            }
            } ;
            
    /*----------------------------- setup user-feat. face */
        for (auto _iter  = 
             this->_tria._set2.head() ;
                  _iter != 
             this->_tria._set2.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0 &&
                _iter->itag() <= -1 )
            {
                this->_tria._set1[
            _iter->node(0)].itag() =  -1 ;
                
                this->_tria._set1[
            _iter->node(1)].itag() =  -1 ;
            }
        }
            
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
             this->_tria._set1.head() ;
                   _iter !=
             this->_tria._set1.tend() ;
                 ++_iter  )
        {
            if (_iter->mark() >= +0 )
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
  
        real_type  _BTOL[2] ;
        _BTOL[0] =( this->_bmax[ 0] - 
                    this->_bmin[ 0] )
                 * _RTOL ;
        _BTOL[1] =( this->_bmax[ 1] - 
                    this->_bmin[ 1] ) 
                 * _RTOL ;
    
    /*-------------------- find sharp "features" in geom. */    
        find_feat (_opts);
      
    /*-------------------- make aabb-tree and init. bbox. */
        aabb_mesh(this->_tria._set1, 
                  this->_tria._set2, 
                  this->_ebox,
       _BTOL,this->_nbox, edge_pred () 
                 ) ;
    
    }

    /*
    --------------------------------------------------------
     * SEED-FEAT: setup initial node set.
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

    /*------------------------- push set of feature nodes */
        for (auto _iter  = 
             this->_tria._set1.head() ;
                  _iter != 
             this->_tria._set1.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0 )
            {
            if (_opts .feat() && 
                _iter->feat() != null_feat)
            {
    /*----------------------------- push any 'real' feat. */
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
            else
            if (_iter->itag() <= -1 )
            {
    /*----------------------------- push any 'user' feat. */
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
    }
     
    /*
    --------------------------------------------------------
     * SEED-MESH: setup initial node set.
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
    /*------------------------- well-distributed sampling */
        while (_mesh._tria._nset.count() 
                < (std::size_t)_opts.seed() + 3)
        {
            typename geom_type::
                     mesh_type::
                     node_list::_write_it _best;
            
            real_type _dmax  = (real_type) +.0 ;
            char_type _fdim  ;
            
            for (_fdim = 1; _fdim != 3; ++_fdim)
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
                        geometry::lensqr_2d(
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
    class line_line_pred
        {
    /*------------------ line-line intersection predicate */
        public  :
        real_type      _ipos[2] ;
        real_type      _jpos[2] ;
        
        geom_type&     _geom ;
        
        hits_func&     _hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        bool_type      _exct ;
        
        geometry::
        hits_type      _hits ;
        
        public  :

        line_line_pred operator = (
            line_line_pred & 
            )                   =     delete ;
        line_line_pred operator = (
            line_line_pred&& 
            )                   =     delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call line_line_pred (
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
        
            for (auto _idim = 2; _idim-- != 0; )
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
            
            /*--------------- line-line intersection test */
                iptr_type  _epos = 
                    _iptr->_data.ipos() ;
    
                iptr_type  _enod[2];
                _enod[0] =_geom.
                    _tria._set2[_epos].node(0) ;
                _enod[1] =_geom.
                    _tria._set2[_epos].node(1) ;
                
                real_type  _xpos[2];
                _hits = geometry::line_line_2d (
                   & this->_ipos[0], 
                   & this->_jpos[0],
                   &_geom ._tria.
                    _set1 [_enod[0]].pval(0),
                   &_geom ._tria.
                    _set1 [_enod[1]].pval(0),
                    _xpos, true, 2 );

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
                    
                if (_hits != geometry::edge_hits)
                    this->_exct = true ;
                }
            }
        }
        
        } ;
     
    template <
    typename   hits_func
             >
    class ball_line_pred
        {
    /*------------------ ball-line intersection predicate */
        public  :
        real_type      _ball[2] ;
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
            this->_hits    =
                geometry::null_hits ;
        
            this->_find    = false ;
            this->_hnum    = + 0   ;
        
            this->_ball[0] = _cmid[0] ;
            this->_ball[1] = _cmid[1] ;
            
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
                
                real_type _ipos[2];
                real_type _jpos[2];
                size_t    _nhit = 
                    geometry::ball_line_2d (
                     this->_ball, 
                     this->_rsiz, 
                   &_geom ._tria.
                    _set1 [_enod[0]].pval(0) ,
                   &_geom ._tria.
                    _set1 [_enod[1]].pval(0) , 
                    _ipos, _jpos ) ;

                switch (_nhit)
                {
                case +2 :
                    {
            /*--------------- call output function on hit */
                this->_hfun (_jpos, _hits ,
                    _geom._tria .
                    _set2[_epos].feat() ,
                    _geom._tria .
                    _set2[_epos].topo() ,
                    _geom._tria .
                    _set2[_epos].itag() ) ;
                   
                this->_hnum += +1;
                    }
                
                case +1 :
                    {
            /*--------------- call output function on hit */
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
             +2        >    tree_pred ; 
            
    /*------------------ ball-line intersection predicate */
        typedef 
        ball_line_pred <
             hits_func >    hits_pred ;

    /*------------------ call actual intersection testing */      
        real_type _bmin[ 2] = {
            _cmid[ 0] - _rsiz ,
            _cmid[ 1] - _rsiz
            } ;
        real_type _bmax[ 2] = {
            _cmid[ 0] + _rsiz ,
            _cmid[ 1] + _rsiz
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
     * INTERSECT: compute line-edge intersections.
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
             +2        >    tree_pred ; 

    /*------------------ tria-line intersection predicate */
        typedef 
        line_line_pred <
             hits_func >    hits_pred ; 

    /*------------------ call actual intersection testing */ 
        tree_pred _pred(_ipos, _jpos) ;
        hits_pred _func(_ipos, _jpos, 
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

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
            _ppos[1] < this->_bmin[1] )
            return ( false ) ;
            
        if (_ppos[0] > this->_bmax[0] ||
            _ppos[1] > this->_bmax[1] )
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

        for(iptr_type _idim = +2; _idim-- != +0; )
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
            real_type _rvec[ +3] ;
            _rvec[0] = (real_type)+.0 ;
            _rvec[1] = (real_type)+.0 ;
            
        /*----------------------- linear search direction */
            _rvec[_vdim] =  _sign ;

            if (_iter > +1)
            {
        /*----------------------- random search direction */
            _rvec[0] = 
        ((real_type)std::rand()) / RAND_MAX ;
            _rvec[1] = 
        ((real_type)std::rand()) / RAND_MAX ;
            
            _rvec[0]-= (real_type)+.5 ;
            _rvec[1]-= (real_type)+.5 ;
            }

            _rvec[2] = 
            geometry::length_2d(_rvec);
            _rvec[0]*= 
                _vlen /_rvec[ +2] ;
            _rvec[1]*= 
                _vlen /_rvec[ +2] ;
            
            real_type  _rpos[ +2] = {
            _ppos[0] + _rvec[ +0] ,
            _ppos[1] + _rvec[ +1] } ;

        /*-------------- tree-line intersection predicate */
            typedef 
            geom_tree::aabb_pred_line_k <
                 real_type, 
                 iptr_type, 
                 +2        >    tree_pred ; 

        /*-------------- tria-line intersection predicate */
            typedef 
            line_line_pred <
                 null_pred >    hits_pred ;

            null_pred _hfun;
            tree_pred _pred(_ppos, _rpos) ;
            hits_pred _func(_ppos, _rpos, 
                            *this, _hfun) ;

            this->_ebox.find(_pred,_func) ;

            if (_func._exct) continue ;
            
            return (_func._hnum%+2) != +0 ;
        }

        return ( false ) ;
    }

    } ;   
    
    
    }
    
#   endif   //__GEOM_MESH_EUCLIDEAN_2__



