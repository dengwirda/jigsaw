
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
     * Last updated: 08 December, 2019
     *
     * Copyright 2013-2019
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
    class geom_mesh_euclidean_3d : 
            public  geom_base_3d<R, I>
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


    class node_type: public mesh_complex_node_3<I,R>
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
        
    class edge_type: public mesh_complex_edge_2<I>
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
        
    class tri3_type: public mesh_complex_tria_3<I>
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
        
    class part_data
        {
    /*------------------------------------ loc. part type */
        public  :
        containers::fixed_array <
                iptr_type,  + 3 >     _ints ;
        
        public  :
    /*------------------------------------ "write" access */
        __inline_call iptr_type&       itag (
            )
        {   return  this->_ints[0] ;
        }
        __inline_call iptr_type&       indx (
            )
        {   return  this->_ints[1] ;
        }
        __inline_call iptr_type&       kind (
            )
        {   return  this->_ints[2] ;
        }
    /*------------------------------------ "const" access */
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_ints[0] ;
        }
        __inline_call iptr_type const& indx (
            ) const
        {   return  this->_ints[1] ;
        }
        __inline_call iptr_type const& kind (
            ) const
        {   return  this->_ints[2] ;
        }
        
        } ;
        
    #define __hashscal \
        sizeof(iptr_type)/sizeof(uint32_t)
        
    class part_hash
        {
        public  :
        __inline_call iptr_type operator() (
            part_data const&_pdat
            ) const
        {
    /*------------------------------- hash for part index */
            return hash::hashword (
            (uint32_t*)&_pdat._ints[+0], 
                    +3 * __hashscal, +137) ;
        }
        } ;
        
    class part_same
        {
        public  :
        __inline_call bool_type operator() (
            part_data const&_idat,
            part_data const&_jdat
            ) const
        {
    /*------------------------------- "equal-to" for part */
            return _idat.itag() == 
                   _jdat.itag() &&
                   _idat.indx() == 
                   _jdat.indx() &&
                   _idat.kind() == 
                   _jdat.kind()  ;
        }
        } ;
        
    #undef __hashscal
   
    iptr_type static 
        constexpr pool_byte_size =  96 * 1024 ;

    typedef allocators::_pool_alloc <
            allocators::basic_alloc ,
            pool_byte_size      >   pool_base ;
    typedef allocators::_wrap_alloc <
                pool_base       >   pool_wrap ;

    typedef containers::hash_table  <
                part_data , 
                part_hash , 
                part_same ,
                pool_wrap       >   part_list ;
    
    typedef typename 
            part_list::item_type    part_item ;
            
    iptr_type static constexpr 
            part_bytes    = sizeof (part_item);
  
    typedef containers::array   <
                iptr_type ,
                allocator       >   iptr_list ;
    
    typedef containers::fixed_array <
                real_type , 3   >   bbox_bnds ;
    
    typedef containers::array   <
                bbox_bnds ,
                allocator       >   bbox_list ;
        
    typedef mesh::tria_complex_2<
                node_type,
                edge_type,
                tri3_type,
                allocator       >   mesh_type ;

    typedef geom_tree::aabb_node_base_k
                                    tree_node ;

    typedef geom_tree::aabb_item_rect_k <
                    float,
                iptr_type,
                node_type::
                    _dims       >   tree_item ;
    typedef geom_tree::aabb_tree<
                tree_item,
                node_type::
                    _dims,
                tree_node,
                allocator       >   tree_type ;

    iptr_type static constexpr _nbox = 4 ;

    public  :
    
    pool_base                      _pool ;
    
    part_list                      _part ;
    iptr_list                      _ptag ;
    
    bbox_list                      _pmin ;
    bbox_list                      _pmax ;
    
    bbox_bnds                      _bmin ;
    bbox_bnds                      _bmax ;

    mesh_type                      _tria ;

    tree_type                      _ebox ;
    tree_type                      _tbox ;

    public  :
    
    /*
    --------------------------------------------------------
     * construct geometry from alloc. etc.
    --------------------------------------------------------
     */
     
    __normal_call geom_mesh_euclidean_3d (
        allocator const&
            _asrc = allocator ()            
        ) : _pool(part_bytes) ,
            _part(part_hash() , 
                  part_same() , 
        .8, (pool_wrap(&_pool))) ,
            _ptag(    _asrc ) ,
            _pmin(    _asrc ) ,
            _pmax(    _asrc ) ,
            _tria(    _asrc ) ,
            _ebox(    _asrc ) ,
            _tbox(    _asrc )
        {
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
        return ( !this->_tbox.empty() ) ;
        else
        if (_fdim == +1)
        return ( !this->_ebox.empty() ) ;
        else
        return ( false ) ;
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
             auto _feat  = 
                this->
            _tria._set2[ *_iter].feat() ;
        
            if (_feat != null_feat)
            {
                _ldst.push_tail(*_iter) ;
            }
        }
    }
     
    template <
        typename  list_type ,
        typename  geom_opts
             >
    __normal_call void_type node_feat (
        iptr_type *_node ,
        list_type &_aset ,
        char_type &_feat ,
        char_type &_topo ,
        geom_opts &_opts
        )
    {
        real_type _DtoR = 
       (real_type) +3.1415926536 / 180.0;

        real_type _ZERO = -1. + 
            std::numeric_limits
                <real_type>::epsilon();
    
        real_type _phi1 = 
       (real_type)+180. - _opts.phi1();
        real_type _eta1 = 
       (real_type)+  0. + _opts.eta1();
    
        real_type _hard = 
            std::cos( _phi1 * _DtoR) ;
        real_type _soft = 
            std::cos( _eta1 * _DtoR) ;
  
        __unreferenced(_node) ;
  
        _feat =  null_feat ;
        _topo = (char_type)_aset.count();

        for (auto _ipos  = _aset.head() ;
                  _ipos != _aset.tend() ;
                ++_ipos  )
        {
        char_type _tbad  = +1 ;
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
                
            if (_acos >= _ZERO)
            {
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
            else
            {
            if (_tbad >= +  1 )
            {
                _topo -= _tbad--;
            }
            } 
        }
        }
        {
            if (_topo != +  0 )
            if (_topo != +  2 )
                _feat = 
            std::max(_feat, soft_feat) ;
        }
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
    __normal_call void_type edge_feat (
        iptr_type *_enod ,
        list_type &_aset ,
        char_type &_feat ,
        char_type &_topo ,
        geom_opts &_opts
        )
    {
        real_type _DtoR = 
       (real_type) +3.1415926536 / 180.0;

        real_type _ZERO = -1. + 
            std::numeric_limits
                <real_type>::epsilon();
    
        real_type _phi2 = 
       (real_type)+180. - _opts.phi2();
        real_type _eta2 = 
       (real_type)+  0. + _opts.eta2();
    
        real_type _hard = 
            std::cos( _phi2 * _DtoR) ;
        real_type _soft = 
            std::cos( _eta2 * _DtoR) ;
            
        _feat =  null_feat ;
        _topo = (char_type)_aset.count();
    
        for (auto _ipos  = _aset.head() ;
                  _ipos != _aset.tend() ;
                ++_ipos  )
        {
        char_type _tbad  = +1 ;
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
           
            if (_acos >= _ZERO)
            {
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
            else
            {
            if (_tbad >= +  1 )
            {
                _topo -= _tbad--;
            }
            }
        }
        }
        {
            if (_topo != +  0 )
            if (_topo != +  2 )
                _feat = 
            std::max(_feat, soft_feat) ;
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
        containers::
            array <iptr_type> _eadj ;
        containers::
            array <iptr_type> _fadj ;
        containers::
            array <iptr_type> _ebnd ;
    
        containers::
            array <iptr_type> _nmrk ;
        containers::
            array <iptr_type> _emrk ;
    
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
        for (auto _iter  = 
             this->_tria._set3.head() ;
                  _iter != 
             this->_tria._set3.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
                _iter->feat () = null_feat ;
                _iter->topo () = +2  ;
            }
        }
  
    /*---------------------------------- find sharp feat. */
        _nmrk.set_count (
         this->_tria._set1.count() , 
            containers::loose_alloc,-1) ;
            
        _emrk.set_count (
         this->_tria._set2.count() , 
            containers::loose_alloc,-1) ;
        
        iptr_type _nnum  = +0 ;
        iptr_type _enum  = +0 ;
    
        for (auto _epos  = 
             this->_tria._set2.head() ;
                  _epos != 
             this->_tria._set2.tend() ;
                ++_epos, ++_enum)
        {
    /*---------------------------------- find sharp 1-dim */
            if (_epos->mark() >= +0)
            {
            if (_epos->itag() <= -1 ||
                _emrk[_enum]  >= +1 ||
                _opts .feat() )
            {
    /*---------------------------------- set geo.-defined */
            _fadj.set_count (0);
            
            this->_tria.edge_tri3 (
               &_epos->node (0), _fadj) ;
                    
            edge_feat (
               &_epos->node (0), 
                _fadj ,
                _epos->feat () ,
                _epos->topo () , 
                _opts ) ;
             
            if (_epos->itag() <= -1)
            {
    /*---------------------------------- set user-defined */
            _epos->feat () = 
                std::max(_epos->feat () , 
                    soft_feat) ;
                    
            _nmrk[_epos->node(0)] = +1;
            _nmrk[_epos->node(1)] = +1;
            }
            }
            }
        }
        
        for (auto _npos  = 
             this->_tria._set1.head() ;
                  _npos != 
             this->_tria._set1.tend() ;
                ++_npos, ++_nnum)
        {
    /*---------------------------------- find sharp 0-dim */
            if (_npos->mark() >= +0 )
            {
            if (_npos->itag() <= -1 ||
                _nmrk[_nnum]  >= +1 ||
                _opts .feat() )
            {
    /*---------------------------------- set geo.-defined */
            _eadj.set_count (0);
            _fadj.set_count (0);
            
            this->_tria.node_edge (
               &_npos->node (0), _eadj) ;
                
            this->_tria.node_tri3 (
               &_npos->node (0), _fadj) ;

            _ebnd.set_count (0);

            feat_list (_eadj, _ebnd);

            node_feat (
               &_npos->node (0), 
                _ebnd ,
                _npos->feat () ,
                _npos->topo () , 
                _opts ) ;
            
            if (_npos->itag () <= -1)
            {
    /*---------------------------------- set user-defined */
            _npos->feat () = 
                std::max(_npos->feat () , 
                    soft_feat) ;
            }
            }
            }
        }
  
        for (auto _iter  = 
             this->_tria._set3.head() ;
                  _iter != 
             this->_tria._set3.tend() ;
                ++_iter  )
        {
    /*----------------------------- assign nodes to trias */
            if (_iter->mark() >= +0)
            {
            this->_tria._set1[
                _iter->node(0)].fdim() = 2;
            this->_tria._set1[
                _iter->node(1)].fdim() = 2;
            this->_tria._set1[
                _iter->node(2)].fdim() = 2;
            }
        }
        
        for (auto _iter  = 
             this->_tria._set2.head() ;
                  _iter != 
             this->_tria._set2.tend() ;
                ++_iter  )
        {
    /*----------------------------- assign nodes to edges */
            if (_iter->mark() >= +0)
            {
            if (_iter->feat() != null_feat)
            {
            this->_tria._set1[
                _iter->node(0)].fdim() = 1;
            this->_tria._set1[
                _iter->node(1)].fdim() = 1;
            }
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
            if (_iter->feat() != null_feat)
            {
    /*----------------------------- assign nodes to feat. */
                _iter->fdim()  = 0 ;
            }
            }
        }
    }
    
    /*
    --------------------------------------------------------
     * INIT-PART: init. enclosed PART within geom.
    --------------------------------------------------------
     */
     
    template <
        typename  geom_opts
             >
    __normal_call void_type init_part (
        geom_opts &_opts
        )
    {
        containers::array <
                typename 
            iptr_list::size_type> _mark ;
    
        __unreferenced (_opts) ;

    /*----------------------------- extrema for PART id's */
        iptr_type _imin = 
            +std::numeric_limits
                <iptr_type>::infinity() ;
        iptr_type _imax = 
            -std::numeric_limits
                <iptr_type>::infinity() ;
    
        for (auto _iter  = 
             this->_part._lptr.head() ; 
                  _iter != 
             this->_part._lptr.tend() ;
                ++_iter  )
        {
            typename part_list
                   ::item_type*_next =*_iter;
            
            for ( ; _next != nullptr ; 
                        _next = _next->_next)
            {
                iptr_type _ppos = 
                    _next->_data.itag() ;
            
                _imin = std::min(_imin, 
                                 _ppos) ;
                _imax = std::max(_imax, 
                                 _ppos) ;
            }
        }
 
        if (_imin < (iptr_type) +0 &&
            _imax < (iptr_type) +0 )
            __assert( _imin >=  +0 && 
        "GEOM::INIT-PART: -ve part index!") ;
        
    /*----------------------------- push unique PART id's */
        auto _flag = 
            std::numeric_limits <
                typename 
            iptr_list::size_type>::max();
        
        _mark.set_count(_imax+1, 
        containers::tight_alloc, _flag) ;
    
        for (auto _iter  = 
             this->_part._lptr.head() ; 
                  _iter != 
             this->_part._lptr.tend() ;
                ++_iter  )
        {
            typename part_list
                   ::item_type*_next =*_iter;
            
            for ( ; _next != nullptr ; 
                        _next = _next->_next)
            {
                iptr_type _ppos  = 
                    _next->_data.itag() ;

                if (_mark[_ppos] == _flag)
                {
                    _mark[_ppos] = 
                        this->
                    _ptag.push_tail(_ppos) ;
                }
            }
        }
        
    /*----------------------------- init. aabb's for PART */
        this->_pmin.set_count (
            this->_ptag.count(), 
        containers::tight_alloc, _bmax) ;
        
        this->_pmax.set_count (
            this->_ptag.count(), 
        containers::tight_alloc, _bmin) ;
        
        for (auto _iter  = 
             this->_part._lptr.head() ; 
                  _iter != 
             this->_part._lptr.tend() ;
                ++_iter  )
        {
            typename part_list
                   ::item_type*_next =*_iter;
            
            for ( ; _next != nullptr ; 
                        _next = _next->_next)
            {
                iptr_type _ppos = 
                    _next->_data.itag () ;
                iptr_type _cell = 
                    _next->_data.indx () ;
                iptr_type _kind = 
                    _next->_data.kind () ;
                    
                if (_kind == TRIA3_tag)
                {
    /*----------------------------- expand aabb via TRIA3 */   
                auto _inod =  this->
                    _tria._set3[_cell].node(0);
                auto _jnod =  this->
                    _tria._set3[_cell].node(1);
                auto _knod =  this->
                    _tria._set3[_cell].node(2);
                    
                auto _iptr = &this->
                    _tria._set1[_inod];
                auto _jptr = &this->
                    _tria._set1[_jnod];
                auto _kptr = &this->
                    _tria._set1[_knod];
        
                real_type _xmin = 
                std::min (_iptr->pval(0) ,
                std::min (_jptr->pval(0) ,
                          _kptr->pval(0) 
                         ) ) ;
                real_type _ymin = 
                std::min (_iptr->pval(1) ,
                std::min (_jptr->pval(1) ,
                          _kptr->pval(1) 
                         ) ) ;
                real_type _zmin = 
                std::min (_iptr->pval(2) ,
                std::min (_jptr->pval(2) ,
                          _kptr->pval(2) 
                         ) ) ;
                         
                real_type _xmax = 
                std::max (_iptr->pval(0) ,
                std::max (_jptr->pval(0) ,
                          _kptr->pval(0) 
                         ) ) ;
                real_type _ymax = 
                std::max (_iptr->pval(1) ,
                std::max (_jptr->pval(1) ,
                          _kptr->pval(1) 
                         ) ) ;
                real_type _zmax = 
                std::max (_iptr->pval(2) ,
                std::max (_jptr->pval(2) ,
                          _kptr->pval(2) 
                         ) ) ;
        
                auto _pmap = _mark[_ppos];

                this-> _pmin[_pmap][ 0] = 
                    std::min(_xmin,
                this-> _pmin[_pmap][ 0]) ;
                this-> _pmin[_pmap][ 1] = 
                    std::min(_ymin,
                this-> _pmin[_pmap][ 1]) ;
                this-> _pmin[_pmap][ 2] = 
                    std::min(_zmin,
                this-> _pmin[_pmap][ 2]) ;
                
                this-> _pmax[_pmap][ 0] = 
                    std::max(_xmax,
                this-> _pmax[_pmap][ 0]) ;
                this-> _pmax[_pmap][ 1] = 
                    std::max(_ymax,
                this-> _pmax[_pmap][ 1]) ;
                this-> _pmax[_pmap][ 2] = 
                    std::max(_zmax,
                this-> _pmax[_pmap][ 2]) ;
                
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
        
        float     static const _RTOL = 
            std::pow (
            std::numeric_limits<float>
            ::epsilon(), (float)+.75) ;
  
        float      _BTOL[3] ;
        _BTOL[0] =
       (float)    ( this->_bmax[ 0] - 
                    this->_bmin[ 0] )
                 * _RTOL ;
        _BTOL[1] =
       (float)    ( this->_bmax[ 1] - 
                    this->_bmin[ 1] ) 
                 * _RTOL ;
        _BTOL[2] =
       (float)    ( this->_bmax[ 2] - 
                    this->_bmin[ 2] ) 
                 * _RTOL ;
      
    /*-------------------- find sharp "features" in geom. */    
        find_feat (_opts);
        
    /*-------------------- indexes for PART def. in geom. */
        init_part (_opts);
      
    /*-------------------- make aabb-tree and init. bbox. */
        aabb_mesh(this->_tria._set1, 
                  this->_tria._set2, 
                  this->_ebox,
       _BTOL,this->_nbox,edge_pred () 
                 ) ;
                  
        aabb_mesh(this->_tria._set1, 
                  this->_tria._set3, 
                  this->_tbox,
       _BTOL,this->_nbox,tri3_pred () 
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
            if (_iter->feat() != null_feat)
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
                        
                _mesh._tria.node
                    (_node)->part() 
                        = _iter->itag() ;
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
                        
                _mesh._tria.node
                    (_node)->part() 
                        = _iter->itag() ;
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
                < (std::size_t)_opts.seed() + 4)
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
                if (_ipos->mark() >= 0 && 
                        _ipos->fdim () == _fdim)
                {
                real_type _dmin  = 
                    +std::numeric_limits
                        <real_type>::infinity();

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
    class flat_line_pred
        {
    /*------------------ flat-line intersection predicate */
        public  :
        real_type      _pmid[3] ;
        real_type      _nvec[3] ;
       
        geom_type     &_geom ;
        
        hits_func     &_hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        geometry::
        hits_type      _hits ;
        
        public  :
        flat_line_pred operator = (
            flat_line_pred & 
            )                   =   delete ;
        flat_line_pred operator = (
            flat_line_pred&& 
            )                   =   delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call flat_line_pred (
            real_type *_qsrc ,
            real_type *_vsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits    = 
                geometry ::null_hits ;

            this->_pmid[0] = _qsrc[0];
            this->_pmid[1] = _qsrc[1];
            this->_pmid[2] = _qsrc[2];
           
            this->_nvec[0] = _vsrc[0];
            this->_nvec[1] = _vsrc[1];
            this->_nvec[2] = _vsrc[2];

            this->_hnum    = + 0   ;
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
                geometry::hits_type 
                    _HITS = geometry::face_hits ;
            
            /*--------------- line-tria intersection test */
                iptr_type  _epos = 
                    _iptr->_data.ipos() ;
                
                iptr_type _enod[2];
                _enod[0] =_geom.
                    _tria._set2[_epos].node(0) ;
                _enod[1] =_geom.
                    _tria._set2[_epos].node(1) ;
            
                real_type _xpos[3];
                bool_type _okay = 
                    geometry::line_flat_3d (
                     this->_pmid, 
                     this->_nvec, 
                    &this->_geom.
                _tria._set1[_enod[0]].pval(0) , 
                    &this->_geom.
                _tria._set1[_enod[1]].pval(0) , 
                    _xpos, true);

                if (_okay)
                {
            /*--------------- call output function on hit */
                    this->_hfun (_xpos, _HITS ,
                        _geom._tria .
                        _set2[_epos].feat() ,
                        _geom._tria .
                        _set2[_epos].topo() ,
                        _geom._tria .
                        _set2[_epos].itag() ) ;
                    
                    this->_hnum+= +1 ;
                    
                    this->_find = true ;
                    this->_hits =_HITS ;
                }
            }
        }
        
        } ;

    template <
    typename   hits_func
             >
    class line_tria_pred
        {
    /*------------------ line-tria intersection predicate */
        public  :
        real_type      _ipos[3] ;
        real_type      _jpos[3] ;
        
        geom_type     &_geom ;
        
        hits_func     &_hfun ;
        
        bool_type      _find ;   
        iptr_type      _hnum ;
        
        bool_type      _exct ;
        
        geometry::
        hits_type      _hits ;
        
        iptr_type      _part ;
        
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
            hits_func &_hsrc ,
            iptr_type  _psrc = -1
            ) : _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits    = 
                geometry ::null_hits ;
     
            this->_ipos[0] = _isrc[0];
            this->_ipos[1] = _isrc[1];
            this->_ipos[2] = _isrc[2];
            
            this->_jpos[0] = _jsrc[0];
            this->_jpos[1] = _jsrc[1];
            this->_jpos[2] = _jsrc[2];
     
            this->_part    = _psrc ;
        
            this->_exct    = false ;
        
            this->_hnum    = + 0   ;
            this->_find    = false ;
        }
    /*----------------------- TRUE if PART is within list */    
        __inline_call  bool_type have_part   (
            iptr_type _tpos
            )
        {
            typename geom_type
                ::part_list::item_type *_same ;
        
            typename geom_type
                ::part_list::data_type  _temp ;
            _temp.itag() = _part ;
            _temp.indx() = _tpos ;
            _temp.kind() = 
                mesh:: TRIA3_tag ;
            
            return  this->
                _geom._part.find(_temp, _same);
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
                    _HITS = geometry::null_hits;
            
            /*--------------- line-tria intersection test */
                iptr_type  _tpos = 
                    _iptr->_data.ipos() ;
    
                if (this->_part >= +0 
                        && !have_part(_tpos) ) 
                    continue ;
    
                iptr_type  _tnod[3];
                _tnod[0] =_geom.
                    _tria._set3[_tpos].node(0) ;
                _tnod[1] =_geom.
                    _tria._set3[_tpos].node(1) ;
                _tnod[2] =_geom.
                    _tria._set3[_tpos].node(2) ;
                
                real_type  _xpos[3];
                _HITS = geometry::line_tria_3d (
                   & this->_ipos[0], 
                   & this->_jpos[0],
                   &_geom ._tria.
                    _set1 [_tnod[0]].pval(0),
                   &_geom ._tria.
                    _set1 [_tnod[1]].pval(0),
                   &_geom ._tria.
                    _set1 [_tnod[2]].pval(0),
                    _xpos, true, 2 );

                if(_HITS != geometry::null_hits)
                {
            /*--------------- call output function on hit */
                    this->_hfun (_xpos, _HITS ,
                        _geom._tria .
                        _set3[_tpos].feat() ,
                        _geom._tria .
                        _set3[_tpos].topo() ,
                        _geom._tria .
                        _set3[_tpos].itag() ) ;
                    
                    this->_hnum+= +1 ;
                    
                    this->_find = true ;
                    this->_hits =_HITS ;
                    
                if (_HITS != geometry::face_hits)
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
        real_type      _ball[3] ;
        real_type      _rsiz ;
        
        geom_type     &_geom ;
        
        hits_func     &_hfun ;
        
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
            real_type  _rsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc)
        {
            this->_hits    =
                geometry ::null_hits ;
        
            this->_find    = false ;
            this->_hnum    = + 0   ;
        
            this->_ball[0] = _cmid[0];
            this->_ball[1] = _cmid[1];
            this->_ball[2] = _cmid[2];
            
            this->_rsiz    = _rsrc ;
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
                geometry::hits_type _HITS = 
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
                case +2 :
                    {
            /*--------------- call output function on hit */
                this->_hfun (_jpos, _HITS ,
                    _geom._tria .
                    _set2[_epos].feat() ,
                    _geom._tria .
                    _set2[_epos].topo() ,
                    _geom._tria .
                    _set2[_epos].itag() ) ;
                   
                this->_hnum += +1;
                    }      // falls through
                case +1 :
                    {
            /*--------------- call output function on hit */
                this->_hfun (_ipos, _HITS ,
                    _geom._tria .
                    _set2[_epos].feat() ,
                    _geom._tria .
                    _set2[_epos].topo() ,
                    _geom._tria .
                    _set2[_epos].itag() ) ;
                   
                this->_hnum += +1;
                
                this->_find = true ;
                this->_hits =_HITS ;
                    }      // falls through
                }
            }
        }
        
        } ;
   
    template <
    typename   hits_func
             >
    class disc_tria_pred
        {
    /*------------------ circ-tria intersection predicate */
        public  :
        real_type      _cmid[3] ;
        real_type      _nvec[3] ;
        real_type      _rsiz ;
        
        geom_type     &_geom ;
        
        hits_func     &_hfun ;
        
        bool_type      _find ;

        public  :
        disc_tria_pred operator = (
            disc_tria_pred & 
            )                   =   delete ;
        disc_tria_pred operator = (
            disc_tria_pred&& 
            )                   =   delete ;

        public  :
    /*------------------------------ construct from _src. */
        __normal_call disc_tria_pred (
            real_type *_csrc ,
            real_type *_nsrc ,
            real_type  _rsrc ,
            geom_type &_gsrc ,
            hits_func &_hsrc
            ) : _geom( _gsrc), 
                _hfun( _hsrc), 
                _find( false)
        {            
            this->_cmid[0] = _csrc[0];
            this->_cmid[1] = _csrc[1];
            this->_cmid[2] = _csrc[2];
            
            this->_nvec[0] = _nsrc[0];
            this->_nvec[1] = _nsrc[1];
            this->_nvec[2] = _nsrc[2];
            
            this->_rsiz    = _rsrc;
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
                geometry::hits_type _HITS = 
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
                case +2 :
                    {
            /*--------------- call output function on hit */
                this->_hfun (_jpos, _HITS ,
                    _geom._tria .
                    _set3[_tpos].feat() ,
                    _geom._tria .
                    _set3[_tpos].topo() ,
                    _geom._tria .
                    _set3[_tpos].itag() ) ;
                    }      // falls through
                case +1 :
                    {
            /*--------------- call output function on hit */
                this->_hfun (_ipos, _HITS ,
                    _geom._tria .
                    _set3[_tpos].feat() ,
                    _geom._tria .
                    _set3[_tpos].topo() ,
                    _geom._tria .
                    _set3[_tpos].itag() ) ;
            
                this->_find = true ;
                    }      // falls through
                }
            }
        }
        
        } ;

    class near_edge_pred
        {
    /*------------------ node-edge "projection" predicate */
        public  :
        real_type              *_ppos ;
        real_type              *_qpos ;

        real_type               _dsqr ;
        real_type               _dtol ;
        
        mesh_type              *_mesh ;
        
        bool_type               _find ;
        iptr_type               _epos ;
        
        public  :
    
    /*------------------------ make a tree-edge predicate */
        __inline_call near_edge_pred(
            real_type*_psrc = nullptr ,
            real_type*_qsrc = nullptr ,
            mesh_type*_msrc = nullptr ,
            real_type _near = 
                (real_type) +0.
            ) : _ppos(_psrc) ,
                _qpos(_qsrc) ,
           _dsqr(+std::numeric_limits
           <real_type>::infinity () ) ,
                _dtol(_near) ,
                _mesh(_msrc) ,
                _find(false) ,
                _epos(   -1) { }

    /*------------------------ call pred. on tree matches */
        __inline_call float operator () (
                typename  
            tree_type::item_data *_iptr
            )
        {
            if (this->_find) return +0. ;
        
            real_type  _qtmp[+3] = {+0.};

            for ( ; _iptr != nullptr; 
                        _iptr = _iptr->_next )
            {
                geometry::hits_type 
                    _HITS = geometry::null_hits;

                iptr_type  _EPOS = 
                    _iptr->_data.ipos() ;

                iptr_type  _enod[+2] = {
                     this->_mesh->
                    _set2 [_EPOS   ].node(0) ,
                     this->_mesh->
                    _set2 [_EPOS   ].node(1) ,
                    } ;
                
                if (geometry::proj_line_3d (
                           _ppos,
                    &this->_mesh->
                    _set1 [_enod[0]].pval(0) ,
                    &this->_mesh->
                    _set1 [_enod[1]].pval(0) , 
                    _qtmp, _HITS) )
                {
                    if (_HITS != 
                        geometry::null_hits)
                    {
        /*-------------------- projected match: keep best */
                    real_type _dtmp = 
                geometry::lensqr_3d(_ppos,_qtmp) ;

                    if (_dtmp<_dsqr )
                    {
                    _qpos[0] = _qtmp[0] ;
                    _qpos[1] = _qtmp[1] ;
                    _qpos[2] = _qtmp[2] ;

                    this->_find =
                    this->_dtol * 
                    this->_dtol > _dtmp ;

                    this->_dsqr = _dtmp ;
                    this->_epos = _EPOS ;
                    }

                    }
                }
            }

            return ( (float)this->_dsqr )  ;
        }
        
        } ;

    class near_tria_pred
        {
    /*------------------ node-tria "projection" predicate */
        public  :
        real_type              *_ppos ;
        real_type              *_qpos ;

        real_type               _dsqr ;
        real_type               _dtol ;
        
        mesh_type              *_mesh ;
        
        bool_type               _find ;
        iptr_type               _tpos ;
        
        public  :
    
    /*------------------------ make a tree-edge predicate */
        __inline_call near_tria_pred(
            real_type*_psrc = nullptr ,
            real_type*_qsrc = nullptr ,
            mesh_type*_msrc = nullptr ,
            real_type _near = 
                (real_type) +0.
            ) : _ppos(_psrc) ,
                _qpos(_qsrc) ,
           _dsqr(+std::numeric_limits
           <real_type>::infinity () ) ,
                _dtol(_near) ,
                _mesh(_msrc) ,
                _find(false) ,
                _tpos(   -1) { }

    /*------------------------ call pred. on tree matches */
        __inline_call float operator () (
                typename  
            tree_type::item_data *_iptr
            )
        {
            if (this->_find) return +0. ;
        
            real_type  _qtmp[+3] = {+0.};

            for ( ; _iptr != nullptr; 
                        _iptr = _iptr->_next )
            {
                geometry::hits_type 
                    _HITS = geometry::null_hits;

                iptr_type  _TPOS = 
                    _iptr->_data.ipos() ;

                iptr_type  _tnod[+3] = {
                     this->_mesh->
                    _set3 [_TPOS   ].node(0) ,
                     this->_mesh->
                    _set3 [_TPOS   ].node(1) ,
                     this->_mesh->
                    _set3 [_TPOS   ].node(2) ,
                    } ;
                
                if (geometry::proj_tria_3d (
                           _ppos,
                    &this->_mesh->
                    _set1 [_tnod[0]].pval(0) ,
                    &this->_mesh->
                    _set1 [_tnod[1]].pval(0) ,
                    &this->_mesh->
                    _set1 [_tnod[2]].pval(0) , 
                    _qtmp, _HITS) )
                {
                    if (_HITS != 
                        geometry::null_hits)
                    {
        /*-------------------- projected match: keep best */
                    real_type _dtmp = 
                geometry::lensqr_3d(_ppos,_qtmp) ;

                    if (_dtmp<_dsqr )
                    {
                    _qpos[0] = _qtmp[0] ;
                    _qpos[1] = _qtmp[1] ;
                    _qpos[2] = _qtmp[2] ;

                    this->_find =
                    this->_dtol * 
                    this->_dtol > _dtmp ;

                    this->_dsqr = _dtmp ;
                    this->_tpos = _TPOS ;
                    }

                    }
                }
            }

            return ( (float)this->_dsqr )  ;
        }
        
        } ;

    
    /*
    --------------------------------------------------------
     * INTERSECT: find FLAT/1-GEOM. intersections.
    --------------------------------------------------------
     */
    
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        flat_type &_flat,
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
        flat_line_pred <
             hits_func >    hits_pred ; 

        float           _PPOS[3] = {
                (float) _flat. _ppos[0] ,
                (float) _flat. _ppos[1] ,
                (float) _flat. _ppos[2] ,
                } ;

        float           _NVEC[3] = {
                (float) _flat. _nvec[0] ,
                (float) _flat. _nvec[1] ,
                (float) _flat. _nvec[2] ,
                } ;

        float           _RMIN[3] = {
                (float) _flat. _rmin[0] ,
                (float) _flat. _rmin[1] ,
                (float) _flat. _rmin[2] ,
                } ;

        float           _RMAX[3] = {
                (float) _flat. _rmax[0] ,
                (float) _flat. _rmax[1] ,
                (float) _flat. _rmax[2] ,
                } ;

    /*------------------ call actual intersection testing */
        tree_pred _pred(_PPOS, _NVEC,
                        _RMIN, _RMAX) ;     
        hits_pred _func(_flat. _ppos, 
                        _flat. _nvec,
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
    typename      hits_func
             >
    __normal_call bool_type intersect (
        line_type &_line,
        hits_func &_hfun
        )
    {
    /*------------------ tree-line intersection predicate */
        typedef 
        geom_tree::aabb_pred_line_3 <
                 float , 
             iptr_type >    tree_pred ;
             
    /*------------------ tria-line intersection predicate */
        typedef 
        line_tria_pred <
             hits_func >    hits_pred ; 

        float           _IPOS[3] = {
                (float) _line. _ipos[0] ,
                (float) _line. _ipos[1] ,
                (float) _line. _ipos[2] ,
                } ;

        float           _JPOS[3] = {
                (float) _line. _jpos[0] ,
                (float) _line. _jpos[1] ,
                (float) _line. _jpos[2] ,
                } ;

    /*------------------ call actual intersection testing */
        tree_pred _pred(_IPOS, _JPOS) ;
        hits_pred _func(_line. _ipos, 
                        _line. _jpos, 
                        *this, _hfun) ;

        this->_tbox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }
    
    /*
    --------------------------------------------------------
     * INTERSECT: find BALL/1-GEOM. intersections.
    --------------------------------------------------------
     */
     
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        ball_type &_ball,
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
        ball_line_pred <
             hits_func >    hits_pred ;

        float           _PMID[3] = {
                (float) _ball. _pmid[0] ,
                (float) _ball. _pmid[1] ,
                (float) _ball. _pmid[2] ,
                } ;

        float           _RRAD  =
                (float) _ball. _rrad;

    /*------------------ call actual intersection testing */          
        tree_pred _pred(_PMID, _RRAD) ;
        hits_pred _func(_ball. _pmid, 
                        _ball. _rrad,
                        *this, _hfun) ;

        this->_ebox.find(_pred,_func) ;

    /*------------------ _TRUE if any intersections found */
        return (   _func._find ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: find DISC/2-GEOM. intersections.
    --------------------------------------------------------
     */
     
    template <
    typename      hits_func
             >
    __normal_call bool_type intersect (
        disc_type &_disc,
        real_type *_sbal,
        hits_func &_hfun
        )
    {
    /*------------------ tree-ball intersection predicate */
        typedef 
        geom_tree::aabb_pred_ball_3 <
                 float , 
             iptr_type >    tree_pred ; 

    /*------------------ circ_tria intersection predicate */
        typedef 
        disc_tria_pred <
             hits_func >    hits_pred ;

        __unreferenced(_sbal) ;

        float           _PMID[3] = {
                (float) _disc. _pmid[0] ,
                (float) _disc. _pmid[1] ,
                (float) _disc. _pmid[2] ,
                } ;

        float           _RRAD  =
                (float) _disc. _rrad;

    /*------------------ call actual intersection testing */ 
        tree_pred _pred(_PMID, _RRAD) ;
        hits_pred _func(_disc. _pmid, 
                        _disc. _nvec, 
                        _disc. _rrad, 
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
    
    __normal_call iptr_type is_inside (
        real_type *_ppos
        )
    {
    /*--------------------------- calc. axis-aligned dir. */  
        iptr_type _vdim = (iptr_type)+0;
        iptr_type _sign = (iptr_type)+0;
        
        real_type _vlen = 
           -std::numeric_limits
                <real_type>::infinity();
        real_type _near = 
           +std::numeric_limits
                <real_type>::infinity();

        for(auto _idim = +3; _idim-- != +0; )
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

        _vlen *= (real_type)+2. ;

    /*--------------------------- calc. "is-inside" state */
        if (this->_ptag.empty() )
        {
        
    /*--------------------------- null PART specification */
        
        if (_ppos[0] < this->_bmin[0] ||
            _ppos[1] < this->_bmin[1] ||
            _ppos[2] < this->_bmin[2] )
            return (iptr_type) -1 ;
            
        if (_ppos[0] > this->_bmax[0] ||
            _ppos[1] > this->_bmax[1] ||
            _ppos[2] > this->_bmax[2] )
            return (iptr_type) -1 ;
        
        for(auto _iter = +0; _iter++ != +8; )
        {
            real_type _rvec[ 4] ;
            _rvec[0] = (real_type)+.0 ;
            _rvec[1] = (real_type)+.0 ;
            _rvec[2] = (real_type)+.0 ;
            
        /*----------------------- linear search direction */
            _rvec[_vdim] =  _sign * 
                            _vlen ;

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
            
            _rvec[3] = 
            geometry::length_3d(_rvec);
            _rvec[0]*= 
                _vlen /_rvec[ +3] ;
            _rvec[1]*= 
                _vlen /_rvec[ +3] ;
            _rvec[2]*= 
                _vlen /_rvec[ +3] ;
            }

            real_type  _rpos[ +3] = {
            _ppos[0] + _rvec[ +0] ,
            _ppos[1] + _rvec[ +1] ,
            _ppos[2] + _rvec[ +2] } ;

        /*-------------- tree-line intersection predicate */
            typedef 
            geom_tree::aabb_pred_line_3 <
                     float , 
                 iptr_type >  tree_pred ; 

        /*-------------- tria-line intersection predicate */
            typedef 
            line_tria_pred <
                 null_pred >    hits_pred ;

            float           _PPOS[3] = {
                    (float) _ppos[0] ,
                    (float) _ppos[1] ,
                    (float) _ppos[2] ,
                    } ;

            float           _RPOS[3] = {
                    (float) _rpos[0] ,
                    (float) _rpos[1] ,
                    (float) _rpos[2] ,
                    } ;

            tree_pred _pred(_PPOS, _RPOS) ;

            null_pred _hfun;
            hits_pred _func(_ppos, _rpos, 
                            *this, _hfun) ;

            this->_tbox.find(_pred,_func) ;

            if (_func._exct) continue ;
            
        /*-------------- done if inside - odd no. crosses */
            if((_func._hnum % 2) != +0)
            {
            return ( (iptr_type) +0 ) ;
            }
        }
        
        }
        else
        {
        
    /*--------------------------- have PART specification */
        
        for(auto _pnum = this->_ptag.count(); 
                 _pnum-- != +0; )
        {        
        
        if (_ppos[0] < 
             this->_pmin[_pnum][0] ||
            _ppos[1] < 
             this->_pmin[_pnum][1] ||
            _ppos[2] < 
             this->_pmin[_pnum][2] )
            continue ;
            
        if (_ppos[0] > 
             this->_pmax[_pnum][0] ||
            _ppos[1] > 
             this->_pmax[_pnum][1] ||
            _ppos[2] > 
             this->_pmax[_pnum][2] )
            continue ;
        
        for(auto _iter = +0; _iter++ != +8; )
        {
            real_type _rvec[ 4] ;
            _rvec[0] = (real_type)+.0 ;
            _rvec[1] = (real_type)+.0 ;
            _rvec[2] = (real_type)+.0 ;
            
        /*----------------------- linear search direction */
            _rvec[_vdim] =  _sign * 
                            _vlen ;

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
            
            _rvec[3] = 
            geometry::length_3d(_rvec);
            _rvec[0]*= 
                _vlen /_rvec[ +3] ;
            _rvec[1]*= 
                _vlen /_rvec[ +3] ;
            _rvec[2]*= 
                _vlen /_rvec[ +3] ;
            }

            real_type  _rpos[ +3] = {
            _ppos[0] + _rvec[ +0] ,
            _ppos[1] + _rvec[ +1] ,
            _ppos[2] + _rvec[ +2] } ;

        /*-------------- tree-line intersection predicate */
            typedef 
            geom_tree::aabb_pred_line_3 <
                     float , 
                 iptr_type >  tree_pred ; 

        /*-------------- tria-line intersection predicate */
            typedef 
            line_tria_pred <
                 null_pred >    hits_pred ;

            float           _PPOS[3] = {
                    (float) _ppos[0] ,
                    (float) _ppos[1] ,
                    (float) _ppos[2] ,
                    } ;

            float           _RPOS[3] = {
                    (float) _rpos[0] ,
                    (float) _rpos[1] ,
                    (float) _rpos[2] ,
                    } ;

            tree_pred _pred(_PPOS, _RPOS) ;

            null_pred _hfun;
            hits_pred _func(_ppos, _rpos, 
                            *this, _hfun) ;

            this->_tbox.find(_pred,_func) ;

            if (_func._exct) continue ;
            
        /*-------------- done if inside - odd no. crosses */
            if((_func._hnum % 2) != +0)
            {
            return this->_ptag[_pnum] ;
            }
        }
        }
        
        }

        return ( (iptr_type) -1) ;
    }
    
    /*
    --------------------------------------------------------
     * PROJECTOR: project a point on to the geometry.
    --------------------------------------------------------
     */

    __normal_call void_type projector (
        real_type *_ppos ,
        iptr_type  _idim ,
        real_type *_proj
        )
    {
        if (_idim == +1)
        {
    /*------------------ project to closest 1-dim feature */
            real_type  _PROJ[3] = {
           (real_type) _ppos[0] ,
           (real_type) _ppos[1] ,
           (real_type) _ppos[2] } ;

            float      _PPOS[3] = {
           (float)     _ppos[0] ,
           (float)     _ppos[1] ,
           (float)     _ppos[2] } ;

            real_type static const  _RTOL= 
                std::pow (
            std::numeric_limits<real_type>
            ::epsilon(), (real_type)+.75);
      
            real_type  _BTOL =  (
                 this->_bmax[0] - 
                 this->_bmin[0] + 
                 this->_bmax[1] - 
                 this->_bmin[1] + 
                 this->_bmax[2] - 
                 this->_bmin[2] ) * _RTOL; 
            
            near_edge_pred  _func ( _ppos,
                     _PROJ, &this-> _tria,
                     _BTOL) ;
        
            this->_ebox.near(_PPOS, _func) ;

            _proj[0] = _PROJ[0] ;
            _proj[1] = _PROJ[1] ;
            _proj[2] = _PROJ[2] ;
        }
        else
        if (_idim == +2)
        {
    /*------------------ project to closest 2-dim feature */
            real_type  _PROJ[3] = {
           (real_type) _ppos[0] ,
           (real_type) _ppos[1] ,
           (real_type) _ppos[2] } ;

            float      _PPOS[3] = {
           (float)     _ppos[0] ,
           (float)     _ppos[1] ,
           (float)     _ppos[2] } ;

            real_type static const  _RTOL= 
                std::pow (
            std::numeric_limits<real_type>
            ::epsilon(), (real_type)+.75);
      
            real_type  _BTOL =  (
                 this->_bmax[0] - 
                 this->_bmin[0] +
                 this->_bmax[1] - 
                 this->_bmin[1] + 
                 this->_bmax[2] - 
                 this->_bmin[2] ) * _RTOL; 
            
            near_tria_pred  _func ( _ppos,
                     _PROJ, &this-> _tria,
                     _BTOL) ;
        
            this->_tbox.near(_PPOS, _func) ;

            _proj[0] = _PROJ[0] ;
            _proj[1] = _PROJ[1] ;
            _proj[2] = _PROJ[2] ;
        }
        else
        {
    /*------------------ NULL projection -- return inputs */
            _proj[0] = _ppos[0] ;
            _proj[1] = _ppos[1] ;
            _proj[2] = _ppos[2] ;
        }
    }

    } ;
    
    
    }
    
#   endif   //__GEOM_MESH_EUCLIDEAN_3__



