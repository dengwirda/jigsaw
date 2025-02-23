
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
     * Last updated: 21 Apr., 2024
     *
     * Copyright 2013-2024
     * Darren Engwirda
     * d.engwirda@gmail.com
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


    class node_type: public mesh_complex_node_3<I,R>
        {
    /*------------------------------------ loc. node type */
        public  :
        iptr_type                     _itag = +0 ;

        char_type                     _fdim = +0 ;
        char_type                     _feat = +0 ;

        char_type                     _topo [ 2] = {+2} ;

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
        __inline_call char_type      & topo (
            char_type _kind = filt_topo
            )
        {   return ( _kind == filt_topo ) ?
                    this->_topo[ +0] :
                    this->_topo[ +1] ;
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
            char_type _kind = filt_topo
            ) const
        {   return ( _kind == filt_topo ) ?
                    this->_topo[ +0] :
                    this->_topo[ +1] ;
        }

        } ;

    class seed_type: public mesh_complex_node_3<I,R>
        {
    /*------------------------------------ loc. node type */
        public  :
        iptr_type                     _itag = +0 ;

        public  :
    /*------------------------------------ "write" access */
        __inline_call iptr_type&       itag (
            )
        {   return  this->_itag ;
        }
    /*------------------------------------ "const" access */
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }

        } ;

    class edge_type: public mesh_complex_edge_2<I>
        {
    /*------------------------------------ loc. edge type */
        public  :
        iptr_type                     _itag = +0 ;

        char_type                     _feat = +0 ;

        char_type                     _topo [ 2] = {+2} ;

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
        __inline_call char_type      & topo (
            char_type _kind = filt_topo
            )
        {   return ( _kind == filt_topo ) ?
                    this->_topo[ +0] :
                    this->_topo[ +1] ;
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
        __inline_call char_type const& topo (
            char_type _kind = filt_topo
            ) const
        {   return ( _kind == filt_topo ) ?
                    this->_topo[ +0] :
                    this->_topo[ +1] ;
        }

        } ;

    typedef mesh::tria_complex_1<
                node_type,
                edge_type,
                allocator       >   mesh_type ;

    typedef containers::array   <
                seed_type ,
                allocator       >   seed_list ;

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

    seed_list                      _seed ;

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
     * construct geometry from alloc. etc.
    --------------------------------------------------------
     */

    __normal_call geom_mesh_ellipsoid_3d (
        allocator const&
            _asrc = allocator ()
        ) : _seed(    _asrc ) ,
            _mesh(    _asrc ) ,
            _ebox(    _asrc )
        {
        }

    /*
    --------------------------------------------------------
     * NODE-FEAT: calc. node feature type.
    --------------------------------------------------------
     */

    template <
        typename  list_type ,
        typename  user_opts
             >
    __normal_call void_type node_feat (
        iptr_type *_node ,
        list_type &_aset ,
        char_type &_feat ,
        char_type &_topo ,
        char_type &_full ,
        user_opts &_opts
        )
    {
    /*------------ "sharp" geometry//topology about node? */
        real_type constexpr _DtoR =
       (real_type)+3.141592653589793 / 180. ;

        real_type _ZERO = -1. +
            std::pow(std::numeric_limits
                <real_type>::epsilon(),.50) ;

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
        _topo = (char_type)_aset.count () ;
        _full = (char_type)_aset.count () ; // unfiltered

        for (auto _ipos  = _aset.head() ;
                  _ipos != _aset.tend() ;
                ++_ipos  )
        {
        char_type _tbad  = +1 ;
        for (auto _jpos  = _ipos+1 ;
                  _jpos != _aset.tend() ;
                ++_jpos  )
        {
             auto _iedg  = _ipos->_cell ;
             auto _jedg  = _jpos->_cell ;

    /*------------ tag as soft feature if dissimilar tags */
             auto _itag  = 
            this->_mesh.edge(_iedg).itag () ;
             auto _jtag  = 
            this->_mesh.edge(_jedg).itag () ;

            if (_itag != _jtag)
            _feat = std::max(_feat, soft_feat) ;

    /*------------ find signed angle between edge vectors */
            iptr_type _inod[2] = {
            this->_mesh.edge(_iedg).node(0) ,
            this->_mesh.edge(_iedg).node(1) ,
                } ;

            iptr_type _jnod[2] = {
            this->_mesh.edge(_jedg).node(0) ,
            this->_mesh.edge(_jedg).node(1) ,
                } ;

            real_type _ivec[3] ;
            geometry::vector_3d(
               & this->_mesh.
                 node(_inod[0]).pval(0) ,
               & this->_mesh.
                 node(_inod[1]).pval(0) ,
                _ivec) ;

            real_type _jvec[3] ;
            geometry::vector_3d(
               & this->_mesh.
                 node(_jnod[0]).pval(0) ,
               & this->_mesh.
                 node(_jnod[1]).pval(0) ,
                _jvec) ;

            real_type _acos = geometry::
                cosine_3d(_ivec, _jvec) ;

            if (_inod[0] == _jnod[1] ||
                _inod[1] == _jnod[0] )
                _acos *= (real_type)+1. ;
            else
                _acos *= (real_type)-1. ;

            if (_acos >= _ZERO)
            {
    /*------------ tag as "feature" if angle sharp enough */
            if (_acos <= _hard)
            {
                _feat  =
            std::max (_feat, hard_feat) ;
            }
            else
            if (_acos <= _soft)
            {
                _feat  =
            std::max (_feat, soft_feat) ;
            }
            }
            else
            {
            if (_tbad >= 1) _topo -= _tbad-- ;
            }
        }
        }
        {
    /*------------ tag as "feature" if topo. is irregular */
            if (_topo != +  0 )
            if (_topo != +  2 )
                _feat =
            std::max (_feat, soft_feat) ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-FEAT: scan geometry and find features.
    --------------------------------------------------------
     */

    template <
        typename  user_opts
             >
    __normal_call void_type find_feat (
        user_opts &_opts
        )
    {
        typename
            mesh_type::connector _eadj ;

    /*---------------------------------- find sharp feat. */
        for (auto _iter  =
             this->_mesh.node().head() ;
                  _iter !=
             this->_mesh.node().tend() ;
                ++_iter  )
        {
    /*---------------------------------- find sharp 0-dim */
            if (_iter->mark() >= +0 )
            {
            if (_iter->itag() <= -1 ||
                _opts .feat() )
            {
    /*---------------------------------- set geo.-defined */
            _eadj.set_count (0);

            this->_mesh.connect_1(
               &_iter->node (0), POINT_tag, _eadj) ;

            node_feat (
               &_iter->node (0),
                _eadj ,
                _iter->feat () ,
                _iter->topo (0),
                _iter->topo (1),
                _opts ) ;

            if (_iter->itag() <= -1)
            {
    /*---------------------------------- set user-defined */
            _iter->feat () =
                std::max(_iter->feat(), soft_feat) ;
            }
            }
            }
        }

        for (auto _iter  =
             this->_mesh.edge().head() ;
                  _iter !=
             this->_mesh.edge().tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
    /*----------------------------- assign nodes to edges */
            this->_mesh.node()[
                _iter->node(0)].fdim() = 1;
            this->_mesh.node()[
                _iter->node(1)].fdim() = 1;
           }
        }

        for (auto _iter  =
             this->_mesh.node().head() ;
                  _iter !=
             this->_mesh.node().tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
            if (_iter->feat() != null_feat)
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
        typename  user_opts
             >
    __normal_call void_type init_geom (
        user_opts &_opts
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
             this->_mesh.node().head();     // POINT
                  _iter !=
             this->_mesh.node().tend();
                ++_iter  )
        {
            if (_iter->mark() < 0) continue ;

            real_type _apos[2];
            real_type _ppos[3];

            _apos[0] = _iter->pval(0) ;
            _apos[1] = _iter->pval(1) ;

            toR3(_apos, _ppos);

            _iter->pval(0) = _ppos[0] ;
            _iter->pval(1) = _ppos[1] ;
            _iter->pval(2) = _ppos[2] ;
        }

    /*--------------------------- convert to R^3 coord.'s */
        for (auto _iter  =
                   this->_seed. head();     // SEEDS
                  _iter !=
                   this->_seed. tend();
                ++_iter  )
        {
            if (_iter->mark() < 0) continue ;

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

        real_type _rBAR = (real_type)0. ;
        _rBAR += this->_radA ;
        _rBAR += this->_radB ;
        _rBAR += this->_radC ;
        _rBAR  = _rBAR  / (real_type)3. ;

        this->_rEPS *= _rBAR ;

    /*--------------------------- sharp feat. in geometry */
        find_feat (_opts);

    /*--------------------------- init. AABB for arc-seg. */
        containers::
        block_array<tree_item> _bbox;

        iptr_type  _inum = +0;

        for (auto _iter  =
             this->_mesh.edge().head();
                  _iter !=
             this->_mesh.edge().tend();
                ++_iter, ++_inum )
        {
            if (_iter->mark() < 0) continue ;

            iptr_type _enod[ +2] ;
            _enod[0] = _iter->node(0) ;
            _enod[1] = _iter->node(1) ;

            tree_item _tdat ;
            _tdat.ipos() = _inum ;

            make_aabb (*_iter ,
           & this->_mesh.
             node(_enod[ 0 ]).pval(0) ,
           & this->_mesh.
             node(_enod[ 1 ]).pval(0) ,
           &_tdat .pmin( 0 ),
           &_tdat .pmax( 0 )) ;

            _bbox.push_tail(_tdat) ;
        }

        this->_ebox.load (
            _bbox.head() ,
            _bbox.tend() , this->_nbox) ;
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
     * SEED-FEAT: init. "seed" vertex set on feat.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  user_opts
             >
    __normal_call void_type seed_root (
        mesh_type &_rdel ,
        user_opts &_opts
        )
    {
        if (_opts.dims() < +3)
        {
    /*------------------------- point at ellipsoid centre */
        real_type  _ppos[4];
        iptr_type  _inod = -1;
        _ppos[0] = (real_type) +0.0E+0 ;
        _ppos[1] = (real_type) +0.0E+0 ;
        _ppos[2] = (real_type) +0.0E+0 ;
        _ppos[3] = (real_type) +0.0E+0 ;
        _rdel.
        _tria.push_node(_ppos, _inod);
        _rdel.
        _tria.node(_inod)->fdim() = +4 ;
        _rdel.
        _tria.node(_inod)->feat() = +0 ;
        _rdel.
        _tria.node(_inod)->topo() = +0 ;
        _rdel.
        _tria.node(_inod)->part() = +0 ;
        }
    }

    template <
        typename  mesh_type ,
        typename  hfun_type ,
        typename  user_opts
             >
    __normal_call void_type seed_feat (
        mesh_type &_rdel ,
        hfun_type &_hfun ,
        user_opts &_opts
        )
    {
        __unreferenced(_opts) ;
        __unreferenced(_hfun) ;

    /*------------------------- push set of feature nodes */
        for (auto _iter  =
             this->_mesh.node().head() ;
                  _iter !=
             this->_mesh.node().tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0 )
            {
            if (_iter->feat() != null_feat)
            {
    /*----------------------------- push any 'real' feat. */
            real_type  _ppos[4] ;
            _ppos[0] = _iter->pval(0) ;
            _ppos[1] = _iter->pval(1) ;
            _ppos[2] = _iter->pval(2) ;
            _ppos[3] = (real_type)+0. ;

            iptr_type _node = -1 ;
            if (_rdel._tria.push_node (
                    &_ppos[ 0], _node))
            {
                _rdel._tria.node
                    (_node)->fdim()
                        = _iter->fdim() ;

                _rdel._tria.node
                    (_node)->feat()
                        = _iter->feat() ;

                _rdel._tria.node
                    (_node)->topo(0)
                        = _iter->topo(0);

                _rdel._tria.node
                    (_node)->topo(1)
                        = _iter->topo(1);

                _rdel._tria.node
                    (_node)->part()
                        = _iter->itag() ;
            }
            }
            else
            if (_iter->itag() <= -1 )
            {
    /*----------------------------- push any 'user' feat. */
            real_type  _ppos[4] ;
            _ppos[0] = _iter->pval(0) ;
            _ppos[1] = _iter->pval(1) ;
            _ppos[2] = _iter->pval(2) ;
            _ppos[3] = (real_type)+0. ;

            iptr_type _node = -1 ;
            if (_rdel._tria.push_node (
                    &_ppos[ 0], _node))
            {
                _rdel._tria.node
                    (_node)->fdim()
                        = _iter->fdim() ;

                _rdel._tria.node
                    (_node)->feat()
                        = _iter->feat() ;

                _rdel._tria.node
                    (_node)->topo(0)
                        = _iter->topo(0);

                _rdel._tria.node
                    (_node)->topo(1)
                        = _iter->topo(1);

                _rdel._tria.node
                    (_node)->part()
                        = _iter->itag() ;
            }
            }
            }
        }
    }

    /*
    --------------------------------------------------------
     * SEED-BASE: init. icosahedral "seed" vertex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call void_type seed_base (
        mesh_type &_rdel
        )
    {
    /*--------------------------- init. reg.-icosahedron */
        real_type  _pi =
       (real_type)std::atan(1.0) * 4. ;

        real_type  _la =
       (real_type)std::atan(0.5) * 1. ;

        real_type  _lo = 2.*_pi / 10. ;

        real_type  _ppos[4];
        iptr_type  _inod = -1;
        _ppos[0] = this->_radA *
        std::cos(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)+0.5) ;
        _ppos[1] = this->_radB *
        std::sin(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)+0.5) ;
        _ppos[2] = this->_radC *
        std::sin(_pi*(real_type)+0.5) ;
        _ppos[3] =   (real_type)+0.0;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)-0.5) ;
        _ppos[1] = this->_radB *
        std::sin(_pi*(real_type)+0.0) *
        std::cos(_pi*(real_type)-0.5) ;
        _ppos[2] = this->_radC *
        std::sin(_pi*(real_type)-0.5) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+0.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+0.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)+1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+1.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+1.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)-1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+2.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+2.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)+1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+3.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+3.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)-1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+4.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+4.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)+1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+5.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+5.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)-1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+6.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+6.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)+1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+7.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+7.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)-1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+8.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+8.0) *
        std::cos(_la*(real_type)+1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)+1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;

        _ppos[0] = this->_radA *
        std::cos(_lo*(real_type)+9.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[1] = this->_radB *
        std::sin(_lo*(real_type)+9.0) *
        std::cos(_la*(real_type)-1.0) ;
        _ppos[2] = this->_radC *
        std::sin(_la*(real_type)-1.0) ;
        _ppos[3] =   (real_type)+0.0;
        _inod = -1 ;
        _rdel.
        _tria.push_node(_ppos, _inod) ;
        _rdel.
        _tria.node(_inod)->fdim() = +0;
        _rdel.
        _tria.node(_inod)->feat() = +0;
        _rdel.
        _tria.node(_inod)->topo() = +2;
        _rdel.
        _tria.node(_inod)->part() = +0;
    }

    /*
    --------------------------------------------------------
     * SEED-MESH: init. "seed" vertex set on geom.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  hfun_type ,
        typename  user_opts
             >
    __normal_call void_type seed_mesh (
        mesh_type &_rdel ,
        hfun_type &_hfun ,
        user_opts &_opts
        )
    {
    /*------------------------- eval. h(x) func. on nodes */  
        containers::array<real_type> _spac(
             this->_mesh.node().count(), 0) ;
        iptr_type _npos  = 0 ;
        for (auto _node  =
             this->_mesh.node().head() ;
                  _node !=
             this->_mesh.node().tend() ;
                ++_node, ++_npos)
        {
            _spac[_npos] = 
                _hfun.eval(&_node->pval(+0));
        }
    
    /*------------------------- well-distributed sampling */
        while (_rdel._tria._nset.count()
                < (std::size_t)_opts.seed() + 5)
        {
            typename geom_type::
                     mesh_type::
                node_list::_write_it _best ;

            real_type _dmax  = (real_type) +.0 ;
            iptr_type _inum  = 0 ;
            for (auto _ipos  =
                 this->_mesh.node().head() ;
                      _ipos !=
                 this->_mesh.node().tend() ;
                    ++_ipos, ++_inum)
            {
    /*------------------------- get current furthest node */
                if (_ipos->mark() >= 0)
                {
                real_type _dmin  =
                   +std::numeric_limits
                        <real_type>::infinity();

                iptr_type _jnum  = 0 ;
                for (auto _jpos  =
                    _rdel._tria._nset.head() ;
                          _jpos !=
                    _rdel._tria._nset.tend() ;
                        ++_jpos, ++_jnum)
                {
                    real_type _dist =
                        geometry::lensqr_3d(
                        &_ipos->pval(+0),
                        &_jpos->pval(+0)) ;
                        
                    real_type _isiz = _spac[_inum];
                    real_type _jsiz = _spac[_jnum];
                        
                    _dist/= _isiz + _jsiz ;

                    _dmin = std::min(_dmin, _dist);
                }

                if (_dmax < _dmin)
                {
                    _dmax = _dmin;
                    _best = _ipos;
                }
                }
            }

            if (_dmax > (real_type)0.)
            {
    /*------------------------- add current furthest node */
                real_type _ppos[4] ;
                _ppos[0] = _best->pval(0) ;
                _ppos[1] = _best->pval(1) ;
                _ppos[2] = _best->pval(2) ;
                _ppos[3] = (real_type)+0. ;

                iptr_type _node = -1;
                if (_rdel._tria.push_node(
                        &_ppos[ 0], _node) )
                {
                    _rdel._tria.node
                        (_node)->fdim()
                            = _best->fdim() ;

                    _rdel._tria.node
                        (_node)->feat()
                            = _best->feat() ;

                    _rdel._tria.node
                        (_node)->topo(0)
                            = _best->topo(0);

                    _rdel._tria.node
                        (_node)->topo(1)
                            = _best->topo(1);

                    _rdel._tria.node
                        (_node)->part()
                            = _best->itag() ;
                }
            }
            else break  ;
        }

        if (_rdel._tria._nset.count() <= +8 )
        {
            seed_base(_rdel) ;
        }
    }

    /*
    --------------------------------------------------------
     * HELPERS: predicates for intersection tests.
    --------------------------------------------------------
     */

    __normal_call void_type make_aabb (
        edge_type  _edge,
        real_type *_apos,
        real_type *_bpos,
        float     *_rmin,
        float     *_rmax
        )
    {
    /*- build an AABB that encloses a spheroidal arc-seg. */
        real_type _rTOL = this->_rEPS ;

        __unreferenced(_edge);

        _rmin[0] = std::min(
        (float)_apos[0], (float)_bpos[0]) ;
        _rmin[1] = std::min(
        (float)_apos[1], (float)_bpos[1]) ;
        _rmin[2] = std::min(
        (float)_apos[2], (float)_bpos[2]) ;

        _rmax[0] = std::max(
        (float)_apos[0], (float)_bpos[0]) ;
        _rmax[1] = std::max(
        (float)_apos[1], (float)_bpos[1]) ;
        _rmax[2] = std::max(
        (float)_apos[2], (float)_bpos[2]) ;

        real_type _rlen =
        geometry::length_3d(_apos, _bpos) ;

        real_type  _rmid[3] = {
        _apos[0] * (real_type)+.5 +
        _bpos[0] * (real_type)+.5 ,
        _apos[1] * (real_type)+.5 +
        _bpos[1] * (real_type)+.5 ,
        _apos[2] * (real_type)+.5 +
        _bpos[2] * (real_type)+.5 ,
            } ;

        _rmin[0] = std::min(
        _rmin[0], (float)(_rmid[0]-_rlen));
        _rmin[1] = std::min(
        _rmin[1], (float)(_rmid[1]-_rlen));
        _rmin[2] = std::min(
        _rmin[2], (float)(_rmid[2]-_rlen));

        _rmax[0] = std::max(
        _rmax[0], (float)(_rmid[0]+_rlen));
        _rmax[1] = std::max(
        _rmax[1], (float)(_rmid[1]+_rlen));
        _rmax[2] = std::max(
        _rmax[2], (float)(_rmid[2]+_rlen));

        _rmin[0] = _rmin[0]-(float)_rTOL;
        _rmin[1] = _rmin[1]-(float)_rTOL;
        _rmin[2] = _rmin[2]-(float)_rTOL;

        _rmax[0] = _rmax[0]+(float)_rTOL;
        _rmax[1] = _rmax[1]+(float)_rTOL;
        _rmax[2] = _rmax[2]+(float)_rTOL;
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

        _pprj[0] =_psrc[0] ;
        _pprj[1] =_psrc[1] ;
        _pprj[2] =_psrc[2] ;

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

        if (_ttaa > (real_type) -1.)
        {
        _pprj[0] =
        _pmid[0] + _ttaa*_pdel[0] ;
        _pprj[1] =
        _pmid[1] + _ttaa*_pdel[1] ;
        _pprj[2] =
        _pmid[2] + _ttaa*_pdel[2] ;
        }
        else
        if (_ttbb > (real_type) -1.)
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
                   &_edge.topo(0),
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
                   &_edge.topo(0),
                    _edge.itag() )  ;

                _hfun(&_qprj[0], _htmp ,
                    _edge.feat() ,
                   &_edge.topo(0),
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

            real_type _xtmp[3] ;
            geometry::proj_flat_3d(
                _xprj ,
                _flat._ppos,
                _flat._nvec, _xtmp) ;

            real_type _xlen =
                geometry::
            lensqr_3d(_xtmp, _xprj) ;

            if (_xlen < this->_rEPS *
                        this->_rEPS )
            {
                _hfun(&_xprj[0], _htmp ,
                    _edge.feat() ,
                   &_edge.topo(0),
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

            iptr_type  _enod[2] ;
            _enod[0] = _geom.
                _mesh.edge(_epos).node(0) ;
            _enod[1] = _geom.
                _mesh.edge(_epos).node(1) ;

        /*--------------- call output function on hit */
            _geom .flat_kern  (
                 this->_flat,
                 this->_geom.
            _mesh.edge(_epos) ,
                &this->_geom.
            _mesh.node(_enod[ 0]).pval(0) ,
                &this->_geom.
            _mesh.node(_enod[ 1]).pval(0) ,
                    + 0 ,
                 this->_hfun,
                 this->_hnum) ;

            this->_find =
            this->_find | (this->_hnum != 0) ;
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

            iptr_type  _enod[2] ;
            _enod[0] = _geom.
                _mesh.edge(_epos).node(0) ;
            _enod[1] = _geom.
                _mesh.edge(_epos).node(1) ;

        /*--------------- call output function on hit */
            _geom. ball_test  (
                 this->_ball,
                 this->_geom.
            _mesh.edge(_epos) ,
                &this->_geom.
            _mesh.node(_enod[ 0]).pval(0) ,
                &this->_geom.
            _mesh.node(_enod[ 1]).pval(0) ,
                 this->_hfun,
                 this->_hnum) ;

            this->_find =
            this->_find | (this->_hnum != 0) ;
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

        real_type _llen = +.5 *
        geometry::length_3d(_flat._nvec) ;

        float           _RMIN[3] = {
                (float)(_flat. _ppos[0] -
                               _llen) ,
                (float)(_flat. _ppos[1] -
                               _llen) ,
                (float)(_flat. _ppos[2] -
                               _llen) ,
                } ;

        float           _RMAX[3] = {
                (float)(_flat. _ppos[0] +
                               _llen) ,
                (float)(_flat. _ppos[1] +
                               _llen) ,
                (float)(_flat. _ppos[2] +
                               _llen) ,
                } ;

        for ( auto
             _iter  = _flat._bnds.head() ;
             _iter != _flat._bnds.tend() ;
           ++_iter  )
        {
            real_type _xprj[3] ;
            projector(*_iter, +2, _xprj) ;

            _RMIN[0] = std::min(
            _RMIN[0] , (float) _xprj[0]) ;
            _RMIN[1] = std::min(
            _RMIN[1] , (float) _xprj[1]) ;
            _RMIN[2] = std::min(
            _RMIN[2] , (float) _xprj[2]) ;

            _RMAX[0] = std::max(
            _RMAX[0] , (float) _xprj[0]) ;
            _RMAX[1] = std::max(
            _RMAX[1] , (float) _xprj[1]) ;
            _RMAX[2] = std::max(
            _RMAX[2] , (float) _xprj[2]) ;
        }

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

        float     static const _RTOL =
            std::pow (
            std::numeric_limits<float>
            ::epsilon(), (float).675) ;

        float           _PMID[3] = {
                (float) _ball. _pmid[0] ,
                (float) _ball. _pmid[1] ,
                (float) _ball. _pmid[2] ,
                } ;

        float           _RRAD  =
                (float) _ball. _rrad;

        _RRAD *= (1.f + _RTOL) ;

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

    /*------------------ compute line-ellipsoid intersect */
        real_type  _ipos[3] ;
        _ipos[0] = _line._ipos[ 0];
        _ipos[1] = _line._ipos[ 1];
        _ipos[2] = _line._ipos[ 2];

        real_type  _jpos[3] ;
        _jpos[0] = _line._jpos[ 0];
        _jpos[1] = _line._jpos[ 1];
        _jpos[2] = _line._jpos[ 2];

        real_type _ttaa, _ttbb;
        if(line_surf(_ipos, _jpos, _ttaa, _ttbb))
        {
        dd_flt     _pmid[3] ;
        _pmid[0].
            from_add(_jpos[0], _ipos[0]);
        _pmid[0] *= +0.5 ;
        _pmid[1].
            from_add(_jpos[1], _ipos[1]);
        _pmid[1] *= +0.5 ;
        _pmid[2].
            from_add(_jpos[2], _ipos[2]);
        _pmid[2] *= +0.5 ;

        dd_flt     _pdel[3] ;
        _pdel[0].
            from_sub(_jpos[0], _ipos[0]);
        _pdel[0] *= +0.5 ;
        _pdel[1].
            from_sub(_jpos[1], _ipos[1]);
        _pdel[1] *= +0.5 ;
        _pdel[2].
            from_sub(_jpos[2], _ipos[2]);
        _pdel[2] *= +0.5 ;

        if (_ttaa >= (real_type)-1.)
        if (_ttaa <= (real_type)+1.)
        {
    /*------------------ push surf.-proj. to output func. */
            dd_flt    _APOS[3] = {
            _pmid[0] + _ttaa * _pdel[0] ,
            _pmid[1] + _ttaa * _pdel[1] ,
            _pmid[2] + _ttaa * _pdel[2]
                } ;

            char_type _hits =
                geometry::face_hits ;
            char_type _feat = +2;
            iptr_type _itag = +0;
            char_type _topo[2] = {2};

            real_type  _apos[3] ;
            _apos[0] = _APOS[0] ;
            _apos[1] = _APOS[1] ;
            _apos[2] = _APOS[2] ;

            _find  =  true  ;
            _hfun( _apos,
            _hits, _feat, _topo, _itag) ;
        }

        if (_ttbb >= (real_type)-1.)
        if (_ttbb <= (real_type)+1.)
        {
    /*------------------ push surf.-proj. to output func. */
            dd_flt    _BPOS[3] = {
            _pmid[0] + _ttbb * _pdel[0] ,
            _pmid[1] + _ttbb * _pdel[1] ,
            _pmid[2] + _ttbb * _pdel[2]
                } ;

            char_type _hits =
                geometry::face_hits ;
            char_type _feat = +2;
            iptr_type _itag = +0;
            char_type _topo[2] = {2};

            real_type  _bpos[3] ;
            _bpos[0] = _BPOS[0] ;
            _bpos[1] = _BPOS[1] ;
            _bpos[2] = _BPOS[2] ;

            _find  =  true  ;
            _hfun( _bpos,
            _hits, _feat, _topo, _itag) ;
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
        bool_type _find = true  ;

    /*------------------ calc. initial dir. to surf.-ball */
        real_type _proj[3], _vdir[3] ;
        proj_surf(_sbal, _proj) ;

        iptr_type _iter;
        for(_iter = 8; _iter-- != 0; )
        {
    /*------------------ iter. to improve dir.-to-surface */
        geometry::vector_3d(
            _disc._pmid, _proj, _vdir) ;

        geometry::normalise_3d( _vdir) ;

        real_type _xpos[3] = {
        _disc._pmid[0] +
            _disc._rrad *_vdir[0] ,
        _disc._pmid[1] +
            _disc._rrad *_vdir[1] ,
        _disc._pmid[2] +
            _disc._rrad *_vdir[2]
            } ;

        proj_surf(_xpos, _proj) ;

        real_type _dsqr =
            geometry::
        lensqr_3d(_xpos, _proj) ;

        if (_dsqr < this->_rEPS *
                    this->_rEPS )
        {
    /*------------------ push surf.-proj. to output func. */
            char_type _hits =
                geometry::face_hits ;
            char_type _feat = +2;
            iptr_type _itag = +0;
            char_type _topo[2] = {2};

            _hfun( _proj ,
            _hits, _feat , _topo, _itag) ;
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
    /*------------------ (x/a)^2 + (y/b)^2 + (z/c)^2 < 1? */
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
        if (_idim == +2)
        {
    /*------------------ project to closest 2-dim feature */
            proj_surf(_ppos, _proj) ;
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

#   endif   //__GEOM_MESH_ELLIPSOID_3__



