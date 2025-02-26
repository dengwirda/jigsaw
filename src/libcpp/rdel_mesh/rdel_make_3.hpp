
    /*
    --------------------------------------------------------
     * RDEL-MAKE-3: restricted delaunay tria. in R^3.
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
     * Last updated: 12 Nov., 2024
     *
     * Copyright 2013-2024
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     *
     * This class defines the "restricted" delaunay
     * tessellation algorithm for domains in R^3. A
     * bounding DT is built for the points, and the rDT
     * constructed by evaluating the dual predicate for
     * each node, edge, face and cell in the tessellation.
     *
     * My implementation is described here:
     *
     * D. Engwirda, (2014): "Locally-optimal Delaunay-
     * refinement and optimisation-based mesh generation",
     * Ph.D. Thesis, School of Mathematics and Statistics,
     * Univ. of Sydney.
     * http://hdl.handle.net/2123/13148
     *
     * building on various previous works on rDT methods,
     * including (primarily):
     *
     * J.D. Boissonnat, S. Oudot, (2005): "Provably Good
     * Sampling and Meshing of Surfaces", Graphical Models,
     * 67, pp. 405-451,
     * https://doi.org/10.1016/j.gmod.2005.01.004
     *
     * C. Jamin, P. Alliez, M. Yvinec, and J.D. Boissonnat,
     * (2015): "CGALmesh: a generic framework for Delaunay
     * mesh generation", ACM Transactions on Mathematical
     * Software (TOMS), 41, pp. 23
     * https://doi.org/10.1145/2699463
     *
     * L. Rineau, M. Yvinec, (2008): "Meshing 3D Domains
     * Bounded by Piecewise Smooth Surfaces", Proc. of the
     * 16th International Meshing Roundtable, pp. 443-460,
     * https://doi.org/10.1007/978-3-540-75103-8_25
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_MAKE_3__
#   define __RDEL_MAKE_3__

    namespace mesh {

    template <
    typename M ,
    typename G ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class rdel_make_3d
    {
    public  :

    /*----------- restricted delaunay tessellation in R^3 */

    typedef M                           mesh_type ;
    typedef G                           geom_type ;
    typedef I                           init_type ;
    typedef A                           allocator ;

    typedef typename
            mesh_type::real_type        real_type ;
    typedef typename
            mesh_type::iptr_type        iptr_type ;

    typedef typename
            allocator::size_type        uint_type ;

    typedef typename
            mesh_type::node_data        node_data ;
    typedef typename
            mesh_type::edge_data        edge_data ;
    typedef typename
            mesh_type::face_data        face_data ;
    typedef typename
            mesh_type::tria_data        tria_data ;

    typedef mesh::rdel_timers           rdel_stat ;

    typedef containers::array       <
                iptr_type           >   iptr_list ;

    typedef containers::array       <
                edge_data           >   edat_list ;
    typedef containers::array       <
                face_data           >   fdat_list ;
    typedef containers::array       <
                tria_data           >   tdat_list ;

    typedef mesh::rdel_pred_base_3  <
                geom_type,
                mesh_type           >   rdel_pred ;

    typedef mesh::mesh_params       <
                real_type,
                iptr_type           >   rdel_opts ;

    /*
    --------------------------------------------------------
     * PUSH-EDGE: add new edge to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_edge (
        mesh_type &_mesh ,
        geom_type &_geom ,
        iptr_type  _tpos ,
        typename
    mesh_type::edge_list & _edge_test ,
        iptr_type &_nedg ,
        iptr_type  _pass ,
        rdel_opts &_opts
        )
    {
    /*-------------------------------- correct node dims? */
        iptr_type _fdim =+0;
        for (auto _node =+4; _node-- != +0; )
        {
        if (_mesh._tria.node (
            _mesh._tria.tria (
            _tpos)->node(_node))->fdim() <= +1)
            _fdim += +1 ;
        }
    /*-------------------------------- quick break if not */
        if (_fdim  < +2 ) return ;

    /*-------------------------------- check "restricted" */
        for (auto _fpos =+6; _fpos-- != +0; )
        {
        /*---------------------------- extract edge nodes */
            iptr_type _tnod[ +4] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_tnod, _fpos, +3, +1) ;
            _tnod[0] = _mesh._tria.
            tria(_tpos)->node(_tnod[0]);
            _tnod[1] = _mesh._tria.
            tria(_tpos)->node(_tnod[1]);

            if (_mesh._tria.node(
                _tnod[0])->fdim() < 2 &&
                _mesh._tria.node(
                _tnod[1])->fdim() < 2 )
            {
        /*--------------- face contains correct k^d nodes */
            algorithms::isort (
                &_tnod[0], &_tnod[2],
                    std::less<iptr_type>()) ;

            edge_data _edat;
            _edat._node[0] = _tnod[ 0] ;
            _edat._node[1] = _tnod[ 1] ;

            bool_type _test = true;
            typename mesh_type::
                     edge_list::
                item_type *_mptr = nullptr  ;
        //  pragma omp critical(update_face_edg2)
            if(_edge_test.
                find( _edat, _mptr) )
            {
        /*--------------------------- don't test repeats! */
                _test = false;
            }

            _edat._pass    = _pass;

            _edat._tadj    = _tpos;
            _edat._eadj    =
                (char_type)  _fpos;

            if (_test)
            {
        /*--------------------------- call edge predicate */
            char_type _topo[ 2], _feat;
            real_type _ebal[ 4];
            real_type _sbal[ 4];

            __unreferenced(_opts);

            bool_type _rBND    =
                rdel_pred::edge_ball (
                    _geom,_mesh,
                    _edat._tadj,
                    _edat._eadj,
                    _ebal,_sbal,
                    _feat,_topo,
                    _edat._part) ;
            
        /*--------------------------- push edge onto mesh */
        //  pragma omp critical(update_face_edg2)
            {
            if (_rBND) _nedg += +1 ;

            if (_rBND)
            _mesh.push_edge (_edat) ;

            _edge_test.push (_edat) ;
            }
            }

            }
        } // for (auto _fpos = +6; _fpos-- != +0; )
    }

    /*
    --------------------------------------------------------
     * PUSH-FACE: add new face to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_face (
        mesh_type &_mesh ,
        geom_type &_geom ,
        iptr_type  _tpos ,
        typename
    mesh_type::face_list & _face_test ,
        iptr_type &_nfac ,
        iptr_type &_ndup ,
        iptr_type  _pass ,
        rdel_opts &_opts
        )
    {
    /*-------------------------------- correct node dims? */
        iptr_type _fdim =+0;
        for (auto _node =+4; _node-- != +0; )
        {
        if (_mesh._tria.node (
            _mesh._tria.tria (
            _tpos)->node(_node))->fdim() <= +2)
            _fdim += +1 ;
        }
    /*-------------------------------- quick break if not */
        if (_fdim  < +3 ) return ;

    /*-------------------------------- check "restricted" */
        for (auto _fpos =+4; _fpos-- != +0; )
        {
        /*---------------------------- extract face nodes */
            iptr_type _tnod[ +4] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_tnod, _fpos, +3, +2) ;
            _tnod[0] = _mesh._tria.
            tria(_tpos)->node(_tnod[0]);
            _tnod[1] = _mesh._tria.
            tria(_tpos)->node(_tnod[1]);
            _tnod[2] = _mesh._tria.
            tria(_tpos)->node(_tnod[2]);

            if (_mesh._tria.node(
                _tnod[0])->fdim() < 3 &&
                _mesh._tria.node(
                _tnod[1])->fdim() < 3 &&
                _mesh._tria.node(
                _tnod[2])->fdim() < 3 )
            {
        /*--------------- face contains correct k^d nodes */
            algorithms::isort (
                &_tnod[0], &_tnod[3],
                    std::less<iptr_type>()) ;

            face_data _fdat;
            _fdat._node[0] = _tnod[ 0] ;
            _fdat._node[1] = _tnod[ 1] ;
            _fdat._node[2] = _tnod[ 2] ;

            bool_type _test = true;
            typename mesh_type::
                     face_list::
                item_type *_mptr = nullptr  ;
        //  pragma omp critical(update_face_tri3)
            if(_face_test.
                find( _fdat, _mptr) )
            {
        /*--------------------------- count bnd. repeats! */
                _ndup +=
                _mptr->_data._dups;

        /*--------------------------- don't test repeats! */
                _test = false ;
            }

            _fdat._pass    = _pass;

            _fdat._tadj    = _tpos;
            _fdat._fadj    =
                (char_type)  _fpos;
            _fdat._dups    = +0; // count num. dup's
                                 // only in hash-set

            if (_test)
            {
        /*--------------------------- call face predicate */
            char_type _topo[ 2], _feat;
            real_type _fbal[ 4];
            real_type _sbal[ 4];

            __unreferenced(_opts);

            bool_type _rBND    =
                rdel_pred::face_ball (
                    _geom,_mesh,
                    _fdat._tadj,
                    _fdat._fadj,
                    _fbal,_sbal,
                    _feat,_topo,
                    _fdat._part) ;
            
        /*--------------------------- push face onto mesh */
        //  pragma omp critical(update_face_tri3)
            {
            if (_rBND) _nfac += +1 ;

            if (_rBND)
                _fdat._dups   = +1 ;

            if (_rBND)
            _mesh.push_face(_fdat) ;

            _face_test.push(_fdat) ;
            }
            }

            }
        } // for (auto _fpos = +4; _fpos-- != +0; )
    }

    /*
    --------------------------------------------------------
     * PUSH-TRIA: add new tria to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_tria (
        mesh_type &_mesh ,
        geom_type &_geom ,
        iptr_type  _tpos ,
        iptr_type &_sign ,
        /* typename
    mesh_type::tria_list & _tria_test , */
        iptr_type &_ntri ,
        iptr_type  _pass ,
        rdel_opts &_opts
        )
    {
        {
        /*---------------------------- check "restricted" */
            iptr_type  _tnod[ +4] ;
            _tnod[0] = _mesh.
            _tria.tria(_tpos)->node(0);
            _tnod[1] = _mesh.
            _tria.tria(_tpos)->node(1);
            _tnod[2] = _mesh.
            _tria.tria(_tpos)->node(2);
            _tnod[3] = _mesh.
            _tria.tria(_tpos)->node(3);

            if (_mesh._tria.node(
                _tnod[0])->fdim() < 4 &&
                _mesh._tria.node(
                _tnod[1])->fdim() < 4 &&
                _mesh._tria.node(
                _tnod[2])->fdim() < 4 &&
                _mesh._tria.node(
                _tnod[3])->fdim() < 4 )
            {
        /*--------------- face contains correct k^d nodes */
            tria_data _tdat;
            _tdat._node[0] = _tnod[ 0] ;
            _tdat._node[1] = _tnod[ 1] ;
            _tdat._node[2] = _tnod[ 2] ;
            _tdat._node[3] = _tnod[ 3] ;

            _tdat._tadj =  _tpos ;

        //  _tria_test.push( _tdat) ; guarantee no repeats!

        /*--------------------------- calc tria cost/kind */
            _tdat._part =  _sign ;
            _tdat._pass =  _pass ;

            __unreferenced(_opts);

            real_type _tbal[ +3] ;
            bool_type _rBND   =
                rdel_pred::tria_ball (
                    _geom,_mesh,
                    _tdat._tadj,
                    _tbal,
                    _tdat._part) ;

            _sign = _tdat. _part ;

        /*--------------------------- push tria onto mesh */
        //  pragma omp critical(update_face_tri4)
            {
            if (_rBND) _ntri += +1 ;

            if (_rBND)
            _mesh.push_tria(_tdat) ;
            }

            }
        }
    }

    /*
    --------------------------------------------------------
     * PUSH-TRIA: add new tria to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_tria (
        mesh_type &_mesh ,
        iptr_type  _tpos ,
        rdel_opts &_opts
        )
    {
    /*-------------------------------- check "restricted" */
        {
            tria_data _tdat;
            _tdat._node[0] = _mesh.
            _tria.tria(_tpos)->node(0) ;
            _tdat._node[1] = _mesh.
            _tria.tria(_tpos)->node(1) ;
            _tdat._node[2] = _mesh.
            _tria.tria(_tpos)->node(2) ;
            _tdat._node[3] = _mesh.
            _tria.tria(_tpos)->node(3) ;

            if (_mesh._tria.node(
            _tdat._node[0])->fdim() < 4 &&
                _mesh._tria.node(
            _tdat._node[1])->fdim() < 4 &&
                _mesh._tria.node(
            _tdat._node[2])->fdim() < 4 &&
                _mesh._tria.node(
            _tdat._node[3])->fdim() < 4 )
            {
        /*--------------------------- call tria predicate */
            _tdat._tadj    = _tpos ;            
            _tdat._part =    +0  ;
            _tdat._pass =    +0  ;

            __unreferenced (_opts) ;

        /*--------------------------- push tria onto mesh */
            _mesh.push_tria(_tdat) ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * TRIA-CIRC: calc. circumball for tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type tria_circ ( // get "exact" tbal?
        mesh_type &_mesh ,
        iptr_type  _tpos
        )
    {
        iptr_type _tnod[4] = {
        _mesh.
        _tria.tria( _tpos)->node(0) ,
        _mesh.
        _tria.tria( _tpos)->node(1) ,
        _mesh.
        _tria.tria( _tpos)->node(2) ,
        _mesh.
        _tria.tria( _tpos)->node(3)
            } ;

        algorithms::isort(
            &_tnod[0], &_tnod[4],
                std::less<iptr_type>()) ;

    /*---------------------- calc. ball in floating-point */
        real_type _tbal[4] ;
        geometry::perp_ball_3d (
            _tbal ,
       &_mesh._tria.
            node(_tnod[0])->pval(0) ,
       &_mesh._tria.
            node(_tnod[1])->pval(0) ,
       &_mesh._tria.
            node(_tnod[2])->pval(0) ,
       &_mesh._tria.
            node(_tnod[3])->pval(0)
            ) ;

        _mesh._tria.tria(
        _tpos)->circ(0) =  _tbal[0] ;
        _mesh._tria.tria(
        _tpos)->circ(1) =  _tbal[1] ;
        _mesh._tria.tria(
        _tpos)->circ(2) =  _tbal[2] ;
    }

    /*
    --------------------------------------------------------
     * INIT-MESH: init. the bounding DT.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_sort (
        init_type &_init,
        iptr_list &_iset
        )
    {
        typedef geom_tree::aabb_node_base_k
                           tree_node ;

        typedef geom_tree::
                    aabb_item_node_k <
            real_type,
            iptr_type, 3>  tree_item ;

        typedef geom_tree::aabb_tree <
            tree_item, 3,
            tree_node,
            allocator   >  tree_type ;

        containers::array<tree_item> _bbox;

    /*------------------------------ initialise aabb-tree */
        iptr_type _npos  = 0 ;
        tree_type _tree  ;
        for (auto _node  =
            _init._mesh.node().head();
                  _node !=
            _init._mesh.node().tend();
                ++_node, ++_npos)
        {
            if (_node->mark() >= +0)
            {
            _bbox.push_tail() ;
            _bbox.tail()->
               pval(0) = _node->pval(0) ;
            _bbox.tail()->
               pval(1) = _node->pval(1) ;
            _bbox.tail()->
               pval(2) = _node->pval(2) ;

            _bbox.tail()->
               ipos () = _npos ;
            }
        }

        iptr_type _NBOX =
            (iptr_type) std::pow (8, 3) ;       // 8^ndim

        _tree.load(_bbox.head(),
                   _bbox.tend(), _NBOX) ;

    /*------------------------------ randomised tree sort */
        _tree.brio(_iset) ;
    }

    __static_call
    __normal_call void_type init_init (
        init_type &_init,
        geom_type &_geom,
        mesh_type &_mesh
        )
    {
    /*------------------------------ form insertion order */
        iptr_type _hint  = -1 ;
        iptr_list _iset  ;
        init_sort(_init, _iset) ;

    /*------------------------------ find "central" point */
        iptr_type _imid  = -1 ;
        real_type _dmin  =
            std::numeric_limits
                <real_type>::infinity();

        real_type _pmid[3] ;
        _pmid[0] = (real_type) +0. ;
        _pmid[1] = (real_type) +0. ;
        _pmid[2] = (real_type) +0. ;

        for (auto _iter  = _iset.head();
                  _iter != _iset.tend();
                ++_iter  )
        {
             auto _node =
           &_init._mesh. node (*_iter) ;

            _pmid[0] +=
                _node->pval(0);
            _pmid[1] +=
                _node->pval(1);
            _pmid[2] +=
                _node->pval(2);
        }

        _pmid[0] /= _iset.count () ;
        _pmid[1] /= _iset.count () ;
        _pmid[2] /= _iset.count () ;

        for (auto _iter  = _iset.head();
                  _iter != _iset.tend();
                ++_iter  )
        {
             auto _node =
           &_init._mesh. node (*_iter) ;

            real_type _dsqr =
            geometry::lensqr_3d(
               &_node->pval(0), _pmid) ;

            if (_dsqr < _dmin)
            {
                _dmin = _dsqr ;

                _imid =*_iter ;
            }
        }

    /*------------------------------ seed mesh from init. */
        iptr_type _base =
            (iptr_type)
            _mesh._tria._nset.count () ;

        if (_imid > -1)
        {
             auto _node =
           &_init._mesh. node ( _imid) ;

            if (_node->itag () >= +0 )
                _geom.projector(
               &_node->pval(0) ,
                _node->itag () ,
               &_node->pval(0) ) ;

            iptr_type _npos =     // undo BRIO reindexing
                _base + _imid - 0;

            if (_mesh._tria.push_node(
               &_node->pval(0) ,
                _npos, _hint ) )
            {
            _mesh._tria.node
                (_npos)->fdim() = 0 ;

            _mesh._tria.node
                (_npos)->feat()
                    = _node->feat() ;

            _mesh._tria.node
                (_npos)->topo() = 2 ;

            _mesh._tria.node
                (_npos)->part()
                    = _node->itag() ;

            _hint = _mesh._tria.
                node(_npos)->next() ;
            }
        }

    /*------------------------------ seed mesh from init. */
        for (auto _iter  = _iset.head();
                  _iter != _iset.tend();
                ++_iter  )
        {
            if (*_iter == _imid) continue ;

             auto _node =
           &_init._mesh. node (*_iter) ;

            if (_node->itag () >= +0 )
                _geom.projector(
               &_node->pval(0) ,
                _node->itag () ,
               &_node->pval(0) ) ;

            iptr_type _npos =     // undo BRIO reindexing
                _base + *_iter -  0 ;

            if (_mesh._tria.push_node(
               &_node->pval(0) ,
                _npos, _hint ) )
            {
            _mesh._tria.node
                (_npos)->fdim() = 0 ;

            _mesh._tria.node
                (_npos)->feat()
                    = _node->feat() ;

            _mesh._tria.node
                (_npos)->topo() = 2 ;

            _mesh._tria.node
                (_npos)->part()
                    = _node->itag() ;

            _hint = _mesh._tria.
                node(_npos)->next() ;
            }
        }
    }

    __static_call
    __normal_call void_type init_mesh (
        geom_type &_geom,
        init_type &_init,
        mesh_type &_mesh,
        rdel_opts &_opts
        )
    {
        __unreferenced(_opts) ;

    /*------------------------------ initialise mesh bbox */
        real_type _pmin[ 3] ;
        real_type _pmax[ 3] ;
        _pmin[ 0] = _geom._bmin[ 0] ;
        _pmin[ 1] = _geom._bmin[ 1] ;
        _pmin[ 2] = _geom._bmin[ 2] ;
        _pmax[ 0] = _geom._bmax[ 0] ;
        _pmax[ 1] = _geom._bmax[ 1] ;
        _pmax[ 2] = _geom._bmax[ 2] ;

        for (auto _node  =
            _init._mesh.node().head() ;
                  _node !=
            _init._mesh.node().tend() ;
                ++_node  )
        {
        if (_node->mark() >= + 0 )
        {
        _pmin[ 0] = std::min(
        _pmin[ 0], _node->pval(0)) ;
        _pmax[ 0] = std::max(
        _pmax[ 0], _node->pval(0)) ;

        _pmin[ 1] = std::min(
        _pmin[ 1], _node->pval(1)) ;
        _pmax[ 1] = std::max(
        _pmax[ 1], _node->pval(1)) ;

        _pmin[ 2] = std::min(
        _pmin[ 2], _node->pval(2)) ;
        _pmax[ 2] = std::max(
        _pmax[ 2], _node->pval(2)) ;
        }
        }

        real_type _plen[ 3] = {
        _pmax[ 0] - _pmin[ 0] ,
        _pmax[ 1] - _pmin[ 1] ,
        _pmax[ 2] - _pmin[ 2] ,
            } ;
        _plen[ 0]*= (real_type)+8.0 ;
        _plen[ 1]*= (real_type)+8.0 ;
        _plen[ 2]*= (real_type)+8.0 ;

        _pmin[ 0]-= _plen[ 0] ;
        _pmin[ 1]-= _plen[ 1] ;
        _pmin[ 2]-= _plen[ 2] ;
        _pmax[ 0]+= _plen[ 0] ;
        _pmax[ 1]+= _plen[ 1] ;
        _pmax[ 2]+= _plen[ 2] ;

        _mesh.
        _tria.push_root (_pmin, _pmax);

    /*------------------------------ initialise mesh root */
        _mesh.
        _tria.node(+0)->fdim() = +4 ;
        _mesh.
        _tria.node(+1)->fdim() = +4 ;
        _mesh.
        _tria.node(+2)->fdim() = +4 ;
        _mesh.
        _tria.node(+3)->fdim() = +4 ;

        _mesh.
        _tria.node(+0)->feat() = +0 ;
        _mesh.
        _tria.node(+1)->feat() = +0 ;
        _mesh.
        _tria.node(+2)->feat() = +0 ;
        _mesh.
        _tria.node(+3)->feat() = +0 ;

        _mesh.
        _tria.node(+0)->topo() = +0 ;
        _mesh.
        _tria.node(+1)->topo() = +0 ;
        _mesh.
        _tria.node(+2)->topo() = +0 ;
        _mesh.
        _tria.node(+3)->topo() = +0 ;

        _geom.seed_root (_mesh, _opts);

    /*------------------------------ seed mesh from init. */
        init_init(_init, _geom, _mesh);
    }

    /*
    --------------------------------------------------------
     * RDEL-MAKE: build an rDT in R^3 .
    --------------------------------------------------------
     */

    template <
    typename      jlog_file
             >
    __static_call
    __normal_call void_type rdel_make (
        geom_type &_geom ,
        init_type &_init ,
        mesh_type &_mesh ,
        rdel_opts &_args ,
        jlog_file &_dump
        )
    {
    /*------------------------------ ensure deterministic */
        std::srand( +1 ) ;

        rdel_stat  _tcpu ;

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic ;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

        iptr_type _ncpu = +1 ;  //!! sequential for now...

    #   ifdef  __use_openmp
        omp_set_num_threads(_ncpu);
    #   else
        __unreferenced (_ncpu) ;
    #   endif//__use_openmp

    /*------------------------------ initialise mesh obj. */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        init_mesh (
            _geom , _init, _mesh, _args) ;

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._mesh_seed +=
            _tcpu.time_span(_ttic,_ttoc) ;
    #   endif//__use_timers

        iptr_type _pass = +0 ;
        iptr_type _sign = -1 ;

        iptr_type _nedg = +0 ;
        iptr_type _nfac = +0 ;
        iptr_type _ntri = +0 ;
        iptr_type _ndup = +0 ;

    /*------------------------- DT cells to check for rDT */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        iptr_list _tnew ;
        _tnew.set_alloc (
            _mesh._tria._tset.count()) ;

        iptr_type _tpos  = +0 ;
        for (auto _iter  =
            _mesh._tria._tset.head() ;
                  _iter !=
            _mesh._tria._tset.tend() ;
                ++_iter, ++_tpos)
        {
            if (_iter->mark() >= +0)
            {
                _tnew. push_tail( _tpos) ;
            }
        }

        if (_geom.have_feat(1) ||
            _geom.have_feat(2) )
        {
    /*------------------------- calc. voronoi-dual points */
    //  pragma omp parallel default(shared)
        {
    //  pragma omp for schedule(static)
        for( auto _iter  = _tnew.head();
                  _iter != _tnew.tend();
                ++_iter  )
        {
            tria_circ(_mesh, *_iter) ;
        }
        }
        }
        else
        {
    /*------------------------- prune scaffolding tria.'s */
        _mesh._tset.set_slots(_tnew.count()* 1) ;

        for( auto _iter  = _tnew.head();
                  _iter != _tnew.tend();
                ++_iter  )
        {
            push_tria(_mesh, *_iter, _args
                ) ;
        }
        }
    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._mesh_seed +=
            _tcpu.time_span(_ttic,_ttoc) ;
    #   endif//__use_timers

    /*------------------------- test for restricted edges */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        if (_args.dims() >= 1 && _geom.have_feat( 1) )
        {
            _mesh._eset.set_slots(_tnew.count() * 6) ;

            _mesh._etwo.clear();
            _mesh._etwo.set_slots(_tnew.count() * 6) ;

        //  pragma omp parallel default(shared)
            {
        //  pragma omp for schedule(static)
            for( auto _iter  = _tnew.head();
                      _iter != _tnew.tend();
                    ++_iter  )
            {
                push_edge( _mesh, _geom,*_iter, 
                           _mesh. _etwo,
                           _nedg, _pass, _args) ;
            }
            }
        }
    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._edge_init +=
            _tcpu.time_span(_ttic,_ttoc) ;
    #   endif//__use_timers

    /*------------------------- test for restricted faces */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        if (_args.dims() >= 2 && _geom.have_feat( 2) )
        {
            _mesh._fset.set_slots(_tnew.count() * 4) ;

            _mesh._ftwo.clear();
            _mesh._ftwo.set_slots(_tnew.count() * 4) ;

        //  pragma omp parallel default(shared)
            {
        //  pragma omp for schedule(static)
            for( auto _iter  = _tnew.head();
                      _iter != _tnew.tend();
                    ++_iter  )
            {
                push_face( _mesh, _geom,*_iter, 
                           _mesh. _ftwo,
                    _nfac, _ndup, _pass, _args) ;
            }
            }
        }
    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._face_init +=
            _tcpu.time_span(_ttic,_ttoc) ;
    #   endif//__use_timers

    /*------------------------- test for restricted cells */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        if (_args.dims() >= 3 && _geom.have_feat( 2) )
        {
        //  _mesh._tset.set_slots(_tnew.count() * 1) ;

        //  _mesh._ttwo.clear();
        //  _mesh._ttwo.set_slots(_tnew.count() * 1) ;

            bool_type _safe = true ;
        //  if (_dim0 >  +0) _safe = false ;
        //  if (_nfac >= +1) _safe = false ;
            if (_ndup >= +1) _safe = false ;

        //  pragma omp parallel default(shared)
            {
        //  pragma omp for schedule(static)
            for( auto _iter  = _tnew.head();
                      _iter != _tnew.tend();
                    ++_iter )
            {
               _sign =   (!_safe) ? -1 : _sign;

                push_tria( _mesh, _geom,*_iter, 
                    _sign, _ntri, _pass, _args) ;
            }
            }
        }
    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._tria_init +=
            _tcpu.time_span(_ttic,_ttoc) ;
    #   endif//__use_timers

        if (_args.verb() >= +2 )
        {
    /*------------------------- push rDT scheme's metrics */
        _dump.push("\n")  ;
        _dump.push("**TIMING statistics...\n") ;

        _dump.push(" *mesh-seed = ") ;
        _dump.push(std::to_string (_tcpu._mesh_seed));
        _dump.push("\n")  ;
        _dump.push(" *edge-init = ") ;
        _dump.push(std::to_string (_tcpu._edge_init));
        _dump.push("\n")  ;
        _dump.push(" *face-init = ") ;
        _dump.push(std::to_string (_tcpu._face_init));
        _dump.push("\n")  ;
        _dump.push(" *tria-init = ") ;
        _dump.push(std::to_string (_tcpu._tria_init));
        _dump.push("\n")  ;

        _dump.push("\n")  ;
        _dump.push("**rDT(x) statistics...\n") ;

        _dump.push(" *rDEL-edge = ");
        _dump.push(std::to_string (_mesh._eset.count())) ;
        _dump.push("\n")  ;
        _dump.push(" *rDEL-face = ");
        _dump.push(std::to_string (_mesh._fset.count())) ;
        _dump.push("\n")  ;
        _dump.push(" *rDEL-tria = ");
        _dump.push(std::to_string (_mesh._tset.count())) ;
        _dump.push("\n")  ;

        _dump.push("\n")  ;
        }

        if (_args.verb() >= +2 )
        {
    /*------------------------- more rDT scheme's metrics */
        _dump.push("\n")  ;
        _dump.push("**MEMORY statistics...\n") ;

        _dump.push("  xDEL-type:\n") ;
        _dump.push(" *node-byte = ") ;
        _dump.push(std::to_string(sizeof(
            typename mesh_type::tria_type::node_type)));
        _dump.push("\n")  ;
        _dump.push(" *nset-size = ") ;
        _dump.push(
            std::to_string(_mesh._tria._nset.alloc())) ;
        _dump.push("\n")  ;
        _dump.push(" *tria-byte = ") ;
        _dump.push(std::to_string(sizeof(
            typename mesh_type::tria_type::tria_type)));
        _dump.push("\n")  ;
        _dump.push(" *tset-size = ") ;
        _dump.push(
            std::to_string(_mesh._tria._tset.alloc())) ;
        _dump.push("\n")  ;
        _dump.push(" *pool-byte = ") ;
        _dump.push(
            std::to_string(_mesh._tria._fpol.bytes())) ;
        _dump.push("\n")  ;

        _dump.push("\n")  ;
        _dump.push("  rDEL-type:\n") ;
        _dump.push(" *edge-byte = ") ;
        _dump.push(std::to_string(
            sizeof(typename mesh_type::edge_item))) ;
        _dump.push("\n")  ;
        _dump.push(" *eset-size = ") ;
        _dump.push(
            std::to_string(_mesh._eset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push(" *pool-byte = ") ;
        _dump.push(std::to_string(_mesh._epol.bytes()));
        _dump.push("\n")  ;
        _dump.push(" *face-byte = ") ;
        _dump.push(std::to_string(
            sizeof(typename mesh_type::face_item))) ;
        _dump.push("\n")  ;
        _dump.push(" *fset-size = ") ;
        _dump.push(
            std::to_string(_mesh._fset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push(" *pool-byte = ") ;
        _dump.push(std::to_string(_mesh._fpol.bytes()));
        _dump.push("\n")  ;
        _dump.push(" *tria-byte = ") ;
        _dump.push(std::to_string(
            sizeof(typename mesh_type::tria_item))) ;
        _dump.push("\n")  ;
        _dump.push(" *tset-size = ") ;
        _dump.push(
            std::to_string(_mesh._tset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push(" *pool-byte = ") ;
        _dump.push(std::to_string(_mesh._tpol.bytes()));
        _dump.push("\n")  ;

        _dump.push("\n")  ;
        }
    }

    } ;


    }

#   endif   //__RDEL_MAKE_3__



