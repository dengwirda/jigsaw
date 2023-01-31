
    /*
    --------------------------------------------------------
     * ITER-MESH-2: mesh-optimisation for 2-complexes.
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
     * Last updated: 12 Dec., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda,
     * Marc Tunnell
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_MESH_2__
#   define __ITER_MESH_2__

    namespace  mesh {

    /*
    --------------------------------------------------------
     * ITER-MESH-2: optimise mixed 2-complex meshes
    --------------------------------------------------------
     */

    template <
    typename G ,
    typename M ,
    typename H ,
    typename P
             >
    class iter_mesh_2
    {
    public  :
    typedef  M                          mesh_type ;
    typedef  G                          geom_type ;
    typedef  H                          hfun_type ;
    typedef  P                          pred_type ;

    typedef typename
            mesh_type::real_type        real_type ;
    typedef typename
            mesh_type::iptr_type        iptr_type ;

    iptr_type static constexpr min_subit = +2 ;
    iptr_type static constexpr max_subit = +8 ;

    iptr_type static constexpr
            topo_dims =      pred_type::topo_dims ;
    iptr_type static constexpr
            geom_dims =      pred_type::geom_dims ;
    iptr_type static constexpr
            real_dims =      pred_type::real_dims ;

    typedef char_type                   kern_kind ;

    char_type static                // optm. kern. selector
        constexpr _odt_optimise = +1 ;
    char_type static
        constexpr _cvt_optimise = +2 ;
    char_type static
        constexpr _h95_optimise = +3 ;
    char_type static
        constexpr dqdx_optimise = +5 ;

    typedef char_type                   flip_kind ;

    flip_kind static                // optm. kern. selector
        constexpr laguerre_flip = +1 ;
    flip_kind static
        constexpr delaunay_flip = +2 ;
    flip_kind static
        constexpr topology_flip = +3 ;

    class   cell_kind{} ;           // dummies for overload
    class   dual_kind{} ;

    typedef mesh::iter_params  <
            real_type ,
            iptr_type          >        iter_opts ;

    typedef mesh::iter_timers           iter_stat ;

    typedef containers
            ::array< iptr_type >        iptr_list ;
    typedef containers
            ::array< size_t    >        uint_list ;
    typedef containers
            ::array< real_type >        real_list ;

    typedef containers::fixed_array <
            real_type ,
            geom_dims          >        line_data ;

    typedef containers
            ::array< line_data >        line_list ;

    class mark_list                 // integer cell markers
        {
        public  :
        iptr_list                      _node;
        iptr_list                      _edge;
        iptr_list                      _tri3;
        iptr_list                      _quad;
        } ;

    class part_sets                 // mesh parallel decomp.
        {
        public  :
        uint_list                      _lptr;
        iptr_list                      _list;
        iptr_list                      _part;
        iptr_list                      _seqs;
        } ;

    class part_data
        {
        public  :
        part_sets                      _full;
        part_sets                      _bnds;
        } ;

    typedef typename
            mesh_type::connector        conn_list ;

    class conn_sets                 // cache neighbourhoods
        {
        public  :
        conn_list                      _adj1;
        uint_list                      _idx1;
        conn_list                      _adj2;
        uint_list                      _idx2;
        } ;

    public  :

    /*
    --------------------------------------------------------
     * FLIP-SIGN: flip cells for +ve iter. cost fn.
    --------------------------------------------------------
     */

    #include "_bfs_mesh_2.inc"


    /*
    --------------------------------------------------------
     * COST-MESH: 1-neighbourhood cost vector util.
    --------------------------------------------------------
     */

    #include "cost_mesh_2.inc"


    /*
    --------------------------------------------------------
     * PART-MESH: paritions for parallel decomposition.
    --------------------------------------------------------
     */

    #include "part_mesh_2.inc"


    /*
    --------------------------------------------------------
     * SORT-MESH: greedy ordering for optimisation ops.
    --------------------------------------------------------
     */

    #include "sort_mesh_2.inc"


    /*
    --------------------------------------------------------
     * MOVE-MESH: do single cell/dual smoothing passes.
    --------------------------------------------------------
     */

    #include "move_mesh_2.inc"


    /*
    --------------------------------------------------------
     * FLIP-MESH: "flip" mesh topology until Delaunay.
    --------------------------------------------------------
     */

    #include "flip_mesh_2.inc"


    /*
    --------------------------------------------------------
     * _ZIP-MESH: merge/split operations on cells.
    --------------------------------------------------------
     */

    #include "_zip_mesh_2.inc"


    /*------------------------------ helper: init. marker */

    __static_call
    __normal_call void_type init_mark (
        mesh_type &_mesh ,
        mark_list &_mark ,
        iptr_type  _flag = +0
        )
    {
        iptr_type _nnN1 = std::max(
           (iptr_type) _mark. _node.count() ,
           (iptr_type) _mesh.node().count()
                ) ;
        iptr_type _nnE2 = std::max(
           (iptr_type) _mark. _edge.count() ,
           (iptr_type) _mesh.edge().count()
                ) ;
        iptr_type _nnT3 = std::max(
           (iptr_type) _mark. _tri3.count() ,
           (iptr_type) _mesh.tri3().count()
                ) ;
        iptr_type _nnQ4 = std::max(
           (iptr_type) _mark. _quad.count() ,
           (iptr_type) _mesh.quad().count()
                ) ;

        _mark._node.set_count(_nnN1,
            containers::loose_alloc, _flag) ;

        _mark._edge.set_count(_nnE2,
            containers::loose_alloc, _flag) ;

        _mark._tri3.set_count(_nnT3,
            containers::loose_alloc, _flag) ;

        _mark._quad.set_count(_nnQ4,
            containers::loose_alloc, _flag) ;
    }

    /*------------------------------ helper: init. bounds */

    __static_call
    __normal_call void_type init_bnds (
        mesh_type &_mesh ,
        mark_list &_mark
        )
    {
    /*------------------------------ mark bnd. item <= -1 */
        conn_list _conn ;
        iptr_type _nnN1  = +0 ;
        for (auto _node  = _mesh.node().head() ;
                  _node != _mesh.node().tend() ;
                ++_node, ++_nnN1 )
        {
            if (_node->mark() >= +0)
            {
            if (_node->feat()
                    != mesh::null_feat)
            {
                _mark._node[_nnN1]  = -1 ;
            }
            }
        }

        iptr_type _nnE2  = +0 ;
        for (auto _edge  = _mesh.edge().head() ;
                  _edge != _mesh.edge().tend() ;
                ++_edge, ++_nnE2 )
        {
            if (_edge->mark() >= +0)
            {
            if (_edge->self() >= +1)
            {
                _mark._node[
                    _edge->node(0)] = -1 ;
                _mark._node[
                    _edge->node(1)] = -1 ;
            }
            else
            {
                _conn.set_count(0) ;

                _mesh.connect_2(
            &_edge->node(0), EDGE2_tag, _conn) ;

                if (_conn.count () != +2)
                {
                _mark._node[
                    _edge->node(0)] = -1 ;
                _mark._node[
                    _edge->node(1)] = -1 ;
                }
            }
            }
        }
    }

    /*------------------------------ helper: init. scores */

    __static_call
    __normal_call real_type init_cost (
        mesh_type &_mesh,
        iter_opts &_opts
        )
    {
    #   ifdef  __use_openmp
        omp_set_num_threads(_opts.nprt()) ;
    #   endif//__use_openmp

        real_type _QMIN = (real_type) +1. ;

    #   pragma omp parallel default(none) \
           shared(_mesh, _opts, _QMIN)
        {
        real_type _qmin = (real_type) +1. ;
    #   pragma omp for nowait  schedule(static)
        for (auto _cpos  = +0u ;
                  _cpos  < _mesh.tri3().count() ;
                ++_cpos  )
        {
             auto _cell  =
                _mesh. tri3().head() + _cpos ;

            if (_cell->mark() >= +0 )
            {
        /*--------------------- test initial cell quality */
            real_type _cost;
            _cost = pred_type::tri3_cost  (
               &_mesh .node(
                _cell->node(0)).pval(0),
               &_mesh .node(
                _cell->node(1)).pval(0),
               &_mesh .node(
                _cell->node(2)).pval(0),
                typename
                pred_type::cell_kind ())  ;

            _qmin = std::min (_qmin, _cost) ;
            _qmin = std::max (_qmin,
                              _opts.qlim()) ;
            }
        }

    #   pragma omp for nowait  schedule(static)
        for (auto _cpos  = +0u ;
                  _cpos  < _mesh.quad().count() ;
                ++_cpos  )
        {
             auto _cell  =
                _mesh. quad().head() + _cpos ;

            if (_cell->mark() >= +0 )
            {
        /*--------------------- test initial cell quality */
            real_type _cost;
            _cost = pred_type::quad_cost  (
               &_mesh .node(
                _cell->node(0)).pval(0),
               &_mesh .node(
                _cell->node(1)).pval(0),
               &_mesh .node(
                _cell->node(2)).pval(0),
               &_mesh .node(
                _cell->node(3)).pval(0),
                typename
                pred_type::cell_kind ())  ;

            _qmin = std::min (_qmin, _cost) ;
            _qmin = std::max (_qmin,
                              _opts.qlim()) ;
            }
        }

    #   ifdef  __use_openmp
        auto _nprt = omp_get_num_threads();

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank = 0 ; _rank <_nprt;
                ++_rank )
        {
    #       pragma omp ordered
            {
            _QMIN = std::min (_QMIN, _qmin) ;
            }
        }
    #   endif//__use_openmp
        }

        return _QMIN ;
    }

    /*
    --------------------------------------------------------
     * ITER-MESH: "hill-climbing" type mesh optimisation.
    --------------------------------------------------------
     */

    template <
        typename  text_dump
             >
    __static_call
    __normal_call void_type iter_mesh (
        geom_type &_geom ,  // geometry object
        hfun_type &_hfun ,  // spacing object
        mesh_type &_mesh ,  // mesh object
        kern_kind  _kern ,  // optim. kernel selector
        iter_opts &_opts ,  // user options
        text_dump &_dump    // log-file object
        )
    {
        iter_stat  _tcpu ;

    /*------------------------------ push log-file header */
        if (_opts.verb() >= 0 )
        {
            _dump.push(
    "#------------------------------------------------------------\n"
    "#    |MOVE.|      |FLIP.|      |MERGE|      |SPLIT| \n"
    "#------------------------------------------------------------\n"
                    ) ;
        }

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic ;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*------------------------------ ensure deterministic */
        std::srand( +1 ) ;

    /*------------------------------ push boundary marker */
        iptr_list _nset ;
        mark_list _mark ;

        init_mark(_mesh, _mark) ;
        init_bnds(_mesh, _mark) ;

        flip_sign(_mesh);

    /*------------------------------ do optimisation loop */
        bool_type
        static constexpr ITER_FLIP = true;

        iptr_type
        static constexpr ITER_MIN_ = min_subit ;
        iptr_type
        static constexpr ITER_MAX_ = max_subit ;

        real_type _QMIN = init_cost (_mesh, _opts) ;
        real_type _QMOV = +0. ;
        real_type _DMOV = +0. ;

        real_list _hval ;         // cache h(x) node val.
        _hval.set_count(
            _mesh. node().count(),
        containers::loose_alloc, (real_type) -1.0) ;

        for (auto _iter = +1;
            _iter <= _opts.iter(); ++_iter)
        {
    /*------------------------------ set-up current iter. */
            init_mark(_mesh ,
                _mark,std::max(_iter - 1, +0)) ;

            _nset.set_count(  +0);

            size_t    _nmov = +0 ;
            size_t    _nflp = +0 ;
            size_t    _nzip = +0 ;
            size_t    _ndiv = +0 ;

    /*------------------------------ scale quality limits */
            iptr_type _nsub = _iter +1 ;

            _nsub =
            std::min(ITER_MAX_, _nsub) ;
            _nsub =
            std::max(ITER_MIN_, _nsub) ;

            real_type _RMIN = _QMIN *
           (real_type)(0.90 + 1./20 * (_iter - 1)) ;
            real_type _RMAX = _QMIN *
           (real_type)(0.90 + 1./30 * (_iter - 1)) ;

            real_type _DLIM =
           (real_type)+1.-_opts.qtol() ;

            _DMOV = std::max(_DMOV, _DLIM) ;
            _DMOV = std::min(_DMOV, _RMAX) ;

            real_type _QLIM = std::min (
                _opts.qlim(),_RMIN) ;

            _QMOV = std::max(_QMOV, _QLIM) ;
            _QMOV = std::min(_QMOV, _RMAX) ;

    /*------------------------------ 1. CELL GEOM. PASSES */

            if (_opts.tria())
            {
    /*------------------------------ update mesh geometry */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            line_data _zero ;
            _zero.fill((real_type) +0.) ;

            line_list _last ;     // for momentum updates
            _last.set_count(
                _mesh.node().count() ,
                    containers::tight_alloc, _zero) ;

            iptr_list _aset ;
            _aset.set_alloc(      // set of "active" node
                _mesh.node().count());

            part_data _part ;
            conn_sets _conn ;
            for (auto _isub = 0; _isub != _nsub; ++_isub)
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL MOVE-NODE...\n" ) ;

                move_node( _geom, _mesh , _conn ,
                    _hfun, _kern, _hval , _last ,
                    _nset, _aset, _mark ,
                    _part, _iter, _isub ,
                    _opts,
                    _QLIM, _QMOV, _DLIM , _DMOV ,
                    _tcpu) ;
            }

            _nmov = _nset.count() ;

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._move_node +=
                  _tcpu.time_span(_ttic , _ttoc);
    #       endif//__use_timers

    /*------------------------------ update mesh topology */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            if (ITER_FLIP)
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL FLIP-MESH...\n" ) ;

                size_t     _nloc;
                flip_mesh( _geom, _mesh , _hfun ,
                    _conn, _nset, _mark ,
                +3 * _iter - 2  , _nloc , _tcpu);

                _nflp +=   _nloc;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._topo_flip +=
                  _tcpu.time_span(_ttic , _ttoc);
    #       endif//__use_timers
            }

    /*------------------------------ 2. DUAL GEOM. PASSES */

            if (_opts.dual())
            {
    /*------------------------------ update mesh geometry */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            real_list _last ;     // for momentum updates
            _last.set_count(
                _mesh.node().count() ,
                    containers::tight_alloc,  +0. ) ;

            iptr_list _aset ;
            _aset.set_alloc(      // set of "active" node
                _mesh.node().count());

            _nsub = std::max(_nsub/2, +1) ;

            part_data _part ;
            conn_sets _conn ;
            for (auto _isub = 0; _isub != _nsub; ++_isub)
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL MOVE-DUAL...\n" ) ;

                move_dual( _geom, _mesh , _conn ,
                    _hfun, _hval, _last ,
                    _nset, _aset, _mark ,
                    _part, _iter, _isub ,
                    _opts,
                    _QLIM, _QMOV, _DLIM , _DMOV ,
                    _tcpu) ;
            }

            _nmov = _nset.count() ;

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._move_dual +=
                  _tcpu.time_span(_ttic , _ttoc);
    #       endif//__use_timers

    /*------------------------------ update mesh topology */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            if (ITER_FLIP)
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL FLIP-MESH...\n" ) ;

                size_t     _nloc;
                flip_mesh( _geom, _mesh , _hfun ,
                    _conn, _nset, _mark ,
                +3 * _iter - 1  , _nloc , _tcpu);

                _nflp +=   _nloc;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._topo_flip +=
                  _tcpu.time_span(_ttic , _ttoc);
    #       endif//__use_timers

            }

    /*------------------------------ 3. ZIP + DIV SUBFACE */

            if (_iter < _opts.iter())
            {
    /*------------------------------ change mesh topology */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            _nset.set_count(+0) ;    // don't flip twice!

            if (_opts.zip_ () ||
                _opts.div_ () )
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL _ZIP-MESH...\n" ) ;

                _zip_mesh( _geom, _mesh , _hfun ,
                    _kern, _hval, _nset ,
                    _mark, _iter, _opts ,
                    _QLIM, _DLIM,
                    _nzip, _ndiv, _tcpu ) ;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._topo_zips +=
                  _tcpu.time_span(_ttic , _ttoc);
    #       endif//__use_timers

    /*------------------------------ update mesh topology */

    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            if (ITER_FLIP)
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL FLIP-MESH...\n" ) ;

                size_t     _nloc;
                flip_mesh( _geom, _mesh , _hfun ,
                    _nset, _mark,
                +3 * _iter - 0  , _nloc , _tcpu);

                _nflp +=   _nloc;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._topo_flip +=
                  _tcpu.time_span(_ttic , _ttoc);
    #       endif//__use_timers

            } // if (_iter < _opts.iter())

    /*------------------------------ dump optim. progress */
            if (_opts.verb() >= 0)
            {
            std::stringstream _sstr ;
            _sstr << std::setw(11) << _nmov
                  << std::setw(13) << _nflp
                  << std::setw(13) << _nzip
                  << std::setw(13) << _ndiv
                  <<   "\n" ;
            _dump.push(_sstr.str()) ;
            }

    /*------------------------------ has iter. converged? */
            _QMOV = std::max(_QLIM, 1. -
                ((real_type) _nmov) /
                    _mesh.node().count()) ;

            _DMOV = std::max(_DLIM, 1. -
                ((real_type) _nmov) /
                    _mesh.node().count()) ;

        //  if (_nset.count() == 0) break ;
            if (_nmov == +0 &&
                _nzip == +0 &&
                _ndiv == +0 &&
                _nflp == +0 )       break ;
        }

        if (_opts.verb() >= +2)
        {
    /*------------------------------ print method metrics */
            _dump.push("\n");
            _dump.push("**TIMING statistics...\n") ;

            _dump.push("  MOVE-NODE: ");
            _dump.push(
            std::to_string(_tcpu._move_node)) ;
            _dump.push("\n");
            _dump.push(" *init-node: ");
            _dump.push(
            std::to_string(_tcpu._init_node)) ;
            _dump.push("\n");
            _dump.push(" *part-node: ");
            _dump.push(
            std::to_string(_tcpu._part_node)) ;
            _dump.push("\n");
            _dump.push(" *core-node: ");
            _dump.push(
            std::to_string(_tcpu._core_node)) ;
            _dump.push("\n");
            _dump.push(" *seqs-node: ");
            _dump.push(
            std::to_string(_tcpu._seqs_node)) ;
            _dump.push("\n");
            _dump.push(" *para-node: ");
            _dump.push(
            std::to_string(_tcpu._para_node)) ;
        //  _dump.push("\n");
        //  _dump.push(" *xDIR-node: ");
        //  _dump.push(
        //  std::to_string(_tcpu._ldir_node)) ;
        //  _dump.push("\n");
        //  _dump.push(" *xOPT-node: ");
        //  _dump.push(
        //  std::to_string(_tcpu._lopt_node)) ;
            _dump.push("\n\n");

            _dump.push("  MOVE-DUAL: ");
            _dump.push(
            std::to_string(_tcpu._move_dual)) ;
            _dump.push("\n");
            _dump.push(" *init-dual: ");
            _dump.push(
            std::to_string(_tcpu._init_dual)) ;
            _dump.push("\n");
            _dump.push(" *part-dual: ");
            _dump.push(
            std::to_string(_tcpu._part_dual)) ;
            _dump.push("\n");
            _dump.push(" *core-dual: ");
            _dump.push(
            std::to_string(_tcpu._core_dual)) ;
            _dump.push("\n");
            _dump.push(" *seqs-dual: ");
            _dump.push(
            std::to_string(_tcpu._seqs_dual)) ;
            _dump.push("\n");
            _dump.push(" *para-dual: ");
            _dump.push(
            std::to_string(_tcpu._para_dual)) ;
        //  _dump.push("\n");
        //  _dump.push(" *xDIR-dual: ");
        //  _dump.push(
        //  std::to_string(_tcpu._ldir_dual)) ;
        //  _dump.push("\n");
        //  _dump.push(" *xOPT-dual: ");
        //  _dump.push(
        //  std::to_string(_tcpu._lopt_dual)) ;
            _dump.push("\n\n");

            _dump.push("  TOPO-FLIP: ");
            _dump.push(
            std::to_string(_tcpu._topo_flip)) ;
            _dump.push("\n");
            _dump.push(" *init-flip: ");
            _dump.push(
            std::to_string(_tcpu._init_flip)) ;
            _dump.push("\n");
            _dump.push(" *core-flip: ");
            _dump.push(
            std::to_string(_tcpu._core_flip)) ;
            _dump.push("\n\n");

            _dump.push("  TOPO-ZIPS: ");
            _dump.push(
            std::to_string(_tcpu._topo_zips)) ;
            _dump.push("\n");
            _dump.push(" *init-zips: ");
            _dump.push(
            std::to_string(_tcpu._init_zips)) ;
            _dump.push("\n");
            _dump.push(" *core-divs: ");
            _dump.push(
            std::to_string(_tcpu._core_divs)) ;
            _dump.push("\n");
            _dump.push(" *core-zips: ");
            _dump.push(
            std::to_string(_tcpu._core_zips)) ;
            _dump.push("\n\n");

            _dump.push("  PARTITION: ");
            _dump.push(
            std::to_string(_tcpu._full_part)) ;
            _dump.push("\n");
            _dump.push(" *tree-part: ");
            _dump.push(
            std::to_string(_tcpu._tree_part)) ;
        //  _dump.push("\n");
        //  _dump.push(" *redo-part: ");
        //  _dump.push(
        //  std::to_string(_tcpu._redo_part)) ;
            _dump.push("\n");
            _dump.push(" *part-part: ");
            _dump.push(
            std::to_string(_tcpu._part_part)) ;
            _dump.push("\n");
            _dump.push(" *seqs-part: ");
            _dump.push(
            std::to_string(_tcpu._seqs_part)) ;
            _dump.push("\n");

            _dump.push("\n");
        }
        else
        {
            _dump.push("\n");
        }
    }

    } ;


    }

#   endif   //__ITER_MESH_2__



