
    /*
    --------------------------------------------------------
     * SORT-MESH-2: greedy ordering for mesh optim.
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
     * Last updated: 18 Dec., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda,
     * Marc Tunnell
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    /*
    --------------------------------------------------------
     * PULL-CONN: cache node-to-cell adjacency lists.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type pull_conn (
        mesh_type &_mesh ,  // mesh object
        iptr_list &_aset ,  // list of nodes active
        iptr_list &_amrk ,  // >= 0 if nodes active
        conn_sets &_conn    // cached node spacing array
        )
    {
        _conn._adj1.set_count( +0 ) ;
        _conn._idx1.set_count(
            _mesh.node().count() + 1,
                containers::tight_alloc , +0 ) ;

        _conn._adj2.set_count( +0 ) ;
        _conn._idx2.set_count(
            _mesh.node().count() + 1,
                containers::tight_alloc , +0 ) ;

         auto _scal =
            ((double)_aset.count()) /
                     _mesh.node().count() ;

        _conn._adj1.set_alloc(
            (size_t)(_scal * (
            _mesh.edge().count() * 2) ) ) ;
        _conn._adj2.set_alloc(
            (size_t)(_scal * (
            _mesh.tri3().count() * 3 +
            _mesh.quad().count() * 4) ) ) ;

        iptr_type _npos  = +0 ;
        for (auto _node  = _mesh.node().head() ;
                  _node != _mesh.node().tend() ;
                ++_node, ++_npos)
        {
            if (_node->mark() >= 0 &&
                _amrk[_npos ] >= 0)
            {
    /*-------------------- append adj. lists per node */
            _mesh.connect_1(
            _npos, POINT_tag, _conn._adj1) ;

            _conn._idx1[_npos+1] =
                _conn._adj1.count() ;

            _mesh.connect_2(
            _npos, POINT_tag, _conn._adj2) ;

            _conn._idx2[_npos+1] =
                _conn._adj2.count() ;
            }
            else
            {
    /*-------------------- empty lists for null nodes */
            _conn._idx1[_npos+1] =
                _conn._idx1[ _npos] ;

            _conn._idx2[_npos+1] =
                _conn._idx2[ _npos] ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * SORT-NODE: permutations for node optimisation.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type sort_node (
        mesh_type &_mesh ,  // mesh object
        conn_sets &_CONN ,  // cached node adj. indexing
        iptr_list &_aset ,  // list of nodes to try next
        iptr_list &_nmrk ,  // sets of node flags
        iptr_type  _iout ,  // outer iteration marker
        iptr_type  _isub ,  // inner iteration marker
        real_type  _QLIM ,  // cell quality threshold
        real_type  _QMOV ,
        real_type  _DLIM ,  // dual quality threshold
        real_type  _DMOV ,
        iter_opts &_opts    // user options
        )
    {
    #   define PUSHCOST(_NODE , _COST)      \
            _qbar[_NODE] += _COST ;         \
            _qmin[_NODE]  = \
        std::min (_qmin[_NODE], _COST) ;

        class cost_pair
            {
            public  :
    /*------------------------ tuple for node re-ordering */
            iptr_type           _node ;
            float               _cost ;

    /*------------------------ construct inline from src. */
            __inline_call cost_pair (
                iptr_type const&_nsrc ,
                real_type const&_csrc
                )   :
                   _node (_nsrc),
                   _cost((float)_csrc) {}
            } ;

        class cost_pred
            {
            public  :
    /*------------------------ less-than op. for cost-tup */
            __inline_call
                bool_type operator () (
                cost_pair const&_idat ,
                cost_pair const&_jdat
                ) const
            { return _idat._cost > _jdat._cost ;
            }
            } ;

        typedef
        containers::array<cost_pair> cost_list ;

        real_list _qbar, _qmin, _cost ;
        iptr_list _amrk ;
        cost_list _sset ;

        __unreferenced ( _QMOV ) ;
        __unreferenced ( _DMOV ) ;
        __unreferenced ( _DLIM ) ;
        __unreferenced ( _opts ) ;
        __unreferenced ( _isub ) ;

    /*-------------------- 1ST SUB-ITER: build full init. */
        _qbar.set_count(
            _mesh.node().count(),
                containers::tight_alloc, +0.0) ;
        _qmin.set_count(
            _mesh.node().count(),
                containers::tight_alloc, +1.0) ;

        _amrk.set_count(
            _mesh.node().count(),
                containers::tight_alloc, -1) ;

    #   ifdef  __use_openmp
        omp_set_num_threads(_opts.nprt()) ;
    #   endif//__use_openmp

        _cost.set_count(
            _mesh.tri3().count(),
                containers::loose_alloc, +0.0) ;

    #   pragma omp parallel for schedule(static) \
           default(none) \
           shared(_mesh, _cost)
        for (auto _tpos = +0u ;
                  _tpos < _mesh.tri3().count() ;
                ++_tpos )
        {
             auto _tria =
                  _mesh.tri3().head() + _tpos;

            if (_tria->mark() >= +0)
            {
        /*-------------------- calc. min. scores at cells */
                _cost[_tpos] =
                    pred_type::tri3_cost (
               &_mesh .node(
                _tria->node(0)).pval(0),
               &_mesh .node(
                _tria->node(1)).pval(0),
               &_mesh .node(
                _tria->node(2)).pval(0),
                typename
                pred_type::cell_kind ()) ;
            }
        }

        iptr_type _tpos  = 0;
        for (auto _tria  = _mesh.tri3().head() ;
                  _tria != _mesh.tri3().tend() ;
                ++_tria, ++_tpos)
        {
            if (_tria->mark() >= +0)
            {
        /*-------------------- dist. min. scores to nodes */
                iptr_type _inod, _jnod, _knod;
                _inod = _tria->node( 0);
                _jnod = _tria->node( 1);
                _knod = _tria->node( 2);

                real_type _cval = _cost[_tpos] ;

                PUSHCOST( _inod, _cval )
                PUSHCOST( _jnod, _cval )
                PUSHCOST( _knod, _cval )
            }
        }

        _cost.set_count(
            _mesh.quad().count(),
                containers::loose_alloc, +0.0) ;

    #   pragma omp parallel for schedule(static) \
           default(none) \
           shared(_mesh, _cost)
        for (auto _qpos = +0u ;
                  _qpos < _mesh.quad().count() ;
                ++_qpos )
        {
             auto _quad =
                  _mesh.quad().head() + _qpos;

            if (_quad->mark() >= +0)
            {
        /*-------------------- calc. min. scores at cells */
                _cost[_qpos] =
                    pred_type::quad_cost (
               &_mesh .node(
                _quad->node(0)).pval(0),
               &_mesh .node(
                _quad->node(1)).pval(0),
               &_mesh .node(
                _quad->node(2)).pval(0),
               &_mesh .node(
                _quad->node(3)).pval(0),
                typename
                pred_type::cell_kind ()) ;
            }
        }

        iptr_type _qpos  = 0;
        for (auto _quad  = _mesh.quad().head() ;
                  _quad != _mesh.quad().tend() ;
                ++_quad, ++_qpos)
        {
            if (_quad->mark() >= +0)
            {
        /*-------------------- dist. min. scores to nodes */
                iptr_type _inod, _jnod, _knod,
                          _lnod;
                _inod = _quad->node( 0);
                _jnod = _quad->node( 1);
                _knod = _quad->node( 2);
                _lnod = _quad->node( 3);

                real_type _cval = _cost[_qpos] ;

                PUSHCOST( _inod, _cval )
                PUSHCOST( _jnod, _cval )
                PUSHCOST( _knod, _cval )
                PUSHCOST( _lnod, _cval )
            }
        }

        iptr_type _FLAG  = _iout - 2 ;  // only "recent"

        _sset.set_alloc(
            _mesh.node().count())  ;

        iptr_type _ipos  = +0 ;
        for (auto _node  = _mesh.node().head() ;
                  _node != _mesh.node().tend() ;
                ++_node, ++_ipos)
        {
    /*------------------------ assign score for each node */
            if (_node->mark() >= +0)
            {
                if (std::abs(
                    _nmrk[_ipos]) >= _FLAG || // recent
                    _qmin[_ipos]  <= _QLIM )  // skewed
                {
                    auto  _nadj =
                        _node->ndeg(TRIA3_tag) +
                        _node->ndeg(QUAD4_tag) ;

                    _nadj = std::max(1, _nadj) ;

                    float _cval =
                        (float)_qbar[_ipos]/_nadj -
                        (float)_qmin[_ipos] ;

                    _sset.push_tail (
                        cost_pair(_ipos , _cval)) ;
                }
            }
        }

        algorithms::qsort( _sset.head() ,
                           _sset.tend() ,
            cost_pred () ) ;

        for (auto _iter  = _sset.head() ;
                  _iter != _sset.tend() ;
                ++_iter  )
        {
    /*------------------------ push sorted wrt. min.-cost */
             auto _spos  = _iter->_node ;
            _amrk[_spos] = _isub ;
            _aset.push_tail(_spos) ;
        }

    /*------------------------ cache adj. for active node */
        pull_conn(
            _mesh, _aset, _amrk, _CONN) ;

    #   undef   PUSHCOST
    }



