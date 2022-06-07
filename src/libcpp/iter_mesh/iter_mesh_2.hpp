
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 05 Jun., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
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

    typedef mesh::iter_timers  <
            real_type ,
            iptr_type          >        iter_stat ;

    typedef containers
            ::array< iptr_type >        iptr_list ;
    typedef containers
            ::array< real_type >        real_list ;

    class mark_list                 // integer cell markers
        {
        public  :
        iptr_list                      _node;
        iptr_list                      _edge;
        iptr_list                      _tri3;
        iptr_list                      _quad;
        } ;

    typedef typename
        mesh_type::connector            conn_list ;

    class conn_sets
        {
        public  :
        conn_list                      _adj1;
        iptr_list                      _idx1;
        conn_list                      _adj2;
        iptr_list                      _idx2;
        } ;

    public  :

    /*
    --------------------------------------------------------
     * FLIP-SIGN: flip cells for +ve iter. cost fn.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type flip_next (
        mesh_type &_mesh,
        conn_list &_list,
        iptr_type  _icel,
        iptr_type *_inod,
        mark_list &_seen,
        conn_list &_conn
        )
    {
    /*------------------------- flip to match cell "sign" */
        _conn.set_count(0) ;
        _mesh.connect_2(
            _inod, EDGE2_tag, _conn);

        if (_conn.count() != 2) return;

        iptr_type _jcel, _jfac;
        char_type _kind;
        if ((iptr_type)
            _conn[0]._cell == _icel)
        {
            _jcel = _conn[ 1]._cell ;
            _kind = _conn[ 1]._kind ;
        }
        else
        {
            _jcel = _conn[ 0]._cell ;
            _kind = _conn[ 0]._kind ;
        }

    /*------------------------- match indexing along edge */
        if (_kind == TRIA3_tag)
        {
        if (_seen._tri3[_jcel] == +0)
        {
            _list.push_tail(typename
                conn_list::data_type(_jcel, _kind)) ;

            _seen._tri3[_jcel] =  +1;

            iptr_type _jnod[3] ;
            for(_jfac = 3 ; _jfac-- != 0; )
            {
            mesh_type::tri3_type::
            face_node(_jnod, _jfac, 2, 1) ;
            _jnod[ 0] = _mesh.
             tri3(_jcel).node(_jnod[0]) ;
            _jnod[ 1] = _mesh.
             tri3(_jcel).node(_jnod[1]) ;
            _jnod[ 2] = _mesh.
             tri3(_jcel).node(_jnod[2]) ;

            if (_jnod[ 2] != _inod[ 0])
            if (_jnod[ 2] != _inod[ 1])
                break  ;
            }

            if (_jnod[ 0] == _inod[ 0] &&
                _jnod[ 1] == _inod[ 1])
            {
                std::swap (
            _mesh.tri3(_jcel).node( 0),
            _mesh.tri3(_jcel).node( 1)) ;
            }
        }
        }
        else
        if (_kind == QUAD4_tag)
        {
        if (_seen._quad[_jcel] == +0)
        {
            _list.push_tail(typename
                conn_list::data_type(_jcel, _kind)) ;

            _seen._quad[_jcel] =  +1;

            // TODO: reverse quad indexes
        }
        }
    }

    __static_call
    __normal_call void_type flip_bfs_ (
        mesh_type &_mesh,
        conn_list &_list,
        mark_list &_seen
        )
    {
        conn_list  _conn;
        for ( ; !_list.empty() ; )
        {
    /*------------------------- flip via BFS on cell adj. */
            typename conn_list::data_type _cell ;
            _list._pop_tail(_cell) ;

            if (_cell._kind == TRIA3_tag)
            {
            for(auto _enum = +3; _enum-- != +0; )
            {
                auto _cpos = _cell._cell;

                iptr_type _enod[3] ;
                mesh_type::tri3_type::
                face_node(_enod, _enum, +2, +1) ;
                _enod[ 0] = _mesh.tri3(
                    _cpos).node( _enod[ 0]) ;
                _enod[ 1] = _mesh.tri3(
                    _cpos).node( _enod[ 1]) ;

                flip_next( _mesh, _list,
                    _cpos, _enod, _seen, _conn) ;
            }
            }
            else
            if (_cell._kind == QUAD4_tag)
            {
            for(auto _enum = +4; _enum-- != +0; )
            {
                auto _cpos = _cell._cell;

                iptr_type _enod[4] ;
                mesh_type::quad_type::
                face_node(_enod, _enum, +2, +1) ;
                _enod[ 0] = _mesh.quad(
                    _cpos).node( _enod[ 0]) ;
                _enod[ 1] = _mesh.quad(
                    _cpos).node( _enod[ 1]) ;

                flip_next( _mesh, _list,
                    _cpos, _enod, _seen, _conn) ;
            }
            }
        }
    }

    __static_call
    __normal_call void_type flip_sign (
        mesh_type &_mesh
        )
    {
    #   define SEENTRI3( _CELL) _seen._tri3[_CELL]
    #   define SEENQUAD( _CELL) _seen._quad[_CELL]

        conn_list _list;
        mark_list _seen;
        init_mark(_mesh, _seen, +0) ;

    /*------------------------- flip to correct mesh sign */
        iptr_type _tnum  = +0 ;
        for (auto _tria  = _mesh.tri3().head();
                  _tria != _mesh.tri3().tend();
                ++_tria, ++_tnum )
        {
            if (_tria->mark  () < +0) continue;
            if (SEENTRI3(_tnum) > +0) continue;

        /*--------------------- flip seed for +ve quality */
            real_type _cost =
                 pred_type::tri3_cost (
               &_mesh. node(
                _tria->node(0)).pval(0),
               &_mesh. node(
                _tria->node(1)).pval(0),
               &_mesh. node(
                _tria->node(2)).pval(0),
                typename
                pred_type::cell_kind()) ;

            if (_cost < (real_type) +0.)
            {
                std::swap (
                    _tria->node(0),
                        _tria->node(1)) ;
            }

        /*--------------------- a BFS from seed via topo. */
            _list.push_tail(typename conn_list
                ::data_type (_tnum, TRIA3_tag)) ;

            SEENTRI3( _tnum ) = 1 ;

            flip_bfs_(_mesh, _list, _seen) ;
        }

        iptr_type _qnum  = +0 ;
        for (auto _quad  = _mesh.quad().head();
                  _quad != _mesh.quad().tend();
                ++_quad, ++_qnum )
        {
            if (_quad->mark  () < +0) continue;
            if (SEENQUAD(_qnum) > +0) continue;

        /*--------------------- flip seed for +ve quality */
            real_type _cost =
                 pred_type::quad_cost (
               &_mesh. node(
                _quad->node(0)).pval(0),
               &_mesh. node(
                _quad->node(1)).pval(0),
               &_mesh. node(
                _quad->node(2)).pval(0),
               &_mesh. node(
                _quad->node(3)).pval(0),
                typename
                pred_type::cell_kind()) ;

            if (_cost < (real_type) +0.)
            {
                // TODO: reverse quad indexes
            }

        /*--------------------- a BFS from seed via topo. */
            _list.push_tail(typename conn_list
                ::data_type (_qnum, QUAD4_tag)) ;

            SEENQUAD( _qnum ) = 1 ;

            flip_bfs_(_mesh, _list, _seen) ;
        }
    #   undef   SEENQUAD
    #   undef   SEENTRI3
    }

    /*
    --------------------------------------------------------
     * MOVE-OKAY: TRUE if state sufficiently good.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type move_okay (
        real_list &_cdst ,
        real_list &_csrc ,
        iptr_type &_move ,
        real_type  _good = +9.25E-01,
        real_type  _qtol = +1.00E-04,
        real_type  _xdel = +0.00E+00,
        real_type  _xtol =
    std::numeric_limits<real_type>::infinity()
        )
    {
        _move = (iptr_type) -1;

        if (_cdst.empty()) return ;
        if (_csrc.empty()) return ;

    /*--------------------- calc. min. + mean metrics */
        real_type _zero =
            +std::numeric_limits
                <real_type>::epsilon ();
        real_type _0src =
            +std::numeric_limits
                <real_type>::infinity();
        real_type _0dst =
            +std::numeric_limits
                <real_type>::infinity();

    /*--------------------- compute harmonic averages */
        _good = std::pow(_good, +7./8.);

        real_type _msrc, _mdst;
        _msrc = (real_type) +0. ;
        _mdst = (real_type) +0. ;

        for (auto _iter  = _csrc.head(),
                  _tend  = _csrc.tend();
                  _iter != _tend;
                ++_iter  )
        {
            _0src  =
            std::min(_0src, *_iter) ;

            _msrc += std::pow(
            (real_type)1. / *_iter, +7);
        }
        for (auto _iter  = _cdst.head(),
                  _tend  = _cdst.tend();
                  _iter != _tend;
                ++_iter  )
        {
            _0dst  =
            std::min(_0dst, *_iter) ;

            _mdst += std::pow (
            (real_type)1. / *_iter, +7);
        }

        _msrc  = std::pow(
        _csrc.count() / _msrc, +1./7.0);
        _mdst  = std::pow(
        _cdst.count() / _mdst, +1./7.0);

        _qtol /=
        std::pow(_csrc.count(), 1./7.0);
        _qtol /=
        std::pow(_cdst.count(), 1./7.0);

    /*---------------------------- test move = "okay" */
        if (_0dst >= _good)
        {
    /*--------------------- okay if moves unconverged */
            if (_xdel > _xtol)
                _move = +1;
            else
                _move = +0;     // no iter; too short

            if (_move > +0) return ;
        }

        if (_0dst >= _zero)
        {
    /*--------------------- okay if mean is improving */
            if (_mdst > _msrc+_qtol)
                _move = +1;

            if (_move > +0) return ;
        }

    }

    /*
    --------------------------------------------------------
     * LOOP-COST: cost vector for 1-neighbourhood.
    --------------------------------------------------------
     */

    __static_call
    __normal_call real_type loop_cost (
        mesh_type &_mesh ,
        conn_list &_cset ,
        real_list &_cost ,
        cell_kind const&        // cell costs on CSET
        )
    {
        real_type _qmin = (real_type)1. ;

        for (auto _next  = _cset.head() ,
                  _tend  = _cset.tend() ;
                  _next != _tend;
                ++_next  )
        {
            real_type _cscr = (real_type)+1. ;

            iptr_type _cell = _next->_cell ;
            iptr_type _kind = _next->_kind ;

            switch (_kind)
            {
            case TRIA3_tag:
                {
           _cscr = pred_type::tri3_cost (
               &_mesh .node(
                _mesh .tri3(
                _cell).node(0)).pval(0),
               &_mesh .node(
                _mesh .tri3(
                _cell).node(1)).pval(0),
               &_mesh .node(
                _mesh .tri3(
                _cell).node(2)).pval(0),
                typename
                pred_type::cell_kind ()) ;
                break ;
                }

            case QUAD4_tag:
                {
           _cscr = pred_type::quad_cost (
               &_mesh .node(
                _mesh .quad(
                _cell).node(0)).pval(0),
               &_mesh .node(
                _mesh .quad(
                _cell).node(1)).pval(0),
               &_mesh .node(
                _mesh .quad(
                _cell).node(2)).pval(0),
               &_mesh .node(
                _mesh .quad(
                _cell).node(3)).pval(0),
                typename
                pred_type::cell_kind ()) ;
                break ;
                }
            }

            _qmin =
            std::min (_qmin, _cscr) ;

            _cost.push_tail (_cscr) ;
        }

        return ( _qmin )  ;
    }

    __static_call
    __normal_call real_type loop_cost (
        mesh_type &_mesh ,
        conn_list &_cset ,
        real_list &_cost ,
        dual_kind const&        // dual costs on CSET
        )
    {
        real_type _qmin = (real_type)1. ;

        for (auto _next  = _cset.head() ,
                  _tend  = _cset.tend() ;
                  _next != _tend;
                ++_next  )
        {
            real_type _cscr = (real_type)+1. ;

            iptr_type _cell = _next->_cell ;
            iptr_type _kind = _next->_kind ;

            switch (_kind)
            {
            case TRIA3_tag:
                {
           _cscr = pred_type::tri3_cost (
               &_mesh .node(
                _mesh .tri3(
                _cell).node(0)).pval(0),
               &_mesh .node(
                _mesh .tri3(
                _cell).node(1)).pval(0),
               &_mesh .node(
                _mesh .tri3(
                _cell).node(2)).pval(0),
                typename
                pred_type::dual_kind ()) ;
                break ;
                }

            case QUAD4_tag:
                {
           _cscr = pred_type::quad_cost (
               &_mesh .node(
                _mesh .quad(
                _cell).node(0)).pval(0),
               &_mesh .node(
                _mesh .quad(
                _cell).node(1)).pval(0),
               &_mesh .node(
                _mesh .quad(
                _cell).node(2)).pval(0),
               &_mesh .node(
                _mesh .quad(
                _cell).node(3)).pval(0),
                typename
                pred_type::dual_kind ()) ;
                break ;
                }
            }

            _qmin =
            std::min (_qmin, _cscr) ;

            _cost.push_tail (_cscr) ;
        }

        return ( _qmin )  ;
    }

    /*
    --------------------------------------------------------
     * MOVE-NODE: "limited" single node coord. update.
    --------------------------------------------------------
     */

    #include "iter_node_2.inc"

    template <
        typename  node_iter
             >
    __static_call
    __inline_call void_type move_node (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iter_opts &_opts ,
        node_iter  _node ,
        kern_kind  _kern ,
        iptr_type &_move ,
        conn_list &_conn ,
        real_list &_qold ,
        real_list &_qnew ,
        real_type  _QMIN ,
        real_type  _QLIM ,
        iter_stat &_tcpu
        )
    {
    /*---------------- try variational; fallback on dQ/dx */
        move_kern( _geom, _mesh, _hfun,
            _hval, _opts, _node,
            _kern, _move, _conn,
            _qold, _qnew,
            _QMIN, _QLIM, _tcpu) ;  // variational

        if (_move >= +0 ) return ;

        move_kern( _geom, _mesh, _hfun,
            _hval, _opts, _node,
             dqdx_optimise,
            _move, _conn,
            _qold, _qnew,
            _QMIN, _QLIM, _tcpu) ;  // local dQ/dx

        if (_move >= +0 ) return ;
    }

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type move_kern (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iter_opts &_opts ,
        node_iter  _node ,
        kern_kind  _kern ,
        iptr_type &_move ,
        conn_list &_conn ,
        real_list &_qold ,
        real_list &_qnew ,
        real_type  _QMIN ,
        real_type  _QLIM ,
        iter_stat &_tcpu
        )
    {
    /*---------------- optimise single node's coordinates */
        iptr_type static
        constexpr _ITER = (iptr_type) +4  ;

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic ;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

        _move = (iptr_type)-1 ;

        real_type _line[geom_dims] =
            {(real_type)+0.0} ;

        real_type _save[geom_dims] =
            {(real_type)+0.0} ;

        real_type _proj[geom_dims] =
            {(real_type)+0.0} ;

    /*---------------- calc. line search direction vector */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        real_type _ladj = (real_type) +0. ;

        if (_kern == _odt_optimise)
        {
    /*--------------------------- ODT-style update vector */
            _odt_move_2 (
                _geom, _mesh, _hfun, _hval,
                _conn, _node,
                _line, _ladj) ;
        }
        else
        if (_kern == _cvt_optimise)
        {
    /*--------------------------- CVT-style update vector */
            _cvt_move_2 (
                _geom, _mesh, _hfun, _hval,
                _conn, _node,
                _line, _ladj) ;
        }
        else
        if (_kern == dqdx_optimise)
        {
            if (_QMIN<=_QLIM)
            {
    /*--------------------------- d./dx Q^T update vector */
            dqdx_move_2 (
                _mesh, _conn, _node, _qold,
                _line, _ladj) ;
            }
            else { return ; }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._ldir_node += _tcpu.nano_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*---------------- scale line search direction vector */
        real_type _xeps =           // delta_x ~= 0.0
       (real_type)+0.01*_opts.qtol() ;

        real_type _xtol =           // delta_x reltol
       +std::sqrt(_opts.qtol()) / +10.0  ;

        if (_kern == dqdx_optimise) // test cost-only
        {
            _QLIM  =
       +std::numeric_limits<real_type>::infinity() ;
        }

        auto _ppos = &_node->pval(0) ;

        real_type _lsqr ;
        _lsqr = std::pow(_ladj, 2) ;
        _xeps = std::pow(_xeps, 2) ;
        _xtol = std::pow(_xtol, 2) ;

        real_type _scal =           // overrelaxation
            (real_type) std::sqrt( 2.0 ) ;

    /*---------------- do backtracking line search iter's */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _idim =
        pred_type::geom_dims; _idim-- != +0; )
        {
            _save[_idim] = _ppos [_idim] ;
        }

        for (auto _iter = +0 ;
                _iter != _ITER; ++_iter )
        {
    /*---------------- push update along search direction */
            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _proj[_idim] =
                    _save[_idim] +
                        _scal* _line[_idim] ;
            }

            pred_type::
            proj_node (_geom, _save, _proj) ;

            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _ppos[_idim] = _proj[_idim] ;
            }

            real_type _XEPS = _xeps * _scal ;
            real_type _XTOL = _xtol * _scal ;

            real_type _lmov =
            pred_type::length_sq(_save, _proj) ;

            if (_lmov <= _XEPS * _lsqr) break;

          //_move  = +1 ; return ;

            _scal *= (real_type)+0.5 ;

    /*---------------- test quasi-monotonicity w.r.t. Q^T */
            _qnew.set_count(0) ;

            loop_cost( _mesh, _conn, _qnew,
                 cell_kind ()) ;

            move_okay( _qnew, _qold, _move,
                _QLIM, _opts.qtol(),
                _lmov, _XTOL* _lsqr) ;

            if (_move >= +0) break ;
        }

    /*---------------- swap with the saved coord. if fail */
        if (_move <= +0)
        {
            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _ppos[_idim] = _save[_idim] ;
            }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._lopt_node += _tcpu.nano_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    }

    /*
    --------------------------------------------------------
     * MOVE-DUAL: "limited" single node weight update.
    --------------------------------------------------------
     */

   #include "iter_dual_2.inc"

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type move_dual (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iter_opts &_opts ,
        node_iter  _node ,
        iptr_type &_move ,
        conn_list &_conn ,
        real_list &_dold ,
        real_list &_dnew ,
        real_type  _DMIN ,
        real_type  _DLIM ,
        iter_stat &_tcpu
        )
    {
    /*---------------- optimise single node's coordinates */
        iptr_type static
        constexpr _ITER = (iptr_type) +4  ;

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic ;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

        __unreferenced(_geom);
        __unreferenced(_hfun);
        __unreferenced(_hval);

        _move = (iptr_type)-1;

        real_type _wadj, _step, _save;

    /*---------------- calc. line search direction vector */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        if (_DMIN <= _DLIM)
        {
            dqdw_move_2 (
                _mesh, _conn, _node, _dold,
                _step, _wadj) ;
        }
        else { return ; }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._ldir_dual += _tcpu.nano_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*---------------- scale line search direction vector */
        real_type _weps =           // delta_w ~= 0.0
            (real_type) .01*_opts.qtol() ;

        _save = _node->pval(
            pred_type::real_dims - 1);

        real_type _scal =           // overrelaxation
            (real_type) std::sqrt( 2.0 ) ;

    /*---------------- do backtracking line search iter's */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _iter = +0 ;
                _iter != _ITER; ++_iter )
        {
            _node->pval(real_dims-1) =
                _save + ( _scal * _step ) ;

            _node->pval(real_dims-1) =
                std::max(-_wadj ,
            _node->pval(real_dims-1));

            _node->pval(real_dims-1) =
                std::min(+_wadj ,
            _node->pval(real_dims-1));

            real_type _wmov =
                std::abs (_save -
            _node->pval(real_dims-1));

            if (_wmov <=
            _weps * _scal * _wadj) break;

            _scal *= (real_type)+0.5 ;

    /*---------------- test quasi-monotonicity w.r.t. Q^D */
            _dnew.set_count(0) ;

            loop_cost( _mesh, _conn, _dnew,
                 dual_kind ()) ;

            move_okay( _dnew, _dold, _move,
                _DLIM, _opts.qtol()) ;

            if (_move >= +0) break ;
        }

    /*---------------- swap with the saved coord. if fail */
        if (_move <= +0)
        {
            _node->pval(real_dims-1) = _save ;
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._lopt_dual += _tcpu.nano_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    }

    /*
    --------------------------------------------------------
     * PULL-CONN: cache node-to-cell adjacency lists.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type pull_conn (
        mesh_type &_mesh ,
        conn_sets &_conn
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

        _conn._adj1.set_alloc(
            _mesh.edge().count() * 2) ;
        _conn._adj2.set_alloc(
            _mesh.tri3().count() * 3 +
            _mesh.quad().count() * 4) ;

        iptr_type _npos  = +0 ;
        for (auto _node  = _mesh.node().head() ;
                  _node != _mesh.node().tend() ;
                ++_node, ++_npos)
        {
        if (_node->mark() >= + 0)
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
        mesh_type &_mesh ,
        conn_sets &_CONN ,
        iptr_list &_lset ,
        iptr_list &_aset ,
        iptr_list &_nmrk ,
        iptr_list &_amrk ,
        iptr_type  _iout ,
        iptr_type  _isub ,
        real_type  _QLIM ,
        real_type  _DLIM ,
        iter_opts &_opts
        )
    {
    #   define PUSHSORT(_NODE)  \
            if (_amrk[_NODE] != _isub)      \
            {               \
                _amrk[_NODE]  = _isub;      \
                _sset.push_tail(            \
                 cost_pair(_NODE, 1.)) ;    \
            }

    #   define PUSHCONN(_NODE)  \
            if (_amrk[_NODE] != _isub)      \
            {               \
                _amrk[_NODE]  = _isub;      \
                _aset.push_tail(_NODE) ;    \
            }

    #   define PUSHCOST(_NODE , _COST)      \
            _qbar[_NODE] += _COST ;         \
            _nadj[_NODE] += +1 ;\
            _qmin[_NODE]  =     \
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

        real_list _qbar, _qmin ;
        iptr_list _nadj ;
        cost_list _sset ;

        __unreferenced ( _DLIM ) ;
        __unreferenced ( _opts ) ;

        if (_isub == (iptr_type) +0)
        {
    /*-------------------- 1ST SUB-ITER: build full init. */
        _qbar.set_count(
            _mesh.node().count(),
                containers::tight_alloc, +0.0) ;
        _qmin.set_count(
            _mesh.node().count(),
                containers::tight_alloc, +1.0) ;
        _nadj.set_count(
            _mesh.node().count(),
                containers::tight_alloc, + 0 ) ;

        _sset.set_alloc(
            _mesh.node().count())  ;

        pull_conn(_mesh, _CONN) ;

        for (auto _tria  = _mesh.tri3().head() ;
                  _tria != _mesh.tri3().tend() ;
                ++_tria  )
        {
            if (_tria->mark() >= +0)
            {
        /*-------------------- calc. min. scores at nodes */
                iptr_type _inod, _jnod, _knod;
                _inod = _tria->node( 0);
                _jnod = _tria->node( 1);
                _knod = _tria->node( 2);

                real_type _cost;
                _cost = pred_type::tri3_cost (
               &_mesh .node(
                _tria->node(0)).pval(0),
               &_mesh .node(
                _tria->node(1)).pval(0),
               &_mesh .node(
                _tria->node(2)).pval(0),
                typename
                pred_type::cell_kind ()) ;

                PUSHSORT( _inod )
                PUSHCOST( _inod, _cost )
                PUSHSORT( _jnod )
                PUSHCOST( _jnod, _cost )
                PUSHSORT( _knod )
                PUSHCOST( _knod, _cost )
            }
        }

        for (auto _quad  = _mesh.quad().head() ;
                  _quad != _mesh.quad().tend() ;
                ++_quad  )
        {
            if (_quad->mark() >= +0)
            {
        /*-------------------- calc. min. scores at nodes */
                iptr_type _inod, _jnod, _knod,
                          _lnod;
                _inod = _quad->node( 0);
                _jnod = _quad->node( 1);
                _knod = _quad->node( 2);
                _lnod = _quad->node( 3);

                real_type _cost;
                _cost = pred_type::quad_cost (
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

                PUSHSORT( _inod )
                PUSHCOST( _inod, _cost )
                PUSHSORT( _jnod )
                PUSHCOST( _jnod, _cost )
                PUSHSORT( _knod )
                PUSHCOST( _knod, _cost )
                PUSHSORT( _lnod )
                PUSHCOST( _lnod, _cost )
            }
        }

        for (auto _iter  =
            _qbar.count(); _iter-- != +0; )
        {
    /*------------------------ assign score for each node */
            _qbar[_iter] /=
                std::max(+1, _nadj [_iter]) ;
        }

        for (auto _iter  = _sset.head() ;
                  _iter != _sset.tend() ;
                ++_iter  )
        {
            _iter->_cost =
            (float)_qbar [_iter->_node] -
            (float)_qmin [_iter->_node] ;
        }

        algorithms::qsort( _sset.head() ,
                           _sset.tend() ,
            cost_pred () ) ;

        iptr_type _FLAG  = _iout - 2 ; // only "recent"

        for (auto _iter  = _sset.head() ;
                  _iter != _sset.tend() ;
                ++_iter  )
        {
    /*------------------------ push sorted wrt. min.-cost */
            if (std::abs(
            _nmrk[_iter->_node]) >= _FLAG || // recent
            _qmin[_iter->_node]  <= _QLIM )  // skewed
            {
            _amrk[_iter->_node]   = _isub ;
            _aset.push_tail( _iter->_node ) ;
            }
        }

        }
        else
        {
    /*-------------------- NTH sub-iter: init. from prev. */
        for (auto _iter  = _lset.head() ;
                  _iter != _lset.tend() ;
                ++_iter  )
        {
                PUSHCONN( *_iter )
        }

    /*-------------------- add adj.: 1-ring neighbourhood */
        for (auto _iter =
            _aset.count() ; _iter-- != 0; )
        {
             auto _apos =   _aset [_iter] ;

        /*-------------------- push any 1-cell neighbours */
            for (auto _next  =
            _CONN._adj1.head()+_CONN._idx1[_apos+0] ;
                      _next !=
            _CONN._adj1.head()+_CONN._idx1[_apos+1] ;
                    ++_next  )
            {
            if (_next->_kind == EDGE2_tag)
            {
                 auto _eptr =
                _mesh. edge().head()+_next->_cell ;

                PUSHCONN( _eptr->node(0) )
                PUSHCONN( _eptr->node(1) )
            }
            }

        /*-------------------- push any 2-cell neighbours */
            for (auto _next  =
            _CONN._adj2.head()+_CONN._idx2[_apos+0] ;
                      _next !=
            _CONN._adj2.head()+_CONN._idx2[_apos+1] ;
                    ++_next  )
            {
            if (_next->_kind == TRIA3_tag)
            {
                 auto _tptr =
                _mesh. tri3().head()+_next->_cell ;

                PUSHCONN( _tptr->node(0) )
                PUSHCONN( _tptr->node(1) )
                PUSHCONN( _tptr->node(2) )
            }
            else
            if (_next->_kind == QUAD4_tag)
            {
                 auto _qptr =
                _mesh. quad().head()+_next->_cell ;

                PUSHCONN( _qptr->node(0) )
                PUSHCONN( _qptr->node(1) )
                PUSHCONN( _qptr->node(2) )
                PUSHCONN( _qptr->node(3) )
            }
            }
        }

        }
    #   undef   PUSHCOST
    #   undef   PUSHSORT
    #   undef   PUSHCONN
    }

    /*
    --------------------------------------------------------
     * MOVE-NODE: do a single node smoothing pass.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type move_node (
        geom_type &_geom ,
        mesh_type &_mesh ,
        conn_sets &_CONN ,
        hfun_type &_hfun ,
        char_type  _kern ,
        real_list &_hval ,
        iptr_list &_nset ,
        iptr_list &_lset ,
        iptr_list &_aset ,
        iptr_list &_amrk ,
        mark_list &_mark ,
        iptr_type  _iout ,
        iptr_type  _isub ,
        iter_opts &_opts ,
        iptr_type &_nmov ,
        real_type  _QLIM ,
        real_type  _DLIM ,
        iter_stat &_tcpu
        )
    {
    #   define MARK(_NODE) _mark._node[_NODE]

    #   define HEAD(_PASS)(_PASS == +0 ?    \
                           _aset.head() :   \
                           _aset.tail())

    #   define STOP(_PASS)(_PASS == +0 ?    \
                           _aset.tend() :   \
                           _aset.hend())

    #   define ITER(_PASS)(_PASS == +0 ? +1 : -1)

        conn_list _conn;
        real_list _qold, _qnew, _dold, _dnew;

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*-------------------- permute nodes for optimisation */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        _aset.set_count(+0) ;

        sort_node(_mesh, _CONN, _lset, _aset,
            _mark._node, _amrk, _iout, _isub,
            _QLIM,_DLIM, _opts) ;

        _nmov = (iptr_type)0;

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*-------------------- SYMMETRIC GAUSS-SEIDEL on CELL */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        _lset.set_count(+0) ;

        if (_opts .tria())  // 1=>GS, 2=>SGS
        {
        for (auto _pass  = 0; _pass < 2; ++_pass)
        for (auto _apos  = HEAD(_pass) ;
                  _apos != STOP(_pass) ;
                  _apos += ITER(_pass) )
        {
             auto _node  =
            _mesh.node().head()+ *_apos ;

        /*---------------- assemble a local tria. stencil */
            _conn.set_count( +0) ;
            _conn.push_tail(
            _CONN._adj2.head()+_CONN._idx2[*_apos+0],
            _CONN._adj2.head()+_CONN._idx2[*_apos+1]
                ) ;

            if (_conn.empty()) continue ;

            if (_mark._node[*_apos] >= +0)
            {
        /*---------------- attempt to optimise CELL geom. */
            _qold.set_count( +0) ;
            _qnew.set_count( +0) ;

            real_type  _QMIN =
            loop_cost( _mesh,
                _conn, _qold, cell_kind());

            iptr_type _move = -1 ;

            if (_move < +0 )
            {
        /*---------------- do optimisation of node coord. */
                move_node( _geom, _mesh,
                    _hfun, _hval,
                    _opts, _node, _kern,
                    _move, _conn,
                    _qold, _qnew,
                    _QMIN, _QLIM, _tcpu) ;
            }

            if (_move > +0 )
            {
        /*---------------- update when state is improving */
            _hval[*_apos] = (real_type)-1;

            if (std::abs(
                MARK( *_apos )) !=_iout)
            {
                if (MARK( *_apos ) >= 0)
                MARK( *_apos ) = +_iout;
                else
                MARK( *_apos ) = -_iout;

                _nset.push_tail (*_apos) ;
            }

            _lset.push_tail (*_apos) ;

            _nmov += +1 ;
            }
            }
        }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    #   undef   MARK
    #   undef   HEAD
    #   undef   STOP
    #   undef   ITER
    }

    /*
    --------------------------------------------------------
     * MOVE-DUAL: do a single dual smoothing pass.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type move_dual (
        geom_type &_geom ,
        mesh_type &_mesh ,
        conn_sets &_CONN ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iptr_list &_nset ,
        iptr_list &_lset ,
        iptr_list &_aset ,
        iptr_list &_amrk ,
        mark_list &_mark ,
        iptr_type  _iout ,
        iptr_type  _isub ,
        iter_opts &_opts ,
        iptr_type &_nmov ,
        real_type  _QLIM ,
        real_type  _DLIM ,
        iter_stat &_tcpu
        )
    {
    #   define MARK(_NODE) _mark._node[_NODE]

    #   define HEAD(_PASS)(_PASS == +0 ?    \
                           _aset.head() :   \
                           _aset.tail())

    #   define STOP(_PASS)(_PASS == +0 ?    \
                           _aset.tend() :   \
                           _aset.hend())

    #   define ITER(_PASS)(_PASS == +0 ? +1 : -1)

        conn_list _conn;
        real_list _qold, _qnew, _dold, _dnew;

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic ;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*-------------------- permute nodes for optimisation */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        _aset.set_count(+0) ;

        sort_node(_mesh, _CONN, _lset, _aset,
            _mark._node, _amrk, _iout, _isub,
            _QLIM,_DLIM, _opts) ;

        _nmov = (iptr_type)0;

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*-------------------- SYMMETRIC GAUSS-SEIDEL on DUAL */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        _lset.set_count(+0) ;

        if (_opts .dual())  // 1=>GS, 2=>SGS
        {
        for (auto _pass  = 0; _pass < 2; ++_pass)
        for (auto _apos  = HEAD(_pass) ;
                  _apos != STOP(_pass) ;
                  _apos += ITER(_pass) )
        {
             auto _node  =
            _mesh.node().head()+ *_apos ;

        /*---------------- assemble a local tria. stencil */
            _conn.set_count( +0) ;
            _conn.push_tail(
            _CONN._adj2.head()+_CONN._idx2[*_apos+0],
            _CONN._adj2.head()+_CONN._idx2[*_apos+1]
                ) ;

            if (_conn.empty()) continue ;

        /*---------------- attempt to optimise DUAL geom. */
            _dold.set_count( +0) ;
            _dnew.set_count( +0) ;

            real_type  _DMIN =
            loop_cost( _mesh,
                _conn, _dold, dual_kind());

            iptr_type _move = -1 ;

            if (_move < +0 )
            {
        /*---------------- do optimisation of node weight */
                move_dual( _geom, _mesh,
                    _hfun, _hval,
                    _opts, _node,
                    _move, _conn,
                    _dold, _dnew,
                    _DMIN, _DLIM, _tcpu) ;
            }

            if (_move > +0 )
            {
        /*---------------- update when state is improving */
            if (std::abs(
                MARK( *_apos )) !=_iout)
            {
                if (MARK( *_apos ) >= 0)
                MARK( *_apos ) = +_iout;
                else
                MARK( *_apos ) = -_iout;

                _nset.push_tail (*_apos) ;
            }

            _lset.push_tail(*_apos) ;

            _nmov += +1 ;
            }
        }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    #   undef   MARK
    #   undef   HEAD
    #   undef   STOP
    #   undef   ITER
    }

    /*
    --------------------------------------------------------
     * FLIP-MESH: "flip" mesh topology.
    --------------------------------------------------------
     */

    #include "iter_flip_2.inc"

    __static_call
    __normal_call void_type flip_mesh (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        iptr_list &_nset ,
        mark_list &_mark ,
        iptr_type  _imrk ,
        iptr_type &_nflp
        )
    {
    #   define MARKTRI3(_CELL) \
                _mark._tri3[_CELL->_cell]

    #   define MARKQUAD(_CELL) \
                _mark._quad[_CELL->_cell]

        init_mark(_mesh, _mark ,
            std::max(+0, _imrk - 1) ) ;

        _nflp = +0 ;

        conn_list _flip, _conn ;

    /*--------------------- init. flip stack as ADJ(NSET) */
        for (auto _iter  = _nset.head();
                  _iter != _nset.tend();
                ++_iter  )
        {
            if (_mesh.node(*_iter).mark() >= 0)
            {
                _conn.set_count(0);
                _mesh.connect_2(
                  &*_iter, POINT_tag, _conn) ;

                for (auto _cell  = _conn.head();
                          _cell != _conn.tend();
                        ++_cell  )
                {
                if (_cell->_kind == TRIA3_tag)
                {
                if (MARKTRI3( _cell ) != _imrk)
                {
                    MARKTRI3( _cell )  = _imrk ;

                    _flip.push_tail( *_cell) ;
                }
                }
                else
                if (_cell->_kind == QUAD4_tag)
                {
                if (MARKQUAD( _cell ) != _imrk)
                {
                    MARKQUAD( _cell )  = _imrk ;

                    _flip.push_tail( *_cell) ;
                }
                }
                }
            }
        }

        flip_core(_geom, _mesh, _hfun, _flip, _mark,
            _nflp) ;

    #   undef   MARKQUAD
    #   undef   MARKTRI3
    }

    __static_call
    __normal_call void_type flip_mesh (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        conn_sets &_CONN ,
        iptr_list &_nset ,
        mark_list &_mark ,
        iptr_type  _imrk ,
        iptr_type &_nflp
        )
    {
    #   define MARKTRI3(_CELL) \
                _mark._tri3[_CELL->_cell]

    #   define MARKQUAD(_CELL) \
                _mark._quad[_CELL->_cell]

        init_mark(_mesh, _mark ,
            std::max(+0, _imrk - 1) ) ;

        _nflp = +0 ;

        conn_list _flip  ;

    /*--------------------- init. flip stack as ADJ(NSET) */
        for (auto _iter  = _nset.head();
                  _iter != _nset.tend();
                ++_iter  )
        {
            if (_mesh.node(*_iter).mark() >= 0)
            {
                for (auto _cell  =
            _CONN._adj2.head()+_CONN._idx2[*_iter+0] ;
                          _cell !=
            _CONN._adj2.head()+_CONN._idx2[*_iter+1] ;
                        ++_cell  )
                {
                if (_cell->_kind == TRIA3_tag)
                {
                if (MARKTRI3( _cell ) != _imrk)
                {
                    MARKTRI3( _cell )  = _imrk ;

                    _flip.push_tail( *_cell) ;
                }
                }
                else
                if (_cell->_kind == QUAD4_tag)
                {
                if (MARKQUAD( _cell ) != _imrk)
                {
                    MARKQUAD( _cell )  = _imrk ;

                    _flip.push_tail( *_cell) ;
                }
                }
                }
            }
        }

        flip_core(_geom, _mesh, _hfun, _flip, _mark,
            _nflp) ;

    #   undef   MARKQUAD
    #   undef   MARKTRI3
    }

    __static_call
    __normal_call void_type flip_core (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        conn_list &_flip ,
        mark_list &_mark ,
        iptr_type &_nflp
        )
    {
        __unreferenced ( _hfun ) ;

    /*--------------------- exhaustive, incremental flips */
        conn_list _next, _conn, _CONN;
        real_list _qold, _qnew ;
        for ( ; !_flip.empty() ; )
        {
        for (auto _cell  = _flip.head();
                  _cell != _flip.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
            if (  _mesh.
                tri3(_cell->_cell).mark() >= +0)
            {
                bool_type  _okay = false ;
                flip_tri3( _geom, _mesh,
                    _mark. _node,
                    _cell->_cell,
                    _okay, _conn, _CONN,
                    _qold, _qnew,
               (flip_kind)laguerre_flip) ;

                if (_okay) _nflp += +1 ;

                _next.push_tail(
                    _CONN.head(), _CONN.tend());
            }
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
            if (  _mesh.
                quad(_cell->_cell).mark() >= +0)
            {
                bool_type  _okay = false ;
              //flip_quad( _geom, _mesh,
              //    _mark. _node,
              //    _cell->_cell,
              //    _okay, _conn, _CONN,
              //    _qold, _qnew,
            // (flip_kind)laguerre_flip) ;

                if (_okay) _nflp += +1 ;

                _next.push_tail(
                    _CONN.head(), _CONN.tend());
            }
            }
        }
            _flip = std::move (_next) ;
        }
    }

    /*
    --------------------------------------------------------
     * _ZIP-MESH: edge merge/split operations.
    --------------------------------------------------------
     */

    #include "iter_zips_2.inc"
    #include "iter_divs_2.inc"

    __static_call
    __normal_call void_type _zip_mesh (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        char_type  _kern ,
        real_list &_hval ,
        iptr_list &_nset ,
        mark_list &_mark ,
        iptr_type  _imrk ,
        iter_opts &_opts ,
        real_type  _QLIM ,
        real_type  _DLIM ,
        iptr_type &_nzip ,
        iptr_type &_ndiv ,
        iter_stat &_tcpu
        )
    {
        class sort_pair
            {
            public  :
    /*------------------------ tuple for edge re-ordering */
            iptr_type           _inod ;
            iptr_type           _jnod ;
            float               _cost ;
            public  :
            __inline_call sort_pair (
                iptr_type _isrc ,
                iptr_type _jsrc ,
                float     _csrc
                ) : _inod(_isrc), _jnod(_jsrc),
                    _cost(_csrc)   {}
            } ;
        class sort_less
            {
            public  :
    /*------------------------ less-than op. for cost-tup */
            __inline_call
                bool_type operator () (
                sort_pair const&_idat ,
                sort_pair const&_jdat
                ) const
            {  return  _idat._cost < _jdat._cost;
            }
            } ;

    #   define MARKNODE(_NODE) _mark._node[_NODE]

    #   define PUSHMARK              \
            init_mark( _mesh, _mark,        \
                std::max(_imrk - 1, +0));   \
            if (std::abs(           \
                _mark._node[_nnew])!= _imrk+1)  \
            {                       \
            if (_mark._node[_nnew] >= +0)   \
            {                       \
                _mark._node[_nnew] = +_imrk+1;  \
            }                       \
            else                    \
            {                       \
                _mark._node[_nnew] = -_imrk-1;  \
            }                       \
                _nset.push_tail(_nnew) ;    \
            }                       \

    /*------------------------ add/pop nodes to fix topo. */
        typedef containers::
            array<sort_pair> sort_list ;

        __unreferenced ( _DLIM );

        sort_list _sort ;
        conn_list _aset , _bset, _cset ;
        conn_list _conn , _iset, _jset ;
        real_list _qold , _qnew, _qtmp ;

        _nzip = +0 ; _ndiv = +0 ;

    //  assemble list of edges attached to "recent" nodes

        _sort.set_alloc(
            _mesh.edge().count()) ;

        for (auto _iter  = _mesh.edge().head() ;
                  _iter != _mesh.edge().tend() ;
                ++_iter  )
        {
             auto _inod =  _iter->node(0) ;
             auto _jnod =  _iter->node(1) ;

             auto _iptr = _mesh.
             node().head()+_iter->node(0) ;
             auto _jptr = _mesh.
             node().head()+_iter->node(1) ;

             auto _flag = _imrk - 2 ; // recent

            if (_iter->mark() >= +0 &&
               (    std::abs(
                _mark._node[_inod]) > _flag ||
                    std::abs(
                _mark._node[_jnod]) > _flag ))
            {
                float _lsqr  =
               (float)pred_type::length_sq (
                    & _iptr->pval(0) ,
                    & _jptr->pval(0) ) ;

                _sort.push_tail(
                 sort_pair(_inod, _jnod, _lsqr)) ;
            }
        }

        if (_sort.empty()) return ;

        algorithms::qsort(      // sort edge list by lsqr
            _sort.head() ,
            _sort.tend() , sort_less());

    //  scan edges longest-to-shortest and try to div any
    //  unvisited edges

        for (auto _iter  = _sort.tail();
                  _iter != _sort.hend();
                --_iter  )
        {
    /*--------------------------- try to "div" local edge */
            iptr_type  _eadj, _enod[2] ;
            _enod[0] = _iter->_inod;
            _enod[1] = _iter->_jnod;

            if (MARKNODE(_enod[0]) < +0 &&
                MARKNODE(_enod[1]) < +0 ) continue ;

            if (MARKNODE(_enod[0])>_imrk||
                MARKNODE(_enod[1])>_imrk) continue ;

            if(!_mesh.find_edge(
                          _enod, _eadj) ) continue ;

            if (_opts.div_())
            {
    /*--------------------------- "div" for topo. + score */
                iptr_type  _nnew = -1;

                bool_type  _move = false ;
                _div_edge( _geom, _mesh,
                    _hfun, _hval, _opts,
                    _imrk, _eadj,
                    _kern, _move, _nnew,
                    _iset, _jset,
                    _aset, _bset,
                    _qold, _qnew,
                    _qtmp, _QLIM, _tcpu) ;

                if (_move)
                {
                    PUSHMARK; _ndiv += +1;
                }
            }
        }

    //  scan edges shortest-to-longest and try to zip any
    //  unvisited edges

        for (auto _iter  = _sort.head();
                  _iter != _sort.tend();
                ++_iter  )
        {
    /*--------------------------- try to "zip" local edge */
            iptr_type  _eadj, _enod[2] ;
            _enod[0] = _iter->_inod;
            _enod[1] = _iter->_jnod;

            if (MARKNODE(_enod[0]) < +0 ||
                MARKNODE(_enod[1]) < +0 ) continue ;

            if (MARKNODE(_enod[0])>_imrk||
                MARKNODE(_enod[1])>_imrk) continue ;

            if(!_mesh.find_edge(
                          _enod, _eadj) ) continue ;

            if (_opts.zip_())
            {
    /*--------------------------- "zip" for topo. + score */
                iptr_type  _nnew = -1;

                bool_type  _move = false ;
                _zip_edge( _geom, _mesh,
                    _hfun, _hval, _opts,
                    _mark. _node,
                    _imrk, _eadj,
                    _kern, _move, _nnew,
                    _iset, _jset,
                    _aset, _bset, _cset,
                    _qold, _qnew,
                    _qtmp, _QLIM, _tcpu) ;

                if (_move)
                {
                    PUSHMARK; _nzip += +1;
                }
            }
        }

        for (auto _iter  = _mark._node.head() ;
                  _iter != _mark._node.tend() ;
                ++_iter  )
        {
    /*--------------------- undo local inc. on node flags */
            if (*_iter > + _imrk)
                *_iter = + _imrk;
            else
            if (*_iter < - _imrk)
                *_iter = - _imrk;
        }

    #   undef   PUSHMARK
    #   undef   MARKNODE
    }

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
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        kern_kind  _kern ,
        iter_opts &_opts ,
        text_dump &_dump
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

        {
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

        flip_sign(_mesh) ;

    /*------------------------------ do optimisation loop */
        bool_type
        static constexpr ITER_FLIP = true;

        iptr_type
        static constexpr ITER_MIN_ = min_subit ;
        iptr_type
        static constexpr ITER_MAX_ = max_subit ;

        real_type _QMIN = (real_type) +1.;

        for (auto _cell  = _mesh.tri3().head() ;
                  _cell != _mesh.tri3().tend() ;
                ++_cell  )
        {
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

            _QMIN = std::min (_QMIN, _cost) ;
            _QMIN = std::max (_QMIN,
                              _opts.qlim()) ;
            }
        }

        for (auto _cell  = _mesh.quad().head() ;
                  _cell != _mesh.quad().tend() ;
                ++_cell  )
        {
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

            _QMIN = std::min (_QMIN, _cost) ;
            _QMIN = std::max (_QMIN,
                              _opts.qlim()) ;
            }
        }

        for (auto _iter = +1 ;
            _iter <= _opts.iter(); ++_iter)
        {
    /*------------------------------ set-up current iter. */
            init_mark(_mesh, _mark,
                std::max(_iter-1, +0)) ;

            real_list _hval;
            _hval.set_count(
                _mesh. node().count(),
        containers::tight_alloc, (real_type)-1.);

            _nset.set_count(  +0);

            iptr_type _nmov = +0 ;
            iptr_type _nflp = +0 ;
            iptr_type _nzip = +0 ;
            iptr_type _ndiv = +0 ;

    /*------------------------------ scale quality thresh */
            iptr_type _nsub = _iter +0 ;

            _nsub =
            std::min(ITER_MAX_, _nsub) ;
            _nsub =
            std::max(ITER_MIN_, _nsub) ;

            real_type _DLIM =
           (real_type)+1.-_opts.qtol() ;

            real_type _QLIM =
           (real_type)+.750 * _QMIN +
           (real_type)+.075 * _iter ;

            _QLIM = std::min(
                _opts.qlim () , _QLIM) ;

    /*------------------------------ 1. CELL GEOM. PASSES */

            if (_opts.tria())
            {
    /*------------------------------ update mesh geometry */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            iptr_list _amrk, _aset, _lset ;
            _amrk.set_count(
                _mesh.node().count() ,
                    containers::tight_alloc,-1) ;

            _lset.set_alloc(
                _mesh.node().count()) ;
            _aset.set_alloc(
                _mesh.node().count()) ;

            conn_sets _conn ;
            for (auto _isub = + 0 ;
                _isub != _nsub; ++_isub )
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL MOVE-NODE...\n" ) ;

                iptr_type  _nloc;
                move_node( _geom, _mesh , _conn ,
                    _hfun, _kern, _hval ,
                    _nset, _lset, _aset ,
                    _amrk, _mark,
                    _iter, _isub, _opts ,
                    _nloc, _QLIM, _DLIM , _tcpu);

                _nloc = _nloc / 2 ;

                _nmov = std::max (_nmov , _nloc);
            }

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

                iptr_type  _nloc;
                flip_mesh( _geom, _mesh , _hfun ,
                    _conn, _nset, _mark ,
                +3 * _iter - 2  , _nloc ) ;

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

            iptr_list _amrk, _lset, _aset ;
            _amrk.set_count(
                _mesh.node().count() ,
                    containers::tight_alloc,-1) ;

            _lset.set_alloc(
                _mesh.node().count()) ;
            _aset.set_alloc(
                _mesh.node().count()) ;

            _nsub = std::max(_nsub/2, +1) ;

            conn_sets _conn ;
            for (auto _isub = + 0 ;
                _isub != _nsub; ++_isub )
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                "**CALL MOVE-DUAL...\n" ) ;

                iptr_type  _nloc;
                move_dual( _geom, _mesh , _conn ,
                    _hfun, _hval,
                    _nset, _lset, _aset ,
                    _amrk, _mark,
                    _iter, _isub, _opts ,
                    _nloc, _QLIM, _DLIM , _tcpu);

                _nloc = _nloc / 2 ;

                _nmov = std::max (_nmov , _nloc);
            }

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

                iptr_type  _nloc;
                flip_mesh( _geom, _mesh , _hfun ,
                    _conn, _nset, _mark ,
                +3 * _iter - 1  , _nloc ) ;

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

                iptr_type  _nloc;
                flip_mesh( _geom, _mesh , _hfun ,
                    _nset, _mark,
                +3 * _iter - 0  , _nloc ) ;

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
            _dump.push(" *core-node: ");
            _dump.push(
            std::to_string(_tcpu._core_node)) ;
            _dump.push("\n");
            _dump.push(" *xDIR-node: ");
            _dump.push(
            std::to_string(_tcpu._ldir_node)) ;
            _dump.push("\n");
            _dump.push(" *xOPT-node: ");
            _dump.push(
            std::to_string(_tcpu._lopt_node)) ;
            _dump.push("\n\n");

            _dump.push("  MOVE-DUAL: ");
            _dump.push(
            std::to_string(_tcpu._move_dual)) ;
            _dump.push("\n");
            _dump.push(" *init-dual: ");
            _dump.push(
            std::to_string(_tcpu._init_dual)) ;
            _dump.push("\n");
            _dump.push(" *core-dual: ");
            _dump.push(
            std::to_string(_tcpu._core_dual)) ;
            _dump.push("\n");
            _dump.push(" *xDIR-dual: ");
            _dump.push(
            std::to_string(_tcpu._ldir_dual)) ;
            _dump.push("\n");
            _dump.push(" *xOPT-dual: ");
            _dump.push(
            std::to_string(_tcpu._lopt_dual)) ;
            _dump.push("\n\n");

            _dump.push("  TOPO-FLIP: ");
            _dump.push(
            std::to_string(_tcpu._topo_flip)) ;
            _dump.push("\n");

            _dump.push("  TOPO-ZIPS: ");
            _dump.push(
            std::to_string(_tcpu._topo_zips)) ;
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



