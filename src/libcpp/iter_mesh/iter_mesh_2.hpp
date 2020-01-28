
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
     * Last updated: 25 January, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * de2363@columbia.edu
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
     * ITER-MESH-2: hill-climbing surf. iter.
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
    typedef  H                          size_type ;
    typedef  P                          pred_type ;

    typedef typename
            mesh_type::real_type        real_type ;
    typedef typename
            mesh_type::iptr_type        iptr_type ;

    iptr_type static
        constexpr _dims = pred_type::_dims ;

    char_type static
        constexpr _odt_kern = +1 ;
    char_type static
        constexpr _cvt_kern = +2 ;
    char_type static
        constexpr dqdx_kern = +5 ;

    class   tria_kind {} ;           // dummy for overloads
    class   dual_kind {} ;

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

    public  :

    /*
    --------------------------------------------------------
     * FLIP-SIGN: flip tria for +ve iter. cost fn.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type flip_next (
        mesh_type &_mesh ,
        iptr_type  _tpos ,
        iptr_type  _epos ,
        iptr_type &_tadj ,
        iptr_type &_eadj ,
        iptr_list &_tset
        )
    {
        iptr_type _inod[3] ;
        iptr_type _jnod[3] ;
        mesh_type::tri3_type::
        face_node(_inod, _epos, 2, 1) ;
        _inod[ 0] = _mesh.
        _set3[_tpos].node(_inod[0]) ;
        _inod[ 1] = _mesh.
        _set3[_tpos].node(_inod[1]) ;
        _inod[ 2] = _mesh.
        _set3[_tpos].node(_inod[2]) ;

        _tset.set_count(0) ;

        _mesh.edge_tri3(_inod, _tset) ;

        if (_tset.count()!=+2) return ;

        if (_tset[0] == _tpos)
            _tadj = _tset[1] ;
        else
            _tadj = _tset[0] ;

        for(_eadj = 3 ; _eadj-- != 0; )
        {
        mesh_type::tri3_type::
        face_node(_jnod, _eadj, 2, 1) ;
        _jnod[ 0] = _mesh.
        _set3[_tadj].node(_jnod[0]) ;
        _jnod[ 1] = _mesh.
        _set3[_tadj].node(_jnod[1]) ;
        _jnod[ 2] = _mesh.
        _set3[_tadj].node(_jnod[2]) ;

        if (_jnod[ 2] != _inod[ 0])
        if (_jnod[ 2] != _inod[ 1])
            break  ;
        }

        if (_jnod[ 0] == _inod[ 0] &&
            _jnod[ 1] == _inod[ 1] )
        {
            std::swap (
        _mesh._set3[_tadj].node(0),
        _mesh._set3[_tadj].node(1)) ;
        }
    }

    /*------------------------- main sign-flip driver fn. */

    __static_call
    __normal_call void_type flip_sign (
        mesh_type &_mesh ,
        pred_type &_pred
        )
    {
        iptr_list _tset, _list, _seen ;

        __unreferenced(_pred); // for MSVC...

        _seen.set_count( _mesh.
            _set3.count(),
        containers::tight_alloc , +0) ;

    /*----------- an incremental BFS to correct mesh sign */
        iptr_type _tnum  = +0 ;
        iptr_type _epos  = +0 ;

        for (auto _tria  = _mesh._set3.head() ;
                  _tria != _mesh._set3.tend() ;
                ++_tria, ++_tnum )
        {
            if (_tria->mark() <  +0) continue ;
            if (_seen[_tnum ] >  +0) continue ;

        /*--------------------- flip seed for +ve quality */
            real_type _cost =
                    _pred.cost_tria (
               &_mesh._set1[
                _tria->node(0)].pval(0),
               &_mesh._set1[
                _tria->node(1)].pval(0),
               &_mesh._set1[
                _tria->node(2)].pval(0)) ;

            if (_cost < (real_type) +0.)
            {
                std::swap (
                    _tria->node(0) ,
                        _tria->node(1));
            }

        /*--------------------- a BFS from seed via topo. */
            _list.push_tail(_tnum) ;
            _seen [_tnum] =  +1;

            for ( ; !_list.empty() ; )
            {
                iptr_type _tpos;
                _list._pop_tail(_tpos) ;

                for (_epos = +3; _epos-- != +0; )
                {
                    iptr_type _tadj = -1 ;
                    iptr_type _eadj = -1 ;

                    flip_next( _mesh, _tpos ,
                        _epos, _tadj, _eadj ,
                        _tset) ;

                    if (_tadj == -1) continue ;

                    if (_seen[_tadj] == +0 )
                    {
                        _seen[_tadj]  = +1 ;
                        _list.push_tail(_tadj);
                    }
                }
            }
        }

    }

    /*
    --------------------------------------------------------
     * MOVE-OKAY: TRUE if state sufficiently good.
    --------------------------------------------------------
     */

    __static_call
    __inline_call void_type move_okay (
        real_list &_cdst ,
        real_list &_csrc ,
        bool_type &_okay ,
        real_type  _good = +9.25E-01,
        real_type  _qtol = +1.00E-04,
        real_type  _xdel = +0.00E+00,
        real_type  _xtol =
    std::numeric_limits<real_type>::infinity()
        )
    {
        iptr_type _move;
        move_okay(_cdst, _csrc, _move ,
            _good,_qtol,
            _xdel,_xtol) ;

        _okay = ( _move>(iptr_type)0 );
    }

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
        real_type _0src =
            +std::numeric_limits
                <real_type>::infinity();

        real_type _0dst =
            +std::numeric_limits
                <real_type>::infinity();

        real_type _msrc, _mdst;
        _msrc = (real_type) +0. ;
        _mdst = (real_type) +0. ;

        for (auto _iter  = _csrc.head(),
                  _tend  = _csrc.tend();
                  _iter != _tend;
                ++_iter  )
        {
            _0src  =
            std::min(_0src, *_iter);

          //_msrc += *_iter ;

            _msrc +=
            (real_type)1. / *_iter ;
        }
        for (auto _iter  = _cdst.head(),
                  _tend  = _cdst.tend();
                  _iter != _tend;
                ++_iter  )
        {
            _0dst  =
            std::min(_0dst, *_iter);

          //_mdst += *_iter ;

            _mdst +=
            (real_type)1. / *_iter ;
        }

      //_msrc /= _csrc.count() ;
      //_mdst /= _cdst.count() ;

        _msrc  =
            _csrc.count()/_msrc;
        _mdst  =
            _cdst.count()/_mdst;

    /*--------------------- calc. min. + mean rel-tol */
        _qtol *= std::max(
            _0src, (real_type) +0.);

        real_type _mtol = _qtol;
        real_type _0tol = _qtol;

        _mtol /= _csrc.count() ;
        _mtol /= _cdst.count() ;

    /*--------------------- calc. min. + mean delta's */
        real_type _0del , _mdel;
        _0del  =  _0dst - _0src;
        _mdel  =  _mdst - _msrc;

        _0del /= std::min(
            _cdst.count(),
                _csrc.count()) ;

        _mdel *= std::min(
            _cdst.count(),
                _csrc.count()) ;

        _0del  =
        std::max(_0del,-_0tol) ;
        _mdel  =
        std::max(_mdel,-_mtol) ;

    /*---------------------------- test move = 'okay' */
        if (true)
        {
    /*--------------------- okay if min. is improving */
            if (_0dst > _0src+_0tol)
            if (_mdst > _msrc-_0del)
                _move = +1;

            if (_move > +0) return ;
        }

        if (_0dst >=
            std::pow(_good, +2))
        {
    /*--------------------- okay if mean is improving */
            if (_mdst > _msrc+_mtol)
            if (_0dst > _0src-_mdel)
                _move = +1;

            if (_move > +0) return ;
        }

        if (_0dst >= _good)
        {
    /*--------------------- okay if moves unconverged */
            if (_xdel > _xtol)
            if (_0dst > _0src-_xdel)
                _move = +1;

            real_type _XDEL =
                _xdel/_cdst.count();

            if (_xdel > _xtol)
            if (_mdst > _msrc-_XDEL)
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
        pred_type &_pred ,
        iptr_list &_tset ,
        real_list &_cost ,
        tria_kind const&
        )
    {
        real_type _qmin =
            +std::numeric_limits
                <real_type>::infinity();

        __unreferenced(_pred) ; // for MSVC...

        for (auto _tria  = _tset.head(),
                  _tend  = _tset.tend();
                  _tria != _tend;
                ++_tria  )
        {
            real_type _tscr =
                    _pred.cost_tria (
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(0)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(1)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(2)].pval(0)) ;

            _qmin =
            std::min (_qmin, _tscr) ;

            _cost.push_tail (_tscr) ;
        }

        return ( _qmin )  ;
    }

    __static_call
    __normal_call real_type loop_cost (
        mesh_type &_mesh ,
        pred_type &_pred ,
        iptr_list &_tset ,
        real_list &_cost ,
        dual_kind const&
        )
    {
        real_type _qmin =
            +std::numeric_limits
                <real_type>::infinity();

        __unreferenced(_pred) ; // for MSVC...

        for (auto _tria  = _tset.head(),
                  _tend  = _tset.tend();
                  _tria != _tend;
                ++_tria  )
        {
            real_type _tscr =
                    _pred.cost_dual (
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(0)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(1)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(2)].pval(0)) ;

            _qmin =
            std::min (_qmin, _tscr) ;

            _cost.push_tail (_tscr) ;
        }

        return ( _qmin )  ;
    }

    /*
    --------------------------------------------------------
     * MOVE-NODE: "smart" coord. update for single node.
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
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iter_opts &_opts ,
        node_iter  _node ,
        char_type  _kern ,
        iptr_type &_move ,
        iptr_list &_tset ,
        real_list &_told ,
        real_list &_tnew ,
        real_type  _TMIN ,
        real_type  _TLIM
        )
    {
        move_kern( _geom, _mesh, _hfun,
            _pred, _hval, _opts, _node,
            _kern, _move, _tset, _told,
            _tnew, _TMIN, _TLIM) ;

        if (_move >= +0 ) return ;

        move_kern( _geom, _mesh, _hfun,
            _pred, _hval, _opts, _node,
        dqdx_kern, _move, _tset, _told,
            _tnew, _TMIN, _TLIM) ;

        if (_move >= +0 ) return ;
    }

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type move_kern (
        geom_type &_geom ,
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iter_opts &_opts ,
        node_iter  _node ,
        char_type  _kern ,
        iptr_type &_move ,
        iptr_list &_tset ,
        real_list &_told ,
        real_list &_tnew ,
        real_type  _TMIN ,
        real_type  _TLIM
        )
    {
        iptr_type static
            constexpr _ITER = (iptr_type)+5 ;

        _move = (iptr_type)-1 ;

    /*---------------- calc. line search direction vector */
        real_type _line [_dims] = {
       (real_type) +0.0 } ;
        real_type _save [_dims] = {
       (real_type) +0.0 } ;
        real_type _proj [_dims] = {
       (real_type) +0.0 } ;

        real_type _ladj = (real_type) +0.00 ;

        if (_kern == _odt_kern)
        {
    /*--------------------------- ODT-style update vector */
            _odt_move_2 (
                _mesh, _hfun, _pred ,
                _hval, _tset, _node ,
                _line, _ladj) ;
        }
        else
        if (_kern == _cvt_kern)
        {
    /*--------------------------- CVT-style update vector */
            _cvt_move_2 (
                _mesh, _hfun, _pred ,
                _hval, _tset, _node ,
                _line, _ladj) ;
        }
        else
        if (_kern == dqdx_kern)
        {
            if (_TMIN<=_TLIM)
            {
    /*--------------------------- d./dx Q^T update vector */
            dqdx_move_2 (
                _mesh, _hfun, _pred ,
                _tset, _node, _told ,
                _line, _ladj) ;
            }
            else { return ; }
        }

    /*---------------- scale line search direction vector */
        real_type _xtol =           // delta_x reltol
            (real_type) +5.00E-004 ;

        real_type _xeps =           // delta_x ~= 0.0
       (real_type)+.01 *_opts.qtol() ;

        if (_kern == dqdx_kern)
        {
            _xtol = (real_type)+1. ;
        }

        real_type _lsqr ;
        _lsqr = std::pow(_ladj, 2) ;
        _xtol = std::pow(_xtol, 2) ;
        _xeps = std::pow(_xeps, 2) ;

        real_type _scal =           // overrelaxation
            (real_type) +5.0 / 4.0 ;

    /*---------------- do backtracking line search iter's */

        for (auto _idim = _dims; _idim-- != +0; )
        {
            _save[_idim]  =
                _node->pval(_idim) ;

          //_line[_idim] /= _llen  ;
        }

        for (auto _iter = +0 ;
                _iter != _ITER ; ++_iter)
        {
            _tnew.set_count(0) ;

    /*---------------- push update along search direction */
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _proj[_idim] =
                    _save[_idim] +
                        _scal * _line[_idim];
            }

            _pred.proj_node (
                _geom, _save, _proj) ;

            for (auto _idim = _dims; _idim-- != +0; )
            {
                _node->pval(_idim)
                    = _proj[_idim] ;
            }

            real_type _XTOL = _xtol * _scal ;
            real_type _XEPS = _xeps * _scal ;

            real_type _lmov = _pred.
            length_sq(_save , _proj)/ _lsqr ;

            if (_lmov <= _XEPS) break ;

          //_move  = +1 ; return ;

            _scal *= (real_type).5 ;

    /*---------------- test quasi-monotonicity w.r.t. Q^T */
            loop_cost( _mesh, _pred,
                       _tset,
                       _tnew,
                 tria_kind()) ;

            move_okay( _tnew,
                _told, _move,
                _TLIM, _opts.qtol(),
                _lmov, _XTOL) ;

            if (_move > 0)   break ;
        }

        if (_move <= (iptr_type)0)
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _node->pval(_idim)
                    = _save[_idim] ;
            }
        }

    }

    /*
    --------------------------------------------------------
     * MOVE-DUAL: "smart" weight update for single node.
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
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iter_opts &_opts ,
        node_iter  _node ,
        iptr_type &_move ,
        iptr_list &_tset ,
        real_list &_dold ,
        real_list &_dnew ,
        real_type  _DMIN ,
        real_type  _DLIM
        )
    {
        iptr_type static
            constexpr _ITER = (iptr_type)+5 ;

        __unreferenced(_hval);
        __unreferenced(_DMIN);

        _move = (iptr_type)-1;

        real_type _wadj, _step, _save;

    /*---------------- calc. line search direction vector */
        if (true) // (_DMIN < _DLIM)
        {
            dqdw_move_2 (
                _geom, _mesh, _hfun,
                _pred, _tset, _node,
                _dold, _step, _wadj) ;
        }
        else { return ; }

    /*---------------- scale line search direction vector */
        real_type _xeps =           // delta_w ~= 0.0
       (real_type)+.01 *_opts.qtol() ;

        real_type _scal =           // overrelaxation
            (real_type) +5.0 / 4.0 ;

        _save = _node->pval(_dims) ;

    /*---------------- do backtracking line search iter's */

        for (auto _iter = +0 ;
                _iter != _ITER; ++_iter )
        {
            _node->pval(_dims) =
                _save   + _scal * _step ;

            _node->pval(_dims) =    // not too large!
                std::max(-_wadj,
                    _node->pval(_dims));
            _node->pval(_dims) =
                std::min(+_wadj,
                    _node->pval(_dims));

            real_type _wmov =
                std::abs (_save
                  - _node->pval(_dims));

            if (_wmov <=
            _xeps * _scal * _wadj) break;

            _dnew.set_count(0);

            _scal *= (real_type).5 ;

    /*---------------- test quasi-monotonicity w.r.t. Q^D */
            loop_cost( _mesh, _pred,
                       _tset,
                       _dnew,
                 dual_kind()) ;

            move_okay( _dnew,
                _dold, _move,
                _DLIM, _opts.qtol()) ;

            if (_move > 0)   break ;
        }

        if (_move <= (iptr_type) +0)
        {
            _node->pval(_dims) = _save ;
        }

    }

    /*
    --------------------------------------------------------
     * SORT-NODE: permutations for node optimisation.
    --------------------------------------------------------
     */

    __static_call void_type sort_node (
        geom_type &_geom ,
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        real_list &_qscr ,
        iptr_list &_nset ,
        iptr_list &_aset ,
        iptr_list &_amrk ,
        iptr_list &_nmrk ,
        iptr_type  _iout ,
        iptr_type  _isub ,
        iter_opts &_opts ,
        real_type  _TLIM ,
        real_type  _DLIM
        )
    {
        class cost_pair
            {
            public  :
    /*------------------------ tuple for node re-ordering */
            iptr_type           _node ;
            float               _cost ;
            } ;

        class cost_less
            {
            public  :
    /*------------------------ less-than op. for cost-tup */
            __inline_call
                bool_type operator () (
                cost_pair const&_idat ,
                cost_pair const&_jdat
                ) const
            {    return
            _idat._cost < _jdat._cost ;
            }
            } ;

        typedef containers::
            array<cost_pair>cost_list ;

        iptr_list _eset ;
        cost_list _sset ;

        __unreferenced(_geom) ;
        __unreferenced(_hfun) ;
        __unreferenced(_pred) ;
        __unreferenced(_hval) ;
        __unreferenced(_opts) ;
        __unreferenced(_DLIM) ;

        if (_isub == (iptr_type) +0 )
        {
    /*-------------------- 1ST SUB-ITER: build full init. */
        for (auto _iter  = _qscr.head();
                  _iter != _qscr.tend();
                ++_iter  )
        {
            *_iter = (real_type) +1.;
        }

        for (auto _tria  = _mesh._set3.head();
                  _tria != _mesh._set3.tend();
                ++_tria  )
        {
            if (_tria->mark() >= +0 )
            {
        /*-------------------- calc. min. scores at nodes */
                iptr_type _inod, _jnod, _knod;
                _inod = _tria->node(0);
                _jnod = _tria->node(1);
                _knod = _tria->node(2);

                real_type _cost =
                    _pred.cost_tria (
               &_mesh._set1[_inod].pval(0),
               &_mesh._set1[_jnod].pval(0),
               &_mesh._set1[_knod].pval(0)
                    ) ;

                iptr_type _flag = _iout-2 ;

                real_type _TURN =
                    std::pow(_TLIM, 1./4.);

                if (_cost <= _TURN ||
                        std::abs(
                    _nmrk[_inod])>= _flag ||
                        std::abs(
                    _nmrk[_jnod])>= _flag ||
                        std::abs(
                    _nmrk[_knod])>= _flag )
                {

                if (_amrk[_inod] != _isub )
                {
                    _amrk[_inod]  = _isub ;
                    _sset.push_tail() ;
                    _sset.tail()->
                        _node = _inod ;
                }

                if (_amrk[_jnod] != _isub )
                {
                    _amrk[_jnod]  = _isub ;
                    _sset.push_tail() ;
                    _sset.tail()->
                        _node = _jnod ;
                }

                if (_amrk[_knod] != _isub )
                {
                    _amrk[_knod]  = _isub ;
                    _sset.push_tail() ;
                    _sset.tail()->
                        _node = _knod ;
                }

                }

                _qscr[_inod] = std::min (
                    _cost,_qscr[_inod]) ;

                _qscr[_jnod] = std::min (
                    _cost,_qscr[_jnod]) ;

                _qscr[_knod] = std::min (
                    _cost,_qscr[_knod]) ;
            }
        }

        for (auto _iter  = _sset.head() ;
                  _iter != _sset.tend() ;
                ++_iter  )
        {
    /*------------------------ assign min.-cost for nodes */
            _iter->_cost =
            (float)_qscr[_iter->_node];
        }

        algorithms::qsort( _sset.head() ,
                           _sset.tend() ,
            cost_less () ) ;

        for (auto _iter  = _sset.head() ;
                  _iter != _sset.tend() ;
                ++_iter  )
        {
    /*------------------------ push sorted wrt. min.-cost */
            _aset.push_tail(_iter->_node) ;
        }

        }
        else
        {
    /*-------------------- NTH sub-iter: init. from prev. */
        for (auto _iter  = _nset.head() ;
                  _iter != _nset.tend() ;
                ++_iter  )
        {
            _amrk[*_iter]  = _isub;
            _aset.push_tail(*_iter) ;
        }

        for (auto _iter  = _nset.head() ;
                  _iter != _nset.tend() ;
                ++_iter  )
        {
            _eset.set_count(0) ;

            _mesh.node_edge (
               &*_iter, _eset) ;

            for (auto _edge  = _eset.head();
                      _edge != _eset.tend();
                    ++_edge  )
            {
                 auto _eptr =
                _mesh._set2.head() + *_edge;

                iptr_type _inod, _jnod;
                _inod = _eptr->node(0);
                _jnod = _eptr->node(1);

                if (_amrk[_inod] != _isub)
                {
                    _amrk[_inod]  = _isub;
                    _aset.push_tail(_inod) ;
                }

                if (_amrk[_jnod] != _isub)
                {
                    _amrk[_jnod]  = _isub;
                    _aset.push_tail(_jnod) ;
                }
            }
        }

        }

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
        size_type &_hfun ,
        pred_type &_pred ,
        char_type  _kern ,
        real_list &_hval ,
        real_list &_qscr ,
        iptr_list &_nset ,
        iptr_list &_amrk ,
        iptr_list &_nmrk ,
        iptr_list &_emrk ,
        iptr_list &_tmrk ,
        iptr_type  _iout ,
        iptr_type  _isub ,
        iter_opts &_opts ,
        iptr_type &_nmov ,
        real_type  _TLIM ,
        real_type  _DLIM
        )
    {
        iptr_list _aset, _tset;
        real_list _told, _tnew, _dold, _dnew;

        __unreferenced ( _emrk) ;
        __unreferenced ( _tmrk) ;

        _nmov = (iptr_type)  +0 ;

    /*-------------------- permute nodes for optimisation */
        sort_node( _geom, _mesh, _hfun,
            _pred, _hval, _qscr,
            _nset, _aset,
            _amrk, _nmrk, _iout, _isub,
            _opts, _TLIM, _DLIM ) ;

    /*-------------------- GAUSS-SEIDEL iteration on TRIA */
        if (_opts .tria())
        {
        for (auto _apos  = _aset.head() ;
                  _apos != _aset.tend() ;
                ++_apos  )
        {
             auto _node  =
            _mesh._set1.head() + *_apos ;

            _tset.set_count( +0);

        /*---------------- assemble a local tria. stencil */
            _mesh.node_tri3(
                &_node->node(+0), _tset);

            if (_tset.empty()) continue ;

            if (_nmrk[*_apos] >= +0)
            {
        /*---------------- attempt to optimise TRIA geom. */
            _told.set_count( +0);
            _tnew.set_count( +0);

            real_type _TMIN =
                loop_cost( _mesh,
                    _pred, _tset,
                    _told,  tria_kind());

            iptr_type _move = -1;

            if(_move < +0)
            {
        /*---------------- do optimisation of node coord. */
                move_node( _geom, _mesh ,
                    _hfun, _pred, _hval ,
                    _opts, _node, _kern ,
                    _move, _tset,
                    _told, _tnew,
                    _TMIN, _TLIM ) ;
            }

            if (_move > +0)
            {
        /*---------------- update when state is improving */
            _hval[*_apos] = (real_type)-1. ;

            if (std::abs(
               _nmrk[*_apos]) != _iout)
            {
                if (_nmrk[*_apos] >= 0)
                _nmrk[*_apos] = +_iout;
                else
                _nmrk[*_apos] = -_iout;

                _nset.push_tail(*_apos) ;
            }

            _nmov += +1 ;
            }
            }
        }
        }

    /*-------------------- GAUSS-SEIDEL iteration on DUAL */
        if (_opts .dual())
        {
        for (auto _apos  = _aset.head() ;
                  _apos != _aset.tend() ;
                ++_apos  )
        {
             auto _node  =
            _mesh._set1.head() + *_apos ;

            _tset.set_count( +0);

        /*---------------- assemble a local tria. stencil */
            _mesh.node_tri3(
                &_node->node(+0), _tset);

            if (_tset.empty()) continue ;

        /*---------------- attempt to optimise DUAL geom. */
            _dold.set_count( +0);
            _dnew.set_count( +0);

            real_type _DMIN =
                loop_cost( _mesh,
                    _pred, _tset,
                    _dold,  dual_kind());

            iptr_type _move = -1;

            if(_move < +0)
            {
        /*---------------- do optimisation of node weight */
                move_dual( _geom, _mesh ,
                    _hfun, _pred, _hval ,
                    _opts, _node,
                    _move, _tset,
                    _dold, _dnew,
                    _DMIN, _DLIM ) ;
            }

            if (_move > +0)
            {
        /*---------------- update when state is improving */
            if (std::abs(
               _nmrk[*_apos]) != _iout)
            {
                if (_nmrk[*_apos] >= 0)
                _nmrk[*_apos] = +_iout;
                else
                _nmrk[*_apos] = -_iout;

                _nset.push_tail(*_apos) ;
            }

            _nmov += +1 ;
            }
        }
        }

    }

    /*
    --------------------------------------------------------
     * FLIP-MESH: "flip" mesh topology.
    --------------------------------------------------------
     */

    #include "iter_flip_2.inc"

    __static_call
    __inline_call void_type flip_tria (
        geom_type &_geom ,
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        iptr_type  _tria ,
        bool_type &_flip ,
        iptr_list &_told ,
        iptr_list &_tnew ,
        real_list &_qold ,
        real_list &_qnew
        )
    {
        _flip  = false ;

         auto
        _coin  = std::rand() % +3 ;

        if (_coin == +0)
        {
    /*--------------------------------- flip edges: 0,1,2 */
            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +0 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;

            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +1 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;

            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +2 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;
        }
        else
        if (_coin == +1)
        {
    /*--------------------------------- flip edges: 1,2,0 */
            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +1 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;

            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +2 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;

            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +0 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;
        }
        else
        {
    /*--------------------------------- flip edges: 2,0,1 */
            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +2 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;

            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +0 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;

            flip_t2t2( _geom, _mesh ,
                _hfun, _pred,
                _tria,   +1 ,
                _told, _tnew, _flip ,
                _qold, _qnew) ;
            if (_flip) return ;
        }

    }

    __static_call
    __normal_call void_type flip_mesh (
        geom_type &_geom ,
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        iptr_list &_nset ,
        iptr_list &_nmrk ,
        iptr_list &_emrk ,
        iptr_list &_tmrk ,
        iptr_type  _imrk ,
        iptr_type &_nflp
        )
    {
        init_mark(_mesh, _nmrk, _emrk, _tmrk,
            std::max(+0, _imrk - 1)) ;

    /*--------------------- init. flip stack as ADJ(NSET) */
        iptr_list _tset, _next;
        iptr_list _told, _tnew;
        real_list _qold, _qnew;

        for (auto _iter  = _nset.head();
                  _iter != _nset.tend();
                ++_iter  )
        {
            if (  _mesh.
                _set1[*_iter].mark() >= +0)
            {
                _tnew.set_count(+0);

                _mesh.node_tri3(
                    &*_iter, _tnew);

                for (auto _tadj  = _tnew.head();
                          _tadj != _tnew.tend();
                        ++_tadj  )
                {
                if (_tmrk[*_tadj] != _imrk)
                {
                    _tmrk[*_tadj]  = _imrk;
                    _tset.push_tail(*_tadj) ;
                }
                }
            }
        }

    /*--------------------- exhaustive, incremental flips */
        _nflp = +0 ;

        for ( ; !_tset.empty() ; )
        {
        for (auto _tria  = _tset.head();
                  _tria != _tset.tend();
                ++_tria  )
        {
            if (  _mesh.
                _set3[*_tria].mark() >= +0)
            {
                bool_type  _flip = false ;
                flip_tria( _geom, _mesh,
                    _hfun, _pred,
                   *_tria, _flip,
                    _told, _tnew,
                    _qold, _qnew );
                if (_flip) _nflp += +1 ;

                for (auto _iter  = _tnew.head();
                          _iter != _tnew.tend();
                        ++_iter  )
                {
                    _next.push_tail(*_iter) ;
                }
            }
        }
            _tset = std::move(_next) ;
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
        size_type &_hfun ,
        pred_type &_pred ,
        char_type  _kern ,
        real_list &_hval ,
        iptr_list &_nset ,
        iptr_list &_nmrk ,
        iptr_list &_emrk ,
        iptr_list &_tmrk ,
        iptr_type  _imrk ,
        iter_opts &_opts ,
        real_type  _TLIM ,
        real_type  _DLIM ,
        iptr_type &_nzip ,
        iptr_type &_ndiv
        )
    {
        class sort_pair
            {
            public  :
    /*------------------------ tuple for edge re-ordering */
            iptr_type           _edge ;
            float               _cost ;
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
            {    return
            _idat._cost < _jdat._cost ;
            }
            } ;

    #   define __marknode                   \
            init_mark( _mesh, _nmrk, _emrk, \
                _tmrk, std::max(_imrk-1, +0)) ; \
            if (std::abs(               \
                _nmrk[_nnew])!= _imrk+1)    \
            {                           \
            if (_nmrk[_nnew] >= +0)     \
            {                           \
                _nmrk[_nnew] = +_imrk+1;    \
            }                           \
            else                        \
            {                           \
                _nmrk[_nnew] = -_imrk-1;    \
            }                           \
                _nset.push_tail(_nnew) ;    \
            }                           \

        iptr_type static
            constexpr _DEG_MIN = (iptr_type) +5 ;
        iptr_type static
            constexpr _DEG_MAX = (iptr_type) +7 ;

        typedef containers::
            array<sort_pair>sort_list ;

        __unreferenced  (_DLIM) ;

        _nzip = +0 ; _ndiv = +0 ;

        sort_list _sort;

        iptr_list _aset, _bset, _cset ;
        iptr_list _eset, _done;
        iptr_list _iset, _jset;
        real_list _told, _tnew, _ttmp ;

        for (auto _node =
            _mesh._set1.count() ; _node-- != +0; )
        {
    /*--------------------- scan nodes and zip//div edges */
            if (  _mesh.
                _set1[_node].mark () >= +0 &&
                   std::abs (
                _nmrk[_node]) >= _imrk - 2 )
            {
                _eset.set_count(+0) ;
                _sort.set_count(+0) ;
                _mesh.node_edge(
                (iptr_type) _node, _eset) ;

                for (auto _eadj  = _eset.head();
                          _eadj != _eset.tend();
                        ++_eadj  )
                {
                     auto _eptr =
                    _mesh._set2.head() + *_eadj;

                     auto _iptr = _mesh.
                    _set1.head()+ _eptr->node(0) ;
                     auto _jptr = _mesh.
                    _set1.head()+ _eptr->node(1) ;

                    real_type _line[_dims] ;
                    iptr_type _idim=_dims;

                    for ( ; _idim-- != +0; )
                    {
                        _line[_idim] =
                            _jptr->pval(_idim) -
                            _iptr->pval(_idim) ;
                    }

                    float _lsqr  = (float)
                    _pred.length_sq(_line) ;

                    sort_pair _pair;
                    _pair._cost = _lsqr ;
                    _pair._edge =*_eadj ;

                    _sort.push_tail(_pair) ;
                }

                if (_sort.empty()) continue;

            /*------------------- scan local edges by len */
                algorithms::isort (
                    _sort.head(),
                    _sort.tend(), sort_less()) ;

                bool_type _move = false ;

                for (auto _iter  = _sort.tail();
                          _iter != _sort.hend();
                        --_iter  )
                {
                     auto _eadj  = _iter->_edge;

                     auto _eptr  =
                    _mesh._set2.head()  + _eadj;

                    iptr_type  _enod[2] ;
                    _enod[0] = _eptr->node(0);
                    _enod[1] = _eptr->node(1);

                    if (std::abs(
                    _nmrk[_enod[0]]) > _imrk ||
                        std::abs(
                    _nmrk[_enod[1]]) > _imrk )
                    continue ;

            /*------------------- try to "div" local edge */
                    if (_nmrk[_enod[0]] >= 0 &&
                        _nmrk[_enod[1]] >= 0 )
                    {
                    if (_eset.count() > _DEG_MAX)
                    {
                        real_type _qinc =
                            (real_type) -1./9. ;
                        real_type _ltol =
                            (real_type) +8./9. ;

                        _qinc *=
                            _eset.count() - 6;

                        _qinc /= // don't ping-pong w zip
                            std::sqrt (_imrk);

                        iptr_type  _nnew = -1;

                        if (_opts.div_())
                        _div_edge( _geom, _mesh,
                            _hfun, _pred,
                            _hval, _opts, _eadj,
                            _kern, _move, _nnew,
                            _iset, _jset,
                            _told, _tnew,
                            _ttmp, _TLIM,
                            _ltol, _qinc) ;

                        if (_move)
                        {
                        __marknode; _ndiv += +1; break ;
                        }
                    }
                    else
                    {
                        iptr_type  _nnew = -1;

                        if (_opts.div_())
                        _div_edge( _geom, _mesh,
                            _hfun, _pred,
                            _hval, _opts, _eadj,
                            _kern, _move, _nnew,
                            _iset, _jset,
                            _told, _tnew,
                            _ttmp, _TLIM) ;

                        if (_move)
                        {
                        __marknode; _ndiv += +1; break ;
                        }
                    }
                    }
                }

                if (_move) continue ;

            /*------------------- scan local edges by len */
                for (auto _iter  = _sort.head();
                          _iter != _sort.tend();
                        ++_iter  )
                {
                     auto _eadj  = _iter->_edge;

                     auto _eptr  =
                    _mesh._set2.head()  + _eadj;

                    iptr_type  _enod[2] ;
                    _enod[0] = _eptr->node(0);
                    _enod[1] = _eptr->node(1);

                    if (std::abs(
                    _nmrk[_enod[0]]) > _imrk ||
                        std::abs(
                    _nmrk[_enod[1]]) > _imrk )
                    continue ;

            /*------------------- try to "zip" local edge */
                    if (_nmrk[_enod[0]] >= 0 &&
                        _nmrk[_enod[1]] >= 0 )
                    {
                    if (_eset.count() < _DEG_MIN)
                    {
                        real_type _qinc =
                            (real_type) -1./9. ;
                        real_type _ltol =
                            (real_type) +9./8. ;

                        _qinc /= // don't ping-pong w div
                            std::sqrt (_imrk);

                        iptr_type  _nnew = -1;

                        if (_opts.zip_())
                        _zip_edge( _geom, _mesh,
                            _hfun, _pred,
                            _hval, _opts, _eadj,
                            _kern, _move, _nnew,
                            _iset, _jset,
                            _aset, _bset, _cset,
                            _told, _tnew,
                            _ttmp, _TLIM,
                            _ltol, _qinc) ;

                        if (_move)
                        {
                        __marknode; _nzip += +1; break ;
                        }
                    }
                    else
                    {
                        iptr_type  _nnew = -1;

                        if (_opts.zip_())
                        _zip_edge( _geom, _mesh,
                            _hfun, _pred,
                            _hval, _opts, _eadj,
                            _kern, _move, _nnew,
                            _iset, _jset,
                            _aset, _bset, _cset,
                            _told, _tnew,
                            _ttmp, _TLIM) ;

                        if (_move)
                        {
                        __marknode; _nzip += +1; break ;
                        }
                    }
                    }
                }

                if (_move) continue ;

            }
        }

        for (auto _iter  = _nmrk.head() ;
                  _iter != _nmrk.tend() ;
                ++_iter  )
        {
    /*--------------------- undo local inc. on node flags */
            if (*_iter > + _imrk)
                *_iter = + _imrk;
            else
            if (*_iter < - _imrk)
                *_iter = - _imrk;
        }

    #   undef  __marknode

    }

    /*------------------------------ helper: init. marker */

    __static_call
    __normal_call void_type init_mark (
        mesh_type &_mesh ,
        iptr_list &_nmrk ,
        iptr_list &_emrk ,
        iptr_list &_tmrk ,
        iptr_type  _flag = +0
        )
    {
        iptr_type _nmax =
       (iptr_type)std::max( _nmrk.count() ,
                      _mesh._set1.count()
                    ) ;
        iptr_type _emax =
       (iptr_type)std::max( _emrk.count() ,
                      _mesh._set2.count()
                    ) ;
        iptr_type _tmax =
       (iptr_type)std::max( _tmrk.count() ,
                      _mesh._set3.count()
                    ) ;

        _nmrk.set_count(_nmax,
            containers::
                loose_alloc, _flag) ;
        _emrk.set_count(_emax,
            containers::
                loose_alloc, _flag) ;
        _tmrk.set_count(_tmax,
            containers::
                loose_alloc, _flag) ;
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
        size_type &_hfun ,
        mesh_type &_mesh ,
        char_type  _kern ,
        pred_type &_pred ,
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
        high_resolution_clock::
            time_point _ttic ;
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*------------------------------ ensure deterministic */
        std::srand( +1 ) ;

    /*------------------------------ push boundary marker */
        iptr_list _nmrk, _emrk, _tmrk,
                  _nset, _tset;

        init_mark(_mesh, _nmrk, _emrk, _tmrk) ;

        iptr_type _nnum = +0  ;
        iptr_type _enum = +0  ;

        for (auto _node  = _mesh._set1.head() ;
                  _node != _mesh._set1.tend() ;
                ++_node, ++_nnum )
        {
            if (_node->mark() >= +0)
            {
            if (_node->feat()
                    != mesh::null_feat)
            {
                _nmrk[_nnum]   = -1;
            }
            }
        }

        for (auto _edge  = _mesh._set2.head() ;
                  _edge != _mesh._set2.tend() ;
                ++_edge, ++_enum )
        {
            if (_edge->mark() >= +0)
            {
            if (_edge->self() >= +1)
            {
                _nmrk[_edge->node(0)] = -1 ;
                _nmrk[_edge->node(1)] = -1 ;
            }
            else
            {
                _tset.set_count(0) ;

                _mesh.edge_tri3 (
                   &_edge->node(0), _tset) ;

                if (_tset.count() != +2)
                {
                _nmrk[_edge->node(0)] = -1 ;
                _nmrk[_edge->node(1)] = -1 ;
                }
            }
            }
        }

        flip_sign (_mesh, _pred) ;

    /*------------------------------ do optimisation loop */
        iptr_type static constexpr
            ITER_MIN_ = +  4 ;
        iptr_type static constexpr
            ITER_MAX_ = +  8 ;

        bool_type static constexpr
            ITER_FLIP = true ;

        real_type _QMIN = (real_type) +1. ;

        for (auto _tria  = _mesh._set3.head() ;
                  _tria != _mesh._set3.tend() ;
                ++_tria )
        {
            if (_tria->mark() >= +0 )
            {
        /*--------------------- test initial cell quality */
            real_type _cost =
                    _pred.cost_tria (
               &_mesh._set1[
                _tria->node(0)].pval(0),
               &_mesh._set1[
                _tria->node(1)].pval(0),
               &_mesh._set1[
                _tria->node(2)].pval(0)
                    ) ;

            _QMIN = std::min (_QMIN, _cost) ;
            }
        }

        for (auto _iter = +1 ;
            _iter <= _opts.iter(); ++_iter)
        {
        /*-------------------------- set-up current iter. */
            init_mark(_mesh, _nmrk,
            _emrk, _tmrk, std::max(_iter-1, +0));

            real_list _hval, _qmin;
            _hval.set_count(
                _mesh._set1.count(),
        containers::tight_alloc, (real_type)-1.);

            _qmin.set_count(
                _mesh._set1.count(),
        containers::tight_alloc, (real_type)+1.);

            iptr_list _amrk;
            _amrk.set_count(
                _mesh._set1.count(),
        containers::tight_alloc, (iptr_type)-1 );

            _nset.set_count(  +0);

            iptr_type _nmov = +0 ;
            iptr_type _nflp = +0 ;
            iptr_type _nzip = +0 ;
            iptr_type _ndiv = +0 ;

    /*------------------------------ scale quality thresh */
            iptr_type _nsub = _iter + 0 ;

            _nsub = std::min(
                ITER_MAX_, _nsub) ;
            _nsub = std::max(
                ITER_MIN_, _nsub) ;

            real_type _TLIM =
           (real_type).750*_opts.qlim() +
           (real_type).075*_iter;

            _TLIM = std::max( _TLIM,
                              _QMIN) ;
            _TLIM = std::min(
                _opts.qlim(), _TLIM) ;

            real_type _DLIM =
                (real_type)(1. -
            1. * std::pow(1.-_TLIM, 2)) ;

    /*------------------------------ update mesh geometry */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            for (auto _isub = + 0 ;
                _isub != _nsub; ++_isub )
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                " CALL MOVE-NODE...\n") ;

                iptr_type  _nloc;
                move_node( _geom, _mesh ,
                    _hfun, _pred, _kern ,
                    _hval, _qmin,
                    _nset, _amrk,
                    _nmrk, _emrk, _tmrk ,
                    _iter, _isub,
                    _opts, _nloc,
                    _TLIM, _DLIM) ;

                _nmov = std::max (_nmov ,
                                  _nloc ) ;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;

            _tcpu._move_full +=
                _tcpu.time_span(_ttic, _ttoc);
    #       endif//__use_timers

    /*------------------------------ update mesh topology */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            if (ITER_FLIP)
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                " CALL FLIP-MESH...\n") ;

                flip_mesh( _geom, _mesh ,
                    _hfun, _pred, _nset ,
                    _nmrk, _emrk, _tmrk ,
                    _iter, _nflp) ;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;

            _tcpu._topo_full +=
                _tcpu.time_span(_ttic, _ttoc);
    #       endif//__use_timers

    /*------------------------------ zip/div mesh subface */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            if (_iter  < _opts.iter() )
            if (_opts.zip_() ||
                _opts.div_() )
            {
                if (_opts.verb() >= +3)
                    _dump.push(
                " CALL _ZIP-MESH...\n") ;

                _zip_mesh( _geom, _mesh ,
                    _hfun, _pred, _kern ,
                    _hval, _nset,
                    _nmrk, _emrk, _tmrk ,
                    _iter, _opts,
                    _TLIM, _DLIM,
                    _nzip, _ndiv) ;
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;

            _tcpu._zips_full +=
                _tcpu.time_span(_ttic, _ttoc);
    #       endif//__use_timers

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
            if (_nset.count() == 0) break ;
            if (_nmov == +0 &&
                _nzip == +0 &&
                _ndiv == +0 &&
                _nflp == +0 )       break ;
        }

        if (_opts.verb() >= +2)
        {
    /*------------------------------ print method metrics */
            _dump.push("\n");

            _dump.push(" MOVE-FULL: ");
            _dump.push(
            std::to_string(_tcpu._move_full)) ;
            _dump.push("\n");

            _dump.push(" TOPO-FULL: ");
            _dump.push(
            std::to_string(_tcpu._topo_full)) ;
            _dump.push("\n");

            _dump.push(" ZIPS-FULL: ");
            _dump.push(
            std::to_string(_tcpu._zips_full)) ;
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



