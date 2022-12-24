
    /*
    --------------------------------------------------------
     * MOVE-MESH-2: update cell/dual geometries.
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
     * Last updated: 11 Dec., 2022
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
     * MOVE-NODE: "limited" single node coord. update.
    --------------------------------------------------------
     */

    #include "iter_node_2.inc"

    template <
        typename  node_iter
             >
    __static_call
    __inline_call void_type move_node (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        iter_opts &_opts ,  // user options
        node_iter  _node ,  // iterator to node to iter.
        real_type *_last ,  // lagged line direction        
        kern_kind  _kern ,  // optim. kernel selector
        iptr_type &_move ,  // > 0 if move successful
        conn_list &_conn ,  // list of adj. cells
        real_list &_qold ,  // list of old adj. costs
        real_list &_qnew ,  // list of new adj. costs
        real_type  _QMIN ,  // minimum adj. old costs
        real_type  _QLIM ,  // cell quality threshold
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    /*---------------- try variational; fallback on dQ/dx */
        move_kern( _geom, _mesh, _hfun,
            _hval, _opts, _node, _last,
            _kern, _move, _conn,
            _qold, _qnew,
            _QMIN, _QLIM, _tcpu) ;  // variational

        if (_move >= +0 ) return ;

        move_kern( _geom, _mesh, _hfun,
            _hval, _opts, _node, _last,
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
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        iter_opts &_opts ,  // user options
        node_iter  _node ,  // iterator to node to iter.
        real_type *_last ,  // lagged line direction
        kern_kind  _kern ,  // optim. kernel selector
        iptr_type &_move ,  // > 0 if move successful
        conn_list &_conn ,  // list of adj. cells
        real_list &_qold ,  // list of old adj. costs
        real_list &_qnew ,  // list of new adj. costs
        real_type  _QMIN ,  // minimum adj. old costs
        real_type  _QLIM ,  // cell quality threshold
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    /*---------------- optimise single node's coordinates */
        iptr_type static
        constexpr _ITER = (iptr_type) +5  ;

    //  ifdef  __use_timers
    //  typename std ::chrono::
    //  high_resolution_clock::time_point  _ttic ;
    //  typename std ::chrono::
    //  high_resolution_clock::time_point  _ttoc ;
    //  typename std ::chrono::
    //  high_resolution_clock _time ;

    //  __unreferenced(_time) ; // why does MSVC need this??
    //  endif//__use_timers

        __unreferenced(_tcpu) ;

        _move = (iptr_type)-1 ;

        real_type _line[geom_dims] =
            {(real_type)+0.0} ;
        real_type _step[geom_dims] =
            {(real_type)+0.0} ;
        real_type _save[geom_dims] =
            {(real_type)+0.0} ;
        real_type _proj[geom_dims] =
            {(real_type)+0.0} ;

    /*---------------- calc. line search direction vector */

    //  ifdef  __use_timers
    //  _ttic = _time.now() ;
    //  endif//__use_timers

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

    //  ifdef  __use_timers
    //  _ttoc = _time.now() ;
    //  _tcpu._ldir_node += _tcpu.nano_span(_ttic, _ttoc) ;
    //  endif//__use_timers

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

        real_type _SCAL[_ITER] = {  // overrelaxation
            (real_type) std::sqrt( 2.0 ) ,
            (real_type) +1.00,
            (real_type) +0.50,
            (real_type) +0.25,
            (real_type) +.125  } ;

    /*---------------- do backtracking line search iter's */

    //  ifdef  __use_timers
    //  _ttic = _time.now() ;
    //  endif//__use_timers

        for (auto _idim =
        pred_type::geom_dims; _idim-- != +0; )
        {
            _save[_idim] = _ppos [_idim] ;
        }

        for (auto _iter = +0 ;
                _iter != _ITER; ++_iter )
        {
    /*---------------- push update along search direction */
            real_type _scal = _SCAL[_iter];

            real_type _zeta = _opts._zeta * 
                std::min ((real_type)1.,_scal) ;

            real_type _beta = _opts._beta ;

            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
    /*---------------- momentum-type bias for search dir. */
                _step[_idim] = 
            _last[_idim]*((real_type)0.+_beta) +
            _line[_idim]*((real_type)1.-_beta) ;
                
                _proj[_idim] =
                _save[_idim] + _scal * (
            _step[_idim]*((real_type)0.+_zeta) +
            _line[_idim]*((real_type)1.-_zeta));
            }

    /*---------------- ensure projection is on to surface */
            pred_type::
            proj_node (_geom, _save, _proj) ;

            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _ppos[_idim] = _proj[_idim] ;
            }

            real_type _XTOL = _xtol *   _scal;  // armijo
            real_type _XEPS = _xeps * 
                std::min ((real_type)1.,_scal) ;
            
            real_type _lmov =
            pred_type::length_sq(_save, _proj) ;

            if (_lmov <= _XEPS * _lsqr) break;

          //_move  = +1 ; return ;

    /*---------------- test quasi-monotonicity w.r.t. Q^T */
            _qnew.set_count(0) ;

            loop_cost( _mesh, _conn, _qnew,
                 cell_kind ()) ;

            move_okay( _qnew, _qold, _move,
                _QLIM, _opts.qtol(),
                _lmov, _XTOL* _lsqr) ;

            if (_move >= +0) break ;
        }

        if (_move <= +0)
        {
    /*---------------- swap with the saved coord. if fail */
            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _ppos[_idim] = _save[_idim] ;

                _last[_idim]/= (real_type)4.;   // reset!
            }
        }
        else
        {
    /*---------------- save line dir. for momentum update */
            for (auto _idim =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _last[_idim] = _step[_idim] ;
            }
        }

    //  ifdef  __use_timers
    //  _ttoc = _time.now() ;
    //  _tcpu._lopt_node += _tcpu.nano_span(_ttic, _ttoc) ;
    //  endif//__use_timers

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
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        iter_opts &_opts ,  // user options
        node_iter  _node ,  // iterator to node to iter.
        real_type  _last ,  // lagged line direction
        iptr_type &_move ,  // > 0 if move successful
        conn_list &_conn ,  // list of adj. cells
        real_list &_dold ,  // list of old adj. costs
        real_list &_dnew ,  // list of new adj. costs
        real_type  _DMIN ,  // minimum adj. old costs
        real_type  _DLIM ,  // dual quality threshold
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    /*---------------- optimise single node's coordinates */
        iptr_type static
        constexpr _ITER = (iptr_type) +5  ;

    //  ifdef  __use_timers
    //  typename std ::chrono::
    //  high_resolution_clock::time_point  _ttic ;
    //  typename std ::chrono::
    //  high_resolution_clock::time_point  _ttoc ;
    //  typename std ::chrono::
    //  high_resolution_clock _time ;

    //  __unreferenced(_time) ; // why does MSVC need this??
    //  endif//__use_timers

        __unreferenced(_geom) ;
        __unreferenced(_hfun) ;
        __unreferenced(_hval) ;
        __unreferenced(_tcpu) ;

        _move = (iptr_type)-1 ;

        real_type _wadj, _line, _step;

    /*---------------- calc. line search direction vector */

    //  ifdef  __use_timers
    //  _ttic = _time.now() ;
    //  endif//__use_timers

        if (_DMIN <= _DLIM)
        {
            dqdw_move_2 (
                _mesh, _conn, _node, _dold,
                _line, _wadj) ;
        }
        else { return ; }

    //  ifdef  __use_timers
    //  _ttoc = _time.now() ;
    //  _tcpu._ldir_dual += _tcpu.nano_span(_ttic, _ttoc) ;
    //  endif//__use_timers

    /*---------------- scale line search direction vector */
        real_type _weps =           // delta_w ~= 0.0
            (real_type)+.01*_opts.qtol() ;

        real_type _save = 
            _node->pval(
                pred_type::real_dims - 1);

        real_type _SCAL[_ITER] = {  // overrelaxation
            (real_type) std::sqrt( 2.0 ) ,
            (real_type) +1.00,
            (real_type) +0.50,
            (real_type) +0.25,
            (real_type) +.125  } ;
            
    /*---------------- do backtracking line search iter's */

    //  ifdef  __use_timers
    //  _ttic = _time.now() ;
    //  endif//__use_timers

        for (auto _iter = +0 ;
                _iter != _ITER; ++_iter )
        {
            real_type _scal = _SCAL[_iter];

            real_type _zeta = _opts._zeta * 
                (real_type)3./4. * 
            std::min((real_type)1.,_scal) ;

            real_type _beta = _opts._beta ;

    /*---------------- momentum-type bias for search dir. */
            _step = 
            _last * ((real_type)0.+_beta) +
            _line * ((real_type)1.-_beta) ;
                
            _node->pval(real_dims-1) =
                _save + _scal * (
            _step * ((real_type)0.+_zeta) +
            _line * ((real_type)1.-_zeta));

            _node->pval(real_dims-1) =
                std::max(-_wadj ,
            _node->pval(real_dims-1));

            _node->pval(real_dims-1) =
                std::min(+_wadj ,
            _node->pval(real_dims-1));

            real_type _wmov =
                std::abs (_save -
            _node->pval(real_dims-1));

            real_type _wtol = _weps * _wadj *
            std::min ((real_type)1.0, _scal);

            if (_wmov <= _wtol) break;

    /*---------------- test quasi-monotonicity w.r.t. Q^D */
            _dnew.set_count(0) ;

            loop_cost( _mesh, _conn, _dnew ,
                 dual_kind ()) ;

            move_okay( _dnew, _dold, _move ,
                _DLIM, _opts.qtol()) ;

            if (_move >= +0) break ;
        }

    /*---------------- swap with the saved coord. if fail */
        if (_move <= +0)
        {
            _node->pval(real_dims-1) = _save ;

            _last/= (real_type) 4. ;  // reset dir., fail
        }
        else
        {
            _last = _step ;
        }

    //  ifdef  __use_timers
    //  _ttoc = _time.now() ;
    //  _tcpu._lopt_dual += _tcpu.nano_span(_ttic, _ttoc) ;
    //  endif//__use_timers

    }

    /*
    --------------------------------------------------------
     * MOVE-NODE: do a single node smoothing pass.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type move_node (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        conn_sets &_CONN ,  // cached node adj. indexing
        hfun_type &_hfun ,  // spacing object
        char_type  _kern ,  // optim. kernel selector
        real_list &_hval ,  // cached node spacing array
        line_list &_last ,  // lagged line direction
        iptr_list &_nset ,  // list of nodes active
        iptr_list &_aset ,  // space for new active list
        mark_list &_mark ,  // sets of node/cell flags
        part_data &_part ,  // mesh parallel partition
        iptr_type  _iout ,  // outer iteration marker
        iptr_type  _isub ,  // inner iteration marker
        iter_opts &_opts ,  // user options
        real_type  _QLIM ,  // cell quality threshold
        real_type  _DLIM ,  // dual quality threshold
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    #   define MARK(_NODE) _mark._node[_NODE]

        class move_task     // inner loop: MOVE-NODE
            {
        public  :
        __inline_call void operator() (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        iptr_type  _apos ,  // index of node to update
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        conn_sets &_CONN ,  // cached node adj. indexing
        char_type  _kern ,  // optim. kernel selector
        line_list &_last ,  // lagged line direction
        iptr_list &_nset ,  // list of nodes active
        mark_list &_mark ,  // sets of node/cell flags
        iptr_type  _iout ,  // outer iteration marker
        iter_opts &_opts ,  // user options
        conn_list &_conn ,  // space for connectivity
        real_list &_qold ,  // space for adj. costs
        real_list &_qnew ,  // space for adj. costs
        real_type  _QLIM ,  // cell quality threshold
        iter_stat &_tcpu    // CPU timers/info
                ) const
            {
        /*---------------- apply single move-node updates */
             auto _node  =
            _mesh.node().head()+ _apos ;

        /*---------------- assemble a local tria. stencil */
            _conn.set_count( +0) ;
            _conn.push_tail(
            _CONN._adj2.head()+_CONN._idx2[_apos+0],
            _CONN._adj2.head()+_CONN._idx2[_apos+1]
                ) ;

            if (_conn.empty() ) return ;

            if (_mark._node[_apos] >= +0)
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
                    _opts, _node,
                   &_last[ _apos][ +0 ], 
                    _kern, _move, _conn,
                    _qold, _qnew,
                    _QMIN, _QLIM, _tcpu) ;
            }

            if (_move > +0 )
            {
        /*---------------- update when state is improving */
            _hval[_apos] = (real_type)-1 ;

            if (std::abs(
                MARK( _apos )) !=_iout)
            {
                if (MARK( _apos ) >= 0)
                MARK( _apos ) = +_iout;
                else
                MARK( _apos ) = -_iout;

            //  pragma omp critical(NS)
                {
                _nset.push_tail (_apos)  ;
                }
            }

            }
            }
            }
            } ;

        iptr_type constexpr PART_SCAL = +32  ;
        iptr_type constexpr BNDS_SCAL = 4096 ;

    #   define SBEG(_PASS)(_PASS == +0 ?    \
            _aset.head() : _aset.tail())

    #   define SEND(_PASS)(_PASS == +0 ?    \
            _aset.tend() : _aset.hend())

    #   define PBEG(_PART, _PASS, _RANK)(   \
            _PASS == +0 ? \
            _PART._lptr[_RANK + 0] + 0 :    \
            _PART._lptr[_RANK + 1] - 1)

    #   define PEND(_PART, _PASS, _RANK)(   \
            _PASS == +0 ? \
            _PART._lptr[_RANK + 1] + 0 :    \
            _PART._lptr[_RANK + 0] - 1)

    #   define ITER(_PASS)(_PASS == +0 ? +1 : -1)

        iptr_type _pful, _pbnd;
        conn_list _conn;
        real_list _qold, _qnew;

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

        if (_isub == +0)
        sort_node(_mesh, _CONN, _aset,
            _mark._node, 
            _iout,_isub, _QLIM, _DLIM, _opts) ;

        _pful = std::sqrt(_aset.count ()) / 
                PART_SCAL ;
        _pful = std::min(
                _opts. nprt(), _pful + 1) ;
    
    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

        if (_pful <= +1 )
        {
    /*-------------------- single part: sequential scheme */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

    /*-------------------- SYMMETRIC GAUSS-SEIDEL on CELL */

        if (_opts .tria())  // 1=>GS, 2=>SGS
        {
        for (auto _pass  = 0; _pass < 2; ++_pass)
        for (auto _apos  = SBEG(_pass) ;
                  _apos != SEND(_pass) ;
                  _apos += ITER(_pass) )
        {
            move_task()(_geom, _mesh,*_apos,
                _hfun , _hval,
                _CONN , _kern, _last,
                _nset , _mark,
                _iout , _opts,
                _conn , _qold, _qnew, _QLIM, 
                _tcpu ) ;
        }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

        }
        else
        {
    /*-------------------- multiple-part: thread-parallel */    

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        if (_isub == +0)
        part_mesh(_mesh, _part._full, _pful,
                  _aset, 
                  _isub, _tcpu) ;

        _pbnd = _part._full._seqs.count() /
                BNDS_SCAL ;
        _pbnd = std::min(
                _opts. nprt(), _pbnd + 1) ;

        if (_isub == +0)
        part_mesh(_mesh, _part._bnds, _pbnd,
                  _part. _full._seqs,
                  _isub, _tcpu) ;

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._part_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers
    
    /*-------------------- SYMMETRIC GAUSS-SEIDEL on CELL */

        if (_opts .tria())  // 1=>GS, 2=>SGS
        {
    /*--------------------------------- purely sequential */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _apos  = 
                  _part._bnds._seqs.head();
                  _apos != 
                  _part._bnds._seqs.tend();
                  _apos += +1)
        {
            move_task()(_geom, _mesh,*_apos,
                _hfun , _hval,
                _CONN , _kern, _last,
                _nset , _mark,
                _iout , _opts,
                _conn , _qold, _qnew, _QLIM, 
                _tcpu ) ;
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._seqs_node += _tcpu.time_span(_ttic, _ttoc) ;
        _tcpu._core_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*--------------------------------- parallel on outer */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

    #   ifdef  __use_openmp
        omp_set_num_threads(
            _part._bnds._lptr.count()-1 ) ;
    #   endif//__use_openmp

    #   pragma omp parallel default(none) \
           shared(_part, _geom, _mesh, \
                  _hfun, _hval, \
                  _CONN, _kern, _last, \
                  _nset, _mark, \
                  _iout, _isub, _opts, \
                  _QLIM, _DLIM, _tcpu)
        {
        conn_list _CSET; // thread-local arrays
        real_list _QOLD, _QNEW;
        iptr_list _NSET;
        iptr_type _pass  = 0u ;
        iptr_type _anum; 
    #   pragma omp for nowait  schedule(static,1)
        for (auto _rank  = 0u ; 
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
        for (auto _apos  = PBEG(_part. _bnds, 
                                _pass, _rank) ;
                  _apos != PEND(_part. _bnds, 
                                _pass, _rank) ;
                  _apos += ITER(_pass) )
        {
            _anum = _part._bnds._list [_apos] ;

            if (_part._bnds._part[_anum] < 0) 
                continue ; // skip if interface

            move_task()(_geom, _mesh, _anum,
                _hfun , _hval,
                _CONN , _kern, _last,
                _NSET , _mark,
                _iout , _opts,
                _CSET , _QOLD, _QNEW, _QLIM, 
                _tcpu ) ;
        }
        }

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank  = 0u ;
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
    #       pragma omp ordered
            {
            _nset.push_tail(
                _NSET.head(), _NSET.tend());
            }
        }
        }

    /*--------------------------------- parallel on inner */
    #   ifdef  __use_openmp
        omp_set_num_threads(
            _part._full._lptr.count()-1 ) ;
    #   endif//__use_openmp

    #   pragma omp parallel default(none) \
           shared(_part, _geom, _mesh, \
                  _hfun, _hval, \
                  _CONN, _kern, _last, \
                  _nset, _mark, \
                  _iout, _isub, _opts, \
                  _QLIM, _DLIM, _tcpu)
        {
        conn_list _CSET; // thread-local arrays
        real_list _QOLD, _QNEW;
        iptr_list _NSET;
        iptr_type _anum;
    #   pragma omp for nowait  schedule(static,1)
        for (auto _rank  = 0u ; 
                  _rank  <
            _part._full._lptr.count()-1; ++_rank)
        {
        for (auto _pass  = 0; _pass < 2; ++_pass)
        for (auto _apos  = PBEG(_part. _full, 
                                _pass, _rank) ;
                  _apos != PEND(_part. _full, 
                                _pass, _rank) ;
                  _apos += ITER(_pass) )
        {
            _anum = _part._full._list [_apos] ;

            if (_part._full._part[_anum] < 0) 
                continue ; // skip if interface

            move_task()(_geom, _mesh, _anum,
                _hfun , _hval,
                _CONN , _kern, _last,
                _NSET , _mark,
                _iout , _opts,
                _CSET , _QOLD, _QNEW, _QLIM, 
                _tcpu ) ;
        }
        }

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank  = 0u ;
                  _rank  <
            _part._full._lptr.count()-1; ++_rank)
        {
    #       pragma omp ordered
            {
            _nset.push_tail(
                _NSET.head(), _NSET.tend());
            }
        }
        }

    /*--------------------------------- parallel on outer */
    #   ifdef  __use_openmp
        omp_set_num_threads(
            _part._bnds._lptr.count()-1 ) ;
    #   endif//__use_openmp

    #   pragma omp parallel default(none) \
           shared(_part, _geom, _mesh, \
                  _hfun, _hval, \
                  _CONN, _kern, _last, \
                  _nset, _mark, \
                  _iout, _isub, _opts, \
                  _QLIM, _DLIM, _tcpu)
        {
        conn_list _CSET; // thread-local arrays
        real_list _QOLD, _QNEW;
        iptr_list _NSET;
        iptr_type _pass  = 1u ;
        iptr_type _anum; 
    #   pragma omp for nowait  schedule(static,1)
        for (auto _rank  = 0u ; 
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
        for (auto _apos  = PBEG(_part. _bnds, 
                                _pass, _rank) ;
                  _apos != PEND(_part. _bnds, 
                                _pass, _rank) ;
                  _apos += ITER(_pass) )
        {
            _anum = _part._bnds._list [_apos] ;

            if (_part._bnds._part[_anum] < 0) 
                continue ; // skip if interface

            move_task()(_geom, _mesh, _anum,
                _hfun , _hval,
                _CONN , _kern, _last,
                _NSET , _mark,
                _iout , _opts,
                _CSET , _QOLD, _QNEW, _QLIM, 
                _tcpu ) ;
        }
        }

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank  = 0u ;
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
    #       pragma omp ordered
            {
            _nset.push_tail(
                _NSET.head(), _NSET.tend());
            }
        }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._para_node += _tcpu.time_span(_ttic, _ttoc) ;
        _tcpu._core_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*--------------------------------- purely sequential */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _apos  = 
                  _part._bnds._seqs.tail();
                  _apos != 
                  _part._bnds._seqs.hend();
                  _apos += -1)
        {
            move_task()(_geom, _mesh,*_apos,
                _hfun , _hval,
                _CONN , _kern, _last,
                _nset , _mark,
                _iout , _opts,
                _conn , _qold, _qnew, _QLIM, 
                _tcpu ) ;
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._seqs_node += _tcpu.time_span(_ttic, _ttoc) ;
        _tcpu._core_node += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

        }    //if opts.tria
    
        }    //if _nptr > 1

    #   undef   MARK
    #   undef   SBEG
    #   undef   SEND
    #   undef   PBEG
    #   undef   PEND
    #   undef   ITER
    }

    /*
    --------------------------------------------------------
     * MOVE-DUAL: do a single dual smoothing pass.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type move_dual (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        conn_sets &_CONN ,  // cached node adj. indexing
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        real_list &_last ,  // lagged line direction
        iptr_list &_nset ,  // list of nodes active
        iptr_list &_aset ,  // space for new active list
        mark_list &_mark ,  // sets of node/cell flags
        part_data &_part ,  // mesh parallel partition        
        iptr_type  _iout ,  // outer iteration marker
        iptr_type  _isub ,  // inner iteration marker
        iter_opts &_opts ,  // user options
        real_type  _QLIM ,  // cell quality threshold
        real_type  _DLIM ,  // dual quality threshold
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    #   define MARK(_NODE) _mark._node[_NODE]

        class move_task     // inner loop: MOVE-DUAL
            {
        public  :
        __inline_call void operator() (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        iptr_type  _apos ,  // index of node to update
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        conn_sets &_CONN ,  // cached node adj. indexing
        real_list &_last ,  // lagged line direction
        iptr_list &_nset ,  // list of nodes active
        mark_list &_mark ,  // sets of node/cell flags
        iptr_type  _iout ,  // outer iteration marker
        iter_opts &_opts ,  // user options
        conn_list &_conn ,  // space for connectivity
        real_list &_dold ,  // space for adj. costs
        real_list &_dnew ,  // space for adj. costs
        real_type  _DLIM ,  // dual quality threshold
        iter_stat &_tcpu    // CPU timers/info
                ) const
            {
        /*---------------- apply single move-dual updates */
             auto _node  =
            _mesh.node().head()+ _apos ;

        /*---------------- assemble a local tria. stencil */
            _conn.set_count( +0) ;
            _conn.push_tail(
            _CONN._adj2.head()+_CONN._idx2[_apos+0],
            _CONN._adj2.head()+_CONN._idx2[_apos+1]
                ) ;

            if (_conn.empty() ) return ;

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
                    _last[ _apos] ,
                    _move, _conn,
                    _dold, _dnew,
                    _DMIN, _DLIM, _tcpu) ;
            }

            if (_move > +0 )
            {
        /*---------------- update when state is improving */
            if (std::abs(
                MARK( _apos )) !=_iout)
            {
                if (MARK( _apos ) >= 0)
                MARK( _apos ) = +_iout;
                else
                MARK( _apos ) = -_iout;

            //  pragma omp critical(NS)
                {
                _nset.push_tail (_apos)  ;
                }
            }

            }
            }
            } ;

        iptr_type constexpr PART_SCAL = +32 ;
        iptr_type constexpr BNDS_SCAL = 4096;

    #   define SBEG(_PASS)(_PASS == +0 ?    \
            _aset.head() : _aset.tail())

    #   define SEND(_PASS)(_PASS == +0 ?    \
            _aset.tend() : _aset.hend())

    #   define PBEG(_PART, _PASS, _RANK)(   \
            _PASS == +0 ? \
            _PART._lptr[_RANK + 0] + 0 :    \
            _PART._lptr[_RANK + 1] - 1)

    #   define PEND(_PART, _PASS, _RANK)(   \
            _PASS == +0 ? \
            _PART._lptr[_RANK + 1] + 0 :    \
            _PART._lptr[_RANK + 0] - 1)

    #   define ITER(_PASS)(_PASS == +0 ? +1 : -1)

        iptr_type _pful, _pbnd;
        conn_list _conn;
        real_list _dold, _dnew;

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

        if (_isub == +0)
        sort_node(_mesh, _CONN, _aset,
            _mark._node, 
            _iout,_isub, _QLIM, _DLIM, _opts) ;

        _pful = std::sqrt(_aset.count ()) / 
                PART_SCAL ;
        _pful = std::min(
                _opts. nprt(), _pful + 1) ;

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

        if (_pful <= +1 )
        {
    /*-------------------- single part: sequential scheme */

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

    /*-------------------- SYMMETRIC GAUSS-SEIDEL on DUAL */

        if (_opts .dual())  // 1=>GS, 2=>SGS
        {
        for (auto _pass  = 0; _pass < 2; ++_pass)
        for (auto _apos  = SBEG(_pass) ;
                  _apos != SEND(_pass) ;
                  _apos += ITER(_pass) )
        {
            move_task()(_geom, _mesh,*_apos,
                _hfun , _hval,
                _CONN , _last,
                _nset , _mark,
                _iout , _opts,
                _conn , _dold, _dnew, _DLIM, 
                _tcpu ) ;
        }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

        }
        else
        {
    /*-------------------- multiple-part: thread-parallel */    

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        if (_isub == +0)
        part_mesh(_mesh, _part._full, _pful,
                  _aset, 
                  _isub, _tcpu) ;

        _pbnd = _part._full._seqs.count() /
                BNDS_SCAL ;
        _pbnd = std::min(
                _opts. nprt(), _pbnd + 1) ;

        if (_isub == +0)
        part_mesh(_mesh, _part._bnds, _pbnd,
                  _part. _full._seqs,
                  _isub, _tcpu) ;

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._part_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*-------------------- SYMMETRIC GAUSS-SEIDEL on DUAL */

        if (_opts .dual())  // 1=>GS, 2=>SGS
        {
    /*--------------------------------- purely sequential */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _apos  = 
                  _part._bnds._seqs.head();
                  _apos != 
                  _part._bnds._seqs.tend();
                  _apos += +1)
        {
            move_task()(_geom, _mesh,*_apos,
                _hfun , _hval,
                _CONN , _last,
                _nset , _mark,
                _iout , _opts,
                _conn , _dold, _dnew, _DLIM, 
                _tcpu ) ;
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._seqs_dual += _tcpu.time_span(_ttic, _ttoc) ;
        _tcpu._core_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*--------------------------------- parallel on outer */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

    #   ifdef  __use_openmp
        omp_set_num_threads(
            _part._bnds._lptr.count()-1 ) ;
    #   endif//__use_openmp

    #   pragma omp parallel default(none) \
           shared(_part, _geom, _mesh, \
                  _hfun, _hval, \
                  _CONN, _last, \
                  _nset, _mark, \
                  _iout, _isub, _opts, \
                  _QLIM, _DLIM, _tcpu)
        {
        conn_list _CSET; // thread-local arrays
        real_list _DOLD, _DNEW;
        iptr_list _NSET;
        iptr_type _pass  = 0u ;
        iptr_type _anum; 
    #   pragma omp for nowait  schedule(static,1)
        for (auto _rank  = 0u ; 
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
        for (auto _apos  = PBEG(_part. _bnds, 
                                _pass, _rank) ;
                  _apos != PEND(_part. _bnds, 
                                _pass, _rank) ;
                  _apos += ITER(_pass) )
        {
            _anum = _part._bnds._list [_apos] ;

            if (_part._bnds._part[_anum] < 0) 
                continue ; // skip if interface

            move_task()(_geom, _mesh, _anum,
                _hfun , _hval,
                _CONN , _last,
                _NSET , _mark,
                _iout , _opts,
                _CSET , _DOLD, _DNEW, _DLIM, 
                _tcpu ) ;
        }
        }

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank  = 0u ;
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
    #       pragma omp ordered
            {
            _nset.push_tail(
                _NSET.head(), _NSET.tend());
            }
        }
        }

    /*--------------------------------- parallel on inner */
    #   ifdef  __use_openmp
        omp_set_num_threads(
            _part._full._lptr.count()-1 ) ;
    #   endif//__use_openmp

    #   pragma omp parallel default(none) \
           shared(_part, _geom, _mesh, \
                  _hfun, _hval, \
                  _CONN, _last, \
                  _nset, _mark, \
                  _iout, _isub, _opts, \
                  _QLIM, _DLIM, _tcpu)
        {
        conn_list _CSET; // thread-local arrays
        real_list _DOLD, _DNEW;
        iptr_list _NSET;
        iptr_type _anum;
    #   pragma omp for nowait  schedule(static,1)
        for (auto _rank  = 0u ; 
                  _rank  <
            _part._full._lptr.count()-1; ++_rank)
        {
        for (auto _pass  = 0; _pass < 2; ++_pass)
        for (auto _apos  = PBEG(_part. _full, 
                                _pass, _rank) ;
                  _apos != PEND(_part. _full, 
                                _pass, _rank) ;
                  _apos += ITER(_pass) )
        {
            _anum = _part._full._list [_apos] ;

            if (_part._full._part[_anum] < 0) 
                continue ; // skip if interface

            move_task()(_geom, _mesh, _anum,
                _hfun , _hval,
                _CONN , _last,
                _NSET , _mark,
                _iout , _opts,
                _CSET , _DOLD, _DNEW, _DLIM, 
                _tcpu ) ;
        }
        }

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank  = 0u ;
                  _rank  <
            _part._full._lptr.count()-1; ++_rank)
        {
    #       pragma omp ordered
            {
            _nset.push_tail(
                _NSET.head(), _NSET.tend());
            }
        }
        }

    /*--------------------------------- parallel on outer */
    #   ifdef  __use_openmp
        omp_set_num_threads(
            _part._bnds._lptr.count()-1 ) ;
    #   endif//__use_openmp

    #   pragma omp parallel default(none) \
           shared(_part, _geom, _mesh, \
                  _hfun, _hval, \
                  _CONN, _last, \
                  _nset, _mark, \
                  _iout, _isub, _opts, \
                  _QLIM, _DLIM, _tcpu)
        {
        conn_list _CSET; // thread-local arrays
        real_list _DOLD, _DNEW;
        iptr_list _NSET;
        iptr_type _pass  = 1u ;
        iptr_type _anum; 
    #   pragma omp for nowait  schedule(static,1)
        for (auto _rank  = 0u ; 
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
        for (auto _apos  = PBEG(_part. _bnds, 
                                _pass, _rank) ;
                  _apos != PEND(_part. _bnds, 
                                _pass, _rank) ;
                  _apos += ITER(_pass) )
        {
            _anum = _part._bnds._list [_apos] ;

            if (_part._bnds._part[_anum] < 0) 
                continue ; // skip if interface

            move_task()(_geom, _mesh, _anum,
                _hfun , _hval,
                _CONN , _last,
                _NSET , _mark,
                _iout , _opts,
                _CSET , _DOLD, _DNEW, _DLIM, 
                _tcpu ) ;
        }
        }

    #   pragma omp for ordered schedule(static,1)
        for (auto _rank  = 0u ;
                  _rank  <
            _part._bnds._lptr.count()-1; ++_rank)
        {
    #       pragma omp ordered
            {
            _nset.push_tail(
                _NSET.head(), _NSET.tend());
            }
        }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._para_dual += _tcpu.time_span(_ttic, _ttoc) ;
        _tcpu._core_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    /*--------------------------------- purely sequential */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _apos  = 
                  _part._bnds._seqs.tail();
                  _apos != 
                  _part._bnds._seqs.hend();
                  _apos += -1)
        {
            move_task()(_geom, _mesh,*_apos,
                _hfun , _hval,
                _CONN , _last,
                _nset , _mark,
                _iout , _opts,
                _conn , _dold, _dnew, _DLIM, 
                _tcpu ) ;
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._seqs_dual += _tcpu.time_span(_ttic, _ttoc) ;
        _tcpu._core_dual += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

        }    //if opts.dual

        }    //if _nptr > 1

    #   undef   MARK
    #   undef   SBEG
    #   undef   SEND
    #   undef   PBEG
    #   undef   PEND
    #   undef   ITER
    }



