
    /*
    --------------------------------------------------------
     * ITER-DUAL-2: optim. schemes to reposition dual.
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
     * Last updated: 25 Feb., 2025
     *
     * Copyright 2013-2025
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_k.hpp


    /*
    --------------------------------------------------------
     * DQDW-TRI3: "local-ascent" node movement vector.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdw_tri3_k (
        mesh_type &_mesh ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type  _COST ,
        real_type  _STEP ,
        real_type  _SAVE ,
        real_type &_DQDW ,
        real_type &_dmin ,
        real_type &_dmax
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        real_type static const _WINC =
            std::pow(std::numeric_limits
                <real_type>::epsilon(), .50) ;

    /*------------------------------------- cell indexing */
        auto _ppos = &_node->pval(0);

        auto _tptr =
       _mesh. tri3().head() + _cell ;

    /*------------------ local iter. on finite-diff. step */
        real_type _wdel = _WINC * _STEP;
        real_type _dsup , _dsdn , _dqdw;
        {
    /*------------------ centred finite-diff. for dQ / dw */
            _ppos[_last] = _SAVE  + _wdel;

            real_type _scr1 =
                pred_type::tri3_cost (
               &_mesh. node(
                _tptr->node(0)).pval(0),
               &_mesh. node(
                _tptr->node(1)).pval(0),
               &_mesh. node(
                _tptr->node(2)).pval(0),
                typename
                pred_type::dual_kind ()
                    ) ;

            _dsup = (_scr1-_COST) / _wdel;

            _ppos[_last] = _SAVE  - _wdel;

            real_type _scr0 =
                pred_type::tri3_cost (
               &_mesh. node(
                _tptr->node(0)).pval(0),
               &_mesh. node(
                _tptr->node(1)).pval(0),
               &_mesh. node(
                _tptr->node(2)).pval(0),
                typename
                pred_type::dual_kind ()
                    ) ;

            _dsdn = (_COST-_scr0) / _wdel;
        }

        _ppos [_last]  = _SAVE ;

    /*------------------ finalise gradient and accumulate */
    //  if (_dsup * _dsdn > (real_type)0.)
        {
            _dqdw  = +.5 * (_dsup + _dsdn) ;

            _DQDW += _dqdw ;

            _dmin = std::min(_dmin, _dqdw) ;
            _dmax = std::max(_dmax, _dqdw) ;
        }
    }

    /*
    --------------------------------------------------------
     * DQDW-QUAD: "local-ascent" node movement vector.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdw_quad_k (
        mesh_type &_mesh ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type  _COST ,
        real_type  _STEP ,
        real_type  _SAVE ,
        real_type &_DQDW ,
        real_type &_dmin ,
        real_type &_dmax
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        real_type static const _WINC =
            std::pow(std::numeric_limits
                <real_type>::epsilon(), .50) ;

    /*------------------------------------- cell indexing */
        auto _ppos = &_node->pval(0);

        auto _qptr =
       _mesh. quad().head() + _cell ;

    /*------------------ local iter. on finite-diff. step */
        real_type _wdel = _WINC * _STEP;
        real_type _dsup , _dsdn , _dqdw;
        {
    /*------------------ centred finite-diff. for dQ / dw */
            _ppos[_last] = _SAVE  + _wdel;

            real_type _scr1 =
                pred_type::quad_cost (
               &_mesh. node(
                _qptr->node(0)).pval(0),
               &_mesh. node(
                _qptr->node(1)).pval(0),
               &_mesh. node(
                _qptr->node(2)).pval(0),
               &_mesh. node(
                _qptr->node(3)).pval(0),
                typename
                pred_type::dual_kind ()
                    ) ;

            _dsup = (_scr1-_COST) / _wdel;

            _ppos[_last] = _SAVE  - _wdel;

            real_type _scr0 =
                pred_type::quad_cost (
               &_mesh. node(
                _qptr->node(0)).pval(0),
               &_mesh. node(
                _qptr->node(1)).pval(0),
               &_mesh. node(
                _qptr->node(2)).pval(0),
               &_mesh. node(
                _qptr->node(3)).pval(0),
                typename
                pred_type::dual_kind ()
                    ) ;

            _dsdn = (_COST-_scr0) / _wdel;
        }

        _ppos [_last]  = _SAVE ;

    /*------------------ finalise gradient and accumulate */
    //  if (_dsup * _dsdn > (real_type)0.)
        {
            _dqdw  = +.5 * (_dsup + _dsdn) ;

            _DQDW += _dqdw ;

            _dmin = std::min(_dmin, _dqdw) ;
            _dmax = std::max(_dmax, _dqdw) ;
        }
    }

    /*
    --------------------------------------------------------
     * DQDW-MOVE: "local-ascent" node movement vector.
    --------------------------------------------------------
     */

    // Shift weight via gradients of quality functions

    // wnew <=== wold + beta * BEST[ d./dw Q(x,w)|_i ]

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdw_move_2 (
        mesh_type &_mesh ,
        conn_list &_conn ,
        node_iter  _node ,
        fp32_list &_cost ,
        real_type &_step ,
        real_type &_wadj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        real_type static const _WINC =
            std::pow(std::numeric_limits
                <real_type>::epsilon(), .50) ;

        _step = (real_type) 0. ;
        _wadj = (real_type) 0. ;

    /*------------------ calc. local characteristic scale */
        float     _qmin =
       +std::numeric_limits<float    >::infinity();

        real_type _bmin[_last] = {
       +std::numeric_limits<real_type>::infinity()
            } ;
        real_type _bmax[_last] = {
       -std::numeric_limits<real_type>::infinity()
            } ;

        iptr_type _cnum  = + 0 ;
        for (auto _next  = _conn.head() ,
                  _tend  = _conn.tend() ;
                  _next != _tend ;
                ++_next, ++_cnum )
        {
             auto _cell  = _next->_cell ;
             auto _kind  = _next->_kind ;

            _qmin = std::min(
                    _qmin, _cost [_cnum]) ;

            if (_kind == TRIA3_tag)
            {
                halo_tri3_k(_mesh,
                    _cell , _node, _bmin, _bmax,
                    _wadj ) ;
            }
            else
            if (_kind == QUAD4_tag)
            {
                halo_quad_k(_mesh,
                    _cell , _node, _bmin, _bmax,
                    _wadj ) ;
            }
        }

        _wadj = _wadj / _cnum ;

    /*------------------ adj. gradients wrt. pos: dQ / dx */
        float     _qlim = _qmin + (float)+1.0E-002 ;

        real_type _dmax =
       -std::numeric_limits<real_type>::infinity() ;
        real_type _dmin =
       +std::numeric_limits<real_type>::infinity() ;

        real_type _DQDW = (real_type) 0.;

        real_type _SAVE =
            _node->pval  (_last) ;

        real_type _qbar = (real_type) 1.;
        real_type _qlow = (real_type) 0.;

        iptr_type _lnum  = 0, _hnum = 1, _nnum = 0 ;
        for (auto _next  = _conn.head() ,
                  _tend  = _conn.tend() ;
                  _next != _tend ;
                ++_next, ++_nnum )
        {
             auto _cell  = _next->_cell ;
             auto _kind  = _next->_kind ;

             auto _COST  = _cost [_nnum];

            if (_COST <= _qlim)
            {
                if (_kind == TRIA3_tag)
                {
                dqdw_tri3_k(_mesh, _cell,
                    _node , _COST, _wadj, _SAVE,
                    _DQDW , _dmin, _dmax) ;
                }
                else
                if (_kind == QUAD4_tag)
                {
                dqdw_quad_k(_mesh, _cell,
                    _node , _COST, _wadj, _SAVE,
                    _DQDW , _dmin, _dmax) ;
                }

                _qlow += _COST ; _lnum += +1 ;
            }
            else
            {
                _qbar += _COST ; _hnum += +1 ;
            }
        }

        if (_dmin  *
            _dmax <=(real_type) +0. )  // check monotone?
            _DQDW  =(real_type) +0. ;

        if (_nnum != +0)
        {
            _DQDW /=  _lnum ;
            _qlow /=  _lnum ;
            _qbar /=  _hnum ;
        }

    /*------------------ 1st ord. Taylor-series step est. */
        real_type _scal = std::abs (_DQDW) ;

        if (_scal*_wadj <= _WINC)
        {
            _step = (real_type) +0. ;
        }
        else
        {
            _step = (_qbar-_qlow) / _DQDW  ;

            _step/= (real_type) +2. ;
        }

        if (_step > (real_type) +0. )
        {
            _step = std::min(_step,+_wadj) ;
        }
        else
        if (_step < (real_type) +0. )
        {
            _step = std::max(_step,-_wadj) ;
        }

    }



