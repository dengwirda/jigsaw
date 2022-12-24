
    /*
    --------------------------------------------------------
     * COST-MESH-2: utilities to manip. cell cost vectors.
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
     * Last updated: 10 Jun., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    /*
    --------------------------------------------------------
     * MOVE-OKAY: TRUE if state sufficiently good.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type move_okay (
        real_list &_cdst ,  // list of new adj. costs
        real_list &_csrc ,  // list of old adj. costs
        iptr_type &_move ,  // > +0 if move acceptable
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
            (real_type)1. / *_iter, +8);
        }
        for (auto _iter  = _cdst.head(),
                  _tend  = _cdst.tend();
                  _iter != _tend;
                ++_iter  )
        {
            _0dst  =
            std::min(_0dst, *_iter) ;

            _mdst += std::pow (
            (real_type)1. / *_iter, +8);
        }

        _msrc  = std::pow(
        _csrc.count() / _msrc, +1./8.0);
        _mdst  = std::pow(
        _cdst.count() / _mdst, +1./8.0);

        _qtol /=
        std::pow(_csrc.count(), 1./8.0);
        _qtol /=
        std::pow(_cdst.count(), 1./8.0);

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
        mesh_type &_mesh ,      // mesh object
        conn_list &_cset ,      // list of adj. cells
        real_list &_cost ,      // list of adj. costs
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
        mesh_type &_mesh ,      // mesh object
        conn_list &_cset ,      // list of adj. cells
        real_list &_cost ,      // list of adj. costs
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



