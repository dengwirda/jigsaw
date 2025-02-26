
    /*
    --------------------------------------------------------
     * ITER-NODE-2: optim. schemes to reposition nodes.
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
     * ~ODT-TRI3: optimal delaunay tessellation update.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _odt_tri3_k (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type *_move ,
        real_type &_wsum ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0;

    #   define EVALHFUN( _NODE, _NPTR ) \
        if (_hval[_NODE] < (real_type)+0.)  \
            _hval[_NODE] = _hfun.eval ( \
           &_NPTR->pval(0), _NPTR->hidx())  ;

    #   define LOOPDIMS( _ITER )    \
            for ( auto _ITER =          \
        pred_type::geom_dims; _ITER-- != 0; )

        __unreferenced(_node) ;

    /*------------------------------------- cell indexing */
        auto _tptr =
       _mesh. tri3().head() + _cell ;

        auto _inod = _tptr->node(0) ;
        auto _jnod = _tptr->node(1) ;
        auto _knod = _tptr->node(2) ;

        auto _imid = _node->node(0) ;

        auto _iptr =
       _mesh. node().head() + _inod ;
        auto _jptr =
       _mesh. node().head() + _jnod ;
        auto _kptr =
       _mesh. node().head() + _knod ;

    /*------------------------------------- dual geometry */
        real_type _ball [_last + 1] ;
        pred_type::tri3_ball(_ball,
            &_iptr->pval(+0),
            &_jptr->pval(+0),
            &_kptr->pval(+0), true) ;

        real_type _tmag =
             pred_type::tri3_mass (
            &_iptr->pval(+0),
            &_jptr->pval(+0),
            &_kptr->pval(+0)) ;

        _tmag = std::max(
            (real_type)+0.0, _tmag) ;

    /*------------------------------------- odt weighting */
        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_ball[+0]) ;

        real_type _lsqr = std::max(
            _ball[_last], (real_type)+0.) ;

        EVALHFUN( _inod , _iptr )
        EVALHFUN( _jnod , _jptr )
        EVALHFUN( _knod , _kptr )

        real_type _hbal = _hfun.eval(
            &_ball[0], _node->hidx()) ;

        real_type _irho =
            _hval[_imid] / _hval[_inod] ;
        real_type _jrho =
            _hval[_imid] / _hval[_jnod] ;
        real_type _krho =
            _hval[_imid] / _hval[_knod] ;

        real_type _0rho =
            _hval[_imid] / _hbal        ;

        real_type _hrho =
            _irho * (real_type)2./9. +
            _jrho * (real_type)2./9. +
            _krho * (real_type)2./9. +
            _0rho * (real_type)3./9. ;

        _hrho =
            _tmag * std::pow(_hrho, +3) ;

        LOOPDIMS( _idim )
        {
            _move[_idim]
               += _hrho * _ball [_idim] ;
        }

        _wsum  += _hrho ;
        _ladj  += _lsqr ;

    #   undef   LOOPDIMS
    #   undef   EVALHFUN
    }

    /*
    --------------------------------------------------------
     * ~ODT-QUAD: optimal delaunay tessellation update.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _odt_quad_k (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type *_move ,
        real_type &_wsum ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0;

    #   define EVALHFUN( _NODE, _NPTR ) \
        if (_hval[_NODE] < (real_type)+0.)  \
            _hval[_NODE] = _hfun.eval ( \
           &_NPTR->pval(0), _NPTR->hidx())  ;

    #   define LOOPDIMS( _ITER )    \
            for ( auto _ITER =          \
        pred_type::geom_dims; _ITER-- != 0; )

        __unreferenced(_node) ;

    /*------------------------------------- cell indexing */
        auto _qptr =
       _mesh. quad().head() + _cell ;

        auto _inod = _qptr->node(0) ;
        auto _jnod = _qptr->node(1) ;
        auto _knod = _qptr->node(2) ;
        auto _lnod = _qptr->node(3) ;

        auto _imid = _node->node(0) ;

        auto _iptr =
       _mesh. node().head() + _inod ;
        auto _jptr =
       _mesh. node().head() + _jnod ;
        auto _kptr =
       _mesh. node().head() + _knod ;
        auto _lptr =
       _mesh. node().head() + _lnod ;

    /*------------------------------------- dual geometry */
        real_type _ball [_last + 1] ;
        pred_type::quad_ball(_ball,
            &_iptr->pval(+0),
            &_jptr->pval(+0),
            &_kptr->pval(+0),
            &_lptr->pval(+0), true) ;

        real_type _qmag =
             pred_type::quad_mass (
            &_iptr->pval(+0),
            &_jptr->pval(+0),
            &_kptr->pval(+0),
            &_lptr->pval(+0)) ;

        _qmag = std::max(
            (real_type)+0.0, _qmag) ;

    /*------------------------------------- odt weighting */
        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_ball[+0]) ;

        real_type _lsqr = std::max(
            _ball[_last], (real_type)+0.) ;

        EVALHFUN( _inod , _iptr )
        EVALHFUN( _jnod , _jptr )
        EVALHFUN( _knod , _kptr )
        EVALHFUN( _lnod , _lptr )

        real_type _hbal = _hfun.eval(
            &_ball[0], _node->hidx()) ;

        real_type _irho =
            _hval[_imid] / _hval[_inod] ;
        real_type _jrho =
            _hval[_imid] / _hval[_jnod] ;
        real_type _krho =
            _hval[_imid] / _hval[_knod] ;
        real_type _lrho =
            _hval[_imid] / _hval[_lnod] ;

        real_type _0rho =
            _hval[_imid] / _hbal        ;

        real_type _hrho =
            _irho * (real_type)2./12 +
            _jrho * (real_type)2./12 +
            _krho * (real_type)2./12 +
            _lrho * (real_type)2./12 +
            _0rho * (real_type)4./12 ;

        _hrho =
            _qmag * std::pow(_hrho, +3) ;

        LOOPDIMS( _idim )
        {
            _move[_idim]
               += _hrho * _ball [_idim] ;
        }

        _wsum  += _hrho ;
        _ladj  += _lsqr ;

    #   undef   LOOPDIMS
    #   undef   EVALHFUN
    }

    /*
    --------------------------------------------------------
     * ~ODT-MOVE: optimal delaunay tessellation update.
    --------------------------------------------------------
     */

    // Shift node to weighted mean of adj. dual vertex
    // positions...

    // xnew <=== SUM( |t_i|_r * c_i ) / SUM( |t_i|_r )
    //
    // with |t_i|_r = INT|i rho(x) dA ~ |t_i| / h_i**3

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _odt_move_2 (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        conn_list &_conn ,
        node_iter  _node ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0;

        real_type _move[_last + 0] = {
                (real_type) +0.0};

        _ladj = (real_type) +0.0 ;

        real_type _wsum =
       +std::numeric_limits<real_type>::epsilon();

        iptr_type _cnum = +0  ;

        for (auto _next  = _conn.head() ,
                  _tend  = _conn.tend() ;
                  _next != _tend ;
                ++_next, ++_cnum )
        {
             auto _cell  = _next->_cell ;
             auto _kind  = _next->_kind ;

            if (_kind == TRIA3_tag)
            {
                _odt_tri3_k(_geom, _mesh, _hfun,
                    _hval , _cell, _node, _move,
                    _wsum , _ladj ) ;
            }
            else
            if (_kind == QUAD4_tag)
            {
                _odt_quad_k(_geom, _mesh, _hfun,
                    _hval , _cell, _node, _move,
                    _wsum , _ladj ) ;
            }
        }

        if (_cnum > +0)
        {
            for (auto _idim =
            pred_type::geom_dims ; _idim-- != 0; )
            {
                _line[_idim] =
                    _move[_idim] / _wsum
                        - _node->pval(_idim) ;
            }

            _ladj = std::sqrt(_ladj / _cnum) ;
        }
    }

    /*
    --------------------------------------------------------
     * ~CVT-TRI3: optimal delaunay tessellation update.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _cvt_tri3_k (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type *_move ,
        real_type &_wsum ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0;

    #   define EVALHFUN( _NODE, _NPTR ) \
        if (_hval[_NODE] < (real_type)+0.)  \
            _hval[_NODE] = _hfun.eval ( \
           &_NPTR->pval(0), _NPTR->hidx())  ;

    #   define LOOPDIMS( _ITER )    \
            for ( auto _ITER =          \
        pred_type::geom_dims; _ITER-- != 0; )

    /*------------------------------------- cell indexing */
        auto _tptr =
       _mesh. tri3().head() + _cell ;

        auto _inod = _tptr->node(0) ;
        auto _jnod = _tptr->node(1) ;
        auto _knod = _tptr->node(2) ;

        if ( _jnod ==_node->node(0) )
        {
           _inod   = _tptr->node(1) ;
           _jnod   = _tptr->node(2) ;
           _knod   = _tptr->node(0) ;
        }
        else
        if ( _knod ==_node->node(0) )
        {
           _inod   = _tptr->node(2) ;
           _jnod   = _tptr->node(0) ;
           _knod   = _tptr->node(1) ;
        }

        auto _iptr =
       _mesh. node().head() + _inod ;
        auto _jptr =
       _mesh. node().head() + _jnod ;
        auto _kptr =
       _mesh. node().head() + _knod ;

    /*------------------------------------- dual geometry */
        real_type _0bal [_last + 1] ;
        pred_type::tri3_ball(_0bal,
            &_iptr->pval(+0),
            &_jptr->pval(+0),
            &_kptr->pval(+0), true) ;

        real_type _1bal [_last + 1] ;
        pred_type::edge_ball(_1bal,
            &_iptr->pval(+0),
            &_jptr->pval(+0), true) ;

        real_type _2bal [_last + 1] ;
        pred_type::edge_ball(_2bal,
            &_iptr->pval(+0),
            &_kptr->pval(+0), true) ;

        real_type _1mag = std::max(
        pred_type::tri3_mass(
            &_iptr->pval(+0),
            &_1bal [+0] ,
            &_0bal [+0]), (real_type)+0.0);

        real_type _2mag = std::max(
        pred_type::tri3_mass(
            &_iptr->pval(+0),
            &_0bal [+0] ,
            &_2bal [+0]), (real_type)+0.0);

    /*------------------------------------- cvt weighting */
        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_0bal[ +0]);

        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_1bal[ +0]);

        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_2bal[ +0]);

        real_type _lsqr = std::max(
            _0bal[_last], (real_type) +0.);

        EVALHFUN( _inod , _iptr )        // eval. at node
        EVALHFUN( _jnod , _jptr )
        EVALHFUN( _knod , _kptr )

        real_type _0rho =
            _hval[_inod] * (real_type)1./3. +
            _hval[_jnod] * (real_type)1./3. +
            _hval[_knod] * (real_type)1./3. ;

        real_type _1rho =
            _hval[_inod] * (real_type)1./2. +
            _hval[_jnod] * (real_type)1./2. ;

        real_type _2rho =
            _hval[_inod] * (real_type)1./2. +
            _hval[_knod] * (real_type)1./2. ;

        real_type _0val = _hfun.eval(    // eval. on dual
            &_0bal[0], _node->hidx()) ;
        real_type _1val = _hfun.eval(
            &_1bal[0], _node->hidx()) ;
        real_type _2val = _hfun.eval(
            &_2bal[0], _node->hidx()) ;

        _0rho = +.5 * (_0rho + _0val) ;
        _1rho = +.5 * (_1rho + _1val) ;
        _2rho = +.5 * (_2rho + _2val) ;

        _0rho = std::pow(
            _hval[_inod] / _0rho, +3) ;
        _1rho = std::pow(
            _hval[_inod] / _1rho, +3) ;
        _2rho = std::pow(
            _hval[_inod] / _2rho, +3) ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != 0; )
        {
            _move[_idim] +=
        // 1st sub-tria in voro. cell
                _1mag * (
            _0rho * _0bal[_idim] +
            _1rho * _1bal[_idim] +
            _iptr->  pval(_idim) ) +

        // 2nd sub-tria in voro. cell
                _2mag * (
            _0rho * _0bal[_idim] +
            _2rho * _2bal[_idim] +
            _iptr->  pval(_idim) ) ;
        }

        _wsum += _1mag * (_0rho + _1rho + 1.) ;
        _wsum += _2mag * (_0rho + _2rho + 1.) ;

        _ladj += _lsqr ;

    #   undef   LOOPDIMS
    #   undef   EVALHFUN
    }

    /*
    --------------------------------------------------------
     * ~CVT-QUAD: optimal delaunay tessellation update.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _cvt_quad_k (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type *_move ,
        real_type &_wsum ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0;

    #   define EVALHFUN( _NODE, _NPTR ) \
        if (_hval[_NODE] < (real_type)+0.)  \
            _hval[_NODE] = _hfun.eval ( \
           &_NPTR->pval(0), _NPTR->hidx())  ;

    #   define LOOPDIMS( _ITER )    \
            for ( auto _ITER =          \
        pred_type::geom_dims; _ITER-- != 0; )

    /*------------------------------------- cell indexing */
        auto _qptr =
       _mesh. quad().head() + _cell ;

        auto _inod = _qptr->node(0) ;
        auto _jnod = _qptr->node(1) ;
        auto _knod = _qptr->node(2) ;
        auto _lnod = _qptr->node(3) ;

        if ( _jnod ==_node->node(0) )
        {
           _inod   = _qptr->node(1) ;
           _jnod   = _qptr->node(2) ;
           _knod   = _qptr->node(3) ;
           _lnod   = _qptr->node(0) ;
        }
        else
        if ( _knod ==_node->node(0) )
        {
           _inod   = _qptr->node(2) ;
           _jnod   = _qptr->node(3) ;
           _knod   = _qptr->node(0) ;
           _lnod   = _qptr->node(1) ;
        }
        else
        if ( _lnod ==_node->node(0) )
        {
           _inod   = _qptr->node(3) ;
           _jnod   = _qptr->node(0) ;
           _knod   = _qptr->node(1) ;
           _lnod   = _qptr->node(2) ;
        }

        auto _iptr =
       _mesh. node().head() + _inod ;
        auto _jptr =
       _mesh. node().head() + _jnod ;
        auto _kptr =
       _mesh. node().head() + _knod ;
        auto _lptr =
       _mesh. node().head() + _lnod ;

    /*------------------------------------- dual geometry */
        real_type _0bal [_last + 1] ;
        pred_type::quad_ball(_0bal,
            &_iptr->pval(+0),
            &_jptr->pval(+0),
            &_kptr->pval(+0),
            &_lptr->pval(+0), true) ;

        real_type _1bal [_last + 1] ;
        pred_type::edge_ball(_1bal,
            &_iptr->pval(+0),
            &_jptr->pval(+0), true) ;

        real_type _2bal [_last + 1] ;
        pred_type::edge_ball(_2bal,
            &_iptr->pval(+0),
            &_lptr->pval(+0), true) ;

        real_type _1mag = std::max(
        pred_type::tri3_mass(
            &_iptr->pval(+0),
            &_1bal [+0] ,
            &_0bal [+0]), (real_type)+0.0);

        real_type _2mag = std::max(
        pred_type::tri3_mass(
            &_iptr->pval(+0),
            &_0bal [+0] ,
            &_2bal [+0]), (real_type)+0.0);

    /*------------------------------------- cvt weighting */
        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_0bal[ +0]);

        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_1bal[ +0]);

        pred_type::proj_node(
            _geom ,
           &_node->pval( +0), &_2bal[ +0]);

        real_type _lsqr = std::max(
            _0bal[_last], (real_type) +0.);

        EVALHFUN( _inod , _iptr )        // eval. at node
        EVALHFUN( _jnod , _jptr )
        EVALHFUN( _knod , _kptr )
        EVALHFUN( _lnod , _lptr )

        real_type _0rho =
            _hval[_inod] * (real_type)1./4. +
            _hval[_jnod] * (real_type)1./4. +
            _hval[_knod] * (real_type)1./4. +
            _hval[_lnod] * (real_type)1./4. ;

        real_type _1rho =
            _hval[_inod] * (real_type)1./2. +
            _hval[_jnod] * (real_type)1./2. ;

        real_type _2rho =
            _hval[_inod] * (real_type)1./2. +
            _hval[_lnod] * (real_type)1./2. ;

        real_type _0val = _hfun.eval(    // eval. on dual
            &_0bal[0], _node->hidx()) ;
        real_type _1val = _hfun.eval(
            &_1bal[0], _node->hidx()) ;
        real_type _2val = _hfun.eval(
            &_2bal[0], _node->hidx()) ;

        _0rho = +.5 * (_0rho + _0val) ;
        _1rho = +.5 * (_1rho + _1val) ;
        _2rho = +.5 * (_2rho + _2val) ;

        _0rho = std::pow(
            _hval[_inod] / _0rho, +3) ;
        _1rho = std::pow(
            _hval[_inod] / _1rho, +3) ;
        _2rho = std::pow(
            _hval[_inod] / _2rho, +3) ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != 0; )
        {
            _move[_idim] +=
        // 1st sub-tria in voro. cell
                _1mag * (
            _0rho * _0bal[_idim] +
            _1rho * _1bal[_idim] +
            _iptr->  pval(_idim)
                ) +

        // 2nd sub-tria in voro. cell
                _2mag * (
            _0rho * _0bal[_idim] +
            _2rho * _2bal[_idim] +
            _iptr->  pval(_idim)
                ) ;
        }

        _wsum += _1mag * (_0rho + _1rho + 1.) ;
        _wsum += _2mag * (_0rho + _2rho + 1.) ;

        _ladj += _lsqr ;

    #   undef   LOOPDIMS
    #   undef   EVALHFUN
    }

    /*
    --------------------------------------------------------
     * _CVT-MOVE: centroidal voronoi tessellation move.
    --------------------------------------------------------
     */

    // Shift vertex to weighted mean of adj. dual cell

    // xnew <=== INT|i rho(x) * x dA / INT|i rho(x) dA

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _cvt_move_2 (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        conn_list &_conn ,
        node_iter  _node ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0;

        real_type _move[_last + 0] = {
                (real_type) +0.0};

        _ladj = (real_type) +0.0 ;

        real_type _wsum =
       +std::numeric_limits<real_type>::epsilon();

        iptr_type _cnum = +0  ;

        for (auto _next  = _conn.head() ,
                  _tend  = _conn.tend() ;
                  _next != _tend ;
                ++_next, ++_cnum )
        {
             auto _cell  = _next->_cell ;
             auto _kind  = _next->_kind ;

            if (_kind == TRIA3_tag)
            {
                _cvt_tri3_k(_geom, _mesh, _hfun,
                    _hval , _cell, _node, _move,
                    _wsum , _ladj ) ;
            }
            else
            if (_kind == QUAD4_tag)
            {
                _cvt_quad_k(_geom, _mesh, _hfun,
                    _hval , _cell, _node, _move,
                    _wsum , _ladj ) ;
            }
        }

        if (_cnum > +0)
        {
            for (auto _idim =
            pred_type::geom_dims ; _idim-- != 0; )
            {
                _line[_idim] =
                    _move[_idim] / _wsum
                        - _node->pval(_idim) ;
            }

            _ladj = std::sqrt(_ladj / _cnum) ;
        }
    }

    /*
    --------------------------------------------------------
     * DQDX-TRI3: "local-ascent" node movement vector.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type halo_tri3_k (
        mesh_type &_mesh ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type *_bmin ,
        real_type *_bmax ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

    /*------------------------------------- cell indexing */
        auto _tptr =
       _mesh. tri3().head() + _cell ;

        auto _inod = _tptr->node(0) ;
        auto _jnod = _tptr->node(1) ;
        auto _knod = _tptr->node(2) ;

        auto _iptr =
       _mesh. node().head() + _inod ;
        auto _jptr =
       _mesh. node().head() + _jnod ;
        auto _kptr =
       _mesh. node().head() + _knod ;

    /*------------------------------------- cell centroid */
        real_type _pmid[_last + 0] = {
            (real_type) +0.0} ;
    
        for (auto _idim =
        pred_type::geom_dims; _idim-- != +0; )
        {
            _pmid[_idim] = (real_type) 1./3. * (
                _iptr->pval(_idim) +
                _jptr->pval(_idim) +
                _kptr->pval(_idim) ) ;
                
            _bmin[_idim] = std::min(
            _bmin[_idim], _pmid[_idim]) ;

            _bmax[_idim] = std::max(
            _bmax[_idim], _pmid[_idim]) ;
        }
        
    /*------------------------------------- cell "length" */
        real_type _lsqr =
            pred_type::length_sq (
                _pmid, &_node->pval(0)) ;

        _ladj +=  _lsqr ;
    }

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdx_tri3_k (
        mesh_type &_mesh ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type  _COST ,
        real_type *_bmin ,
        real_type *_bmax ,
        real_type *_SAVE ,
        real_type *_DQDX
        )
    {
        real_type static const _HINC =
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +0.75) ;

    /*------------------------------------- cell indexing */
        auto _ppos = &_node->pval(0);

        auto _tptr =
       _mesh. tri3().head() + _cell ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != +0; )
        {
    /*------------------ local iter. on finite-diff. step */
        real_type _binc = _bmax[_idim] -
                          _bmin[_idim] ;

        real_type _hdel = _HINC*_binc;
        real_type _dsup , _dsdn ;
        {
    /*------------------ centred finite-diff. for dQ / dx */
            _ppos[_idim] =
                _SAVE[_idim] + _hdel;

            real_type _scr1 =
                pred_type::tri3_cost (
               &_mesh. node(
                _tptr->node(0)).pval(0),
               &_mesh. node(
                _tptr->node(1)).pval(0),
               &_mesh. node(
                _tptr->node(2)).pval(0),
                typename
                pred_type::cell_kind ()
                    ) ;

            _dsup = (_scr1 - _COST) / _hdel ;

            _ppos[_idim] =
                _SAVE[_idim] - _hdel;

            real_type _scr0 =
                pred_type::tri3_cost (
               &_mesh. node(
                _tptr->node(0)).pval(0),
               &_mesh. node(
                _tptr->node(1)).pval(0),
               &_mesh. node(
                _tptr->node(2)).pval(0),
                typename
                pred_type::cell_kind ()
                    ) ;

            _dsdn = (_COST - _scr0) / _hdel ;
        }

        _ppos[_idim]  = _SAVE  [_idim] ;

    /*------------------ finalise gradient and accumulate */
      //if (_dsup * _dsdn >= (real_type)+0.)
            _DQDX [ _idim ] +=
        (real_type) +0.50 * (_dsup + _dsdn) ;

        }
    }

    /*
    --------------------------------------------------------
     * DQDX-QUAD: "local-ascent" node movement vector.
    --------------------------------------------------------
     */

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type halo_quad_k (
        mesh_type &_mesh ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type *_bmin ,
        real_type *_bmax ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

    /*------------------------------------- cell indexing */
        auto _qptr =
       _mesh. quad().head() + _cell ;

        auto _inod = _qptr->node(0) ;
        auto _jnod = _qptr->node(1) ;
        auto _knod = _qptr->node(2) ;
        auto _lnod = _qptr->node(3) ;

        auto _iptr =
       _mesh. node().head() + _inod ;
        auto _jptr =
       _mesh. node().head() + _jnod ;
        auto _kptr =
       _mesh. node().head() + _knod ;
        auto _lptr =
       _mesh. node().head() + _lnod ;

    /*------------------------------------- cell centroid */
        real_type _pmid[_last + 0] = {
            (real_type) +0.0} ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != +0; )
        {
            _pmid[_idim] = (real_type) 1./4. * (
                _iptr->pval(_idim) +
                _jptr->pval(_idim) +
                _kptr->pval(_idim) +
                _lptr->pval(_idim) ) ;
                
            _bmin[_idim] = std::min(
            _bmin[_idim], _pmid[_idim]) ;

            _bmax[_idim] = std::max(
            _bmax[_idim], _pmid[_idim]) ;
        }

    /*------------------------------------- cell "length" */
        real_type _lsqr =
            pred_type::length_sq (
                _pmid, &_node->pval(0)) ;

        _ladj +=  _lsqr ;
    }

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdx_quad_k (
        mesh_type &_mesh ,
        iptr_type  _cell ,
        node_iter  _node ,
        real_type  _COST ,
        real_type *_bmin ,
        real_type *_bmax ,
        real_type *_SAVE ,
        real_type *_DQDX
        )
    {
        real_type static const _HINC =
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +0.75) ;

    /*------------------------------------- cell indexing */
        auto _ppos = &_node->pval(0);

        auto _qptr =
       _mesh. quad().head() + _cell ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != +0; )
        {
    /*------------------ local iter. on finite-diff. step */
        real_type _binc = _bmax[_idim] -
                          _bmin[_idim] ;

        real_type _hdel = _HINC*_binc;
        real_type _dsup , _dsdn ;
        {
    /*------------------ centred finite-diff. for dQ / dx */
            _ppos[_idim] =
                _SAVE[_idim] + _hdel ;

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
                pred_type::cell_kind ()
                    ) ;

            _dsup = (_scr1 - _COST) / _hdel ;

            _ppos[_idim] =
                _SAVE[_idim] - _hdel ;

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
                pred_type::cell_kind ()
                    ) ;

            _dsdn = (_COST - _scr0) / _hdel ;
        }

        _ppos[_idim]  = _SAVE  [_idim] ;

    /*------------------ finalise gradient and accumulate */
      //if (_dsup * _dsdn >= (real_type)+0.)
            _DQDX [ _idim ] +=
        (real_type) +0.50 * (_dsup + _dsdn) ;

        }
    }

    /*
    --------------------------------------------------------
     * DQDX-MOVE: "local-ascent" node movement vector.
    --------------------------------------------------------
     */

    // Shift vertex via gradients of quality functions

    // xnew <=== xold + beta * BEST[ d./dx Q(x,w)|_i ]

    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdx_move_2 (
        mesh_type &_mesh ,
        conn_list &_conn ,
        node_iter  _node ,
        fp32_list &_cost ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        real_type static const _HINC =
            std::pow(std::numeric_limits
                <real_type>::epsilon(), .50) ;

        _ladj = (real_type) 0. ;

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
                    _ladj ) ;
            }
            else
            if (_kind == QUAD4_tag)
            {
                halo_quad_k(_mesh,
                    _cell , _node, _bmin, _bmax,
                    _ladj ) ;
            }
        }

        _ladj = std::sqrt (_ladj / _cnum) ;

    /*------------------ adj. gradients wrt. pos: dQ / dx */
        float    _qlim = _qmin +
            (float    ) +1.0E-002  ;

        real_type _DQDX[_last] = {
            (real_type) +0.0E+000} ;

        real_type _qbar = (real_type) 1.;
        real_type _qlow = (real_type) 0.;

        real_type _SAVE[_last] ;
        for (auto _idim  =
        pred_type::geom_dims; _idim-- != +0; )
        {
            _SAVE[_idim] = _node->pval(_idim);
        }

        iptr_type _lnum  = 0, _hnum  = 1, _nnum  = 0;
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
                dqdx_tri3_k(_mesh, _cell,
                    _node , _COST, _bmin, _bmax,
                    _SAVE , _DQDX) ;
                }
                else
                if (_kind == QUAD4_tag)
                {
                dqdx_quad_k(_mesh, _cell,
                    _node , _COST, _bmin, _bmax,
                    _SAVE , _DQDX) ;
                }

                _qlow += _COST ; _lnum += +1 ;
            }
            else
            {
                _qbar += _COST ; _hnum += +1 ;
            }
        }

        if (_lnum > +0)
        {
            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _DQDX[_idim] /= _lnum ;
            }

            _qlow /=  _lnum ;
            _qbar /=  _hnum ;
        }

    /*------------------ 1st ord. Taylor-series step est. */
        real_type _scal =
        std::sqrt (pred_type::length_sq(_DQDX));

        if (_scal*_ladj <= _HINC)
        {
            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _line[_idim] = (real_type)+0.;
            }
        }
        else
        {
            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _line[_idim] =
                    _DQDX[_idim]/ _scal ;
            }

            _scal = (_qbar - _qlow) /
            pred_type::innerprod(_line, _DQDX) ;

            _scal = std::min(
                _scal / (real_type)2. , _ladj) ;

            for (auto _idim  =
            pred_type::geom_dims; _idim-- != +0; )
            {
                _line[_idim] *= _scal ;
            }
        }

    }



