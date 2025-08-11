
    /*
    --------------------------------------------------------
     * ITER-ZIPS-2: optim. schemes to merge nodes.
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
     * Last updated: 04 Aug., 2025
     *
     * Copyright 2013-2025
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    __static_call
    __normal_call void_type _pop_sets (
        mesh_type &_mesh ,  // mesh object
        conn_list &_ladj    // space for adj. cells
        )
    {
        for (auto _cell  = _ladj.head() ;
                  _cell != _ladj.tend() ;
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
                _mesh.
                _pop_tri3(_cell->_cell) ;
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
                _mesh.
                _pop_quad(_cell->_cell) ;
            }
        }
    }

    template <
        typename  cell_type
             >
    __static_call
    __normal_call void_type _adj_sets (
        mesh_type &_mesh ,  // mesh object
        cell_type &_cdat ,  // cell to zip
        iptr_type  _nnod ,  // no. nodes in cell
        conn_list &_list ,  // space for adj. cells
        conn_list &_ladj    // space for adj. cells
        )
    {
        for (auto _cell  = _list.head();
                  _cell != _list.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
    /*--------------------------------- set a = { i \ t } */
             auto _tptr  =
            _mesh. tri3().head()+_cell->_cell;

            iptr_type _nnum, _mnum;
            iptr_type _M = _nnod, _same = +0 ;
            for(_nnum = 3; _nnum-- != 0; )
            for(_mnum =_M; _mnum-- != 0; )
            {
                if (_tptr->node(_nnum) ==
                    _cdat. node(_mnum) )
                {
                    _same += +1 ;
                }
            }

            if (_same > +1 ) continue;
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
    /*--------------------------------- set a = { i \ q } */
             auto _qptr  =
            _mesh. quad().head()+_cell->_cell;

            iptr_type _nnum, _mnum;
            iptr_type _M = _nnod, _same = +0 ;
            for(_nnum = 4; _nnum-- != 0; )
            for(_mnum =_M; _mnum-- != 0; )
            {
                if (_qptr->node(_nnum) ==
                    _cdat. node(_mnum) )
                {
                    _same += +1 ;
                }
            }

            if (_same > +1 ) continue;
            }

    /*----------------- if we got this far, add adj. cell */
            _ladj.push_tail( *_cell );
        }
    }

    template <
        typename  cell_type
             >
    __static_call
    __normal_call void_type _adj_redo (
        mesh_type &_mesh ,  // mesh object
        cell_type &_cdat ,  // cell to zip
        iptr_type  _inod ,  // old node in cell
        iptr_type  _inew ,  // new node to push
        conn_list &_ladj ,  // space for adj. cells
        conn_list &_lcav    // space for adj. cells
        )
    {
        typename mesh_type::tri3_type _tdat;
        typename mesh_type::quad_type _qdat;

        for (auto _cell  = _ladj.head();
                  _cell != _ladj.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
             auto _tptr  =
            _mesh. tri3().head()+_cell->_cell ;

            _tdat.node(0) = _tptr->node(0) ;
            _tdat.node(1) = _tptr->node(1) ;
            _tdat.node(2) = _tptr->node(2) ;

            _tdat.itag () = _tptr->itag () ;

            iptr_type _nnum ;
            for(_nnum = 3; _nnum-- != 0; )
            {
            if (_tdat. node(_nnum) ==
                _cdat. node(_inod) )
                _tdat. node(_nnum) = _inew ;
            }

             auto _tnum =
                _mesh.push_tri3(_tdat) ;

            _lcav.push_tail( typename
        conn_list::data_type(_tnum,TRIA3_tag));
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
             auto _qptr  =
            _mesh. quad().head()+_cell->_cell ;

            _qdat.node(0) = _qptr->node(0) ;
            _qdat.node(1) = _qptr->node(1) ;
            _qdat.node(2) = _qptr->node(2) ;
            _qdat.node(3) = _qptr->node(3) ;

            _qdat.itag () = _qptr->itag () ;

            iptr_type _nnum ;
            for(_nnum = 4; _nnum-- != 0; )
            {
            if (_qdat. node(_nnum) ==
                _cdat. node(_inod) )
                _qdat. node(_nnum) = _inew ;
            }

             auto _qnum =
                _mesh.push_quad(_qdat) ;

            _lcav.push_tail( typename
        conn_list::data_type(_qnum,QUAD4_tag));
            }
        }
    }

    __static_call
    __normal_call void_type _adj_ball (
        mesh_type &_mesh ,  // mesh object
        conn_list &_ladj ,  // space for adj. cells
        real_type *_ball ,  // accumulate cell ball
        real_type &_mass    // accumulate cell mass
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        for (auto _cell  = _ladj.head() ;
                  _cell != _ladj.tend() ;
                ++_cell  )
        {
            real_type _BALL [_last + 1] = {0};
            real_type _MASS = 0 ;
            if (_cell->_kind == TRIA3_tag)
            {
    /*--------------------------------- accum. ball; tri3 */
            auto _tptr =
           _mesh. tri3().head() + _cell->_cell ;

            auto _INOD = _tptr->node(0) ;
            auto _JNOD = _tptr->node(1) ;
            auto _KNOD = _tptr->node(2) ;

            auto _IPTR =
           _mesh. node().head() + _INOD ;
            auto _JPTR =
           _mesh. node().head() + _JNOD ;
            auto _KPTR =
           _mesh. node().head() + _KNOD ;

            _MASS = std::abs(
                  pred_type::tri3_mass(
                &_IPTR->pval(+0),
                &_JPTR->pval(+0),
                &_KPTR->pval(+0)) ) ;

            pred_type::tri3_ball(_BALL,
                &_IPTR->pval(+0),
                &_JPTR->pval(+0),
                &_KPTR->pval(+0), true) ;
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
    /*--------------------------------- accum. ball; quad */
            auto _qptr =
           _mesh. quad().head() + _cell->_cell ;

            auto _INOD = _qptr->node(0) ;
            auto _JNOD = _qptr->node(1) ;
            auto _KNOD = _qptr->node(2) ;
            auto _LNOD = _qptr->node(3) ;

            auto _IPTR =
           _mesh. node().head() + _INOD ;
            auto _JPTR =
           _mesh. node().head() + _JNOD ;
            auto _KPTR =
           _mesh. node().head() + _KNOD ;
            auto _LPTR =
           _mesh. node().head() + _LNOD ;

            _MASS = std::abs(
                  pred_type::quad_mass(
                &_IPTR->pval(+0),
                &_JPTR->pval(+0),
                &_KPTR->pval(+0),
                &_LPTR->pval(+0)) ) ;

            pred_type::quad_ball(_BALL,
                &_IPTR->pval(+0),
                &_JPTR->pval(+0),
                &_KPTR->pval(+0),
                &_LPTR->pval(+0), true) ;
            }

            for (auto _idim = _last;
                      _idim-- != +0; )
            {
            _ball[_idim] += _MASS * _BALL[_idim] ;
            }
            _mass +=  _MASS ;
        }
    }


    /*
    --------------------------------------------------------
     * _ZIP-EDGE: try edge merge to improve adj. cost-fun.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type _zip_edge (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        iter_opts &_opts ,  // user options
        iptr_list &_nmrk ,  // list of node flags
        iptr_type  _iout ,  // outer iteration marker
        iptr_type  _edge ,  // edge number to zip.
        kern_kind  _kern ,  // optim. kernel selector
        bool_type &_okay ,  // TRUE if zip. successful
        iptr_type &_nnew ,  // new node index added
        conn_list &_iset ,  // space for adj. cells
        conn_list &_jset ,  // space for adj. cells
        conn_list &_iadj ,  // space for adj. cells
        conn_list &_jadj ,  // space for adj. cells
        conn_list &_cset ,  // space for adj. cells
        fp32_list &_qsrc ,  // space for adj. costs
        fp32_list &_qdst ,  // space for adj. costs
        fp32_list &_qtmp ,  // space for adj. costs
        real_type  _QLIM ,  // cell quality threshold
        real_type  _QMOV ,  // cell quality tolerance
        iter_stat &_tcpu ,  // cpu timers/info
        real_type  _lmax    // zip. spacing threshold
            = (real_type) +8.91667E-01 ,
        real_type  _qinc    // zip. quality threshold
            = (real_type) +3.33333E-03
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        iptr_type static constexpr
            _DEG_TRIA3 = (iptr_type)+6 ;
    //  iptr_type static constexpr
    //      _DEG_QUAD4 = (iptr_type)+4 ;

    /*--------------------------------- get edge indexing */
         auto _eptr =
        _mesh. edge().head()+_edge  ;

        typename mesh_type::
            edge_type _edat(*_eptr) ;

        iptr_type _inod, _jnod ;
        _inod = _eptr->node(0) ;
        _jnod = _eptr->node(1) ;

         auto _iptr = _mesh.
         node().head() + _eptr->node(0);
         auto _jptr = _mesh.
         node().head() + _eptr->node(1);

        _okay = false ;

        _iset.set_count(
            0, containers::loose_alloc);
        _jset.set_count(
            0, containers::loose_alloc);
        _iadj.set_count(
            0, containers::loose_alloc);
        _jadj.set_count(
            0, containers::loose_alloc);

        _cset.set_count(
            0, containers::loose_alloc);

        _qsrc.set_count(
            0, containers::loose_alloc);
        _qdst.set_count(
            0, containers::loose_alloc);
        _qtmp.set_count(
            0, containers::loose_alloc);

    /*--------------------------------- exit if FEAT node */
        if (_iptr->feat() != null_feat)
            return ;
        if (_jptr->feat() != null_feat)
            return ;

    /*--------------------------------- get adjacent face */
        _mesh.connect_2(&_eptr->node(0),
            EDGE2_tag , _cset) ;

        if (_cset.count() != 2) return ;

        auto _atri =
       _mesh.tri3().head()+_cset[0]._cell ;
        auto _btri =
       _mesh.tri3().head()+_cset[1]._cell ;

        iptr_type _anod = -1;           // find cell apex
        if (_atri->node(0) != _inod &&
            _atri->node(0) != _jnod )
            _anod = _atri->node( 0 ) ;
        if (_atri->node(1) != _inod &&
            _atri->node(1) != _jnod )
            _anod = _atri->node( 1 ) ;
        if (_atri->node(2) != _inod &&
            _atri->node(2) != _jnod )
            _anod = _atri->node( 2 ) ;

        iptr_type _bnod = -1;           // find cell apex
        if (_btri->node(0) != _inod &&
            _btri->node(0) != _jnod )
            _bnod = _btri->node( 0 ) ;
        if (_btri->node(1) != _inod &&
            _btri->node(1) != _jnod )
            _bnod = _btri->node( 1 ) ;
        if (_btri->node(2) != _inod &&
            _btri->node(2) != _jnod )
            _bnod = _btri->node( 2 ) ;

    /*--------------------------------- calc. local topo. */
        auto _ideg =
       _mesh .node(_inod).ndeg(TRIA3_tag) ;
        auto _jdeg =
       _mesh .node(_jnod).ndeg(TRIA3_tag) ;
        auto _adeg =
       _mesh .node(_anod).ndeg(TRIA3_tag) ;
        auto _bdeg =
       _mesh .node(_bnod).ndeg(TRIA3_tag) ;

        if (_nmrk[_inod] < 0) _ideg = 6 ;
        if (_nmrk[_jnod] < 0) _jdeg = 6 ;
        if (_nmrk[_anod] < 0) _adeg = 6 ;
        if (_nmrk[_bnod] < 0) _bdeg = 6 ;

        auto _ierr =
            (iptr_type)(_DEG_TRIA3-_ideg) ;
        auto _jerr =
            (iptr_type)(_DEG_TRIA3-_jdeg) ;
        auto _aerr =
            (iptr_type)(_DEG_TRIA3-_adeg) ;
        auto _berr =
            (iptr_type)(_DEG_TRIA3-_bdeg) ;

        auto _ndeg = _ideg + _jdeg -  4 ;
        auto _ADEG = _adeg - 1 ;
        auto _BDEG = _bdeg - 1 ;

        auto _nerr =
            (iptr_type)(_DEG_TRIA3-_ndeg) ;
        auto _AERR =
            (iptr_type)(_DEG_TRIA3-_ADEG) ;
        auto _BERR =
            (iptr_type)(_DEG_TRIA3-_BDEG) ;

        auto _dnow = std::abs(_ierr) +  // current defect
                     std::abs(_jerr) +
                     std::abs(_aerr) +
                     std::abs(_berr) ;

        auto _dnew = std::abs(_nerr) +  // updated defect
                     std::abs(_AERR) +
                     std::abs(_BERR) ;

    //  bail-out early if the topological defect would be
    //  made worse if the zip is done

        if (_dnew > _dnow + 0) return;

    //  if more regular topo. is constructed via the edge
    //  zip, make it easier to do so!

        real_type _qtol = 
            std::min( +1./24., +1. - _QMOV) ;

        real_type _qerr = (real_type)
           -_qtol *std::max(0, _dnow-_dnew) ;

        real_type _lerr = (real_type)
           +_qtol *std::max(0, _dnow-_dnew) ;

        if (_dnew < _dnow - 0)
        {                            // no oscl. wrt. div
        _qerr/= std::pow ( _iout, +1./4.) ;
        _qinc = std::min ( _qinc,
                           _qinc + _qerr) ;

        _lerr/= std::pow ( _iout, +1./4.) ;
        _lmax = std::max ( _lmax,
                           _lmax + _lerr) ;
        }

        iptr_type static constexpr
            _dMIN  = (iptr_type)+2 ;

        if (_ierr >= _dMIN || _jerr >= _dMIN)
        {
        _qinc = (real_type) -0.250 ; // force a quad. zip
        _lmax = (real_type) +1.250 ;
        }

    /*--------------------------------- get edge h-sizing */
        real_type _pmid[_last + 1] ;
        for (auto _idim =
        pred_type::real_dims; _idim-- != 0; )
        {
            _pmid[_idim]  = (real_type)0. ;
            _pmid[_idim] +=
                _iptr->pval(_idim) ;
            _pmid[_idim] +=
                _jptr->pval(_idim) ;

            _pmid[_idim] /= (real_type)2. ;
        }

        pred_type::
        proj_node(_geom, _pmid, _pmid);

        real_type _isiz, _jsiz, _msiz ;
        _isiz = _hfun.eval(
            &_iptr->pval(0), _iptr->hidx()) ;
        _jsiz = _hfun.eval(
            &_jptr->pval(0), _jptr->hidx()) ;

         auto _hint = _iptr->hidx();

        _msiz = _hfun.eval(_pmid , _hint) ;

        real_type _lsqr =
            pred_type::length_sq (
           &_iptr->pval(0),&_jptr->pval(0)) ;

        real_type _hbar =
            _isiz * (real_type)(1./4.) +
            _jsiz * (real_type)(1./4.) +
            _msiz * (real_type)(2./4.) ;

    /*--------------------------------- exit if too large */
        if (_lsqr >= _hbar * _lmax *
                     _hbar * _lmax )
            return  ;

    /*--------------------------------- get disjoint sets */
        _mesh.connect_2(_eptr->node(+0),
            POINT_tag , _iset) ;

        _adj_sets(
            _mesh, _edat, 2, _iset, _iadj) ;

        _mesh.connect_2(_eptr->node(+1),
            POINT_tag , _jset) ;

        _adj_sets(
            _mesh, _edat, 2, _jset, _jadj) ;

    /*--------------------------------- get adjacent cost */
        real_type  _amin = loop_cost(
            _mesh,
            _iset, _qsrc, cell_kind()) ;

        real_type  _bmin = loop_cost(
            _mesh,
            _jadj, _qsrc, cell_kind()) ;

        real_type  _qmin =
            std::min(_amin, _bmin) ;

    /*--------------------------------- get adjacent ball */
        real_type _ball[_last + 1] =
            {(real_type) 0.} ;
        real_type _prev[_last + 1] =
            {(real_type) 0.} ;
        real_type _mass = (real_type)0.;

        _adj_ball(
            _mesh,_iset, _ball, _mass) ;
        _adj_ball(
            _mesh,_jset, _ball, _mass) ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != 0; )
        {
            _ball[_idim] /= _mass ;
        }

        pred_type::
        proj_node(_geom, _ball, _ball) ;

    /*--------------------------------- try to merge edge */
    #   define NULLHINT hfun_type::null_hint()
    #   define NULLFEAT null_feat

        typename mesh_type::node_type _ndat ;

        for (auto _idim =
        pred_type::real_dims; _idim-- != 0; )
        {
            _ndat.pval(_idim) = _ball[_idim];
        }

        _ndat.fdim() = (iptr_type) +2 ;
        _ndat.feat() = NULLFEAT ;
        _ndat.hidx() = NULLHINT ;
        _nnew = _mesh.push_node(_ndat);

         auto _nptr  =
        _mesh. node().head () + _nnew ;

        _hval.set_count( std::max(
        _nnew+1, (iptr_type)_hval.count())) ;

        _hval[_nnew] = (real_type)-1. ;

    /*--------------------------------- push a new cavity */
        _cset.set_count( +0 ) ;
        _adj_redo( _mesh,
        _edat, +0, _nnew, _iadj, _cset) ;

        _adj_redo( _mesh,
        _edat, +1, _nnew, _jadj, _cset) ;

    /*--------------------------------- optim. node coord */
        iptr_type static
        constexpr _INUM = (iptr_type) +4;

        iptr_type _move = -1;
        for (auto _iter = +0; _iter != _INUM;
                ++_iter )
        {
            _qtmp.set_count(
                +0, containers::loose_alloc);

            real_type  _minC =
                loop_cost( _mesh,
                    _cset, _qtmp, cell_kind()
                ) ;

            move_node( _geom, _mesh,
                _hfun, _hval,
                _opts, _nptr, _prev,
                _kern, _move, _cset,
                _qtmp, _qdst,
                _minC, _QLIM, _QMOV,
                _tcpu) ;

            if (_move <= +0 ) break;
        }

    /*--------------------------------- is cost improved? */
        _qdst.set_count(
            0, containers::loose_alloc) ;

        real_type  _QMIN =
        loop_cost( _mesh,
            _cset, _qdst,  cell_kind()) ;

        move_okay( _qdst, _qsrc, _move,
            +1.0 , _qinc) ;

        if((_okay  =  (_move > 0) && 
            _QMIN >=  (_qmin + _qinc)))
        {
    /*--------------------------------- delete old cavity */
            _pop_sets( _mesh , _iset );
            _pop_sets( _mesh , _jset );

            _mesh._pop_node ( &_inod );
            _mesh._pop_node ( &_jnod );
        }
        else
        {
    /*--------------------------------- delete new cavity */
            _pop_sets( _mesh , _cset );

            _mesh._pop_node ( &_nnew );
        }

    #   undef   NULLHINT
    #   undef   NULLFEAT
    }

    /*
    --------------------------------------------------------
     * _ZIP-TRI3: try tri3 merge to improve adj. cost-fun.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type _zip_tri3 (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        real_list &_hval ,  // cached node spacing array
        iter_opts &_opts ,  // user options
        iptr_list &_nmrk ,  // list of node flags
        iptr_type  _iout ,  // outer iteration marker
        iptr_type  _tri3 ,  // tri3 number to zip.
        kern_kind  _kern ,  // optim. kernel selector
        bool_type &_okay ,  // TRUE if zip. successful
        iptr_type &_nnew ,  // new node index added
        conn_list &_iset ,  // space for adj. cells
        conn_list &_jset ,  // space for adj. cells
        conn_list &_kset ,  // space for adj. cells
        conn_list &_iadj ,  // space for adj. cells
        conn_list &_jadj ,  // space for adj. cells
        conn_list &_kadj ,  // space for adj. cells
        conn_list &_cset ,  // space for adj. cells
        fp32_list &_qsrc ,  // space for adj. costs
        fp32_list &_qdst ,  // space for adj. costs
        fp32_list &_qtmp ,  // space for adj. costs
        real_type  _QLIM ,  // cell quality threshold
        real_type  _QMOV ,  // cell quality tolerance
        iter_stat &_tcpu ,  // cpu timers/info
        real_type  _lmax    // zip. spacing threshold
            = (real_type) +8.91667E-01 ,
        real_type  _qinc    // zip. quality threshold
            = (real_type) +3.33333E-03
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims + 0 ;

        iptr_type static constexpr
            _DEG_TRIA3 = (iptr_type)+6 ;
    //  iptr_type static constexpr
    //      _DEG_QUAD4 = (iptr_type)+4 ;

        __unreferenced (_iout) ;

    /*--------------------------------- get face indexing */
         auto _fptr =
        _mesh. tri3().head()+_tri3  ;

        typename mesh_type::
            tri3_type _fdat(*_fptr) ;

        iptr_type _inod, _jnod, _knod ;
        _inod = _fptr->node(0) ;
        _jnod = _fptr->node(1) ;
        _knod = _fptr->node(2) ;

         auto _iptr = _mesh.
         node().head() + _fptr->node(0);
         auto _jptr = _mesh.
         node().head() + _fptr->node(1);
         auto _kptr = _mesh.
         node().head() + _fptr->node(2);

        _okay = false ;

        _iset.set_count(
            0, containers::loose_alloc);
        _jset.set_count(
            0, containers::loose_alloc);
        _kset.set_count(
            0, containers::loose_alloc);
        _iadj.set_count(
            0, containers::loose_alloc);
        _jadj.set_count(
            0, containers::loose_alloc);
        _kadj.set_count(
            0, containers::loose_alloc);

        _cset.set_count(
            0, containers::loose_alloc);

        _qsrc.set_count(
            0, containers::loose_alloc);
        _qdst.set_count(
            0, containers::loose_alloc);
        _qtmp.set_count(
            0, containers::loose_alloc);

    /*--------------------------------- exit if FEAT node */
        if (_iptr->feat() != null_feat)
            return ;
        if (_jptr->feat() != null_feat)
            return ;
        if (_kptr->feat() != null_feat)
            return ;

    /*--------------------------------- calc. local topo. */
        auto _ideg =
       _mesh .node(_inod).ndeg(TRIA3_tag) ;
        auto _jdeg =
       _mesh .node(_jnod).ndeg(TRIA3_tag) ;
        auto _kdeg =
       _mesh .node(_knod).ndeg(TRIA3_tag) ;

        if (_nmrk[_inod] < 0) _ideg = 6 ;
        if (_nmrk[_jnod] < 0) _jdeg = 6 ;
        if (_nmrk[_knod] < 0) _kdeg = 6 ;

        auto _ierr =
            (iptr_type)(_DEG_TRIA3-_ideg) ;
        auto _jerr =
            (iptr_type)(_DEG_TRIA3-_jdeg) ;
        auto _kerr =
            (iptr_type)(_DEG_TRIA3-_kdeg) ;

        auto _ndeg = _ideg + _jdeg + _kdeg -  3*3 ;

        auto _nerr =
            (iptr_type)(_DEG_TRIA3-_ndeg) ;

        auto _dnew =(std::abs(_nerr) * 3) / 2 ;
        auto _dnow = std::abs(_ierr) +
                     std::abs(_jerr) +
                     std::abs(_kerr) ;

    //  bail-out early if the topological defect would be
    //  made worse if the zip is done

        if (_dnew > _dnow + 0) return ;

    //  if more regular topo. is constructed via the edge
    //  zip, make it easier to do so!

        real_type _qtol = 
            std::min( +1./36., +1. - _QMOV) ;

        real_type _qerr = (real_type)
           -_qtol *std::max(0, _dnow-_dnew) ;

        real_type _lerr = (real_type)
           +_qtol *std::max(0, _dnow-_dnew) ;

        if (_dnew < _dnow - 0)
        {                            // no oscl. wrt. div
        _qerr/= std::pow ( _iout, +1./4.) ;
        _qinc = std::min ( _qinc,
                           _qinc + _qerr) ;

        _lerr/= std::pow ( _iout, +1./4.) ;
        _lmax = std::max ( _lmax,
                           _lmax + _lerr) ;
        }
        
    /*--------------------------------- get face h-sizing */
        real_type _pbal[ _last + 1 ] ;
        pred_type::tri3_ball ( _pbal ,
           &_iptr->pval( +0 ) ,
           &_jptr->pval( +0 ) ,
           &_kptr->pval( +0 ) , true) ;

        pred_type::
        proj_node(_geom, _pbal, _pbal);

        real_type _isiz, _jsiz, _ksiz ;
        real_type _bsiz;
        _isiz = _hfun.eval(
            &_iptr->pval(0), _iptr->hidx()) ;
        _jsiz = _hfun.eval(
            &_jptr->pval(0), _jptr->hidx()) ;
        _ksiz = _hfun.eval(
            &_kptr->pval(0), _kptr->hidx()) ;

         auto _hint = _iptr->hidx() ;

        _bsiz = _hfun.eval(_pbal , _hint) ;

        real_type _lsqr =
            _pbal[_last] * (real_type)+3. ;
        real_type _hbar =
            _isiz * (real_type)(2./9.) +
            _jsiz * (real_type)(2./9.) +
            _ksiz * (real_type)(2./9.) +
            _bsiz * (real_type)(3./9.) ;

    /*--------------------------------- exit if too large */
        if (_lsqr >= _hbar * _lmax *
                     _hbar * _lmax )
            return  ;

    /*--------------------------------- get disjoint sets */
        _mesh.connect_2(_fptr->node(+0),
            POINT_tag , _iset) ;

        _adj_sets(
            _mesh, _fdat, 3, _iset, _iadj) ;

        _mesh.connect_2(_fptr->node(+1),
            POINT_tag , _jset) ;

        _adj_sets(
            _mesh, _fdat, 3, _jset, _jadj) ;

        _mesh.connect_2(_fptr->node(+2),
            POINT_tag , _kset) ;

        _adj_sets(
            _mesh, _fdat, 3, _kset, _kadj) ;

    /*--------------------------------- get adjacent cost */
        real_type  _amin = loop_cost(
            _mesh,
            _iset, _qsrc, cell_kind()) ;

        real_type  _bmin = loop_cost(
            _mesh,
            _jadj, _qsrc, cell_kind()) ;

        real_type  _cmin = loop_cost(
            _mesh,
            _kadj, _qsrc, cell_kind()) ;

        real_type  _qmin =
        std::min( {_amin, _bmin, _cmin} ) ;

    /*--------------------------------- get adjacent ball */
        real_type _ball[_last + 1] =
            {(real_type) 0.} ;
        real_type _prev[_last + 1] =
            {(real_type) 0.} ;
        real_type _mass = (real_type)0.;

        _adj_ball(
            _mesh,_iset, _ball, _mass) ;
        _adj_ball(
            _mesh,_jset, _ball, _mass) ;
        _adj_ball(
            _mesh,_kset, _ball, _mass) ;

        for (auto _idim =
        pred_type::geom_dims; _idim-- != 0; )
        {
            _ball[_idim] /= _mass ;
        }

        pred_type::
        proj_node(_geom, _ball, _ball) ;

    /*--------------------------------- try to merge face */
    #   define NULLHINT hfun_type::null_hint()
    #   define NULLFEAT null_feat

        typename mesh_type::node_type _ndat ;

        typename mesh_type::tri3_type _tdat ;
        typename mesh_type::quad_type _qdat ;

        for (auto _idim =
        pred_type::real_dims; _idim-- != 0; )
        {
            _ndat.pval(_idim) = _ball[_idim];
        }

        _ndat.fdim() = (iptr_type) +2 ;
        _ndat.feat() = NULLFEAT ;
        _ndat.hidx() = NULLHINT ;
        _nnew = _mesh.push_node(_ndat);

         auto _nptr  =
        _mesh. node().head () + _nnew ;

        _hval.set_count( std::max(
        _nnew+1, (iptr_type)_hval.count())) ;

        _hval[_nnew] = (real_type)-1. ;

    /*--------------------------------- push a new cavity */
        _adj_redo( _mesh,
        _fdat, +0, _nnew, _iadj, _cset) ;

        _adj_redo( _mesh,
        _fdat, +1, _nnew, _jadj, _cset) ;

        _adj_redo( _mesh,
        _fdat, +2, _nnew, _kadj, _cset) ;

    /*--------------------------------- optim. node coord */
        iptr_type static
        constexpr _INUM = (iptr_type) +4;

        iptr_type _move = -1;
        for (auto _iter = +0; _iter != _INUM;
                ++_iter )
        {
            _qtmp.set_count(
                +0, containers::loose_alloc);

            real_type  _minC =
                loop_cost( _mesh,
                    _cset, _qtmp, cell_kind()
                ) ;

            move_node( _geom, _mesh,
                _hfun, _hval,
                _opts, _nptr, _prev,
                _kern, _move, _cset,
                _qtmp, _qdst,
                _minC, _QLIM, _QMOV,
                _tcpu) ;

            if (_move <= +0 ) break;
        }

    /*--------------------------------- is cost improved? */
        _qdst.set_count(
            0, containers::loose_alloc) ;

        real_type  _QMIN =
        loop_cost( _mesh,
            _cset, _qdst,  cell_kind()) ;

        move_okay( _qdst, _qsrc, _move,
            +1.0 , _qinc) ;

        if((_okay  =  (_move > 0) && 
            _QMIN >=  (_qmin + _qinc)))
        {
    /*--------------------------------- delete old cavity */
            _pop_sets( _mesh , _iset );
            _pop_sets( _mesh , _jset );
            _pop_sets( _mesh , _kset );

            _mesh._pop_node ( &_inod );
            _mesh._pop_node ( &_jnod );
            _mesh._pop_node ( &_knod );
        }
        else
        {
    /*--------------------------------- delete new cavity */
            _pop_sets( _mesh , _cset );

            _mesh._pop_node ( &_nnew );
        }

    #   undef   NULLHINT
    #   undef   NULLFEAT
    }



