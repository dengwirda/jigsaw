
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 30 April, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    /*
    --------------------------------------------------------
     * _ZIP-EDGE: try edge merge to improve adj. cost-fun.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type _zip_edge (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iter_opts &_opts ,
        iptr_type  _edge ,
        char_type  _kern ,
        bool_type &_okay ,
        iptr_type &_nnew ,
        conn_list &_iset ,
        conn_list &_jset ,
        conn_list &_aset ,
        conn_list &_bset ,
        conn_list &_cset ,
        real_list &_qsrc ,
        real_list &_qdst ,
        real_list &_qtmp ,
        real_type  _QLIM ,
        real_type  _lmax
            = (real_type) +7.50E-01 ,
        real_type  _qinc
            = (real_type) +1.00E-02
        )
    {
        iptr_type static constexpr
            _last = pred_type::geom_dims+0 ;

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

        _iset.set_count(+0) ;
        _jset.set_count(+0) ;
        _aset.set_count(+0) ;
        _bset.set_count(+0) ;
        _cset.set_count(+0) ;

        _qsrc.set_count(+0) ;
        _qdst.set_count(+0) ;
        _qtmp.set_count(+0) ;

    /*--------------------------------- exit if FEAT node */
        if (_iptr->feat() != null_feat)
            return ;
        if (_jptr->feat() != null_feat)
            return ;

    /*--------------------------------- get adjacent face */
        _mesh.connect_2(_eptr->node(+0),
             POINT_tag, _iset) ;
        _mesh.connect_2(_eptr->node(+1),
             POINT_tag, _jset) ;

    /*--------------------------------- calc. local topo. */
        auto _ndeg = _iset.count()
                   + _jset.count() - 4 ;

        auto _ideg = _iset.count() ;
        auto _ierr =
            (iptr_type)(_DEG_TRIA3-_ideg) ;

        auto _jdeg = _jset.count() ;
        auto _jerr =
            (iptr_type)(_DEG_TRIA3-_jdeg) ;

        auto _iabs = std::abs(_ierr) ;
        auto _jabs = std::abs(_jerr) ;

        auto _nerr = std::abs(
            (iptr_type)(_DEG_TRIA3-_ndeg));

    //  bail-out early if the topological defect would be
    //  made worse if the zip is done

        if (_nerr> std::max(_ierr, _jerr))
            return  ;

    //  if more regular topo. is constructed via the edge
    //  zip, make it easier to do so!

        real_type _qerr =
            (real_type)-1./9. * _ierr +
            (real_type)-1./9. * _jerr ;

        real_type _lerr =
            (real_type)+6./5. * _lmax ;

        if (2 * _nerr < _iabs + _jabs)
        {
        _qinc = std::min(_qinc, _qerr);
        _lmax = std::max(_lmax, _lerr);
        }

    /*--------------------------------- get edge h-sizing */
        real_type _isiz, _jsiz ;
        _isiz = _hfun.eval(
            &_iptr->pval(0), _iptr->hidx()) ;
        _jsiz = _hfun.eval(
            &_jptr->pval(0), _jptr->hidx()) ;

        real_type _lsqr =
            pred_type::length_sq (
           &_iptr->pval(0),&_jptr->pval(0)) ;

        real_type _hbar =
            std::min(_isiz , _jsiz);

    /*--------------------------------- exit if too large */
        if (_lsqr >= _hbar * _lmax *
                     _hbar * _lmax )
            return  ;

    /*--------------------------------- get disjoint sets */
        for (auto _cell  = _iset.head();
                  _cell != _iset.tend();
                ++_cell  )
        {
    /*--------------------------------- set a = { i \ e } */
            if (_cell->_kind == TRIA3_tag)
            {
             auto _tptr  =
            _mesh. tri3().head()+_cell->_cell ;

            iptr_type _nnum, _same = +0;
            for(_nnum = 3; _nnum-- != 0; )
            {
                if (_tptr->node(_nnum) ==
                    _eptr->node(  0) )
                    _same += +1 ;

                if (_tptr->node(_nnum) ==
                    _eptr->node(  1) )
                    _same += +1 ;
            }

            if (_same >= +2) continue ;
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
             auto _qptr  =
            _mesh. quad().head()+_cell->_cell ;

            iptr_type _nnum, _same = +0;
            for(_nnum = 4; _nnum-- != 0; )
            {
                if (_qptr->node(_nnum) ==
                    _eptr->node(  0) )
                    _same += +1 ;

                if (_qptr->node(_nnum) ==
                    _eptr->node(  1) )
                    _same += +1 ;
            }

            if (_same >= +2) continue ;
            }

            _aset.push_tail( *_cell ) ;
        }

        for (auto _cell  = _jset.head();
                  _cell != _jset.tend();
                ++_cell  )
        {
    /*--------------------------------- set b = { j \ e } */
            if (_cell->_kind == TRIA3_tag)
            {
             auto _tptr  =
            _mesh. tri3().head()+_cell->_cell ;

            iptr_type _nnum, _same = +0;
            for(_nnum = 3; _nnum-- != 0; )
            {
                if (_tptr->node(_nnum) ==
                    _eptr->node(  0) )
                    _same += +1 ;

                if (_tptr->node(_nnum) ==
                    _eptr->node(  1) )
                    _same += +1 ;
            }

            if (_same >= +2) continue ;
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
             auto _qptr  =
            _mesh. quad().head()+_cell->_cell ;

            iptr_type _nnum, _same = +0;
            for(_nnum = 4; _nnum-- != 0; )
            {
                if (_qptr->node(_nnum) ==
                    _eptr->node(  0) )
                    _same += +1 ;

                if (_qptr->node(_nnum) ==
                    _eptr->node(  1) )
                    _same += +1 ;
            }

            if (_same >= +2) continue ;
            }

            _bset.push_tail( *_cell ) ;
        }

    /*--------------------------------- get adjacent cost */
        real_type  _amin = loop_cost(
            _mesh, _iset, _qsrc, cell_kind()) ;

        real_type  _bmin = loop_cost(
            _mesh, _bset, _qsrc, cell_kind()) ;

        real_type  _qmin =
            std::min(_amin, _bmin) ;

    /*--------------------------------- get adjacent ball */
        real_type _ball[_last + 1] ;
        for (auto _idim =
        pred_type::real_dims; _idim-- != 0; )
        {
            _ball[_idim]  = (real_type) +0. ;
            _ball[_idim] +=
                _iptr->pval(_idim) ;
            _ball[_idim] +=
                _jptr->pval(_idim) ;

            _ball[_idim] /= (real_type) +2. ;
        }

        pred_type::
        proj_node(_geom, _ball, _ball);

    /*--------------------------------- try to merge edge */
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
        _nnew + 1, (iptr_type)_hval.count())) ;

        _hval[_nnew] = (real_type)-1. ;

    /*--------------------------------- push a new cavity */
        for (auto _cell  = _aset.head();
                  _cell != _aset.tend();
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
                _edat. node(0) )
                _tdat. node(_nnum) = _nnew ;
            }

             auto _tnum =
                _mesh.push_tri3(_tdat) ;

            _cset.push_tail( typename
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
                _edat. node(0) )
                _qdat. node(_nnum) = _nnew ;
            }

             auto _qnum =
                _mesh.push_quad(_qdat) ;

            _cset.push_tail( typename
        conn_list::data_type(_qnum,QUAD4_tag));
            }
        }

        for (auto _cell  = _bset.head();
                  _cell != _bset.tend();
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
                _edat. node(1) )
                _tdat. node(_nnum) = _nnew ;
            }

             auto _tnum =
                _mesh.push_tri3(_tdat) ;

            _cset.push_tail( typename
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
                _edat. node(1) )
                _qdat. node(_nnum) = _nnew ;
            }

             auto _qnum =
                _mesh.push_quad(_qdat) ;

            _cset.push_tail( typename
        conn_list::data_type(_qnum,QUAD4_tag));
            }
        }

    /*--------------------------------- optim. node coord */
        iptr_type static
            constexpr _INUM = (iptr_type) +8;

        iptr_type _move = -1;
        for (auto _iter = +0; _iter != _INUM;
                ++_iter )
        {
            _qtmp.set_count(0) ;

            real_type  _minC =
                loop_cost( _mesh,
                    _cset, _qtmp, cell_kind()
                ) ;

            move_node( _geom, _mesh,
                _hfun, _hval,
                _opts, _nptr, _kern,
                _move, _cset,
                _qtmp, _qdst,
                _minC, _QLIM ) ;

            if (_move <= +0) break ;
        }

    /*--------------------------------- is cost improved? */
        _qdst.set_count(0);

        real_type  _QMIN =
        loop_cost( _mesh,
            _cset, _qdst,  cell_kind ()
            ) ;

        move_okay( _qdst, _qsrc, _move,
            +1.0 , _qinc) ;

        if((_okay = _move > 0 &&
                    _QMIN >= _qmin+_qinc))
        {
    /*--------------------------------- delete old cavity */
        for (auto _cell  = _aset.head();
                  _cell != _aset.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
                _mesh.
                _pop_tri3(_cell->_cell);
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
                _mesh.
                _pop_quad(_cell->_cell);
            }
        }
        for (auto _cell  = _jset.head();
                  _cell != _jset.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
                _mesh.
                _pop_tri3(_cell->_cell);
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
                _mesh.
                _pop_quad(_cell->_cell);
            }
        }

            _mesh._pop_node ( &_inod ) ;
            _mesh._pop_node ( &_jnod ) ;
        }
        else
        {
    /*--------------------------------- delete new cavity */
        for (auto _cell  = _cset.head();
                  _cell != _cset.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
                _mesh.
                _pop_tri3(_cell->_cell);
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
                _mesh.
                _pop_quad(_cell->_cell);
            }
        }

            _mesh._pop_node ( &_nnew ) ;
        }

    #   undef   NULLHINT
    #   undef   NULLFEAT
    }



