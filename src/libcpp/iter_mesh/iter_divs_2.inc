
    /*
    --------------------------------------------------------
     * ITER-DIVS-2: optim. schemes to split edges.
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
     * Last updated: 27 Dec., 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    /*
    --------------------------------------------------------
     * _DIV-EDGE: try edge split to improve adj. cost-fun.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type _div_edge (
        geom_type &_geom ,
        mesh_type &_mesh ,
        hfun_type &_hfun ,
        real_list &_hval ,
        iter_opts &_opts ,
        iptr_type  _iout ,
        iptr_type  _enum ,
        char_type  _kern ,
        bool_type &_okay ,
        iptr_type &_inew ,
        conn_list &_cset ,
        conn_list &_cnew ,
        conn_list &_iset ,
        conn_list &_jset ,
        real_list &_qsrc ,
        real_list &_qdst ,
        real_list &_qtmp ,
        real_type  _QLIM ,
        real_type  _lmin
            = (real_type) +1.00E+00 ,
        real_type  _qinc
            = (real_type) +1.00E-04
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
        _mesh. edge().head()+_enum  ;

        iptr_type _enod[2] ;
        _enod[ 0] = _eptr->node(0) ;
        _enod[ 1] = _eptr->node(1) ;

         auto _iptr = _mesh.
         node().head() + _eptr->node(0);
         auto _jptr = _mesh.
         node().head() + _eptr->node(1);

        _cset.set_count(
            0, containers::loose_alloc);
        _cnew.set_count(
            0, containers::loose_alloc);
        _iset.set_count(
            0, containers::loose_alloc);
        _jset.set_count(
            0, containers::loose_alloc);

        _qsrc.set_count(
            0, containers::loose_alloc);
        _qdst.set_count(
            0, containers::loose_alloc);
        _qtmp.set_count(
            0, containers::loose_alloc);

        _mesh.connect_2(
            _enum, EDGE2_tag, _cset) ;

        _okay = false ;

        if (_cset.count()<=1) return ;

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

    /*--------------------------------- exit if too small */
        if (_lsqr <= _hbar * _lmin *
                     _hbar * _lmin )
            return  ;

    /*--------------------------------- split if too long */
        real_type static constexpr
        _lBIG  = (real_type)+1.375 ;

        if (_lsqr >= _hbar * _lBIG *
                     _hbar * _lBIG )
        {
    //  if more regular geom. is constructed via the edge
    //  div, make it easier to do so!

        _qinc  = 0.5 * (_hbar / std::sqrt(_lsqr) - 1.0) ;

        _qinc /= std::abs( _iout ) ; // no oscl. wrt. zip

        }
        else
        {
    /*--------------------------------- get adjacent face */
        _mesh.connect_2(_eptr->node(+0),
            POINT_tag , _iset) ;
        _mesh.connect_2(_eptr->node(+1),
            POINT_tag , _jset) ;

    /*--------------------------------- calc. local topo. */
        auto _ideg = _iset.count() ;
        auto _ierr =
            (iptr_type)(_ideg-_DEG_TRIA3) ;

        auto _jdeg = _jset.count() ;
        auto _jerr =
            (iptr_type)(_jdeg-_DEG_TRIA3) ;

        if (std::max(_ierr, _jerr) < 2)
            return ;

        }   // if (_lBIG)

    /*--------------------------------- get adjacent cost */
        real_type  _qmin = loop_cost(
            _mesh, _cset, _qsrc, cell_kind()) ;

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

    /*--------------------------------- try to split edge */
    #   define NULLFEAT null_feat
    #   define NULLHINT hfun_type::null_hint()

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
        _inew = _mesh.push_node(_ndat);

         auto _nptr  =
        _mesh. node().head () + _inew ;

        _hval.set_count( std::max(
        _inew + 1, (iptr_type)_hval.count())) ;

        _hval[_inew] = (real_type)-1. ;

    /*--------------------------------- push a new cavity */
        for (auto _cell  = _cset.head();
                  _cell != _cset.tend();
                ++_cell  )
        {
            if (_cell->_kind == TRIA3_tag)
            {
    /*--------------------------------- sub-div adj. tria */
             auto _tptr  =
            _mesh. tri3().head()+_cell->_cell ;

            iptr_type _ENOD[ 3], _ENUM, _TNUM ;
            for(_ENUM = 3; _ENUM-- != 0; )
            {
                mesh_type::tri3_type::
                face_node(_ENOD, _ENUM, 2, 1) ;
                _ENOD[0] =
                _tptr->node(_ENOD[0]);
                _ENOD[1] =
                _tptr->node(_ENOD[1]);
                _ENOD[2] =
                _tptr->node(_ENOD[2]);

                if (_ENOD[2] != _enod[0])
                if (_ENOD[2] != _enod[1])
                    break ;
            }

    /*--------------------------------- div into 2 tria's */
            _tdat.node(0) = _inew;
            _tdat.node(1) = _ENOD[2] ;
            _tdat.node(2) = _ENOD[0] ;
            _tdat.itag () = _tptr->itag () ;

            _TNUM = _mesh.push_tri3(_tdat) ;

            _cnew.push_tail(typename
        conn_list::data_type(_TNUM,TRIA3_tag));

            _tdat.node(0) = _inew;
            _tdat.node(1) = _ENOD[1] ;
            _tdat.node(2) = _ENOD[2] ;
            _tdat.itag () = _tptr->itag () ;

            _TNUM = _mesh.push_tri3(_tdat) ;

            _cnew.push_tail(typename
        conn_list::data_type(_TNUM,TRIA3_tag));
            }
            else
            if (_cell->_kind == QUAD4_tag)
            {
    /*--------------------------------- sub-div adj. quad */
             auto _qptr  =
            _mesh. quad().head()+_cell->_cell ;

            iptr_type _ENOD[ 4], _ENUM, _TNUM ;
            for(_ENUM = 4; _ENUM-- != 0; )
            {
                mesh_type::quad_type::
                face_node(_ENOD, _ENUM, 2, 1) ;
                _ENOD[0] =
                _qptr->node(_ENOD[0]);
                _ENOD[1] =
                _qptr->node(_ENOD[1]);
                _ENOD[2] =
                _qptr->node(_ENOD[2]);
                _ENOD[3] =
                _qptr->node(_ENOD[3]);

                if (_ENOD[2] != _enod[0])
                if (_ENOD[2] != _enod[1])
                if (_ENOD[3] != _enod[0])
                if (_ENOD[3] != _enod[1])
                    break ;
            }

    /*--------------------------------- div into 3 tria's */
            _tdat.node(0) = _inew;
            _tdat.node(1) = _ENOD[3] ;
            _tdat.node(2) = _ENOD[0] ;
            _tdat.itag () = _qptr->itag () ;

            _TNUM = _mesh.push_tri3(_tdat) ;

            _cnew.push_tail(typename
        conn_list::data_type(_TNUM,TRIA3_tag));

            _tdat.node(0) = _inew;
            _tdat.node(1) = _ENOD[2] ;
            _tdat.node(2) = _ENOD[3] ;
            _tdat.itag () = _qptr->itag () ;

            _TNUM = _mesh.push_tri3(_tdat) ;

            _cnew.push_tail(typename
        conn_list::data_type(_TNUM,TRIA3_tag));

            _tdat.node(0) = _inew;
            _tdat.node(1) = _ENOD[1] ;
            _tdat.node(2) = _ENOD[2] ;
            _tdat.itag () = _qptr->itag () ;

            _TNUM = _mesh.push_tri3(_tdat) ;

            _cnew.push_tail(typename
        conn_list::data_type(_TNUM,TRIA3_tag));
            }
        }

    /*--------------------------------- optim. node coord */
        iptr_type static
            constexpr _INUM = (iptr_type) +4;

        iptr_type _move = -1;
        for (auto _iloc = +0; _iloc != _INUM;
                ++_iloc )
        {
            _qtmp.set_count(
                +0, containers::loose_alloc);

            real_type  _minC =
                loop_cost( _mesh,
                    _cnew, _qtmp, cell_kind()
                ) ;

            move_node( _geom, _mesh,
                _hfun, _hval,
                _opts, _nptr, _kern,
                _move, _cnew,
                _qtmp, _qdst,
                _minC, _QLIM ) ;

            if (_move <= +0) break ;
        }

    /*--------------------------------- is cost improved? */
        _qdst.set_count(
            0, containers::loose_alloc) ;

        real_type  _QMIN =
        loop_cost( _mesh,
            _cnew, _qdst,  cell_kind()) ;

        move_okay( _qdst, _qsrc, _move,
        std::sqrt( _QLIM) ,
                   _qinc) ;

        if((_okay = _move > 0 &&
                    _QMIN >= _qmin+_qinc))
        {
    /*--------------------------------- delete old cavity */
        for (auto _cell  = _cset.head() ;
                  _cell != _cset.tend() ;
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
        else
        {
    /*--------------------------------- delete new cavity */
        for (auto _cell  = _cnew.head() ;
                  _cell != _cnew.tend() ;
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

            _mesh._pop_node ( &_inew );
        }

    #   undef   NULLHINT
    #   undef   NULLFEAT
    }



