
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
     * Last updated: 01 Sept., 2021
     *
     * Copyright 2013-2021
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
        kern_kind  _kern ,
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
        iter_stat &_tcpu ,
        real_type  _lmin
            = (real_type) +1.00E+00 ,
        real_type  _qinc
            = (real_type) +1.00E-02
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
        _mesh. edge().head()+_enum ;

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

        real_type _isiz, _jsiz, _bsiz ;
        _isiz = _hfun.eval(
            &_iptr->pval(0), _iptr->hidx()) ;
        _jsiz = _hfun.eval(
            &_jptr->pval(0), _jptr->hidx()) ;

         auto _hint = _iptr->hidx();

        _bsiz = _hfun.eval(_ball , _hint) ;

        real_type _lsqr =
            pred_type::length_sq (
           &_iptr->pval(0),&_jptr->pval(0)) ;

        real_type _hbar =
            _isiz * (real_type)(1./4.) +
            _jsiz * (real_type)(1./4.) +
            _bsiz * (real_type)(2./4.) ;

    /*--------------------------------- exit if too small */
        if (_lsqr <= _hbar * _lmin *
                     _hbar * _lmin )
            return  ;

        auto _ideg = _mesh.
            node(_enod[0]).ndeg(TRIA3_tag ) ;
        auto _ierr =
            (iptr_type)(_ideg-_DEG_TRIA3) ;

        auto _jdeg = _mesh.
            node(_enod[1]).ndeg(TRIA3_tag ) ;
        auto _jerr =
            (iptr_type)(_jdeg-_DEG_TRIA3) ;

    /*--------------------------------- split if too long */
        real_type static constexpr
            _lBIG  = (real_type)+1.375 ;

        iptr_type static constexpr
            _dMAX  = (iptr_type)+3 ;

        real_type _lmax =
            _lBIG *std::pow( _iout ,  +1./8.) ;

        if (_lsqr >= _hbar * _lmax * // no oscl. wrt. zip
                     _hbar * _lmax )
        {
    //  if more regular geom. is constructed via the edge
    //  div, make it easier to do so!

        _qinc = (real_type)(+0.50 *
            (_hbar / std::sqrt(_lsqr) - 1.0)) ;

        }
        else
        if (_ierr >= _dMAX || _jerr >= _dMAX)
        {
    //  try to reduce large vertex valence via refinement

        _qinc = (real_type) -0.25 ;
        _qinc/=  std::pow ( _iout ,   +1./8.) ;

        }
        else
        if (_ierr > 0 && _jerr > 0 &&
                    (_ierr + _jerr) >= _dMAX)
        {
    //  try to reduce large vertex valence via refinement

        _qinc = (real_type) -.125 ;
        _qinc/=  std::pow ( _iout ,   +1./8.) ;

        }
        else
        {
    //  don't refine otherwise: bail-out early!
            return  ;
        }

    /*--------------------------------- get adjacent cost */
        real_type  _qmin = loop_cost(
            _mesh, _cset, _qsrc, cell_kind()) ;

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
            constexpr _INUM = (iptr_type) +2;

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
                _minC, _QLIM, _tcpu)  ;

            if (_move <= +0 ) break;
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



