
    /*
    --------------------------------------------------------
     * ITER-FLIP-2: optim. schemes to "flip" topology.
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


    __static_call
    __normal_call bool_type need_flip_delaunay (
        geom_type &_geom,
        mesh_type &_mesh,
      __const_ptr(iptr_type) _inod,
      __const_ptr(iptr_type) _jnod
        )
    {
    /*----------------------------- assess flip criterion */
        iptr_type static constexpr
            _last=pred_type::geom_dims+0 ;

        real_type static const _LTOL =
        std::pow(std::numeric_limits
            <real_type>::epsilon(), +.8) ;

        iptr_type  _iloc[3] ;
        _iloc[0] = _inod[0] ;
        _iloc[1] = _inod[1] ;
        _iloc[2] = _inod[2] ;

        algorithms::isort(
            &_iloc[0], &_iloc[3] ,
                std::less<iptr_type> ()) ;

        iptr_type  _jloc[3] ;
        _jloc[0] = _jnod[0] ;
        _jloc[1] = _jnod[1] ;
        _jloc[2] = _jnod[2] ;

        algorithms::isort(
            &_jloc[0], &_jloc[3] ,
                std::less<iptr_type> ()) ;

        real_type _ibal[_last + 1] ;
        pred_type::tri3_ball(_ibal ,
           &_mesh.node(_iloc[0]).pval(0) ,
           &_mesh.node(_iloc[1]).pval(0) ,
           &_mesh.node(_iloc[2]).pval(0)
            ) ;

        real_type _jbal[_last + 1] ;
        pred_type::tri3_ball(_jbal ,
           &_mesh.node(_jloc[0]).pval(0) ,
           &_mesh.node(_jloc[1]).pval(0) ,
           &_mesh.node(_jloc[2]).pval(0)
            ) ;

        real_type _null[_last] = {
            (real_type) +0.0 } ;

        pred_type::proj_node (
            _geom, _null, _ibal) ;

        pred_type::proj_node (
            _geom, _null, _jbal) ;

    /*----------------------------- measure dual distance */
        real_type  _ilen =
        pred_type::length_sq(_jbal,
           &_mesh.node(_inod[2]).pval(0));

        real_type  _jlen =
        pred_type::length_sq(_ibal,
           &_mesh.node(_jnod[2]).pval(0));

        real_type  _btol =
       (real_type)+.5 * _LTOL * (_ilen + _jlen) ;

        real_type  _idel =
            _ilen - _jbal [_last] ;
        real_type  _jdel =
            _jlen - _ibal [_last] ;

        if (_idel >= -_btol || _jdel >= -_btol)
        {
            return false ;
        }
        else
        {
            return  true ;
        }
    }

    __static_call
    __normal_call bool_type need_flip_laguerre (
        geom_type &_geom,
        mesh_type &_mesh,
      __const_ptr(iptr_type) _inod,
      __const_ptr(iptr_type) _jnod
        )
    {
    /*----------------------------- assess flip criterion */
        iptr_type static constexpr
            _last=pred_type::geom_dims+0 ;

        real_type static const _LTOL =
        std::pow(std::numeric_limits
            <real_type>::epsilon(), +.8) ;

        iptr_type  _iloc[3] ;
        _iloc[0] = _inod[0] ;
        _iloc[1] = _inod[1] ;
        _iloc[2] = _inod[2] ;

        algorithms::isort(
            &_iloc[0], &_iloc[3] ,
                std::less<iptr_type> ()) ;

        iptr_type  _jloc[3] ;
        _jloc[0] = _jnod[0] ;
        _jloc[1] = _jnod[1] ;
        _jloc[2] = _jnod[2] ;

        algorithms::isort(
            &_jloc[0], &_jloc[3] ,
                std::less<iptr_type> ()) ;

        real_type _ibal[_last + 1] ;
        pred_type::tri3_ball(_ibal ,
           &_mesh.node(_iloc[0]).pval(0) ,
           &_mesh.node(_iloc[1]).pval(0) ,
           &_mesh.node(_iloc[2]).pval(0)
            ) ;

        real_type _jbal[_last + 1] ;
        pred_type::tri3_ball(_jbal ,
           &_mesh.node(_jloc[0]).pval(0) ,
           &_mesh.node(_jloc[1]).pval(0) ,
           &_mesh.node(_jloc[2]).pval(0)
            ) ;

        real_type _null[_last] = {
            (real_type) +0.0 } ;

        pred_type::proj_node (
            _geom, _null, _ibal) ;

        pred_type::proj_node (
            _geom, _null, _jbal) ;

    /*----------------------------- measure dual distance */
        real_type  _ilen =
        pred_type::length_sq(_jbal,
           &_mesh.node(_inod[2]).pval(0));

        real_type  _jlen =
        pred_type::length_sq(_ibal,
           &_mesh.node(_jnod[2]).pval(0));

        real_type  _ipwr =
       _mesh.node(_inod[2]).pval (_last) ;

        real_type  _jpwr =
       _mesh.node(_jnod[2]).pval (_last) ;

        _ilen -=   _ipwr ;
        _jlen -=   _jpwr ;

        real_type  _btol =
       (real_type)+.5 * _LTOL * (_ilen + _jlen) ;

        real_type  _idel =
            _ilen - _jbal [_last] ;
        real_type  _jdel =
            _jlen - _ibal [_last] ;

        if (_idel >= -_btol || _jdel >= -_btol)
        {
            return false ;
        }
        else
        {
            return  true ;
        }
    }

    __static_call
    __normal_call bool_type need_flip_topology (
        geom_type &/*_geom*/ ,
        mesh_type &_mesh,
        iptr_list &_nmrk,
      __const_ptr(iptr_type) _enod,
      __const_ptr(iptr_type) _inod,
      __const_ptr(iptr_type) _jnod
        )
    {
        iptr_type static constexpr
            _DEG_TRIA3 = (iptr_type)+6 ;
    //  iptr_type static constexpr
    //      _DEG_QUAD4 = (iptr_type)+4 ;

    /*----------------------------- assess flip criterion */
        auto _ideg = _mesh.
            node(_enod[0]).ndeg(TRIA3_tag);
        auto _jdeg = _mesh.
            node(_enod[1]).ndeg(TRIA3_tag);
        auto _adeg = _mesh.
            node(_inod[2]).ndeg(TRIA3_tag);
        auto _bdeg = _mesh.
            node(_jnod[2]).ndeg(TRIA3_tag);

        if(_nmrk[_enod[0]] < 0) return false ;
        if(_nmrk[_enod[1]] < 0) return false ;
        if(_nmrk[_inod[2]] < 0) return false ;
        if(_nmrk[_jnod[2]] < 0) return false ;

        auto _ierr =
            (iptr_type)(_DEG_TRIA3-_ideg) ;
        auto _jerr =
            (iptr_type)(_DEG_TRIA3-_jdeg) ;
        auto _aerr =
            (iptr_type)(_DEG_TRIA3-_adeg) ;
        auto _berr =
            (iptr_type)(_DEG_TRIA3-_bdeg) ;

        auto _IDEG = _ideg - 1;
        auto _IERR =
            (iptr_type)(_DEG_TRIA3-_IDEG) ;
        auto _JDEG = _jdeg - 1;
        auto _JERR =
            (iptr_type)(_DEG_TRIA3-_JDEG) ;
        auto _ADEG = _adeg + 1;
        auto _AERR =
            (iptr_type)(_DEG_TRIA3-_ADEG) ;
        auto _BDEG = _bdeg + 1;
        auto _BERR =
            (iptr_type)(_DEG_TRIA3-_BDEG) ;

        auto _dnow = std::abs(_ierr) +  // current defect
                     std::abs(_jerr) +
                     std::abs(_aerr) +
                     std::abs(_berr) ;

        auto _dnew = std::abs(_IERR) +  // updated defect
                     std::abs(_JERR) +
                     std::abs(_AERR) +
                     std::abs(_BERR) ;

        return _dnew < _dnow - 1 ? true : false ;
    }

    /*
    --------------------------------------------------------
     * FLIP-T2T2: 2-simplex topological flip.
    --------------------------------------------------------
     */

  //__static_call
  //__normal_call void_type flip_t1q1 (

  //__static_call
  //__normal_call void_type flip_q2q2 (

    __static_call
    __normal_call void_type flip_t2t2 (
        geom_type &_geom ,
        mesh_type &_mesh ,
        iptr_list &_nmrk ,
        iptr_type *_enod ,
        conn_list &_conn ,
        conn_list &_next ,
        bool_type &_flip ,
        real_list &_qold ,
        real_list &_qnew ,
        flip_kind  _kind
        )
    {
        iptr_type _epos =  -1 ;
        _mesh.find_edge(_enod, _epos) ;

        if (_epos==-1) return ;

         auto _eptr =
        _mesh. edge().head() + _epos;

        if (_eptr->self() > 0) return ;

         auto _itri = _conn[0]._cell;
         auto _jtri = _conn[1]._cell;

         auto _iptr =
        _mesh. tri3().head() + _itri;
         auto _jptr =
        _mesh. tri3().head() + _jtri;

        if (_iptr->itag () !=
            _jptr->itag () )   return ;

    /*----------------------------- extract cell indexing */
        __unreferenced( _geom ) ;
        __unreferenced( _qold ) ;
        __unreferenced( _qnew ) ;

        iptr_type _inod [3] ;
        iptr_type _jnod [3] ;
        for(auto _inum = 3; _inum-- != 0; )
        {
            mesh_type::tri3_type::
            face_node(_inod, _inum, 2, 1) ;
            _inod[0] =
                _iptr->node(_inod[0]) ;
            _inod[1] =
                _iptr->node(_inod[1]) ;
            _inod[2] =
                _iptr->node(_inod[2]) ;

            if (_inod[2] != _enod[0])
            if (_inod[2] != _enod[1])
                break ;
        }
        for(auto _inum = 3; _inum-- != 0; )
        {
            mesh_type::tri3_type::
            face_node(_jnod, _inum, 2, 1) ;
            _jnod[0] =
                _jptr->node(_jnod[0]) ;
            _jnod[1] =
                _jptr->node(_jnod[1]) ;
            _jnod[2] =
                _jptr->node(_jnod[2]) ;

            if (_jnod[2] != _enod[0])
            if (_jnod[2] != _enod[1])
                break ;
        }

        assert( _inod[0] == _jnod[1] &&
                _inod[1] == _jnod[0] &&
        "ITER.FLIP-T2T2: bad orientation!") ;

    /*----------------------------- flip cells about edge */
        if ( _kind == laguerre_flip )
        {
        if(!need_flip_laguerre (
                _geom, _mesh ,
                _inod, _jnod)) return ;
        }
        else
        if ( _kind == delaunay_flip )
        {
        if(!need_flip_delaunay (
                _geom, _mesh ,
                _inod, _jnod)) return ;
        }
        else
        if ( _kind == topology_flip )
        {
        if(!need_flip_topology (
                _geom, _mesh ,
                _nmrk, _enod ,
                _inod, _jnod)) return ;
        }

    /*----------------------------- break if cells tangle */
        real_type _ivec[ +3 ] ;
        pred_type::tri3_norm (
           &_mesh.node(_inod[0]).pval(0),
           &_mesh.node(_jnod[2]).pval(0),
           &_mesh.node(_inod[2]).pval(0),
           &_ivec[0]  ) ;

        real_type _jvec[ +3 ] ;
        pred_type::tri3_norm (
           &_mesh.node(_jnod[0]).pval(0),
           &_mesh.node(_inod[2]).pval(0),
           &_mesh.node(_jnod[2]).pval(0),
           &_jvec[0]  ) ;

        real_type _sign =
        geometry::dot_3d (_ivec, _jvec) ;

        if (_sign <= (real_type)+0.) return ;

    /*----------------------------- update local topology */
        _mesh._pop_tri3(_itri) ;
        _mesh._pop_tri3(_jtri) ;

        _flip = true  ;

        typename mesh_type::tri3_type _tdat ;
        _tdat.node(0) = _inod[0] ;
        _tdat.node(1) = _jnod[2] ;
        _tdat.node(2) = _inod[2] ;

        _tdat.itag () = _iptr->itag() ;

        _next.push_tail(typename
            conn_list::data_type(
        _mesh.push_tri3(_tdat), TRIA3_tag)) ;

        _tdat.node(0) = _jnod[0] ;
        _tdat.node(1) = _inod[2] ;
        _tdat.node(2) = _jnod[2] ;

        _tdat.itag () = _jptr->itag() ;

        _next.push_tail(typename
            conn_list::data_type(
        _mesh.push_tri3(_tdat), TRIA3_tag)) ;
    }

    /*
    --------------------------------------------------------
     * FLIP-TRI3: "flip" mesh topology.
    --------------------------------------------------------
     */

    __static_call
    __inline_call void_type flip_tri3 (
        geom_type &_geom ,
        mesh_type &_mesh ,
        iptr_list &_nmrk ,
        iptr_type  _cell ,
        bool_type &_flip ,
        conn_list &_conn ,
        conn_list &_next ,
        real_list &_qold ,
        real_list &_qnew ,
        flip_kind  _kind
        )
    {
        auto _coin = (std::rand() % 3) + 1  ;
        for (auto _enum = _coin; _enum-- > 0; )
        {
            _flip = false ;
    /*----------------------------- flip cells about edge */
            iptr_type _enod[3] ;
            mesh_type::tri3_type::
            face_node(_enod, _enum, +2, +1) ;
            _enod[ 0] = _mesh.tri3(
                _cell).node( _enod[ 0]) ;
            _enod[ 1] = _mesh.tri3(
                _cell).node( _enod[ 1]) ;

            _conn.set_count(0) ;
            _next.set_count(0) ;
            _mesh.connect_2(
               &_enod[0], EDGE2_tag, _conn) ;

            if (_conn.count() != +2) continue ;

            if (_conn[0]._kind==TRIA3_tag &&
                _conn[1]._kind==TRIA3_tag )
            {
                flip_t2t2( _geom, _mesh, _nmrk, _enod,
                    _conn, _next, _flip,
                    _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
            else
            {
            //  flip_t1q1( _geom, _mesh, _nmrk, _enod,
            //      _conn, _next, _flip,
            //      _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
        }

        for (auto _enum = 3; _enum-- > _coin; )
        {
            _flip = false ;
    /*----------------------------- flip cells about edge */
            iptr_type _enod[3] ;
            mesh_type::tri3_type::
            face_node(_enod, _enum, +2, +1) ;
            _enod[ 0] = _mesh.tri3(
                _cell).node( _enod[ 0]) ;
            _enod[ 1] = _mesh.tri3(
                _cell).node( _enod[ 1]) ;

            _conn.set_count(0) ;
            _next.set_count(0) ;
            _mesh.connect_2(
               &_enod[0], EDGE2_tag, _conn) ;

            if (_conn.count() != +2) continue ;

            if (_conn[0]._kind==TRIA3_tag &&
                _conn[1]._kind==TRIA3_tag )
            {
                flip_t2t2( _geom, _mesh, _nmrk, _enod,
                    _conn, _next, _flip,
                    _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
            else
            {
            //  flip_t1q1( _geom, _mesh, _nmrk, _enod,
            //      _conn, _next, _flip,
            //      _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
        }

    }

    /*
    --------------------------------------------------------
     * FLIP-QUAD: "flip" mesh topology.
    --------------------------------------------------------
     */

    __static_call
    __inline_call void_type flip_quad (
        geom_type &/*_geom*/ ,
        mesh_type &_mesh ,
        iptr_list &_nmrk ,
        iptr_type  _cell ,
        bool_type &_flip ,
        conn_list &_conn ,
        conn_list &_next ,
        real_list &/*_qold*/ ,
        real_list &/*_qnew*/ ,
        flip_kind  /*_kind*/
        )
    {
        auto _coin = (std::rand() % 4) + 1  ;
        for (auto _enum = _coin; _enum-- > 0; )
        {
            _flip = false ;
    /*----------------------------- flip cells about edge */
            iptr_type _enod[4] ;
            mesh_type::quad_type::
            face_node(_enod, _enum, +2, +1) ;
            _enod[ 0] = _mesh.quad(
                _cell).node( _enod[ 0]) ;
            _enod[ 1] = _mesh.quad(
                _cell).node( _enod[ 1]) ;

            _conn.set_count(0) ;
            _next.set_count(0) ;
            _mesh.connect_2(
               &_enod[0], EDGE2_tag, _conn) ;

            if (_conn.count() != +2) continue ;

            if (_conn[0]._kind==QUAD4_tag &&
                _conn[1]._kind==QUAD4_tag )
            {
            //  flip_q2q2( _geom, _mesh, _enod, _nmrk,
            //      _conn, _next, _flip,
            //      _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
            else
            {
            //  flip_t1q1( _geom, _mesh, _enod, _nmrk,
            //      _conn, _next, _flip,
            //      _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
        }

        for (auto _enum = 4; _enum-- > _coin; )
        {
            _flip = false ;
    /*----------------------------- flip cells about edge */
            iptr_type _enod[4] ;
            mesh_type::quad_type::
            face_node(_enod, _enum, +2, +1) ;
            _enod[ 0] = _mesh.quad(
                _cell).node( _enod[ 0]) ;
            _enod[ 1] = _mesh.quad(
                _cell).node( _enod[ 1]) ;

            _conn.set_count(0) ;
            _next.set_count(0) ;
            _mesh.connect_2(
               &_enod[0], EDGE2_tag, _conn) ;

            if (_conn.count() != +2) continue ;

            if (_conn[0]._kind==QUAD4_tag &&
                _conn[1]._kind==QUAD4_tag )
            {
            //  flip_q2q2( _geom, _mesh, _enod, _nmrk,
            //      _conn, _next, _flip,
            //      _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
            else
            {
            //  flip_t1q1( _geom, _mesh, _enod, _nmrk,
            //      _conn, _next, _flip,
            //      _qold, _qnew, _kind)  ;

                if (_flip) return ;
            }
        }

    }



