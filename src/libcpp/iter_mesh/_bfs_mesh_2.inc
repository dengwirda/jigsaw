
    /*
    --------------------------------------------------------
     * _BFS-MESH-2: breadth-first-search mesh flips.
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
     * FLIP-SIGN: flip cells for +ve iter. cost fn.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type flip_next (
        mesh_type &_mesh,
        conn_list &_list,
        iptr_type  _icel,
        iptr_type *_inod,
        mark_list &_seen,
        conn_list &_conn
        )
    {
    /*------------------------- flip to match cell "sign" */
        _conn.set_count(0) ;
        _mesh.connect_2(
            _inod, EDGE2_tag, _conn);

        if (_conn.count() != 2) return;

        iptr_type _jcel, _jfac;
        char_type _kind;
        if ((iptr_type)
            _conn[0]._cell == _icel)
        {
            _jcel = _conn[ 1]._cell ;
            _kind = _conn[ 1]._kind ;
        }
        else
        {
            _jcel = _conn[ 0]._cell ;
            _kind = _conn[ 0]._kind ;
        }

    /*------------------------- match indexing along edge */
        if (_kind == TRIA3_tag)
        {
        if (_seen._tri3[_jcel] == +0)
        {
            _list.push_tail(typename
                conn_list::data_type(_jcel, _kind)) ;

            _seen._tri3[_jcel] =  +1;

            iptr_type _jnod[3] ;
            for(_jfac = 3 ; _jfac-- != 0; )
            {
            mesh_type::tri3_type::
            face_node(_jnod, _jfac, 2, 1) ;
            _jnod[ 0] = _mesh.
             tri3(_jcel).node(_jnod[0]) ;
            _jnod[ 1] = _mesh.
             tri3(_jcel).node(_jnod[1]) ;
            _jnod[ 2] = _mesh.
             tri3(_jcel).node(_jnod[2]) ;

            if (_jnod[ 2] != _inod[ 0])
            if (_jnod[ 2] != _inod[ 1])
                break  ;
            }

            if (_jnod[ 0] == _inod[ 0] &&
                _jnod[ 1] == _inod[ 1])
            {
                std::swap (
            _mesh.tri3(_jcel).node( 0),
            _mesh.tri3(_jcel).node( 1)) ;
            }
        }
        }
        else
        if (_kind == QUAD4_tag)
        {
        if (_seen._quad[_jcel] == +0)
        {
            _list.push_tail(typename
                conn_list::data_type(_jcel, _kind)) ;

            _seen._quad[_jcel] =  +1;

            // TODO: reverse quad indexes
        }
        }
    }

    __static_call
    __normal_call void_type flip_bfs_ (
        mesh_type &_mesh,
        conn_list &_list,
        mark_list &_seen
        )
    {
        conn_list  _conn;
        for ( ; !_list.empty() ; )
        {
    /*------------------------- flip via BFS on cell adj. */
            typename conn_list::data_type _cell ;
            _list._pop_tail(_cell) ;

            if (_cell._kind == TRIA3_tag)
            {
            for(auto _enum = +3; _enum-- != +0; )
            {
                auto _cpos = _cell._cell;

                iptr_type _enod[3] ;
                mesh_type::tri3_type::
                face_node(_enod, _enum, +2, +1) ;
                _enod[ 0] = _mesh.tri3(
                    _cpos).node( _enod[ 0]) ;
                _enod[ 1] = _mesh.tri3(
                    _cpos).node( _enod[ 1]) ;

                flip_next( _mesh, _list,
                    _cpos, _enod, _seen, _conn) ;
            }
            }
            else
            if (_cell._kind == QUAD4_tag)
            {
            for(auto _enum = +4; _enum-- != +0; )
            {
                auto _cpos = _cell._cell;

                iptr_type _enod[4] ;
                mesh_type::quad_type::
                face_node(_enod, _enum, +2, +1) ;
                _enod[ 0] = _mesh.quad(
                    _cpos).node( _enod[ 0]) ;
                _enod[ 1] = _mesh.quad(
                    _cpos).node( _enod[ 1]) ;

                flip_next( _mesh, _list,
                    _cpos, _enod, _seen, _conn) ;
            }
            }
        }
    }

    __static_call
    __normal_call void_type flip_sign (
        mesh_type &_mesh
        )
    {
    #   define SEENTRI3( _CELL) _seen._tri3[_CELL]
    #   define SEENQUAD( _CELL) _seen._quad[_CELL]

        conn_list _list;
        mark_list _seen;
        init_mark(_mesh, _seen, +0) ;

    /*------------------------- flip to correct mesh sign */
        iptr_type _tnum  = +0 ;
        for (auto _tria  = _mesh.tri3().head();
                  _tria != _mesh.tri3().tend();
                ++_tria, ++_tnum )
        {
            if (_tria->mark  () < +0) continue;
            if (SEENTRI3(_tnum) > +0) continue;

        /*--------------------- flip seed for +ve quality */
            real_type _cost =
                 pred_type::tri3_cost (
               &_mesh. node(
                _tria->node(0)).pval(0),
               &_mesh. node(
                _tria->node(1)).pval(0),
               &_mesh. node(
                _tria->node(2)).pval(0),
                typename
                pred_type::cell_kind()) ;

            if (_cost < (real_type) +0.)
            {
                std::swap (
                    _tria->node(0),
                        _tria->node(1)) ;
            }

        /*--------------------- a BFS from seed via topo. */
            _list.push_tail(typename conn_list
                ::data_type (_tnum, TRIA3_tag)) ;

            SEENTRI3( _tnum ) = 1 ;

            flip_bfs_(_mesh, _list, _seen) ;
        }

        iptr_type _qnum  = +0 ;
        for (auto _quad  = _mesh.quad().head();
                  _quad != _mesh.quad().tend();
                ++_quad, ++_qnum )
        {
            if (_quad->mark  () < +0) continue;
            if (SEENQUAD(_qnum) > +0) continue;

        /*--------------------- flip seed for +ve quality */
            real_type _cost =
                 pred_type::quad_cost (
               &_mesh. node(
                _quad->node(0)).pval(0),
               &_mesh. node(
                _quad->node(1)).pval(0),
               &_mesh. node(
                _quad->node(2)).pval(0),
               &_mesh. node(
                _quad->node(3)).pval(0),
                typename
                pred_type::cell_kind()) ;

            if (_cost < (real_type) +0.)
            {
                // TODO: reverse quad indexes
            }

        /*--------------------- a BFS from seed via topo. */
            _list.push_tail(typename conn_list
                ::data_type (_qnum, QUAD4_tag)) ;

            SEENQUAD( _qnum ) = 1 ;

            flip_bfs_(_mesh, _list, _seen) ;
        }

    #   undef   SEENQUAD
    #   undef   SEENTRI3
    }



