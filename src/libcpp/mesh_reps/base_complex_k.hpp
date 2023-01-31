
    /*
    --------------------------------------------------------
     * BASE-COMPLEX-k: mesh complex utilities in T^k.
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
     * Copyright 2020--2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __BASE_COMPLEX_K__
#   define __BASE_COMPLEX_K__

    namespace mesh {

    /*
    --------------------------------------------------------
     * BASE-COMPLEX-k utilities for k-cell complexes.
    --------------------------------------------------------
     */

    // Essentially an "unrolled" base class for mesh/tria
    // complexes...

    // Concrete tria_complex_k and mesh_complex_k classes
    // are composed of the various bits and pieces
    // implemented here, depending on dimensionality, and
    // whether simplical vs mixed-cell, etc.


    struct tria_cell {} ;       // simplicial vs mixed...
    struct mesh_cell {} ;

#   define TOPO_DIMS mesh_type::topo_dims

#   define IPTR_TYPE typename \
                     mesh_type::iptr_type
#   define INDX_TYPE typename \
                     mesh_type::indx_type
#   define REAL_TYPE typename \
                     mesh_type::real_type

#   define NODE_TYPE typename \
                     mesh_type::node_type
#   define EDGE_TYPE typename \
                     mesh_type::edge_type
#   define TRI3_TYPE typename \
                     mesh_type::tri3_type
#   define QUAD_TYPE typename \
                     mesh_type::quad_type
#   define TRI4_TYPE typename \
                     mesh_type::tri4_type
#   define HEXA_TYPE typename \
                     mesh_type::hexa_type
#   define WEDG_TYPE typename \
                     mesh_type::wedg_type
#   define PYRA_TYPE typename \
                     mesh_type::pyra_type

    /*
    --------------------------------------------------------
     * CONN-LIST: cell/kind pair for adjacency.
    --------------------------------------------------------
     */

    class conn_pair
        {
    /*----------------------- handle adj. mixed cell kind */
        public  :
        uint32_t               _cell ;
        char_type              _kind ;
        public  :
    /*----------------------- build conn _obj. from _src. */
        __inline_call conn_pair (
            uint32_t  _csrc = +0,
            char_type _ksrc = +0
            ) : _cell(_csrc), _kind(_ksrc) {}
        } ;

    /*
    --------------------------------------------------------
     * INIT-LIST: re-size adj. lists on demand.
    --------------------------------------------------------
     */

    template <
        typename  conn_list
             >
    __inline_call void_type init_list (
        conn_list&_list,
            typename
        conn_list::size_type _lpos
        )
    {
        if (_lpos >= _list._lptr.count())
        _list._lptr.set_count (
            _lpos + 1,
        containers::loose_alloc, nullptr) ;
    }

    /*
    --------------------------------------------------------
     * _POP-LIST: delete item from adj. list.
    --------------------------------------------------------
     */

    template <
        typename  conn_list
             >
    __normal_call void_type _pop_list (
        conn_list &_list,
            typename
        conn_list::size_type _lpos ,
            typename
        conn_list::data_type _item
        )
    {
    /*--------------------------- scan list and pop match */
        typename conn_list::_write_it
            _prev = _list.hend(_lpos) ,
            _iter = _list.head(_lpos) ,
            _tend = _list.tend(_lpos) ;

        for ( ; _iter != _tend;
                _prev  = _iter, ++_iter )
        {
        if (_iter->_cell == _item._cell &&
            _iter->_kind == _item._kind )
        {
            _list._pop (
                _prev, _iter, _lpos ) ;

            break  ;
        }
        }
    }

    /*
    --------------------------------------------------------
     * GET-NODE: "create" new node, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_node (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffN1.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffN1._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llN1.push_tail() ;
        }

        _mesh._llN1[_ipos].mark() = 0 ;
        _mesh._llN1[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-EDGE: "create" new edge, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_edge (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffE2.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffE2._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llE2.push_tail() ;
        }

        _mesh._llE2[_ipos].mark() = 0 ;
        _mesh._llE2[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-TRIA: "create" new tria, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_tri3 (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffT3.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffT3._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llT3.push_tail() ;
        }

        _mesh._llT3[_ipos].mark() = 0 ;
        _mesh._llT3[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-QUAD: "create" new quad, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_quad (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffQ4.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffQ4._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llQ4.push_tail() ;
        }

        _mesh._llQ4[_ipos].mark() = 0 ;
        _mesh._llQ4[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-TRIA: "create" new tria, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_tri4 (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffT4.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffT4._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llT4.push_tail() ;
        }

        _mesh._llT4[_ipos].mark() = 0 ;
        _mesh._llT4[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-HEXA: "create" new hexa, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_hexa (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffH8.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffH8._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llH8.push_tail() ;
        }

        _mesh._llH8[_ipos].mark() = 0 ;
        _mesh._llH8[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-WEDG: "create" new wedg, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_wedg (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffW6.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffW6._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llW6.push_tail() ;
        }

        _mesh._llW6[_ipos].mark() = 0 ;
        _mesh._llW6[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }

    /*
    --------------------------------------------------------
     * GET-PYRA: "create" new pyra, push onto active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call INDX_TYPE _get_pyra (
        mesh_type &_mesh
        )
    {
        INDX_TYPE _ipos;
        if (_mesh._ffP5.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            _mesh._ffP5._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos =
            _mesh._llP5.push_tail() ;
        }

        _mesh._llP5[_ipos].mark() = 0 ;
        _mesh._llP5[_ipos].self() = 0 ;

        return ( _ipos ) ;
    }



    /*
    --------------------------------------------------------
     * PUT-ITEM: "delete" old item, _pop from active set.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_node (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffN1.push_tail (_ipos);
        _mesh._llN1[_ipos].mark() = -1 ;
        _mesh._llN1[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_edge (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffE2.push_tail (_ipos);
        _mesh._llE2[_ipos].mark() = -1 ;
        _mesh._llE2[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_tri3 (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffT3.push_tail (_ipos);
        _mesh._llT3[_ipos].mark() = -1 ;
        _mesh._llT3[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_quad (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffQ4.push_tail (_ipos);
        _mesh._llQ4[_ipos].mark() = -1 ;
        _mesh._llQ4[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_tri4 (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffT4.push_tail (_ipos);
        _mesh._llT4[_ipos].mark() = -1 ;
        _mesh._llT4[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_hexa (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffH8.push_tail (_ipos);
        _mesh._llH8[_ipos].mark() = -1 ;
        _mesh._llH8[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_wedg (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffW6.push_tail (_ipos);
        _mesh._llW6[_ipos].mark() = -1 ;
        _mesh._llW6[_ipos].self() = -1 ;
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _put_pyra (
        mesh_type&_mesh ,
        INDX_TYPE _ipos
        )
    {
        _mesh._ffP5.push_tail (_ipos);
        _mesh._llP5[_ipos].mark() = -1 ;
        _mesh._llP5[_ipos].self() = -1 ;
    }

    /*
    --------------------------------------------------------
     * PUSH-NODE: append new 0-node to complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call IPTR_TYPE push_node (
        mesh_type      & _mesh,
        NODE_TYPE const& _ndat,
        bool_type _link = true,
        IPTR_TYPE _itop =  -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _ipos =  -1 ;

        __unreferenced(_kind) ;     // maybe, per TOPODIM

        if (!_link)
        {

        __assert( _itop == -1 &&
        "mesh-complex: non-top node!" ) ;

    /*------------------------ init. external d-face data */
       _ipos = (IPTR_TYPE)_get_node(_mesh) ;

        _mesh._llN1[_ipos]  = _ndat ;
        _mesh._llN1[_ipos].null() ;
        _mesh._llN1[_ipos].mark() = 0 ;

        _mesh._llN1[_ipos].self() = 1 ;
        _mesh.
        _llN1[_ipos].node(0) =_ipos ;

    /*------------------------ init. local adj. index set */
        if constexpr (TOPO_DIMS > +0)
        init_list(_mesh._aaN1, _ipos) ;

        }
        else
        {

        if (_itop == -1)
        {
    /*------------------------ init. external d-face data */
       _ipos = (IPTR_TYPE)_get_node(_mesh) ;

        _mesh._llN1[_ipos]  = _ndat ;
        _mesh._llN1[_ipos].null() ;
        _mesh._llN1[_ipos].mark() = 0 ;
        _mesh._llN1[_ipos].self() = 1 ;

        _mesh.
        _llN1[_ipos].node(0) = _ipos;

    /*------------------------ init. local adj. index set */
        if constexpr (TOPO_DIMS > +0)
        init_list(_mesh._aaN1, _ipos) ;

        }
        else
        {
    /*------------------------ init. internal d-face data */
        __assert(find_node(
            _mesh, &_ndat.node(0), _ipos ) &&
            "mesh-complex: node not here") ;

    /*------------------------ append index to adj. lists */
        if constexpr (TOPO_DIMS > +0)
        _mesh._aaN1.push(conn_pair(
            _itop, _kind), _ndat.node( 0)) ;

        }

        }

        return _ipos ;
    }

    /*
    --------------------------------------------------------
     * PUSH-EDGE: append new 1-edge to complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call IPTR_TYPE push_edge (
        mesh_type      & _mesh,
        EDGE_TYPE const& _edat,
        bool_type _link = true,
        IPTR_TYPE _itop =  -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        typedef typename
            mesh_type::edge_maps EDGE_MAPS;

        __unreferenced(_kind);      // maybe, per TOPODIM

        IPTR_TYPE _ipos =  -1, _npos;

        if (_itop == -1)
        {
        _mesh.node(_edat.node(0)).
                ndeg(EDGE2_tag)  += 1 ;
        _mesh.node(_edat.node(1)).
                ndeg(EDGE2_tag)  += 1 ;
        }

        if (!_link)
        {

        __assert( _itop == -1 &&
        "mesh-complex: non-top cell!" ) ;

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_edge(_mesh) ;

        _mesh._llE2[_ipos]  = _edat ;
        _mesh._llE2[_ipos].mark() = 0 ;
        _mesh._llE2[_ipos].self() = 1 ;

    /*------------------------ init. local adj. index set */
        if constexpr (TOPO_DIMS > +1)
        init_list(_mesh._aaE2, _ipos) ;

        }
        else
        {

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_edge(_mesh) ;

        _mesh._llE2[_ipos]  = _edat ;
        _mesh._llE2[_ipos].mark() = 0 ;

        if (_itop == -1)
        _mesh._llE2[_ipos].self() = 1 ;
        else
        _mesh._llE2[_ipos].self() = 0 ;

        typename
        EDGE_MAPS::_write_it  _same ;
        if (_mesh.
            _mmE2.find(_ipos, _same))
        {
    /*---- existing d-face found: append to existing data */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, copy data */
                _mesh._llE2[*_same] =
                _mesh._llE2[ _ipos] ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +1)
                _mesh._aaE2.push(
            conn_pair(_itop, _kind),*_same) ;
            }

            _put_edge(_mesh, _ipos) ;
        }
        else
        {
    /*---- d-face is new: push and descend to (d-1)-faces */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, keep data */
            if constexpr (TOPO_DIMS > +1)
                init_list(
                    _mesh._aaE2, _ipos) ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +1)
                _mesh._aaE2.push(
            conn_pair(_itop, _kind), _ipos) ;
            }

        /*-------------- descend into (d-1)-face data */
            for (_npos = 2; _npos-- != 0; )
            {
                NODE_TYPE _ndat;
                _ndat.node(0) =
                    _edat.node(_npos) ;

                push_node(_mesh, _ndat,
                    _link,_ipos, EDGE2_tag) ;
            }

        /*-------------- push new face data onto hash */
            _mesh._mmE2.push(_ipos) ;
        }

        }

        return _ipos ;
    }

    /*
    --------------------------------------------------------
     * PUSH-TRIA: append new 2-tria to complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call IPTR_TYPE push_tri3 (
        mesh_type      & _mesh,
        TRI3_TYPE const& _tdat,
        bool_type _link = true,
        IPTR_TYPE _itop =  -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        typedef typename
            mesh_type::tri3_maps TRI3_MAPS;

        __unreferenced(_kind);      // maybe, per TOPODIM

        IPTR_TYPE _ipos =  -1, _epos;

        if (_itop == -1)
        {
        _mesh.node(_tdat.node(0)).
                ndeg(TRIA3_tag)  += 1 ;
        _mesh.node(_tdat.node(1)).
                ndeg(TRIA3_tag)  += 1 ;
        _mesh.node(_tdat.node(2)).
                ndeg(TRIA3_tag)  += 1 ;
        }

        if (!_link)
        {

        __assert( _itop == -1 &&
        "mesh-complex: non-top cell!" ) ;

        __unreferenced(_kind) ;

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_tri3(_mesh) ;

        _mesh._llT3[_ipos]  = _tdat ;
        _mesh._llT3[_ipos].mark() = 0 ;
        _mesh._llT3[_ipos].self() = 1 ;

    /*------------------------ init. local adj. index set */
        if constexpr (TOPO_DIMS > +2)
        init_list(_mesh._aaT3, _ipos) ;

        }
        else
        {

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_tri3(_mesh) ;

        _mesh._llT3[_ipos]  = _tdat ;
        _mesh._llT3[_ipos].mark() = 0 ;

        if (_itop == -1)
        _mesh._llT3[_ipos].self() = 1 ;
        else
        _mesh._llT3[_ipos].self() = 0 ;

        typename
        TRI3_MAPS::_write_it  _same ;
        if (_mesh.
            _mmT3.find(_ipos, _same))
        {
    /*---- existing d-face found: append to existing data */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, copy data */
                _mesh._llT3[*_same] =
                _mesh._llT3[ _ipos] ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +2)
                _mesh._aaT3.push(
            conn_pair(_itop, _kind),*_same) ;
            }

            _put_tri3(_mesh, _ipos) ;
        }
        else
        {
    /*---- d-face is new: push and descend to (d-1)-faces */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, keep data */
            if constexpr (TOPO_DIMS > +2)
                init_list(
                    _mesh._aaT3, _ipos) ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +2)
                _mesh._aaT3.push(
            conn_pair(_itop, _kind), _ipos) ;
            }

        /*-------------- descend into (d-1)-face data */
            for (_epos = +3; _epos-- != 0 ; )
            {
                IPTR_TYPE _enod [3];
                mesh_type::tri3_type::
                face_node(_enod, _epos, 2, 1) ;

                EDGE_TYPE _edat;
                _edat.node(0) =
                    _tdat.node(_enod[0]) ;
                _edat.node(1) =
                    _tdat.node(_enod[1]) ;

                push_edge(_mesh, _edat,
                    _link,_ipos, TRIA3_tag) ;
            }

        /*-------------- push new face data onto hash */
            _mesh._mmT3.push(_ipos) ;
        }

        }

        return _ipos ;
    }

    /*
    --------------------------------------------------------
     * PUSH-QUAD: append new 2-quad to complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call IPTR_TYPE push_quad (
        mesh_type      & _mesh,
        QUAD_TYPE const& _qdat,
        bool_type _link = true,
        IPTR_TYPE _itop =  -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        typedef typename
            mesh_type::quad_maps QUAD_MAPS;

        __unreferenced(_kind);      // maybe, per TOPODIM

        IPTR_TYPE _ipos =  -1, _epos;

        if (_itop == -1)
        {
        _mesh.node(_qdat.node(0)).
                ndeg(QUAD4_tag)  += 1 ;
        _mesh.node(_qdat.node(1)).
                ndeg(QUAD4_tag)  += 1 ;
        _mesh.node(_qdat.node(2)).
                ndeg(QUAD4_tag)  += 1 ;
        _mesh.node(_qdat.node(3)).
                ndeg(QUAD4_tag)  += 1 ;
        }

        if (!_link)
        {

        __assert( _itop == -1 &&
        "mesh-complex: non-top cell!" ) ;

        __unreferenced(_kind) ;

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_quad(_mesh) ;

        _mesh._llQ4[_ipos]  = _qdat ;
        _mesh._llQ4[_ipos].mark() = 0 ;
        _mesh._llQ4[_ipos].self() = 1 ;

    /*------------------------ init. local adj. index set */
        if constexpr (TOPO_DIMS > +2)
        init_list(_mesh._aaQ4, _ipos) ;

        }
        else
        {

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_quad(_mesh) ;

        _mesh._llQ4[_ipos]  = _qdat ;
        _mesh._llQ4[_ipos].mark() = 0 ;

        if (_itop == -1)
        _mesh._llQ4[_ipos].self() = 1 ;
        else
        _mesh._llQ4[_ipos].self() = 0 ;

        typename
        QUAD_MAPS::_write_it  _same ;
        if (_mesh.
            _mmQ4.find(_ipos, _same))
        {
    /*---- existing d-face found: append to existing data */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, copy data */
                _mesh._llQ4[*_same] =
                _mesh._llQ4[ _ipos] ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +2)
                _mesh._aaQ4.push(
            conn_pair(_itop, _kind),*_same) ;
            }

            _put_quad(_mesh, _ipos) ;
        }
        else
        {
    /*---- d-face is new: push and descend to (d-1)-faces */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, keep data */
            if constexpr (TOPO_DIMS > +2)
                init_list(
                    _mesh._aaQ4, _ipos) ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +2)
                _mesh._aaQ4.push(
            conn_pair(_itop, _kind), _ipos) ;
            }

        /*-------------- descend into (d-1)-face data */
            for (_epos = +4; _epos-- != 0 ; )
            {
                IPTR_TYPE _enod [4];
                mesh_type::quad_type::
                face_node(_enod, _epos, 2, 1) ;

                EDGE_TYPE _edat;
                _edat.node(0) =
                    _qdat.node(_enod[0]) ;
                _edat.node(1) =
                    _qdat.node(_enod[1]) ;

                push_edge(_mesh, _edat,
                    _link,_ipos, QUAD4_tag) ;
            }

        /*-------------- push new face data onto hash */
            _mesh._mmQ4.push(_ipos) ;
        }

        }

        return _ipos ;
    }

    /*
    --------------------------------------------------------
     * PUSH-TRIA: append new 3-tria to complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call IPTR_TYPE push_tri4 (
        mesh_type      & _mesh,
        TRI4_TYPE const& _tdat,
        bool_type _link = true,
        IPTR_TYPE _itop =  -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        typedef typename
            mesh_type::tri4_maps TRI4_MAPS;

        __unreferenced(_kind);      // maybe, per TOPODIM

        IPTR_TYPE _ipos =  -1, _fpos;

        if (_itop == -1)
        {
        _mesh.node(_tdat.node(0)).
                ndeg(TRIA4_tag)  += 1 ;
        _mesh.node(_tdat.node(1)).
                ndeg(TRIA4_tag)  += 1 ;
        _mesh.node(_tdat.node(2)).
                ndeg(TRIA4_tag)  += 1 ;
        _mesh.node(_tdat.node(3)).
                ndeg(TRIA4_tag)  += 1 ;
        }

        if (!_link)
        {

        __assert( _itop == -1 &&
        "mesh-complex: non-top cell!" ) ;

        __unreferenced(_kind) ;

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_tri4(_mesh) ;

        _mesh._llT4[_ipos]  = _tdat ;
        _mesh._llT4[_ipos].mark() = 0 ;
        _mesh._llT4[_ipos].self() = 1 ;

    /*------------------------ init. local adj. index set */
        if constexpr (TOPO_DIMS > +3)
        init_list(_mesh._aaT4, _ipos) ;

        }
        else
        {

    /*------------------------ init. external d-face data */
        _ipos = (IPTR_TYPE)_get_tri4(_mesh) ;

        _mesh._llT4[_ipos]  = _tdat ;
        _mesh._llT4[_ipos].mark() = 0 ;

        if (_itop == -1)
        _mesh._llT4[_ipos].self() = 1 ;
        else
        _mesh._llT4[_ipos].self() = 0 ;

        typename
        TRI4_MAPS::_write_it  _same ;
        if (_mesh.
            _mmT4.find(_ipos, _same))
        {
    /*---- existing d-face found: append to existing data */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, copy data */
                _mesh._llT4[*_same] =
                _mesh._llT4[ _ipos] ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +3)
                _mesh._aaT4.push(
            conn_pair(_itop, _kind),*_same) ;
            }

            _put_tri4(_mesh, _ipos) ;
        }
        else
        {
    /*---- d-face is new: push and descend to (d-1)-faces */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, keep data */
            if constexpr (TOPO_DIMS > +3)
                init_list(_mesh,
                    _mesh._aaT4, _ipos) ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
            if constexpr (TOPO_DIMS > +3)
                _mesh._aaT4.push(
            conn_pair(_itop, _kind), _ipos) ;
            }

        /*-------------- descend into (d-1)-face data */
            for (_fpos = +4; _fpos-- != 0 ; )
            {
                IPTR_TYPE _fnod [4];
                mesh_type::tri4_type::
                face_node(_fnod, _fpos, 3, 2) ;

                TRI3_TYPE _fdat;
                _fdat.node(0) =
                    _tdat.node(_fnod[0]) ;
                _fdat.node(1) =
                    _tdat.node(_fnod[1]) ;
                _fdat.node(2) =
                    _tdat.node(_fnod[2]) ;

                push_tri3(_mesh, _fdat,
                    _link,_ipos, TRIA4_tag) ;
            }

        /*-------------- push new face data onto hash */
            _mesh._mmT4.push(_ipos) ;
        }

        }

        return _ipos ;
    }

    /*
    --------------------------------------------------------
     * MAKE-LINK: build item-to-item adj.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call void_type link_edge (
        mesh_type&_mesh
        )
    {
        IPTR_TYPE _epos  = +0 ;
        for (auto _iter  = _mesh._llE2.head() ;
                  _iter != _mesh._llE2.tend() ;
                ++_iter, ++_epos  )
        {
            if (_iter->mark() < 0) continue ;

        /*-------------- push face data onto hash set */
            _mesh._mmE2.push(_epos) ;

        /*-------------- descend into (d-1)-face data */
            IPTR_TYPE _ipos;
            for (_ipos = +2; _ipos-- != 0 ; )
            {
                NODE_TYPE _ndat;
                _ndat.node(0) =
                    _iter->node( _ipos) ;

                push_node(_mesh, _ndat,
                    true, _epos, EDGE2_tag) ;
            }
        }
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type link_tri3 (
        mesh_type&_mesh
        )
    {
        IPTR_TYPE _tpos  = +0 ;
        for (auto _iter  = _mesh._llT3.head() ;
                  _iter != _mesh._llT3.tend() ;
                ++_iter, ++_tpos  )
        {
            if (_iter->mark() < 0) continue ;

        /*-------------- push face data onto hash set */
            _mesh._mmT3.push(_tpos) ;

        /*-------------- descend into (d-1)-face data */
            IPTR_TYPE _ipos;
            for (_ipos = +3; _ipos-- != 0 ; )
            {
                IPTR_TYPE  _enod [3];
                mesh_type::tri3_type::
                face_node(_enod, _ipos, 2, 1) ;

                EDGE_TYPE  _edat ;
                _edat.node(0) =
                   _iter->node(_enod[0]) ;
                _edat.node(1) =
                   _iter->node(_enod[1]) ;

                push_edge(_mesh, _edat,
                    true, _tpos, TRIA3_tag) ;
            }
        }
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type link_quad (
        mesh_type&_mesh
        )
    {
        IPTR_TYPE _qpos = +0 ;
        for (auto _iter  = _mesh._llQ4.head() ;
                  _iter != _mesh._llQ4.tend() ;
                ++_iter, ++_qpos  )
        {
            if (_iter->mark() < 0) continue ;

        /*-------------- push face data onto hash set */
            _mesh._mmQ4.push(_qpos) ;

        /*-------------- descend into (d-1)-face data */
            IPTR_TYPE _ipos;
            for (_ipos = +4; _ipos-- != 0 ; )
            {
                IPTR_TYPE  _enod [4];
                mesh_type::quad_type::
                face_node(_enod, _ipos, 2, 1) ;

                EDGE_TYPE  _edat ;
                _edat.node(0) =
                   _iter->node(_enod[0]) ;
                _edat.node(1) =
                   _iter->node(_enod[1]) ;

                push_edge(_mesh, _edat,
                    true, _qpos, QUAD4_tag) ;
            }
        }
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type link_tri4 (
        mesh_type&_mesh
        )
    {
        IPTR_TYPE _Tpos = +0 ;
        for (auto _iter  = _mesh._llT4.head() ;
                  _iter != _mesh._llT4.tend() ;
                ++_iter, ++_Tpos  )
        {
            if (_iter->mark() < 0) continue ;

        /*-------------- push face data onto hash set */
            _mesh._mmT4.push(_Tpos) ;

        /*-------------- descend into (d-1)-face data */
            IPTR_TYPE _ipos;
            for (_ipos = +4; _ipos-- != 0 ; )
            {
                IPTR_TYPE  _fnod [4];
                mesh_type::tri4_type::
                face_node(_fnod, _ipos, 3, 2) ;

                TRI3_TYPE  _fdat ;
                _fdat.node(0) =
                   _iter->node(_fnod[0]) ;
                _fdat.node(1) =
                   _iter->node(_fnod[1]) ;
                _fdat.node(2) =
                   _iter->node(_fnod[2]) ;

                push_tri3(_mesh, _fdat,
                    true, _Tpos, TRIA4_tag) ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * NULL-ITEM: TRUE if item is delete-able.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_node (
        mesh_type&_mesh ,
        IPTR_TYPE _npos
        )
    {
        if constexpr (TOPO_DIMS >  +0 )
        {
        return _mesh._llN1[_npos].self() == 0
            && _mesh._aaN1.empty(_npos);
        }
        else
        return _mesh._llN1[_npos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_edge (
        mesh_type&_mesh ,
        IPTR_TYPE _epos
        )
    {
        if constexpr (TOPO_DIMS >  +1 )
        {
        return _mesh._llE2[_epos].self() == 0
            && _mesh._aaE2.empty(_epos);
        }
        else
        return _mesh._llE2[_epos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_tri3 (
        mesh_type&_mesh ,
        IPTR_TYPE _tpos
        )
    {
        if constexpr (TOPO_DIMS >  +2 )
        {
        return _mesh._llT3[_tpos].self() == 0
            && _mesh._aaT3.empty(_tpos);
        }
        else
        return _mesh._llT3[_tpos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_quad (
        mesh_type&_mesh ,
        IPTR_TYPE _qpos
        )
    {
        if constexpr (TOPO_DIMS >  +2 )
        {
        return _mesh._llQ4[_qpos].self() == 0
            && _mesh._aaQ4.empty(_qpos);
        }
        else
        return _mesh._llQ4[_qpos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_tri4 (
        mesh_type&_mesh ,
        IPTR_TYPE _tpos
        )
    {
        if constexpr (TOPO_DIMS >  +3 )
        {
        return _mesh._llT4[_tpos].self() == 0
            && _mesh._aaT4.empty(_tpos);
        }
        else
        return _mesh._llT4[_tpos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_hexa (
        mesh_type&_mesh ,
        IPTR_TYPE _hpos
        )
    {
        if constexpr (TOPO_DIMS >  +3 )
        {
        return _mesh._llH8[_hpos].self() == 0
            && _mesh._aaH8.empty(_hpos);
        }
        else
        return _mesh._llH8[_hpos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_wedg (
        mesh_type&_mesh ,
        IPTR_TYPE _wpos
        )
    {
        if constexpr (TOPO_DIMS >  +3 )
        {
        return _mesh._llW6[_wpos].self() == 0
            && _mesh._aaW6.empty(_wpos);
        }
        else
        return _mesh._llW6[_wpos].self() == 0;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type null_pyra (
        mesh_type&_mesh ,
        IPTR_TYPE _ppos
        )
    {
        if constexpr (TOPO_DIMS >  +3 )
        {
        return _mesh._llP5[_ppos].self() == 0
            && _mesh._aaP5.empty(_ppos);
        }
        else
        return _mesh._llP5[_ppos].self() == 0;
    }

    /*
    --------------------------------------------------------
     * _POP-NODE: delete 0-node from complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call void_type _pop_node (
        mesh_type&_mesh ,
        IPTR_TYPE const*_nptr,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _npos = -1 ;

        IPTR_TYPE  _node[1];
        _node[0] = _nptr[0];

    /*-------------------------- find current 0-node pos. */
        if ( !find_node (_mesh, _node, _npos))
        {
            return ;
        }

    /*-------------------------- _pop current 0-node data */
        _pop_node(_mesh, _npos, _itop, _kind);
    }

    template <
        typename  mesh_type
             >
    __inline_call void_type _pop_node (
        mesh_type&_mesh ,
        IPTR_TYPE _npos ,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        __unreferenced(_kind);      // maybe, per TOPODIM

        if (_mesh.
            _llN1[_npos].mark() < +0 )
            return ;                // is already removed

        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
        if constexpr ( TOPO_DIMS > 0 )
        {
            _pop_list(_mesh. _aaN1 ,
                _npos,
            conn_pair(_itop, _kind)) ;
        }
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            _mesh.
            _llN1[_npos].self() = +0 ;
        }

        if (null_node(_mesh, _npos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */
            _put_node(_mesh, _npos);
        }
    }

    /*
    --------------------------------------------------------
     * _POP-EDGE: delete 1-edge from complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call void_type _pop_edge (
        mesh_type&_mesh ,
        IPTR_TYPE const*_nptr,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _epos = -1 ;

        IPTR_TYPE  _node[2];
        _node[0] = _nptr[0];
        _node[1] = _nptr[1];

    /*-------------------------- find current 1-cell pos. */
        if ( !find_edge (_mesh, _node, _epos))
        {
            return ;
        }

    /*-------------------------- _pop current 1-cell data */
        _pop_edge(_mesh, _epos, _itop, _kind);
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type _pop_edge (
        mesh_type&_mesh ,
        IPTR_TYPE _epos ,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _npos = -1 ;

        __unreferenced(_kind);      // maybe, per TOPODIM

        if (_mesh.
            _llE2[_epos].mark() < +0 )
            return ;                // is already removed

        if (_itop == -1)
        {
        _mesh.node(_mesh.edge(_epos).node(0)).
                ndeg(EDGE2_tag) -= 1 ;
        _mesh.node(_mesh.edge(_epos).node(1)).
                ndeg(EDGE2_tag) -= 1 ;
        }

        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
        if constexpr ( TOPO_DIMS > 1 )
            _pop_list(_mesh. _aaE2 ,
                _epos,
            conn_pair(_itop, _kind)) ;
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            _mesh.
            _llE2[_epos].self() = +0 ;
        }

        if (null_edge(_mesh, _epos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */
        IPTR_TYPE _same ;
        _mesh._mmE2._pop(_epos, _same) ;

        IPTR_TYPE _node [ 2] ;
        _node[0] =
            _mesh._llE2[_epos].node(0) ;
        _node[1] =
            _mesh._llE2[_epos].node(1) ;

        _put_edge(_mesh, _epos) ;

        for (_npos = 2; _npos-- != 0; )
        {
            _pop_node(_mesh,
                &_node [_npos],
                    _epos, EDGE2_tag) ;
        }
        }
    }

    /*
    --------------------------------------------------------
     * _POP-TRIA: delete 2-tria from complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call void_type _pop_tri3 (
        mesh_type&_mesh ,
        IPTR_TYPE const*_nptr,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _tpos = -1 ;

        IPTR_TYPE  _node[3];
        _node[0] = _nptr[0];
        _node[1] = _nptr[1];
        _node[2] = _nptr[2];

    /*-------------------------- find current 2-cell pos. */
        if ( !find_tri3 (_mesh, _node, _tpos))
        {
            return ;
        }

    /*-------------------------- _pop current 2-cell data */
        _pop_tri3(_mesh, _tpos, _itop, _kind);
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type _pop_tri3 (
        mesh_type&_mesh ,
        IPTR_TYPE _tpos ,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _epos = -1 ;

        __unreferenced(_kind);      // maybe, per TOPODIM

        if (_mesh.
            _llT3[_tpos].mark() < +0 )
            return ;                // is already removed

        if (_itop == -1)
        {
        _mesh.node(_mesh.tri3(_tpos).node(0)).
                ndeg(TRIA3_tag) -= 1 ;
        _mesh.node(_mesh.tri3(_tpos).node(1)).
                ndeg(TRIA3_tag) -= 1 ;
        _mesh.node(_mesh.tri3(_tpos).node(2)).
                ndeg(TRIA3_tag) -= 1 ;
        }

        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
        if constexpr ( TOPO_DIMS > 2 )
            _pop_list(_mesh. _aaT3 ,
                _tpos,
            conn_pair(_itop, _kind)) ;
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            _mesh.
            _llT3[_tpos].self() = +0 ;
        }

        if (null_tri3(_mesh, _tpos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */

        IPTR_TYPE _same ;
        _mesh._mmT3._pop(_tpos, _same) ;

        IPTR_TYPE _node [ 3] ;
        _node[0] =
            _mesh._llT3[_tpos].node(0) ;
        _node[1] =
            _mesh._llT3[_tpos].node(1) ;
        _node[2] =
            _mesh._llT3[_tpos].node(2) ;

        _put_tri3(_mesh, _tpos);

        for (_epos = 3; _epos-- != 0; )
        {
            IPTR_TYPE _enod [3];
            mesh_type::tri3_type::
                face_node(_enod, _epos, 2, 1) ;

            _enod[0] =
                _node[_enod [0]] ;
            _enod[1] =
                _node[_enod [1]] ;

            _pop_edge(_mesh,
                _enod,_tpos, TRIA3_tag) ;
        }
        }
    }

    /*
    --------------------------------------------------------
     * _POP-QUAD: delete 2-quad from complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call void_type _pop_quad (
        mesh_type&_mesh ,
        IPTR_TYPE const*_nptr,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _qpos = -1 ;

        IPTR_TYPE  _node[4];
        _node[0] = _nptr[0];
        _node[1] = _nptr[1];
        _node[2] = _nptr[2];
        _node[3] = _nptr[3];

    /*-------------------------- find current 2-cell pos. */
        if ( !find_quad (_mesh, _node, _qpos))
        {
            return ;
        }

    /*-------------------------- _pop current 2-cell data */
        _pop_quad(_mesh, _qpos, _itop, _kind);
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type _pop_quad (
        mesh_type&_mesh ,
        IPTR_TYPE _qpos ,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _epos = -1 ;

        __unreferenced(_kind);      // maybe, per TOPODIM

        if (_mesh.
            _llQ4[_qpos].mark() < +0 )
            return ;                // is already removed

        if (_itop == -1)
        {
        _mesh.node(_mesh.quad(_qpos).node(0)).
                ndeg(QUAD4_tag) -= 1 ;
        _mesh.node(_mesh.quad(_qpos).node(1)).
                ndeg(QUAD4_tag) -= 1 ;
        _mesh.node(_mesh.quad(_qpos).node(2)).
                ndeg(QUAD4_tag) -= 1 ;
        _mesh.node(_mesh.quad(_qpos).node(3)).
                ndeg(QUAD4_tag) -= 1 ;
        }

        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
        if constexpr ( TOPO_DIMS > 2 )
            _pop_list(_mesh. _aaQ4 ,
                _qpos,
            conn_pair(_itop, _kind)) ;
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            _mesh.
            _llQ4[_qpos].self() = +0 ;
        }

        if (null_quad(_mesh, _qpos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */

        IPTR_TYPE _same ;
        _mesh._mmQ4._pop(_qpos, _same) ;

        IPTR_TYPE _node [ 4] ;
        _node[0] =
            _mesh._llQ4[_qpos].node(0) ;
        _node[1] =
            _mesh._llQ4[_qpos].node(1) ;
        _node[2] =
            _mesh._llQ4[_qpos].node(2) ;
        _node[3] =
            _mesh._llQ4[_qpos].node(3) ;

        _put_quad (_mesh, _qpos);

        for (_epos = 4; _epos-- != 0; )
        {
            IPTR_TYPE _enod [4];
            mesh_type::quad_type::
                face_node(_enod, _epos, 2, 1) ;

            _enod[0] =
                _node[_enod [0]] ;
            _enod[1] =
                _node[_enod [1]] ;

            _pop_edge(_mesh,
                _enod,_qpos, QUAD4_tag) ;
        }
        }
    }

    /*
    --------------------------------------------------------
     * _POP-TRIA: delete 3-tria from complex.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call void_type _pop_tri4 (
        mesh_type&_mesh ,
        IPTR_TYPE const*_nptr,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _tpos = -1 ;

        IPTR_TYPE  _node[4];
        _node[0] = _nptr[0];
        _node[1] = _nptr[1];
        _node[2] = _nptr[2];
        _node[3] = _nptr[3];

    /*-------------------------- find current 3-cell pos. */
        if ( !find_tri4 (_mesh, _node, _tpos))
        {
            return ;
        }

    /*-------------------------- _pop current 3-cell data */
        _pop_tri4(_mesh, _tpos, _itop, _kind);
    }

    template <
        typename  mesh_type
             >
    __normal_call void_type _pop_tri4 (
        mesh_type&_mesh ,
        IPTR_TYPE _tpos ,
        IPTR_TYPE _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        IPTR_TYPE _fpos = -1 ;

        __unreferenced(_kind);      // maybe, per TOPODIM

        if (_mesh.
            _llT4[_tpos].mark() < +0 )
            return ;                // is already removed

        if (_itop == -1)
        {
        _mesh.node(_mesh.tri4(_tpos).node(0)).
                ndeg(TRIA4_tag) -= 1 ;
        _mesh.node(_mesh.tri4(_tpos).node(1)).
                ndeg(TRIA4_tag) -= 1 ;
        _mesh.node(_mesh.tri4(_tpos).node(2)).
                ndeg(TRIA4_tag) -= 1 ;
        _mesh.node(_mesh.tri4(_tpos).node(3)).
                ndeg(TRIA4_tag) -= 1 ;
        }

        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
        if constexpr ( TOPO_DIMS > 3 )
            _pop_list(_mesh. _aaT4 ,
                _tpos,
            conn_pair(_itop, _kind)) ;
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            _mesh.
            _llT4[_tpos].self() = +0 ;
        }

        if (null_tri4(_mesh, _tpos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */

        IPTR_TYPE _same ;
        _mesh._mmT4._pop(_tpos, _same) ;

        IPTR_TYPE _node [ 4] ;
        _node[0] =
            _mesh._llT4[_tpos].node(0) ;
        _node[1] =
            _mesh._llT4[_tpos].node(1) ;
        _node[2] =
            _mesh._llT4[_tpos].node(2) ;
        _node[3] =
            _mesh._llT4[_tpos].node(3) ;

        _put_tri4 (_mesh, _tpos) ;

        for (_fpos = 4; _fpos-- != 0; )
        {
            IPTR_TYPE _fnod [4];
            mesh_type::tri4_type::
                face_node(_fnod, _fpos, 3, 2) ;

            _fnod[0] =
                _node[_fnod [0]] ;
            _fnod[1] =
                _node[_fnod [1]] ;
            _fnod[2] =
                _node[_fnod [2]] ;

            _pop_tri3(_mesh,
                _fnod,_tpos, TRIA4_tag) ;
        }
        }
    }

    /*
    --------------------------------------------------------
     * FIND-NODE: return index of assoc. 0-node.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call bool_type find_node (
        mesh_type&_mesh ,
        IPTR_TYPE const*_nptr ,
        IPTR_TYPE&_npos
        )
    {
    /*-------------------------- find current 0-node pos. */
        IPTR_TYPE _node =
            (IPTR_TYPE) _nptr[0] ;

        if (_mesh._llN1.bounds(_node) &&
            _mesh._llN1[_node].mark() >= +0)
        {
    /*------------------------------- found matching node */
            _npos =_node ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-EDGE: return index of assoc. 1-edge.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_edge (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_epos
        )
    {
    /*-------------------------- find current 2-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_edge(_mesh);

       _mesh.
       _llE2[_ipos].node(0) = _node[0];
       _mesh.
       _llE2[_ipos].node(1) = _node[1];

        typename
        mesh_type::
        edge_maps::_write_it _same ;
        if (_mesh.
            _mmE2.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_edge(_mesh, _ipos);

            _epos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_edge(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-TRIA: return index of assoc. 2-tria.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_tri3 (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_tpos
        )
    {
    /*-------------------------- find current 2-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_tri3(_mesh);

       _mesh.
       _llT3[_ipos].node(0) = _node[0];
       _mesh.
       _llT3[_ipos].node(1) = _node[1];
       _mesh.
       _llT3[_ipos].node(2) = _node[2];

        typename
        mesh_type::
        tri3_maps::_write_it _same ;
        if (_mesh.
            _mmT3.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_tri3(_mesh, _ipos);

            _tpos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_tri3(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-QUAD: return index of assoc. 2-quad.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_quad (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_qpos
        )
    {
    /*-------------------------- find current 2-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_quad(_mesh);

       _mesh.
       _llQ4[_ipos].node(0) = _node[0];
       _mesh.
       _llQ4[_ipos].node(1) = _node[1];
       _mesh.
       _llQ4[_ipos].node(2) = _node[2];
       _mesh.
       _llQ4[_ipos].node(3) = _node[3];

        typename
        mesh_type::
        quad_maps::_write_it _same ;
        if (_mesh.
            _mmQ4.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_quad(_mesh, _ipos);

            _qpos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_quad(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-TRIA: return index of assoc. 3-tria.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_tri4 (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_tpos
        )
    {
    /*-------------------------- find current 3-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_tri4(_mesh);

       _mesh.
       _llT4[_ipos].node(0) = _node[0];
       _mesh.
       _llT4[_ipos].node(1) = _node[1];
       _mesh.
       _llT4[_ipos].node(2) = _node[2];
       _mesh.
       _llT4[_ipos].node(3) = _node[3];

        typename
        mesh_type::
        tri4_maps::_write_it _same ;
        if (_mesh.
            _mmT4.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_tri4(_mesh, _ipos);

            _tpos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_tri4(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-HEXA: return index of assoc. 3-hexa.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_hexa (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_hpos
        )
    {
    /*-------------------------- find current 3-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_hexa(_mesh);

       _mesh.
       _llH8[_ipos].node(0) = _node[0];
       _mesh.
       _llH8[_ipos].node(1) = _node[1];
       _mesh.
       _llH8[_ipos].node(2) = _node[2];
       _mesh.
       _llH8[_ipos].node(3) = _node[3];
       _mesh.
       _llH8[_ipos].node(4) = _node[4];
       _mesh.
       _llH8[_ipos].node(5) = _node[5];
       _mesh.
       _llH8[_ipos].node(6) = _node[6];
       _mesh.
       _llH8[_ipos].node(7) = _node[7];

        typename
        mesh_type::
        hexa_maps::_write_it _same ;
        if (_mesh.
            _mmH8.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_hexa(_mesh, _ipos);

            _hpos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_hexa(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-WEDG: return index of assoc. 3-wedg.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_wedg (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_wpos
        )
    {
    /*-------------------------- find current 3-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_wedg(_mesh);

       _mesh.
       _llW6[_ipos].node(0) = _node[0];
       _mesh.
       _llW6[_ipos].node(1) = _node[1];
       _mesh.
       _llW6[_ipos].node(2) = _node[2];
       _mesh.
       _llW6[_ipos].node(3) = _node[3];
       _mesh.
       _llW6[_ipos].node(4) = _node[4];
       _mesh.
       _llW6[_ipos].node(5) = _node[5];

        typename
        mesh_type::
        wedg_maps::_write_it _same ;
        if (_mesh.
            _mmW6.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_wedg(_mesh, _ipos);

            _wpos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_wedg(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-PYRA: return index of assoc. 3-pyra.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __normal_call bool_type find_pyra (
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        IPTR_TYPE&_ppos
        )
    {
    /*-------------------------- find current 3-cell pos. */
        IPTR_TYPE _ipos =
       (IPTR_TYPE) _get_pyra(_mesh);

       _mesh.
       _llP5[_ipos].node(0) = _node[0];
       _mesh.
       _llP5[_ipos].node(1) = _node[1];
       _mesh.
       _llP5[_ipos].node(2) = _node[2];
       _mesh.
       _llP5[_ipos].node(3) = _node[3];
       _mesh.
       _llP5[_ipos].node(4) = _node[4];

        typename
        mesh_type::
        pyra_maps::_write_it _same ;
        if (_mesh.
            _mmP5.find(_ipos,_same))
        {
    /*------------------------------- found matching cell */
            _put_pyra(_mesh, _ipos);

            _ppos=*_same ;

            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_pyra(_mesh, _ipos);

            return false ;
        }
    }

    /*
    --------------------------------------------------------
     * FIND-CELL: return index of assoc. k-cell.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type
             >
    __inline_call bool_type find_cell (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        IPTR_TYPE&_ipos ,
        tria_cell const &
        )
    {
        switch ( _kind )
        {
    /*-------------------------- find current k-cell pos. */
        case POINT_tag:
            {
            if constexpr (TOPO_DIMS >= +0)
            {
            if (find_node(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case EDGE2_tag:
            {
            if constexpr (TOPO_DIMS >= +1)
            {
            if (find_edge(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case TRIA3_tag:
            {
            if constexpr (TOPO_DIMS >= +2)
            {
            if (find_tri3(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case TRIA4_tag:
            {
            if constexpr (TOPO_DIMS >= +3)
            {
            if (find_tri4(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }
        }

        return  false ;
    }

    template <
        typename  mesh_type
             >
    __inline_call bool_type find_cell (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        IPTR_TYPE&_ipos ,
        mesh_cell const &
        )
    {
        switch ( _kind )
        {
    /*-------------------------- find current k-cell pos. */
        case POINT_tag:
            {
            if constexpr (TOPO_DIMS >= +0)
            {
            if (find_node(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case EDGE2_tag:
            {
            if constexpr (TOPO_DIMS >= +1)
            {
            if (find_edge(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case TRIA3_tag:
            {
            if constexpr (TOPO_DIMS >= +2)
            {
            if (find_tri3(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case QUAD4_tag:
            {
            if constexpr (TOPO_DIMS >= +2)
            {
            if (find_quad(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case TRIA4_tag:
            {
            if constexpr (TOPO_DIMS >= +3)
            {
            if (find_tri4(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case HEXA8_tag:
            {
            if constexpr (TOPO_DIMS >= +3)
            {
            if (find_hexa(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case WEDG6_tag:
            {
            if constexpr (TOPO_DIMS >= +3)
            {
            if (find_wedg(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }

        case PYRA5_tag:
            {
            if constexpr (TOPO_DIMS >= +3)
            {
            if (find_pyra(_mesh, _node, _ipos))
                return  true ;
            else
                return false ;
            }
            else
                return false ;
            }
        }

        return  false ;
    }

    /*
    --------------------------------------------------------
     * CONNECT-1: form this-to-1-cell adj. list.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  connector
             >
    __inline_call void_type connect_1 (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        connector&_conn ,
        tria_cell _type
        )
    {
    /*-------------------------- get adj. for k-cell pos. */
        IPTR_TYPE _ipos ;
        if (find_cell(_mesh, _node, _kind, _ipos,
                      _type) )
        {
            connect_1(_mesh, _ipos, _kind, _conn,
                      _type) ;
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __inline_call void_type connect_1 (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        connector&_conn ,
        mesh_cell _type
        )
    {
    /*-------------------------- get adj. for k-cell pos. */
        IPTR_TYPE _ipos ;
        if (find_cell(_mesh, _node, _kind, _ipos,
                      _type) )
        {
            connect_1(_mesh, _ipos, _kind, _conn,
                      _type) ;
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __normal_call void_type connect_1 (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE _ipos ,
        char_type _kind ,
        connector&_conn ,
        tria_cell const &
        )
    {
        _mesh._tmp1.set_count(  +0 ) ;

        switch ( _kind )
        {
    /*-------------------------- init. "source" adj. list */
        case POINT_tag:
            {
            _mesh._tmp1.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }
        }

    /*-------------------------- find set of adj. 1-cells */
        for (auto _iter  = _mesh._tmp1.head();
                  _iter != _mesh._tmp1.tend();
                ++_iter  )
        {
            if (_iter->_kind == POINT_tag)
            for (auto _iadj  =
                _mesh._aaN1.head(_iter->_cell);
                      _iadj !=
                _mesh._aaN1.tend(_iter->_cell);
                    ++_iadj  )
            {
            //  1-adj. list is disjoint, guaranteed
            //  if(_iadj->_kind == EDGE2_tag &&
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() == 0)
            //  {
                   _conn.push_tail (*_iadj) ;
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() += 1;
            //  }
            }
        }

    /*-------------------------- flip d-face marker lists */
        //  done!
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __normal_call void_type connect_1 (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE _ipos ,
        char_type _kind ,
        connector&_conn ,
        mesh_cell const &
        )
    {
        _mesh._tmp1.set_count(  +0 ) ;

        switch ( _kind )
        {
    /*-------------------------- init. "source" adj. list */
        case POINT_tag:
            {
            _mesh._tmp1.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }
        }

    /*-------------------------- find set of adj. 1-cells */
        for (auto _iter  = _mesh._tmp1.head();
                  _iter != _mesh._tmp1.tend();
                ++_iter  )
        {
            if (_iter->_kind == POINT_tag)
            for (auto _iadj  =
                _mesh._aaN1.head(_iter->_cell);
                      _iadj !=
                _mesh._aaN1.tend(_iter->_cell);
                    ++_iadj  )
            {
            //  1-adj. list is disjoint, guaranteed
            //  if(_iadj->_kind == EDGE2_tag &&
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() == 0)
            //  {
                   _conn.push_tail (*_iadj) ;
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() += 1;
            //  }
            }
        }

    /*-------------------------- flip d-face marker lists */
        //  done!
    }

    /*
    --------------------------------------------------------
     * CONNECT-2: form this-to-2-cell adj. list.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  connector
             >
    __inline_call void_type connect_2 (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        connector&_conn ,
        tria_cell _type
        )
    {
    /*-------------------------- get adj. for k-cell pos. */
        IPTR_TYPE _ipos ;
        if (find_cell(_mesh, _node, _kind, _ipos,
                      _type) )
        {
            connect_2(_mesh, _ipos, _kind, _conn,
                      _type) ;
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __inline_call void_type connect_2 (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        connector&_conn ,
        mesh_cell _type
        )
    {
    /*-------------------------- get adj. for k-cell pos. */
        IPTR_TYPE _ipos ;
        if (find_cell(_mesh, _node, _kind, _ipos,
                      _type) )
        {
            connect_2(_mesh, _ipos, _kind, _conn,
                      _type) ;
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __normal_call void_type connect_2 (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE _ipos ,
        char_type _kind ,
        connector&_conn ,
        tria_cell const &
        )
    {
        _mesh._tmp1.set_count(  +0 ) ;
        _mesh._tmp2.set_count(  +0 ) ;

         auto _ioff = _conn.count () ;

        switch ( _kind )
        {
    /*-------------------------- init. "source" adj. list */
        case POINT_tag:
            {
            _mesh._tmp1.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }

        case EDGE2_tag:
            {
            _mesh._tmp2.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }
        }

    /*-------------------------- find set of adj. 1-cells */
        for (auto _iter  = _mesh._tmp1.head();
                  _iter != _mesh._tmp1.tend();
                ++_iter  )
        {
            if (_iter->_kind == POINT_tag)
            for (auto _iadj  =
                _mesh._aaN1.head(_iter->_cell);
                      _iadj !=
                _mesh._aaN1.tend(_iter->_cell);
                    ++_iadj  )
            {
            //  1-adj. list is disjoint, guaranteed
            //  if(_iadj->_kind == EDGE2_tag &&
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() == 0)
            //  {
                   _mesh.
                   _tmp2.push_tail (*_iadj);
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() += 1;
            //  }
            }
        }

    /*-------------------------- find set of adj. 2-cells */
        for (auto _iter  = _mesh._tmp2.head();
                  _iter != _mesh._tmp2.tend();
                ++_iter  )
        {
            if (_iter->_kind == EDGE2_tag)
            for (auto _iadj  =
                _mesh._aaE2.head(_iter->_cell);
                      _iadj !=
                _mesh._aaE2.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == TRIA3_tag &&
                   _mesh. _llT3[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh._llT3[
                   _iadj->_cell].mark() += 1;
                }
            }
        }

    /*-------------------------- flip d-face marker lists */
        for (auto _iter  = _conn.head()+_ioff ;
                  _iter != _conn.tend() ;
                ++_iter  )
        {
            if(_iter->_kind == TRIA3_tag)
            {
               _mesh. _llT3[
                    _iter->_cell].mark() = 0 ;
            }
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __normal_call void_type connect_2 (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE _ipos ,
        char_type _kind ,
        connector&_conn ,
        mesh_cell const &
        )
    {
        _mesh._tmp1.set_count(  +0 ) ;
        _mesh._tmp2.set_count(  +0 ) ;

         auto _ioff = _conn.count () ;

        switch ( _kind )
        {
    /*-------------------------- init. "source" adj. list */
        case POINT_tag:
            {
            _mesh._tmp1.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }

        case EDGE2_tag:
            {
            _mesh._tmp2.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }
        }

    /*-------------------------- find set of adj. 1-cells */
        for (auto _iter  = _mesh._tmp1.head();
                  _iter != _mesh._tmp1.tend();
                ++_iter  )
        {
            if (_iter->_kind == POINT_tag)
            for (auto _iadj  =
                _mesh._aaN1.head(_iter->_cell);
                      _iadj !=
                _mesh._aaN1.tend(_iter->_cell);
                    ++_iadj  )
            {
            //  1-adj. list is disjoint, guaranteed
            //  if(_iadj->_kind == EDGE2_tag &&
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() == 0)
            //  {
                   _mesh.
                   _tmp2.push_tail (*_iadj);
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() += 1;
            //  }
            }
        }

    /*-------------------------- find set of adj. 2-cells */
        for (auto _iter  = _mesh._tmp2.head();
                  _iter != _mesh._tmp2.tend();
                ++_iter  )
        {
            if (_iter->_kind == EDGE2_tag)
            for (auto _iadj  =
                _mesh._aaE2.head(_iter->_cell);
                      _iadj !=
                _mesh._aaE2.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == TRIA3_tag &&
                   _mesh. _llT3[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh._llT3[
                   _iadj->_cell].mark() += 1;
                }
                else
                if(_iadj->_kind == QUAD4_tag &&
                   _mesh. _llQ4[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llQ4[
                   _iadj->_cell].mark() += 1;
                }
            }
        }

    /*-------------------------- flip d-face marker lists */
        for (auto _iter  = _conn.head()+_ioff ;
                  _iter != _conn.tend() ;
                ++_iter  )
        {
            if(_iter->_kind == TRIA3_tag)
            {
               _mesh. _llT3[
                    _iter->_cell].mark() = 0 ;
            }
            else
            if(_iter->_kind == QUAD4_tag)
            {
               _mesh. _llQ4[
                    _iter->_cell].mark() = 0 ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * CONNECT-3: form this-to-3-cell adj. list.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  connector
             >
    __inline_call void_type connect_3 (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        connector&_conn ,
        tria_cell _type
        )
    {
    /*-------------------------- get adj. for k-cell pos. */
        IPTR_TYPE _ipos ;
        if (find_cell(_mesh, _node, _kind, _ipos,
                      _type) )
        {
            connect_3(_mesh, _ipos, _kind, _conn,
                      _type) ;
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __inline_call void_type connect_3 (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE const*_node ,
        char_type _kind ,
        connector&_conn ,
        mesh_cell _type
        )
    {
    /*-------------------------- get adj. for k-cell pos. */
        IPTR_TYPE _ipos ;
        if (find_cell(_mesh, _node, _kind, _ipos,
                      _type) )
        {
            connect_3(_mesh, _ipos, _kind, _conn,
                      _type) ;
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __normal_call void_type connect_3 (   // tria-complex
        mesh_type&_mesh ,
        IPTR_TYPE _ipos ,
        char_type _kind ,
        connector&_conn ,
        tria_cell const &
        )
    {
        _mesh._tmp1.set_count(  +0 ) ;
        _mesh._tmp2.set_count(  +0 ) ;
        _mesh._tmp3.set_count(  +0 ) ;

         auto _ioff = _conn.count () ;

        switch ( _kind )
        {
    /*-------------------------- init. "source" adj. list */
        case POINT_tag:
            {
            _mesh._tmp1.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }

        case EDGE2_tag:
            {
            _mesh._tmp2.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }

        case TRIA3_tag:
            {
            _mesh._tmp3.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }
        }

    /*-------------------------- find set of adj. 1-cells */
        for (auto _iter  = _mesh._tmp1.head();
                  _iter != _mesh._tmp1.tend();
                ++_iter  )
        {
            if (_iter->_kind == POINT_tag)
            for (auto _iadj  =
                _mesh._aaN1.head(_iter->_cell);
                      _iadj !=
                _mesh._aaN1.tend(_iter->_cell);
                    ++_iadj  )
            {
            //  1-adj. list is disjoint, guaranteed
            //  if(_iadj->_kind == EDGE2_tag &&
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() == 0)
            //  {
                   _mesh.
                   _tmp2.push_tail (*_iadj);
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() += 1;
            //  }
            }
        }

    /*-------------------------- find set of adj. 2-cells */
        for (auto _iter  = _mesh._tmp2.head();
                  _iter != _mesh._tmp2.tend();
                ++_iter  )
        {
            if (_iter->_kind == EDGE2_tag)
            for (auto _iadj  =
                _mesh._aaE2.head(_iter->_cell);
                      _iadj !=
                _mesh._aaE2.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == TRIA3_tag &&
                   _mesh. _llT3[
                   _iadj->_cell].mark() == 0)
                {
                   _mesh.
                   _tmp3.push_tail (*_iadj);
                   _mesh. _llT3[
                   _iadj->_cell].mark() += 1;
                }
            }
        }

    /*-------------------------- find set of adj. 3-cells */
        for (auto _iter  = _mesh._tmp3.head();
                  _iter != _mesh._tmp3.tend();
                ++_iter  )
        {
            if (_iter->_kind == TRIA3_tag) // from TRIA-3
            for (auto _iadj  =
                _mesh._aaT3.head(_iter->_cell);
                      _iadj !=
                _mesh._aaT3.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == TRIA4_tag &&
                   _mesh. _llT4[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llT4[
                   _iadj->_cell].mark() += 1;
                }
            }
        }

    /*-------------------------- flip d-face marker lists */
        for (auto _iter  = _conn.head()+_ioff ;
                  _iter != _conn.tend() ;
                ++_iter  )
        {
            if(_iter->_kind == TRIA4_tag)
            {
               _mesh. _llT4[
                    _iter->_cell].mark() = 0 ;
            }
        }

        for (auto _iter  = _mesh._tmp3.head();
                  _iter != _mesh._tmp3.tend();
                ++_iter  )
        {
            if(_iter->_kind == TRIA3_tag)
            {
               _mesh. _llT3[
                    _iter->_cell].mark() = 0 ;
            }
        }
    }

    template <
        typename  mesh_type ,
        typename  connector
             >
    __normal_call void_type connect_3 (   // mesh-complex
        mesh_type&_mesh ,
        IPTR_TYPE _ipos ,
        char_type _kind ,
        connector&_conn ,
        mesh_cell const &
        )
    {
        _mesh._tmp1.set_count(  +0 ) ;
        _mesh._tmp2.set_count(  +0 ) ;
        _mesh._tmp3.set_count(  +0 ) ;

         auto _ioff = _conn.count () ;

        switch ( _kind )
        {
    /*-------------------------- init. "source" adj. list */
        case POINT_tag:
            {
            _mesh._tmp1.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }

        case EDGE2_tag:
            {
            _mesh._tmp2.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }

        case TRIA3_tag:     // fall-through
        case QUAD4_tag:
            {
            _mesh._tmp3.push_tail(
                conn_pair (_ipos, _kind)) ;
            break ;
            }
        }

    /*-------------------------- find set of adj. 1-cells */
        for (auto _iter  = _mesh._tmp1.head();
                  _iter != _mesh._tmp1.tend();
                ++_iter  )
        {
            if (_iter->_kind == POINT_tag)
            for (auto _iadj  =
                _mesh._aaN1.head(_iter->_cell);
                      _iadj !=
                _mesh._aaN1.tend(_iter->_cell);
                    ++_iadj  )
            {
            //  1-adj. list is disjoint, guaranteed
            //  if(_iadj->_kind == EDGE2_tag &&
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() == 0)
            //  {
                   _mesh.
                   _tmp2.push_tail (*_iadj);
            //     _mesh. _llE2[
            //     _iadj->_cell].mark() += 1;
            //  }
            }
        }

    /*-------------------------- find set of adj. 2-cells */
        for (auto _iter  = _mesh._tmp2.head();
                  _iter != _mesh._tmp2.tend();
                ++_iter  )
        {
            if (_iter->_kind == EDGE2_tag)
            for (auto _iadj  =
                _mesh._aaE2.head(_iter->_cell);
                      _iadj !=
                _mesh._aaE2.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == TRIA3_tag &&
                   _mesh. _llT3[
                   _iadj->_cell].mark() == 0)
                {
                   _mesh.
                   _tmp3.push_tail (*_iadj);
                   _mesh. _llT3[
                   _iadj->_cell].mark() += 1;
                }
                else
                if(_iadj->_kind == QUAD4_tag &&
                   _mesh. _llQ4[
                   _iadj->_cell].mark() == 0)
                {
                   _mesh.
                   _tmp3.push_tail (*_iadj);
                   _mesh. _llQ4[
                   _iadj->_cell].mark() += 1;
                }
            }
        }

    /*-------------------------- find set of adj. 3-cells */
        for (auto _iter  = _mesh._tmp3.head();
                  _iter != _mesh._tmp3.tend();
                ++_iter  )
        {
            if (_iter->_kind == TRIA3_tag) // from TRIA-3
            for (auto _iadj  =
                _mesh._aaT3.head(_iter->_cell);
                      _iadj !=
                _mesh._aaT3.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == TRIA4_tag &&
                   _mesh. _llT4[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llT4[
                   _iadj->_cell].mark() += 1;
                }
                else
                if(_iadj->_kind == WEDG6_tag &&
                   _mesh. _llW6[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llW6[
                   _iadj->_cell].mark() += 1;
                }
                else
                if(_iadj->_kind == PYRA5_tag &&
                   _mesh. _llP5[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llP5[
                   _iadj->_cell].mark() += 1;
                }
            }

            if (_iter->_kind == QUAD4_tag) // from QUAD-4
            for (auto _iadj  =
                _mesh._aaQ4.head(_iter->_cell);
                      _iadj !=
                _mesh._aaQ4.tend(_iter->_cell);
                    ++_iadj  )
            {
                if(_iadj->_kind == HEXA8_tag &&
                   _mesh. _llH8[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llH8[
                   _iadj->_cell].mark() += 1;
                }
                else
                if(_iadj->_kind == WEDG6_tag &&
                   _mesh. _llW6[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llW6[
                   _iadj->_cell].mark() += 1;
                }
                else
                if(_iadj->_kind == PYRA5_tag &&
                   _mesh. _llP5[
                   _iadj->_cell].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                   _mesh. _llP5[
                   _iadj->_cell].mark() += 1;
                }
            }
        }

    /*-------------------------- flip d-face marker lists */
        for (auto _iter  = _conn.head()+_ioff ;
                  _iter != _conn.tend() ;
                ++_iter  )
        {
            if(_iter->_kind == TRIA4_tag)
            {
               _mesh. _llT4[
                    _iter->_cell].mark() = 0 ;
            }
            else
            if(_iter->_kind == HEXA8_tag)
            {
               _mesh. _llH8[
                    _iter->_cell].mark() = 0 ;
            }
            else
            if(_iter->_kind == WEDG6_tag)
            {
               _mesh. _llW6[
                    _iter->_cell].mark() = 0 ;
            }
            else
            if(_iter->_kind == PYRA5_tag)
            {
               _mesh. _llP5[
                    _iter->_cell].mark() = 0 ;
            }
        }

        for (auto _iter  = _mesh._tmp3.head();
                  _iter != _mesh._tmp3.tend();
                ++_iter  )
        {
            if(_iter->_kind == TRIA3_tag)
            {
               _mesh. _llT3[
                    _iter->_cell].mark() = 0 ;
            }
            else
            if(_iter->_kind == QUAD4_tag)
            {
               _mesh. _llQ4[
                    _iter->_cell].mark() = 0 ;
            }
        }
    }

#   undef  TOPO_DIMS

#   undef  IPTR_TYPE
#   undef  INDX_TYPE
#   undef  REAL_TYPE

#   undef  NODE_TYPE
#   undef  EDGE_TYPE
#   undef  TRI3_TYPE
#   undef  QUAD_TYPE
#   undef  TRI4_TYPE
#   undef  HEXA_TYPE
#   undef  WEDG_TYPE
#   undef  PYRA_TYPE

    }

#   endif  //__BASE_COMPLEX_k__



