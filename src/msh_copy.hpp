
    /*
    --------------------------------------------------------
     * MSH-COPY: copy MESH data into a tria-complex.
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
     * Last updated: 18 Jul., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_COPY__
#   define __MSH_COPY__

    template <
    typename      iptr_list
             >
    __normal_call void_type sort_node (
        mesh_data &_mesh ,
        iptr_list &_perm
        )
    {
        iptr_type static constexpr SIZE = +96 ;

        if (_mesh._ndim == +2 &&
            _mesh._kind ==
            jmsh_kind::euclidean_mesh)
        {
    /*---------------------------------- sort 2-dim. mesh */
        typedef geom_tree::aabb_item_node_k <
            real_type,
            iptr_type, 2>   tree_item ;

        typedef geom_tree::aabb_tree<
            tree_item, 2>   tree_type ;

        tree_type _tree ;

        containers::array<tree_item>  _item ;
        containers::array<  size_t >  _itmp ;

        iptr_type _ipos  = +0 ;
        for (auto _iter  = _mesh.
        _euclidean_rdel_2d._tria._nset.head() ;
                  _iter != _mesh.
        _euclidean_rdel_2d._tria._nset.tend() ;
                ++_iter, ++_ipos)
        {
            if (_iter->mark() >= +0 &&
                _iter->fdim() <= +2 )  // is not isolated
            {
                _item.push_tail();
                _item.tail()->ipos () = _ipos ;
                _item.tail()->pval(0) =
                             _iter->pval(0) ;
                _item.tail()->pval(1) =
                             _iter->pval(1) ;
            }
        }

        _tree.load(_item.head(), _item.tend() ,
            SIZE ,  // pop.=#nodes per leaf
        _tree.median_split ,
            1.0) ;  // set LONG=1, only div. by pop.

        _tree.dcmp(_itmp, _perm) ;
        }
        else
        if (_mesh._ndim == +3 &&
            _mesh._kind ==
            jmsh_kind::euclidean_mesh)
        {
    /*---------------------------------- sort 3-dim. mesh */
        typedef geom_tree::aabb_item_node_k <
            real_type,
            iptr_type, 3>   tree_item ;

        typedef geom_tree::aabb_tree<
            tree_item, 3>   tree_type ;

        tree_type _tree ;

        containers::array<tree_item>  _item ;
        containers::array<  size_t >  _itmp ;

        iptr_type _ipos  = +0 ;
        for (auto _iter  = _mesh.
        _euclidean_rdel_3d._tria._nset.head() ;
                  _iter != _mesh.
        _euclidean_rdel_3d._tria._nset.tend() ;
                ++_iter, ++_ipos)
        {
            if (_iter->mark() >= +0 &&
                _iter->fdim() <= +3 )  // is not isolated
            {
                _item.push_tail();
                _item.tail()->ipos () = _ipos ;
                _item.tail()->pval(0) =
                             _iter->pval(0) ;
                _item.tail()->pval(1) =
                             _iter->pval(1) ;
                _item.tail()->pval(2) =
                             _iter->pval(2) ;
            }
        }

        _tree.load(_item.head(), _item.tend() ,
            SIZE ,  // pop.=#nodes per leaf
        _tree.median_split ,
            1.0) ;  // set LONG=1, only div. by pop.

        _tree.dcmp(_itmp, _perm) ;
        }
    }


    template <
    typename      cell_list ,
    typename      iptr_list ,
    typename      cell_perm
             >
    __normal_call void_type sort_cell (
        cell_list &_cell ,
        size_t     _nsiz ,
        iptr_list &_nmap ,
        cell_perm &_perm
        )
    {
        class sort_pair
        {
            public  :
    /*------------------------ tuple for cell re-ordering */
            typename cell_list::item_type *
                                _cell ;
            float               _cost ;

    /*------------------------ construct inline from src. */
            __inline_call sort_pair (
            typename cell_list::item_type *
                                _psrc ,
                float     const&_csrc
                )   :
                   _cell (_psrc),
                   _cost((float)_csrc) {}
            } ;

        class sort_pred
            {
            public  :
    /*------------------------ less-than op. for cost-tup */
            __inline_call
                bool_type operator () (
                sort_pair const&_idat ,
                sort_pair const&_jdat
                ) const
            { return _idat._cost < _jdat._cost ;
            }
            } ;

        typedef
        containers::array<sort_pair> sort_list ;

        sort_list _sort ;
        for (auto _iter  = _cell._lptr.head() ;
                  _iter != _cell._lptr.tend() ;
                ++_iter  )
        {
            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
                auto _cost =            // min. node indx
                    _nmap[_item->_data._node [0]] ;
                for (auto
                    _inod = _nsiz; _inod-- != 1; )
                {
                    _cost = std::min(_cost,
                _nmap[_item->_data._node[_inod]]) ;
                }

                float _cval = (float)_cost;

                _sort.push_tail(
                    sort_pair(_item, _cval));
            }
        }

        algorithms::qsort( _sort.head() ,
                           _sort.tend() ,
            sort_pred () ) ;

        for (auto _iter  = _sort.head() ;
                  _iter != _sort.tend() ;
                ++_iter  )
        {
            _perm.push_tail( _iter->_cell ) ;
        }
    }


    /*
    --------------------------------------------------------
     * COPY-MESH: copy rdel-complex to tria-complex.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type copy_mesh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_mesh
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;
        __unreferenced(_jlog) ;

        if (_mesh._ndim == +2 &&
            _mesh._kind ==
            jmsh_kind::euclidean_mesh)
        {
    /*---------------------------------- copy 2-dim. mesh */
        containers::array<iptr_type> _nmap, _perm;

        sort_node(_mesh, _perm) ;

        _nmap.set_count(_mesh.
        _euclidean_rdel_2d._tria._nset.count(),
        containers::loose_alloc, -1) ;

        _mesh._euclidean_mesh_2d._mesh.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _perm.head() ;
                  _iter != _perm.tend() ;
                ++_iter  )
        {
             auto _nptr  = _mesh.
        _euclidean_rdel_2d._tria._nset.head() + *_iter ;

            typename mesh_data::
                euclidean_mesh_2d::
            mesh_type::node_type _node ;
            _node.pval(0) = _nptr->pval(0) ;
            _node.pval(1) = _nptr->pval(1) ;

            _node.pval(2) = _nptr->pval(2) ;

            _node.hidx () = _nptr->idxh () ;

            _node.itag () = _nptr->part () ;
            _node.fdim () = _nptr->fdim () ;
            _node.feat () = _nptr->feat () ;

            _nmap[ *_iter ] =
            _mesh._euclidean_mesh_2d.
                _mesh.push_node (_node, false) ;
        }

        _mesh._euclidean_rdel_2d._tria.
            clear(containers::tight_alloc) ;

        containers::array<
            mesh_data::euclidean_rdel_2d::
                edge_list::item_type*> _eprm ;

        sort_cell(_mesh.
            _euclidean_rdel_2d._eset, +2, _nmap, _eprm);

        for (auto _iter  = _eprm.head() ;
                  _iter != _eprm.tend() ;
                ++_iter  )
        {
             auto _item = *_iter ;

            typename mesh_data::
                euclidean_mesh_2d::
            mesh_type::edge_type _face;
            _face.node(0) =
                _nmap[_item->_data._node[0]] ;
            _face.node(1) =
                _nmap[_item->_data._node[1]] ;

            _face.itag () =
                _item->_data._part;

            _mesh._euclidean_mesh_2d.
                _mesh.push_edge (_face, false) ;
        }

        _mesh._euclidean_rdel_2d._eset.
            clear(containers::tight_alloc) ;

        containers::array<
            mesh_data::euclidean_rdel_2d::
                tria_list::item_type*> _tprm ;

        sort_cell(_mesh.
            _euclidean_rdel_2d._tset, +3, _nmap, _tprm);

        for (auto _iter  = _tprm.head() ;
                  _iter != _tprm.tend() ;
                ++_iter  )
        {
             auto _item = *_iter ;

            typename mesh_data::
                euclidean_mesh_2d::
            mesh_type::tri3_type _face;
            _face.node(0) =
                _nmap[_item->_data._node[0]] ;
            _face.node(1) =
                _nmap[_item->_data._node[1]] ;
            _face.node(2) =
                _nmap[_item->_data._node[2]] ;

            _face.itag () =
                _item->_data._part;

            _mesh._euclidean_mesh_2d.
                _mesh.push_tri3 (_face, false) ;
        }

        _mesh._euclidean_rdel_2d._tset.
            clear(containers::tight_alloc) ;

        _mesh._euclidean_rdel_2d.
            clear(containers::tight_alloc) ;

        _mesh._euclidean_mesh_2d.
                        _mesh.make_link () ;

        }
        else
        if (_mesh._ndim == +3 &&
            _mesh._kind ==
            jmsh_kind::euclidean_mesh)
        {
    /*---------------------------------- copy 3-dim. mesh */
        containers::array<iptr_type> _nmap, _perm;

        sort_node(_mesh, _perm) ;

        _nmap.set_count(_mesh.
        _euclidean_rdel_3d._tria._nset.count(),
        containers::loose_alloc, -1) ;

        _mesh._euclidean_mesh_3d._mesh.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _perm.head() ;
                  _iter != _perm.tend() ;
                ++_iter  )
        {
             auto _nptr  = _mesh.
        _euclidean_rdel_3d._tria._nset.head() + *_iter ;

            typename mesh_data::
                euclidean_mesh_3d::
            mesh_type::node_type _node ;
            _node.pval(0) = _nptr->pval(0) ;
            _node.pval(1) = _nptr->pval(1) ;
            _node.pval(2) = _nptr->pval(2) ;

            _node.pval(3) = _nptr->pval(3) ;

            _node.hidx () = _nptr->idxh () ;

            _node.itag () = _nptr->part () ;
            _node.fdim () = _nptr->fdim () ;
            _node.feat () = _nptr->feat () ;

            _nmap[ *_iter ] =
            _mesh._euclidean_mesh_3d.
                _mesh.push_node (_node, false) ;
        }

        _mesh._euclidean_rdel_3d._tria.
            clear(containers::tight_alloc) ;

        containers::array<
            mesh_data::euclidean_rdel_3d::
                edge_list::item_type*> _eprm ;

        sort_cell(_mesh.
            _euclidean_rdel_3d._eset, +2, _nmap, _eprm);

        for (auto _iter  = _eprm.head() ;
                  _iter != _eprm.tend() ;
                ++_iter  )
        {
             auto _item = *_iter ;

            typename mesh_data::
                euclidean_mesh_3d::
            mesh_type::edge_type _face;
            _face.node(0) =
                _nmap[_item->_data._node[0]] ;
            _face.node(1) =
                _nmap[_item->_data._node[1]] ;

            _face.itag () =
                _item->_data._part;

            _mesh._euclidean_mesh_3d.
                _mesh.push_edge (_face, false) ;
        }

        _mesh._euclidean_rdel_3d._eset.
            clear(containers::tight_alloc) ;

        containers::array<
            mesh_data::euclidean_rdel_3d::
                face_list::item_type*> _fprm ;

        sort_cell(_mesh.
            _euclidean_rdel_3d._fset, +3, _nmap, _fprm);

        for (auto _iter  = _fprm.head() ;
                  _iter != _fprm.tend() ;
                ++_iter  )
        {
             auto _item = *_iter ;

            typename mesh_data::
                euclidean_mesh_3d::
            mesh_type::tri3_type _face;
            _face.node(0) =
                _nmap[_item->_data._node[0]] ;
            _face.node(1) =
                _nmap[_item->_data._node[1]] ;
            _face.node(2) =
                _nmap[_item->_data._node[2]] ;

            _face.itag () =
                _item->_data._part;

            _mesh._euclidean_mesh_3d.
                _mesh.push_tri3 (_face, false) ;
        }

        _mesh._euclidean_rdel_3d._fset.
            clear(containers::tight_alloc) ;

        containers::array<
            mesh_data::euclidean_rdel_3d::
                tria_list::item_type*> _tprm ;

        sort_cell(_mesh.
            _euclidean_rdel_3d._tset, +4, _nmap, _tprm);

        for (auto _iter  = _tprm.head() ;
                  _iter != _tprm.tend() ;
                ++_iter  )
        {
             auto _item = *_iter ;

            typename mesh_data::
                euclidean_mesh_3d::
            mesh_type::tri4_type _face;
            _face.node(0) =
                _nmap[_item->_data._node[0]] ;
            _face.node(1) =
                _nmap[_item->_data._node[1]] ;
            _face.node(2) =
                _nmap[_item->_data._node[2]] ;
            _face.node(3) =
                _nmap[_item->_data._node[3]] ;

            _face.itag () =
                _item->_data._part;

            _mesh._euclidean_mesh_3d.
                _mesh.push_tri4 (_face, false) ;
        }

        _mesh._euclidean_rdel_3d._tset.
            clear(containers::tight_alloc) ;

        _mesh._euclidean_rdel_3d.
            clear(containers::tight_alloc) ;

        _mesh._euclidean_mesh_3d.
                        _mesh.make_link () ;

        }

        return ( _errv ) ;
    }


#   endif   //__MSH_COPY__



