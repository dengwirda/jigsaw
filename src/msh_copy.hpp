
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 19 Apr., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_COPY__
#   define __MSH_COPY__

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
        containers::array<iptr_type> _nmap ;

        _nmap.set_count(_mesh.
        _euclidean_rdel_2d._tria._nset.count(),
        containers::loose_alloc, -1) ;

        _mesh._euclidean_mesh_2d._mesh.
            clear(containers::tight_alloc) ;

        iptr_type _npos  = +0 ;
        for (auto _iter  = _mesh.
        _euclidean_rdel_2d._tria._nset.head() ;
                  _iter != _mesh.
        _euclidean_rdel_2d._tria._nset.tend() ;
                ++_iter, ++_npos)
        {
            if (_iter->mark() >= +0)
            {
            typename mesh_data::
                euclidean_mesh_2d::
            mesh_type::node_type _node ;
            _node.pval(0) = _iter->pval(0) ;
            _node.pval(1) = _iter->pval(1) ;

            _node.pval(2) = _iter->pval(2) ;

            _node.hidx () = _iter->idxh () ;

            _node.itag () = _iter->part () ;
            _node.fdim () = _iter->fdim () ;
            _node.feat () = _iter->feat () ;

            _nmap[_npos] = 
            _mesh._euclidean_mesh_2d.
                _mesh.push_node (_node, false) ;
            }
        }

        _mesh._euclidean_rdel_2d._tria.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _mesh.
        _euclidean_rdel_2d._eset._lptr.head() ;
                  _iter != _mesh.
        _euclidean_rdel_2d._eset._lptr.tend() ;
                ++_iter  )
        {
            if(*_iter != nullptr)
            {
            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
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
            }
        }

        _mesh._euclidean_rdel_2d._eset.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _mesh.
        _euclidean_rdel_2d._tset._lptr.head() ;
                  _iter != _mesh.
        _euclidean_rdel_2d._tset._lptr.tend() ;
                ++_iter  )
        {
            if(*_iter != nullptr)
            {
            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
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
            }
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
        containers::array<iptr_type> _nmap ;

        _nmap.set_count(_mesh.
        _euclidean_rdel_3d._tria._nset.count(),
        containers::loose_alloc, -1) ;

        _mesh._euclidean_mesh_3d._mesh.
            clear(containers::tight_alloc) ;

        iptr_type _npos  = +0 ;
        for (auto _iter  = _mesh.
        _euclidean_rdel_3d._tria._nset.head() ;
                  _iter != _mesh.
        _euclidean_rdel_3d._tria._nset.tend() ;
                ++_iter, ++_npos)
        {
            if (_iter->mark() >= +0)
            {
            typename mesh_data::
                euclidean_mesh_3d::
            mesh_type::node_type _node ;
            _node.pval(0) = _iter->pval(0) ;
            _node.pval(1) = _iter->pval(1) ;
            _node.pval(2) = _iter->pval(2) ;

            _node.pval(3) = _iter->pval(3) ;

            _node.hidx () = _iter->idxh () ;

            _node.itag () = _iter->part () ;
            _node.fdim () = _iter->fdim () ;
            _node.feat () = _iter->feat () ;

            _nmap[_npos] = 
            _mesh._euclidean_mesh_3d.
                _mesh.push_node (_node, false) ;
            }
        }

        _mesh._euclidean_rdel_3d._tria.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _mesh.
        _euclidean_rdel_3d._eset._lptr.head() ;
                  _iter != _mesh.
        _euclidean_rdel_3d._eset._lptr.tend() ;
                ++_iter  )
        {
            if(*_iter != nullptr)
            {
            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
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
            }
        }

        _mesh._euclidean_rdel_3d._eset.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _mesh.
        _euclidean_rdel_3d._fset._lptr.head() ;
                  _iter != _mesh.
        _euclidean_rdel_3d._fset._lptr.tend() ;
                ++_iter  )
        {
            if(*_iter != nullptr)
            {
            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
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
            }
        }

        _mesh._euclidean_rdel_3d._fset.
            clear(containers::tight_alloc) ;

        for (auto _iter  = _mesh.
        _euclidean_rdel_3d._tset._lptr.head() ;
                  _iter != _mesh.
        _euclidean_rdel_3d._tset._lptr.tend() ;
                ++_iter  )
        {
            if(*_iter != nullptr)
            {
            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
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
            }
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



