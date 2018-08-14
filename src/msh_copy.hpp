
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
     * Last updated: 21 March, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
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
        rdel_data &_rdel ,
        mesh_data &_mesh
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;
        __unreferenced(_jlog) ;
    
    //!! need to free memory here as well somehow?
    
        if (_rdel._ndim == +2 &&
            _rdel._kind ==
            jmsh_kind::euclidean_mesh)
        {
      
        _mesh._euclidean_mesh_2d._mesh.
            clear(containers::loose_alloc) ;
      
        _mesh._kind = _rdel._kind ;
        _mesh._ndim = _rdel._ndim ;
        
        _mesh._euclidean_mesh_2d._mesh.
            _map1._lptr.set_count (
           (_rdel._euclidean_rdel_2d.
            _tria._nset.count()*5)/4, 
        containers::loose_alloc, nullptr);
                
        _mesh._euclidean_mesh_2d._mesh.
            _map2._lptr.set_count (
           (_rdel._euclidean_rdel_2d.
            _eset.count()*5) / 4 +
           (_rdel._euclidean_rdel_2d.
            _tset.count()*7) / 3 , 
        containers::loose_alloc, nullptr);
        
        _mesh._euclidean_mesh_2d._mesh.
            _map3._lptr.set_count (
           (_rdel._euclidean_rdel_2d.
            _tset.count()*5) / 4 , 
        containers::loose_alloc, nullptr);
        
        
        for (auto _iter  = _rdel.
        _euclidean_rdel_2d._tria._nset.head() ;
                  _iter != _rdel.
        _euclidean_rdel_2d._tria._nset.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
            typename mesh_data::
                euclidean_mesh_2d::
            mesh_type::node_type _node ;
            _node.pval(0) = _iter->pval(0) ;
            _node.pval(1) = _iter->pval(1) ;   
            
            _node.pval(2) = (real_type)+0. ;
                  
            _node.hidx () = _iter->idxh () ;
            
            _node.fdim () = _iter->fdim () ;
            _node.feat () = _iter->feat () ;
            
            _mesh._euclidean_mesh_2d.
                _mesh.push_node (_node) ;
            }
        }
    
        for (auto _iter  = _rdel.
        _euclidean_rdel_2d._eset._lptr.head() ;
                  _iter != _rdel.
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
                _item->_data._node[0] ;
            _face.node(1) = 
                _item->_data._node[1] ;
            
            _mesh._euclidean_mesh_2d.
                _mesh.push_edge (_face) ;
            }
            }
        }
    
        for (auto _iter  = _rdel.
        _euclidean_rdel_2d._tset._lptr.head() ;
                  _iter != _rdel.
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
                _item->_data._node[0] ;
            _face.node(1) = 
                _item->_data._node[1] ;
            _face.node(2) = 
                _item->_data._node[2] ;
            
            _mesh._euclidean_mesh_2d.
                _mesh.push_tri3 (_face) ;
            }
            }
        }
        
        }
        else
        if (_rdel._ndim == +3 &&
            _rdel._kind ==
            jmsh_kind::euclidean_mesh)
        {
        
        _mesh._euclidean_mesh_3d._mesh.
            clear(containers::loose_alloc) ;
        
        _mesh._kind = _rdel._kind ;
        _mesh._ndim = _rdel._ndim ;
        
        _mesh._euclidean_mesh_3d._mesh.
            _map1._lptr.set_count (
           (_rdel._euclidean_rdel_3d.
            _tria._nset.count()*5)/4, 
        containers::loose_alloc, nullptr);
        
        _mesh._euclidean_mesh_3d._mesh.
            _map2._lptr.set_count (
           (_rdel._euclidean_rdel_3d.
            _eset.count()*5) / 4 +
           (_rdel._euclidean_rdel_3d.
            _fset.count()*7) / 3 +
           (_rdel._euclidean_rdel_3d.
            _tset.count()*9) / 2 , 
        containers::loose_alloc, nullptr);
        
        _mesh._euclidean_mesh_3d._mesh.
            _map3._lptr.set_count (
           (_rdel._euclidean_rdel_3d.
            _fset.count()*5) / 4 +
           (_rdel._euclidean_rdel_3d.
            _tset.count()*8) / 3 , 
        containers::loose_alloc, nullptr);
        
        _mesh._euclidean_mesh_3d._mesh.
            _map4._lptr.set_count (
           (_rdel._euclidean_rdel_3d.
            _tset.count()*5) / 4 , 
        containers::loose_alloc, nullptr);
        
        
        for (auto _iter  = _rdel.
        _euclidean_rdel_3d._tria._nset.head() ;
                  _iter != _rdel.
        _euclidean_rdel_3d._tria._nset.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
            typename mesh_data::
                euclidean_mesh_3d::
            mesh_type::node_type _node ;
            _node.pval(0) = _iter->pval(0) ;
            _node.pval(1) = _iter->pval(1) ;
            _node.pval(2) = _iter->pval(2) ;
            
            _node.pval(3) = (real_type)+0. ;
                     
            _node.hidx () = _iter->idxh () ;
            
            _node.fdim () = _iter->fdim () ;
            _node.feat () = _iter->feat () ;
            
            _mesh._euclidean_mesh_3d.
                _mesh.push_node (_node) ;
            }
        }
        
        for (auto _iter  = _rdel.
        _euclidean_rdel_3d._eset._lptr.head() ;
                  _iter != _rdel.
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
                _item->_data._node[0] ;
            _face.node(1) = 
                _item->_data._node[1] ;
            
            _mesh._euclidean_mesh_3d.
                _mesh.push_edge (_face) ;
            }
            }
        }
        
        for (auto _iter  = _rdel.
        _euclidean_rdel_3d._fset._lptr.head() ;
                  _iter != _rdel.
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
                _item->_data._node[0] ;
            _face.node(1) = 
                _item->_data._node[1] ;
            _face.node(2) = 
                _item->_data._node[2] ;
            
            _mesh._euclidean_mesh_3d.
                _mesh.push_tri3 (_face) ;
            }
            }
        }
        
        for (auto _iter  = _rdel.
        _euclidean_rdel_3d._tset._lptr.head() ;
                  _iter != _rdel.
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
            mesh_type::tri4_type _tria;
            _tria.node(0) = 
                _item->_data._node[0] ;
            _tria.node(1) = 
                _item->_data._node[1] ;
            _tria.node(2) = 
                _item->_data._node[2] ;
            _tria.node(3) = 
                _item->_data._node[3] ;
            
            _mesh._euclidean_mesh_3d.
                _mesh.push_tri4 (_tria) ;
            }
            }
        }
        
        }
    
        return ( _errv ) ;
    }
    
    
#   endif   //__MSH_COPY__


    
