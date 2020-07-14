
    /*
    --------------------------------------------------------
     * MESH-COMPLEX-2: piecewise mixed complex in T^2.
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
     * Last updated: 28 April, 2020
     *
     * Copyright 2020--
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MESH_COMPLEX_2__
#   define __MESH_COMPLEX_2__

    namespace mesh {

    /*
    --------------------------------------------------------
     * MESH-COMPLEX-2 generic mixed-cell 2-complex.
    --------------------------------------------------------
     * NODE-TYPE - 0-cell type + coord.'s
     * EDGE-TYPE - 1-cell type
     * TRI3-TYPE - 2-cell type
     * QUAD-TYPE - 2-cell type
     * ALLOCATOR - memory allocator
    --------------------------------------------------------
     */

    template <
    typename NN ,
    typename E2 ,
    typename T3 ,
    typename Q4 ,
    typename AA = allocators::basic_alloc
             >
    class mesh_complex_2
    {
/*----------------------------------- mesh-complex in T^2 */
    public  :
    typedef NN                          node_type ;
    typedef E2                          edge_type ;
    typedef T3                          tri3_type ;
    typedef Q4                          quad_type ;
    typedef AA                          allocator ;

    typedef typename
            node_type::real_type        real_type ;
    typedef typename
            edge_type::iptr_type        iptr_type ;
    typedef typename
            allocator::size_type        indx_type ;

    indx_type static constexpr topo_dims = +2 ;

    typedef containers::block_array <
                node_type,
                allocator           >   node_list ;

    typedef containers::block_array <
                edge_type,
                allocator           >   edge_list ;

    typedef containers::block_array <
                tri3_type,
                allocator           >   tri3_list ;

    typedef containers::block_array <
                quad_type,
                allocator           >   quad_list ;

    typedef containers::array       <
                iptr_type,
                allocator           >   iptr_list ;

    typedef containers::array       <
                indx_type,
                allocator           >   indx_list ;

    typedef mesh::node_hash<node_list>  node_hash ;
    typedef mesh::node_pred<node_list>  node_pred ;

    typedef mesh::edge_hash<edge_list>  edge_hash ;
    typedef mesh::edge_pred<edge_list>  edge_pred ;

    typedef mesh::tri3_hash<tri3_list>  tri3_hash ;
    typedef mesh::tri3_pred<tri3_list>  tri3_pred ;

    typedef mesh::quad_hash<quad_list>  quad_hash ;
    typedef mesh::quad_pred<quad_list>  quad_pred ;

    indx_type static
        constexpr pool_byte_size = 96 * 1024 ;

    typedef allocators::_pool_alloc <
            allocators::basic_alloc ,
                pool_byte_size      >   pool_base ;

    typedef allocators::_wrap_alloc <
                pool_base           >   pool_wrap ;

    typedef containers::hash_table  <
                iptr_type,
                edge_hash,
                edge_pred,
                pool_wrap>              edge_maps ;

    typedef containers::hash_table  <
                iptr_type,
                tri3_hash,
                tri3_pred,
                pool_wrap>              tri3_maps ;

    typedef containers::hash_table  <
                iptr_type,
                quad_hash,
                quad_pred,
                pool_wrap>              quad_maps ;

    typedef containers::array       <
                conn_pair,
                allocator           >   connector ;

    typedef containers::array_list  <
                conn_pair,
                pool_wrap           >   conn_list ;

    public  :

    pool_base                          _hsrc ;
    pool_base                          _csrc ;

    conn_list                          _aaN1 ;
    conn_list                          _aaE2 ;

    edge_maps                          _mmE2 ;
    tri3_maps                          _mmT3 ;
    quad_maps                          _mmQ4 ;

    node_list                          _llN1 ;
    edge_list                          _llE2 ;
    tri3_list                          _llT3 ;
    quad_list                          _llQ4 ;

    indx_list                          _ffN1 ;
    indx_list                          _ffE2 ;
    indx_list                          _ffT3 ;
    indx_list                          _ffQ4 ;

    connector                          _tmp1 ;
    connector                          _tmp2 ;

    public  :

    /*
    --------------------------------------------------------
     * Public access to underlying data for k-CELL types.
    --------------------------------------------------------
     */

    __inline_call node_list const& node (
        ) const
    {   return  this->_llN1 ;
    }

    __inline_call node_list      & node (
        )
    {   return  this->_llN1 ;
    }

    __inline_call node_type const& node (
        indx_type _ipos
        ) const
    {   return  this->_llN1[_ipos] ;
    }

    __inline_call node_type      & node (
        indx_type _ipos
        )
    {   return  this->_llN1[_ipos] ;
    }

    __inline_call edge_list const& edge (
        ) const
    {   return  this->_llE2 ;
    }

    __inline_call edge_list      & edge (
        )
    {   return  this->_llE2 ;
    }

    __inline_call edge_type const& edge (
        indx_type _ipos
        ) const
    {   return  this->_llE2[_ipos] ;
    }

    __inline_call edge_type      & edge (
        indx_type _ipos
        )
    {   return  this->_llE2[_ipos] ;
    }

    __inline_call tri3_list const& tri3 (
        ) const
    {   return  this->_llT3 ;
    }

    __inline_call tri3_list      & tri3 (
        )
    {   return  this->_llT3 ;
    }

    __inline_call tri3_type const& tri3 (
        indx_type _ipos
        ) const
    {   return  this->_llT3[_ipos] ;
    }

    __inline_call tri3_type      & tri3 (
        indx_type _ipos
        )
    {   return  this->_llT3[_ipos] ;
    }

    __inline_call quad_list const& quad (
        ) const
    {   return  this->_llQ4 ;
    }

    __inline_call quad_list      & quad (
        )
    {   return  this->_llQ4 ;
    }

    __inline_call quad_type const& quad (
        indx_type _ipos
        ) const
    {   return  this->_llQ4[_ipos] ;
    }

    __inline_call quad_type      & quad (
        indx_type _ipos
        )
    {   return  this->_llQ4[_ipos] ;
    }

    public  :

    /*
    --------------------------------------------------------
     * construct mesh-complex from alloc. etc...
    --------------------------------------------------------
     */

    __normal_call mesh_complex_2 (
        allocator const& _asrc = allocator()
        ) : _hsrc(sizeof (
        typename edge_maps::item_type)),
            _csrc(sizeof (
        typename conn_list::item_type)),
    /*---------------------------------- init. adj. lists */
        _aaN1(pool_wrap(&_csrc)),
        _aaE2(pool_wrap(&_csrc)),
    /*---------------------------------- init. hash lists */
        _mmE2(
         edge_hash(& this->_llE2) ,
         edge_pred(& this->_llE2) ,
        +.8, (pool_wrap(&_hsrc))) ,
        _mmT3(
         tri3_hash(& this->_llT3) ,
         tri3_pred(& this->_llT3) ,
        +.8, (pool_wrap(&_hsrc))) ,
        _mmQ4(
         quad_hash(& this->_llQ4) ,
         quad_pred(& this->_llQ4) ,
        +.8, (pool_wrap(&_hsrc))) ,
    /*---------------------------------- init. face lists */
        _llN1(_asrc),_llE2(_asrc) ,
        _llT3(_asrc),_llQ4(_asrc) ,
    /*---------------------------------- init. free lists */
        _ffN1(_asrc),_ffE2(_asrc) ,
        _ffT3(_asrc),_ffQ4(_asrc) ,
    /*---------------------------------- init. work lists */
        _tmp1(_asrc),_tmp2(_asrc)
    {
    }

    __inline_call~mesh_complex_2 (
        )
    {
        clear(containers::tight_alloc);
    }

    /*
    --------------------------------------------------------
     * "clear" a mesh-complex
    --------------------------------------------------------
     */

    __normal_call void_type clear (
        containers::alloc_types _kind =
        containers::loose_alloc
        )
    {
        this->_aaN1.clear (_kind) ;
        this->_aaE2.clear (_kind) ;

        this->_mmE2.clear (_kind) ;
        this->_mmT3.clear (_kind) ;
        this->_mmQ4.clear (_kind) ;

        this->_llN1.clear (_kind) ;
        this->_llE2.clear (_kind) ;
        this->_llT3.clear (_kind) ;
        this->_llQ4.clear (_kind) ;

        this->_ffN1.clear (_kind) ;
        this->_ffE2.clear (_kind) ;
        this->_ffT3.clear (_kind) ;
        this->_ffQ4.clear (_kind) ;

        this->_tmp1.clear (_kind) ;
        this->_tmp2.clear (_kind) ;

        this->_hsrc.clear ();
        this->_csrc.clear ();
    }

    /*
    --------------------------------------------------------
     * MAKE-LINK: build item-to-item adj.
    --------------------------------------------------------
     */

    __normal_call void_type make_link (
        )
    {
        this->_mmE2.set_slots(
       (this->_llE2.count()*5) / 4 + 1 +
       (this->_llT3.count()*3) / 2 + 1 +
       (this->_llQ4.count()*4) / 2 + 1
            ) ;

        this->_mmT3.set_slots(
       (this->_llT3.count()*5) / 4 + 1
            ) ;
        this->_mmQ4.set_slots(
       (this->_llQ4.count()*5) / 4 + 1
            ) ;

        this->_aaN1.empty () ;
        this->_aaE2.empty () ;

        mesh::link_edge(*this) ;
        mesh::link_tri3(*this) ;
        mesh::link_quad(*this) ;
    }

    /*
    --------------------------------------------------------
     * PUSH-NODE: append new 0-node to complex.
    --------------------------------------------------------
     */

    __inline_call iptr_type push_node (
        node_type const& _ndat,
        bool_type _link = true,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        return mesh::push_node(
            *this, _ndat, _link, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * PUSH-EDGE: append new 1-edge to complex.
    --------------------------------------------------------
     */

    __inline_call iptr_type push_edge (
        edge_type const& _edat,
        bool_type _link = true,
        iptr_type _itop = -1  ,
        char_type _kind = EMPTY_tag
        )
    {
        return mesh::push_edge(
            *this, _edat, _link, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * PUSH-TRIA: append new 2-tria to complex.
    --------------------------------------------------------
     */

    __inline_call iptr_type push_tri3 (
        tri3_type const& _tdat,
        bool_type _link = true,
        iptr_type _itop = -1  ,
        char_type _kind = EMPTY_tag
        )
    {
        return mesh::push_tri3(
            *this, _tdat, _link, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * PUSH-QUAD: append new 2-quad to complex.
    --------------------------------------------------------
     */

    __inline_call iptr_type push_quad (
        quad_type const& _qdat,
        bool_type _link = true,
        iptr_type _itop = -1  ,
        char_type _kind = EMPTY_tag
        )
    {
        return mesh::push_quad(
            *this, _qdat, _link, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * _POP-NODE: delete 0-node from complex.
    --------------------------------------------------------
     */

    __inline_call void_type _pop_node (
        iptr_type const*_nptr,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_node(*this, _nptr, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * _POP-EDGE: delete 1-edge from complex.
    --------------------------------------------------------
     */

    __inline_call void_type _pop_edge (
        iptr_type const*_nptr,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_edge(*this, _nptr, _itop, _kind) ;
    }

    __inline_call void_type _pop_edge (
        iptr_type _epos ,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_edge(*this, _epos, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * _POP-TRIA: delete 2-tria from complex.
    --------------------------------------------------------
     */

    __inline_call void_type _pop_tri3 (
        iptr_type const*_nptr,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_tri3(*this, _nptr, _itop, _kind) ;
    }

    __inline_call void_type _pop_tri3 (
        iptr_type _tpos ,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_tri3(*this, _tpos, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * _POP-QUAD: delete 2-quad from complex.
    --------------------------------------------------------
     */

    __inline_call void_type _pop_quad (
        iptr_type const*_nptr,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_quad(*this, _nptr, _itop, _kind) ;
    }

    __inline_call void_type _pop_quad (
        iptr_type _qpos ,
        iptr_type _itop = -1 ,
        char_type _kind = EMPTY_tag
        )
    {
        mesh::_pop_quad(*this, _qpos, _itop, _kind) ;
    }

    /*
    --------------------------------------------------------
     * FIND-NODE: return index of assoc. 0-node.
    --------------------------------------------------------
     */

    __inline_call bool_type find_node (
        iptr_type const*_nptr ,
        iptr_type&_npos
        ) const
    {
        return mesh::find_node(*this, _nptr, _npos) ;
    }

    /*
    --------------------------------------------------------
     * FIND-EDGE: return index of assoc. 1-edge.
    --------------------------------------------------------
     */

    __inline_call bool_type find_edge (
        iptr_type const*_node ,
        iptr_type&_epos
        )
    {
        return mesh::find_edge(*this, _node, _epos) ;
    }

    /*
    --------------------------------------------------------
     * FIND-TRIA: return index of assoc. 2-tria.
    --------------------------------------------------------
     */

    __inline_call bool_type find_tri3 (
        iptr_type const*_node ,
        iptr_type&_tpos
        )
    {
        return mesh::find_tri3(*this, _node, _tpos) ;
    }

    /*
    --------------------------------------------------------
     * FIND-QUAD: return index of assoc. 2-quad.
    --------------------------------------------------------
     */

    __inline_call bool_type find_quad (
        iptr_type const*_node ,
        iptr_type&_qpos
        )
    {
        return mesh::find_quad(*this, _node, _qpos) ;
    }

    /*
    --------------------------------------------------------
     * FIND-CELL: return index of assoc. k-cell.
    --------------------------------------------------------
     */

    __inline_call bool_type find_cell (
        iptr_type const*_node ,
        char_type _kind ,
        iptr_type&_ipos
        )
    {
        return mesh::find_cell(
       *this, _node, _kind, _ipos, mesh_cell()) ;
    }

    /*
    --------------------------------------------------------
     * CONNECT-1: form this-to-1-cell adj. list.
    --------------------------------------------------------
     */

    __inline_call void_type connect_1 (
        iptr_type const*_node ,
        char_type _kind ,
        connector&_conn
        )
    {
        mesh::connect_1(
       *this, _node, _kind, _conn, mesh_cell()) ;
    }

    __inline_call void_type connect_1 (
        iptr_type _ipos ,
        char_type _kind ,
        connector&_conn
        )
    {
        mesh::connect_1(
       *this, _ipos, _kind, _conn, mesh_cell()) ;
    }

    /*
    --------------------------------------------------------
     * CONNECT-2: form this-to-2-cell adj. list.
    --------------------------------------------------------
     */

    __inline_call void_type connect_2 (
        iptr_type const*_node ,
        char_type _kind ,
        connector&_conn
        )
    {
        mesh::connect_2(
       *this, _node, _kind, _conn, mesh_cell()) ;
    }

    __inline_call void_type connect_2 (
        iptr_type _ipos ,
        char_type _kind ,
        connector&_conn
        )
    {
        mesh::connect_2(
       *this, _ipos, _kind, _conn, mesh_cell()) ;
    }

    } ;


    }

#   endif  //__MESH_COMPLEX_2__



