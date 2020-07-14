
    /*
    --------------------------------------------------------
     * TRIA-COMPLEX-1: piecewise tria. complex in T^1.
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
     * Copyright 2013-2020
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __TRIA_COMPLEX_1__
#   define __TRIA_COMPLEX_1__

    namespace mesh {

    /*
    --------------------------------------------------------
     * TRIA-COMPLEX-1 generic 1-tria. simplicial complex.
    --------------------------------------------------------
     * NODE-TYPE - 0-face type + coord.'s
     * EDGE-TYPE - 1-face type
     * ALLOCATOR - memory allocator
    --------------------------------------------------------
     */

    template <
    typename NN ,
    typename E2 ,
    typename AA = allocators::basic_alloc
             >
    class tria_complex_1
    {
/*----------------------------------- tria-complex in R^1 */
    public  :
    typedef NN                          node_type ;
    typedef E2                          edge_type ;
    typedef AA                          allocator ;

    typedef typename
            node_type::real_type        real_type ;
    typedef typename
            edge_type::iptr_type        iptr_type ;
    typedef typename
            allocator::size_type        indx_type ;

    indx_type static constexpr topo_dims = +1 ;

    typedef containers::block_array <
                node_type,
                allocator           >   node_list ;

    typedef containers::block_array <
                edge_type,
                allocator           >   edge_list ;

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

    edge_maps                          _mmE2 ;

    node_list                          _llN1 ;
    edge_list                          _llE2 ;

    indx_list                          _ffN1 ;
    indx_list                          _ffE2 ;

    connector                          _tmp1 ;

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

    public  :

    /*
    --------------------------------------------------------
     * construct tria_complex from alloc. etc...
    --------------------------------------------------------
     */

    __normal_call tria_complex_1 (
        allocator const& _asrc = allocator()
        ) : _hsrc(sizeof (
        typename edge_maps::item_type)),
            _csrc(sizeof (
        typename conn_list::item_type)),
    /*---------------------------------- init. adj. lists */
        _aaN1(pool_wrap(&_csrc)),
    /*---------------------------------- init. hash lists */
        _mmE2(
         edge_hash(& this->_llE2) ,
         edge_pred(& this->_llE2) ,
        +.8, (pool_wrap(&_hsrc))) ,
    /*---------------------------------- init. face lists */
        _llN1(_asrc),_llE2(_asrc) ,
    /*---------------------------------- init. free lists */
        _ffN1(_asrc),_ffE2(_asrc) ,
    /*---------------------------------- init. work lists */
        _tmp1(_asrc)
    {
    }

    __inline_call~tria_complex_1 (
        )
    {
        clear(containers::tight_alloc);
    }

    /*
    --------------------------------------------------------
     * "clear" a tria-complex
    --------------------------------------------------------
     */

    __normal_call void_type clear (
        containers::alloc_types _kind =
        containers::loose_alloc
        )
    {
        this->_aaN1.clear (_kind) ;

        this->_mmE2.clear (_kind) ;

        this->_llN1.clear (_kind) ;
        this->_llE2.clear (_kind) ;

        this->_ffN1.clear (_kind) ;
        this->_ffE2.clear (_kind) ;

        this->_tmp1.clear (_kind) ;

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
       (this->_llE2.count()*5) / 4 + 1
            ) ;

        this->_aaN1.empty () ;

        mesh::link_edge(*this) ;
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
       *this, _node, _kind, _ipos, tria_cell()) ;
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
       *this, _node, _kind, _conn, tria_cell()) ;
    }

    __inline_call void_type connect_1 (
        iptr_type _ipos ,
        char_type _kind ,
        connector&_conn
        )
    {
        mesh::connect_1(
       *this, _ipos, _kind, _conn, tria_cell()) ;
    }

    } ;


    }

#   endif  //__TRIA_COMPLEX_1__



