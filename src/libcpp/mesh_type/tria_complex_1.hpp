
    /*
    --------------------------------------------------------
     * TRIA-COMPLEX-1: piecewise tria. complex in R^1. 
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
     * Last updated: 11 April, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    #pragma once

    #ifndef __TRIA_COMPLEX_1___
    #define __TRIA_COMPLEX_1___

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
            allocator::size_type        size_type ;

    iptr_type static constexpr _dims = +1 ;

    typedef containers::block_array <
                node_type, 
                allocator           >   node_list ;

    typedef containers::block_array <
                edge_type, 
                allocator           >   edge_list ;

    typedef containers::array       <
                iptr_type, 
                allocator           >   iptr_list ;
    
    #define __hashscal sizeof(iptr_type)/sizeof(uint32_t)

    class node_hash
        {
        public  :
        __const_ptr(node_list) _nset ;    
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call node_hash  (
            node_list *_nsrc
            ) : _nset( _nsrc) {}
    /*----------------------- hash node indexing for node */
        __inline_call 
            iptr_type operator() (
            iptr_type _npos
            ) const
        {
            return this->_nset->
                operator[](_npos).node( 0) ;
        }
        } ;

    class edge_hash
        {
        public  :
        __const_ptr(edge_list) _eset ;    
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call edge_hash  (
            edge_list *_esrc
            ) : _eset( _esrc) {}
    /*----------------------- hash node indexing for edge */
        __inline_call 
            iptr_type operator() (
            iptr_type _epos
            ) const
        {
            iptr_type _enod[2];
            _enod[0] = this->_eset->
                operator[](_epos).node( 0) ;
            _enod[1] = this->_eset->
                operator[](_epos).node( 1) ;

            algorithms::isort (
                &_enod[0], &_enod[2] , 
                    std::less<iptr_type>());
            
            return hash::hashword (
                (uint32_t*)&_enod[0] ,
                    +2 * __hashscal, +137) ;
        }
        } ;
        
 
    class node_pred
        {
        public  :
        __const_ptr(node_list) _nset ;    
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call node_pred  (
            node_list *_nsrc
            ) : _nset( _nsrc) {}
    /*----------------------- compute "equal-to" for node */
        __inline_call 
            iptr_type operator() (
            iptr_type _ipos,
            iptr_type _jpos
            ) const
        {
            iptr_type _inod[1];
            _inod[0] = this->_nset->
                operator[](_ipos).node( 0) ;
                
            iptr_type _jnod[1];
            _jnod[0] = this->_nset->
                operator[](_jpos).node( 0) ;
            
            return ( _inod[0] ==
                     _jnod[0] );
        }
        } ;
 
    class edge_pred
        {
        public  :
        __const_ptr(edge_list) _eset ;    
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call edge_pred  (
            edge_list *_esrc
            ) : _eset( _esrc) {}
    /*----------------------- compute "equal-to" for edge */
        __inline_call 
            iptr_type operator() (
            iptr_type _ipos,
            iptr_type _jpos
            ) const
        {
            iptr_type _inod[2];
            _inod[0] = this->_eset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_eset->
                operator[](_ipos).node( 1) ;
                
            iptr_type _jnod[2];
            _jnod[0] = this->_eset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_eset->
                operator[](_jpos).node( 1) ;

            algorithms::isort (
                &_inod[0], &_inod[2], 
                    std::less<iptr_type>());

            algorithms::isort (
                &_jnod[0], &_jnod[2], 
                    std::less<iptr_type>());
            
            return ( _inod[0] ==
                     _jnod[0] &&
                     _inod[1] ==
                     _jnod[1] );
        }
        } ;
        
    #undef __hashscal

    iptr_type static const pool_byte_size=96*1024 ;

    typedef allocators::_pool_alloc <
            allocators::basic_alloc ,
                pool_byte_size      >   pool_base ;

    typedef allocators::_wrap_alloc <
                pool_base           >   pool_wrap ;

    typedef containers::hash_table  <
                iptr_type, 
                node_hash, 
                node_pred,
                pool_wrap>              node_maps ;

    typedef containers::hash_table  <
                iptr_type, 
                edge_hash, 
                edge_pred,
                pool_wrap>              edge_maps ;
                
    typedef containers::array_list  <
                iptr_type,
                pool_wrap           >   conn_list ;

    public  :

    pool_base                   _hsrc ;
    pool_base                   _csrc ;

    conn_list                   _adj1 ;
    conn_list                   _adj2 ;

    node_maps                   _map1 ;
    edge_maps                   _map2 ;
    
    node_list                   _set1 ;
    edge_list                   _set2 ;
   
    iptr_list                   _del1 ;
    iptr_list                   _del2 ;
    
    iptr_list                   _tmp1 ;
    iptr_list                   _tmp2 ;


    public  :

    /*
    --------------------------------------------------------
     * GET-NODE: "create" new node, push onto active set.
    --------------------------------------------------------
     */
     
    __inline_call iptr_type _get_node (
        )
    {
        iptr_type _ipos = -1;
        if (this->_del1.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            this->_del1._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos = (iptr_type)
                this->_set1.count() ;
            this->_set1.push_tail() ;
        }

        this->_set1[_ipos].mark() = +0 ;
        this->_set1[_ipos].self() = +0 ;
        
        return ( _ipos ) ;
    }
    
    /*
    --------------------------------------------------------
     * GET-EDGE: "create" new edge, push onto active set.
    --------------------------------------------------------
     */
     
    __inline_call iptr_type _get_edge (
        )
    {
        iptr_type _ipos = -1;
        if (this->_del2.count() != +0 )
        {
    /*---------------------------- recycle from free list */
            this->_del2._pop_tail(_ipos) ;
        }
        else
        {
    /*---------------------------- alloc. from underlying */
            _ipos = (iptr_type)
                this->_set2.count() ;
            this->_set2.push_tail() ;
        }

        this->_set2[_ipos].mark() = +0 ;
        this->_set2[_ipos].self() = +0 ;

        return ( _ipos ) ;
    }
    
    /*
    --------------------------------------------------------
     * PUT-ITEM: "delete" old item, _pop from active set.
    --------------------------------------------------------
     */
     
    __inline_call void_type _put_node (
        iptr_type _ipos
        )
    {
        this->_del1.push_tail (_ipos);
        this->_set1[_ipos].mark() = -1 ;
        this->_set1[_ipos].self() = -1 ;
    }
    
    __inline_call void_type _put_edge (
        iptr_type _ipos
        )
    {
        this->_del2.push_tail (_ipos);
        this->_set2[_ipos].mark() = -1 ;
        this->_set2[_ipos].self() = -1 ;
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
        /*------------------------------ init. adj. lists */
            _adj1(pool_wrap(&_csrc)),
            _adj2(pool_wrap(&_csrc)),
        /*------------------------------ init. hash lists */
            _map1(
             node_hash(& this->_set1) , 
             node_pred(& this->_set1) , 
            +.8, (pool_wrap(&_hsrc))) ,
            _map2(
             edge_hash(& this->_set2) , 
             edge_pred(& this->_set2) , 
            +.8, (pool_wrap(&_hsrc))) ,
        /*------------------------------ init. face lists */
            _set1(_asrc),_set2(_asrc) ,
        /*------------------------------ init. free lists */
            _del1(_asrc),_del2(_asrc) ,
        /*------------------------------ init. work lists */
            _tmp1(_asrc),_tmp2(_asrc)
        {
        }
        
    /*
    --------------------------------------------------------
     * 'clear' a tria-complex
    --------------------------------------------------------
     */
        
    __normal_call void_type clear (
        containers::alloc_types _kind = 
        containers::loose_alloc
        )
    {
            this->_adj1.clear (_kind) ;
            this->_adj2.clear (_kind) ;
            
            this->_map1.clear (_kind) ;
            this->_map2.clear (_kind) ;

            this->_set1.clear (_kind) ;
            this->_set2.clear (_kind) ;
       
            this->_del1.clear (_kind) ;
            this->_del2.clear (_kind) ;
        
            this->_tmp1.clear (_kind) ;
            this->_tmp2.clear (_kind) ;
            
        //!!do something with pool'd-alloc?
    }
    
    /*
    --------------------------------------------------------
     * INIT-LIST: re-size adj. lists on demand.
    --------------------------------------------------------
     */
     
    __inline_call void_type init_list (
        conn_list&_list,
        size_type _lpos
        )
    {
        if (_lpos >= _list._lptr.count())
        _list._lptr.set_count (
            _lpos + 1, 
        containers::loose_alloc, nullptr) ;
    }
    
    /*
    --------------------------------------------------------
     * PUSH-NODE: append new 0-node to complex.
    --------------------------------------------------------
     */
     
    __normal_call iptr_type push_node (
        node_type const&_ndat,
        bool_type _link = true,
        iptr_type _itop = -1
        )
    {
        iptr_type _ipos = _get_node() ;
        
        if (!_link)
        {
  
    /*------------------------ init. external d-face data */
        this->_set1[_ipos]  = _ndat ;
        this->_set1[_ipos].mark() = 0 ;
        
        this->_set1[_ipos].self() = 1 ;
        this->
       _set1 [_ipos].node(0) =_ipos ;
       
    /*------------------------ push face data to hash set */
        this->_map1.push(_ipos) ;
        
    /*------------------------ init. local adj. index set */
        init_list(this->_adj1, _ipos) ;
  
        }
        else
        {
              
    /*------------------------ init. external d-face data */
        this->_set1[_ipos]  = _ndat ;
        this->_set1[_ipos].mark() = 0 ;
        
        if (_itop == -1)
        {
        this->_set1[_ipos].self() = 1 ;
        this->
       _set1 [_ipos].node(0) =_ipos ;
        }
        else
        this->_set1[_ipos].self() = 0 ;
   
        typename
        node_maps::_write_it  _same ;
        if (this->
           _map1.find(_ipos , _same))
        {
    /*---- existing d-face found: append to existing data */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, copy data */
                this->_set1[*_same] =
                this->_set1[ _ipos] ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
                this->
               _adj1.push(_itop,*_same) ;
            }        
                           
            _put_node(_ipos) ; 
        }
        else
        {
    /*---- d-face is new: push and descend to (d-1)-faces */   
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, keep data */
                init_list(
                    this->_adj1, _ipos) ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
                this->
               _adj1.push(_itop, _ipos) ;
            }
            
        /*-------------- push new face data onto hash */
            this->_map1.push(_ipos) ;
        }
        
        }
        
        return _ipos ;
    }
    
    /*
    --------------------------------------------------------
     * PUSH-EDGE: append new 1-edge to complex.
    --------------------------------------------------------
     */
     
    __normal_call iptr_type push_edge (
        edge_type const&_edat,
        bool_type _link = true,
        iptr_type _itop = -1   
        )
    {
        iptr_type _ipos = _get_edge() ;
        iptr_type _npos ;
        
        if (!_link)
        {
        
    /*------------------------ init. external d-face data */
        this->_set2[_ipos]  = _edat ;
        this->_set2[_ipos].mark() = 0 ;
        this->_set2[_ipos].self() = 1 ;
        
    /*------------------------ push face data to hash set */
        this->_map2.push(_ipos) ;
        
    /*------------------------ init. local adj. index set */
        init_list(this->_adj2, _ipos) ;
        
        }
        else
        {
        
    /*------------------------ init. external d-face data */
        this->_set2[_ipos]  = _edat ;
        this->_set2[_ipos].mark() = 0 ;
        
        if (_itop == -1)
        this->_set2[_ipos].self() = 1 ;
        else
        this->_set2[_ipos].self() = 0 ;
    
        typename
        edge_maps::_write_it  _same ;
        if (this->
           _map2.find(_ipos , _same))
        {
    /*---- existing d-face found: append to existing data */
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, copy data */
                this->_set2[*_same] =
                this->_set2[ _ipos] ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
                this->
               _adj2.push(_itop,*_same) ;
            }     
                      
            _put_edge(_ipos) ; 
        }
        else
        {
    /*---- d-face is new: push and descend to (d-1)-faces */   
            if (_itop == -1)
            {
        /*----- if it's a "top"-level face, keep data */
                init_list(
                    this->_adj2, _ipos) ;
            }
            else
            {
        /*----- otherwise, append index to adj. lists */
                this->
               _adj2.push(_itop, _ipos) ;
            }
            
        /*-------------- descend into (d-1)-face data */
            for (_npos = 2; _npos-- != 0; )
            {
                node_type _ndat;                
                _ndat.node(0) = 
                    _edat.node(_npos) ;
                
                push_node( 
                 _ndat, _link, _ipos) ;
            }
            
        /*-------------- push new face data onto hash */
            this->_map2.push(_ipos) ;
        }
        
        }
        
        return _ipos ;
    }
    
    /*
    --------------------------------------------------------
     * MAKE-PTRS: build item-to-item adj.
    --------------------------------------------------------
     */
  
    __normal_call void_type make_ptrs (
        )
    {
        this->_adj1.empty () ;
        this->_adj2.empty () ;
    
        iptr_type _epos = +0 ;
        for (auto _iter  = this->_set2.head();
                  _iter != this->_set2.tend();
                ++_iter, ++_epos  )
        {
        /*-------------- descend into (d-1)-face data */
            iptr_type _ipos;
            for (_ipos = +2; _ipos-- != 0; )
            {
                node_type _ndat;                
                _ndat.node(0) = 
                    _iter->node( _ipos) ;
                
                push_node(
                    _ndat, true, _epos) ;
            }
        }
    }
    
    /*
    --------------------------------------------------------
     * NULL-ITEM: TRUE if item is delete-able.
    --------------------------------------------------------
     */
     
    __inline_call bool_type null_node (
        iptr_type _npos
        )
    {
        return this->_set1[_npos].self()==0
            && this->_adj1.empty(_npos);
    }
    
    __inline_call bool_type null_edge (
        iptr_type _epos
        )
    {
        return this->_set2[_epos].self()==0
            && this->_adj2.empty(_epos);
    }
  
    /*
    --------------------------------------------------------
     * _POP-LIST: delete item from adj. list.
    --------------------------------------------------------
     */
     
    __normal_call void_type _pop_list (
        conn_list &_list,
        iptr_type  _lpos,
        iptr_type  _itop
        )
    {
    /*--------------------------- scan list and pop match */
        typename conn_list::_write_it
           _prev = _list.hend(_lpos),
           _iter = _list.head(_lpos),
           _tend = _list.tend(_lpos);
        
        for ( ; _iter != _tend; 
                _prev  = _iter, ++_iter)
        {
            if(*_iter == _itop)
            {
                _list._pop (
                    _prev, _iter, _lpos) ; 
                
                break  ;
            }
        }
    }
  
    /*
    --------------------------------------------------------
     * _POP-NODE: delete 0-node from complex.
    --------------------------------------------------------
     */
     
    __normal_call void_type _pop_node (
        iptr_type*_nptr ,
        iptr_type _itop = -1
        )
    {
        iptr_type _npos = -1 ;
        
        iptr_type  _node[1];
        _node[0] = _nptr[0];
        
    /*-------------------------- find current 0-node pos. */
        if ((_npos = _node[0]) < +0)
        {
            return ;
        }
    
    /*-------------------------- _pop current 0-node data */   
        _pop_node(_npos, _itop) ;
    }
    
    __normal_call void_type _pop_node (
        iptr_type _npos ,
        iptr_type _itop
        )
    {
        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
            _pop_list(
            this->_adj1, _npos, _itop);
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            this->
           _set1[_npos].self() = 0 ;
        }
       
        if (null_node (_npos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */
            
        iptr_type _same ;
        this->_map1._pop(_npos, _same);
        
        _put_node(_npos); 
        }
    }
  
    /*
    --------------------------------------------------------
     * _POP-EDGE: delete 1-edge from complex.
    --------------------------------------------------------
     */
     
    __normal_call void_type _pop_edge (
        iptr_type*_nptr ,
        iptr_type _itop = -1
        )
    {
        iptr_type _epos = -1 ;
        
        iptr_type  _node[2];
        _node[0] = _nptr[0];
        _node[1] = _nptr[1];
        
    /*-------------------------- find current 1-edge pos. */
        if ( !find_edge (_node, _epos))
        {
            return ;
        }
        
    /*-------------------------- _pop current 1-edge data */   
        _pop_edge(_epos, _itop) ;
    }
     
    __normal_call void_type _pop_edge (
        iptr_type _epos ,
        iptr_type _itop = -1
        )
    {   
        iptr_type _npos = -1 ;
    
        if (_itop != -1)
        {
    /*---- deleteing (d+k)-face: scan adj. and _pop _itop */
            _pop_list(
            this->_adj2, _epos, _itop);
        }
        else
        {
    /*---- deleting (d+0)-face: set "self" marker to null */
            this->
           _set2[_epos].self() = 0 ;
        }
       
        if (null_edge (_epos))
        {
    /*---- ref. count: delete (d+0), (d-1)-faces if empty */
        iptr_type _same ;
        this->_map2._pop(_epos, _same);
        
        iptr_type _node [ 2] ;
        _node[0] = this->
            _set2[_epos].node(0) ;
        _node[1] = this->
            _set2[_epos].node(1) ;
        
        _put_edge(_epos);
        
        for (_npos = 2; _npos-- != 0; )
        {
            _pop_node(
               &_node[_npos], _epos);
        }
        }
    }
    
    /*
    --------------------------------------------------------
     * FIND-EDGE: return index of assoc. 1-edge.
    --------------------------------------------------------
     */
     
    __normal_call bool_type find_edge (
        iptr_type*_node,
        iptr_type&_epos
        )
    {
    /*-------------------------- find current 1-edge pos. */        
        iptr_type _ipos = _get_edge() ;
        
        this->
       _set2[_ipos].node(0)=_node[0];
        this->
       _set2[_ipos].node(1)=_node[1];      
        
        typename edge_maps::_write_it 
           _same ;
        if (this->
           _map2.find(_ipos, _same))
        {
    /*------------------------------- found matching face */
            _put_edge(_ipos) ;

            _epos=*_same ;
        
            return  true ;
        }
        else
        {
    /*------------------------------- couldn't find match */
            _put_edge(_ipos) ;

            return false ;
        }
    }
    
    /*
    --------------------------------------------------------
     * NODE-EDGE: form node-to-edge adj. list.
    --------------------------------------------------------
     */
     
    template <
    typename      list_type
             >
    __normal_call void_type node_edge (
        iptr_type*_node,
        list_type&_conn
        )
    {
        node_edge(_node [ 0], _conn) ;
    }     
        
    template <
    typename      list_type
             >
    __normal_call void_type node_edge (
        iptr_type _npos,
        list_type&_conn
        )
    {   
        this->_tmp1.set_count(  +0 ) ;
        this->_tmp2.set_count(  +0 ) ;
      
        this->_tmp1.push_tail(_npos) ;
       
    /*-------------------------- find set of adj. 1-edges */
        for (auto _iter  = this->_tmp1.head(); 
                  _iter != this->_tmp1.tend();
                ++_iter  )
        {
            for (auto _iadj  = 
                 this->_adj1.head(*_iter) ;
                      _iadj != 
                 this->_adj1.tend(*_iter) ;
                    ++_iadj  )
            {
                if (this->
                   _set2[*_iadj].mark() == 0)
                {
                   _conn.push_tail (*_iadj);
                    this->
                   _set2[*_iadj].mark() += 1;
                }
            }
        }
  
    /*-------------------------- flip d-face marker lists */
        for (auto _iter  = _conn.head() ; 
                  _iter != _conn.tend() ;
                ++_iter  )
        {
            this->_set2 [*_iter].mark() = +0;
        }
    }
            
    } ;


    }

#   endif  //__TRIA_COMPLEX_1__



