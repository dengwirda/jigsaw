
/*
------------------------------------------------------------
 * an array of linked-lists.
------------------------------------------------------------
 *
 * ARRAY-LIST is a linear array of singly-linked list 
 * objects. It's useful as a base-class for HASH-TABLE 
 * like implementations, amongst other things...
 *
------------------------------------------------------------
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
 
#   pragma once

#   ifndef __ARRAY_LIST__
#   define __ARRAY_LIST__

#   include "single_list.hpp"
#   include "block_array.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename T ,
    typename A = allocators::basic_alloc
             >
    class array_list : public 
            allocators::_item_alloc <
                single_item < T>, A >
    {
/*----------------------- a dynamic array of linked-lists */
    public  :     
    
    typedef T                           data_type ;
    typedef A                           allocator ;
    
    typedef typename 
            allocator::size_type        size_type ;
    typedef typename 
            allocator::diff_type        diff_type ;

    typedef containers::array_list  <
                data_type, 
                allocator       >       self_type ;

    typedef containers::single_item <
                data_type       >       item_type ;

    typedef containers::const_single_iterator <
                self_type       >       _const_it ;
                
    typedef containers::write_single_iterator <
                self_type       >       _write_it ;

    typedef allocators::_item_alloc <
                item_type, 
                allocator       >       item_pool ;

    typedef containers::array   <
                item_type * ,
                allocator       >       lptr_list ;

    public  :
    
    lptr_list              _lptr;

    size_type              _size;

    public  :
    
    /*
    --------------------------------------------------------
     MAKE-ITEM (helper): _new/construct a node.
    --------------------------------------------------------
     */
     
    __inline_call void_type make_item (// copy
        data_type const&_data,
        item_type *&    _item
        )
    {
        _item = self_type::allocate(1); 
        self_type::construct(_item,
                nullptr, //link! 
            __copy(data_type,_data) ) ;
    }
    
    __inline_call void_type make_item (// move
        data_type &&    _data,
        item_type *&    _item
        )
    {
        _item = self_type::allocate(1); 
        self_type::construct(_item,
                nullptr, //link! 
            __move(data_type,_data) ) ;
    }
    
    /*
    --------------------------------------------------------
     FREE-ITEM (helper): _destruct/free a node.
    --------------------------------------------------------
     */
     
    __inline_call void_type free_item (
        item_type *_item
        )
    {
        self_type:: _destruct(_item);       
        self_type::deallocate(_item,+1) ;
    }

    /*
    --------------------------------------------------------
     PUSH/_POP-ITEM (helper): push/_pop item onto list.
    --------------------------------------------------------
     */
     
    __inline_call void_type push_item (
        item_type*&_head,
        item_type *_item
        )
    {
    /*--------------------------------- append about head */
        _item->_next = _head;
        _head =_item ;
    }
    
    __inline_call void_type _pop_item (
        item_type*&_head,
        item_type *_item,
        item_type *_prev
        )
    {
        if (_prev != nullptr)
    /*--------------------------------- relink about item */
            _prev->_next = 
            _item->_next ; 
        else
    /*--------------------------------- relink about head */
            _head = 
            _item->_next ;
    }
    
    public  :
    
    /*
    --------------------------------------------------------
     * container c'tor's/d'tor's/assignment op's etc.
    --------------------------------------------------------
     */
     
    __inline_call  array_list (
        allocator const&_asrc = allocator()
    /*----------------------------- c'tor alloc from obj. */
        ) : item_pool(  _asrc),
    /*----------------------------- c'tor other from obj. */
        _lptr(_asrc) ,
        _size(+0) {}
    
    __inline_call ~array_list (
        ) { clear () ; }

/*-------------------------------------------- move c'tor */    
    __inline_call  array_list (
        self_type &&    _src
        ) : item_pool(  
        __move(item_pool,_src))
    {
        this->_lptr = 
        std::move(_src. _lptr);
        
        this->_size = 
        _src. _size ;     
        _src. _size = +0 ;
    }

/*-------------------------------------------- copy c'tor */
    __inline_call  array_list (
        self_type const&_src
        ) : item_pool(  
        __copy(item_pool,_src))
    {
    /*--------------------------------- build a full copy */
        this->_lptr.set_count(
            _src._lptr.count() ,
        containers::tight_alloc,nullptr);
    
        size_type _lpos = +0 ;
    
        typename lptr_list::_const_it
        _iter = _src._lptr.head(),
        _tend = _src._lptr.tend();
        for( ; _iter != _tend; 
                ++_iter, ++_lpos)
        {
    /*--------------------------------- if non-empty list */
        if (*_iter == nullptr) continue ;
        
        item_type *_next =*_iter ;
        item_type *_item = nullptr ;
        item_type *_head = nullptr ;
        item_type *_prev = nullptr ;
    /*--------------------------------- push current list */
        for( ; _next != nullptr; )
        {
            make_item(
              _next->_data, _item) ;
            
            if (_prev == nullptr)
            _head =_item ;
            else
            _prev->_next = _item;
            
            _prev =_item ;                
            _next =_next-> _next;
            
            this->_size += 1 ;
        }
        this->_lptr[_lpos] = _head ;
        }        
    }
    
/*-------------------------------------------- copy a-op. */
    __inline_call self_type& operator = (
        self_type const& _src
        )
    {
        if (this != &_src)
        {
        self_type _copy (
            __copy(self_type, _src)) ;

        using std::swap ;
        swap (
        static_cast<item_pool&>(*this),
        static_cast<item_pool&>(_copy)) ;

        swap(this->_size, 
            _copy. _size) ;
        swap(this->_lptr, 
            _copy. _lptr) ;
        }
        
        return ( *this )  ;
    }

/*-------------------------------------------- move a-op. */
    __inline_call self_type& operator = (
        self_type &&  _src
        )
    {
        if (this != &_src)
        {
        self_type _copy (
            __move(self_type, _src)) ;

        using std::swap ;
        swap (
        static_cast<item_pool&>(*this),
        static_cast<item_pool&>(_copy)) ;

        swap(this->_size, 
            _copy. _size) ;
        swap(this->_lptr, 
            _copy. _lptr) ;
        }
        
        return ( *this )  ;
    }

    /*
    --------------------------------------------------------
     * GET-ALLOC: return container alloc.
    --------------------------------------------------------
     */
     
    __inline_call allocator get_alloc(
        ) const
    {   return static_cast<allocator>(*this);
    }

    /*
    --------------------------------------------------------
     * CLEAR: clear container lists + backing array.
    --------------------------------------------------------
     */
     
    __normal_call void_type clear (
        containers::alloc_types _kind = 
        containers::loose_alloc
        )
    {
    /*------------------------------ clear lists of items */
        this->       empty(_kind) ;
    /*------------------------------ resize backing array */
        this->_lptr .clear(_kind) ;
    }
     
    __normal_call void_type empty (
        containers::alloc_types _kind = 
        containers::loose_alloc
        )
    {
        __unreferenced(_kind);
    
    /*------------------------------ clear lists of items */
        typename lptr_list::_write_it
        _iter = this->_lptr.head(),
        _tend = this->_lptr.tend();
        for( ; _iter != _tend; ++_iter)
        {
        /*----------------------------- if non-empty list */
            if (*_iter == nullptr) continue ;
            
            item_type*_item =*_iter ;
            item_type*_next ;
        /*----------------------------- free current list */
            for(;_item != nullptr; 
                 _item  = _next)
            {
               _next = _item->_next ; 
               
                free_item(_item) ;
                this->_size -= 1 ;
            }
        }
    }
    
    /*
    --------------------------------------------------------
     * COUNT: return container count.
    --------------------------------------------------------
     */
    
    __inline_call size_type count (
        ) const 
    {   return  this->_size ; 
    }
    
    __inline_call bool_type empty (
        size_type _lpos
        ) const
    {   return  this->
            _lptr[_lpos]==nullptr ;
    }
    
    /*
    --------------------------------------------------------
     * HEAD/TEND: return iterator range.
    --------------------------------------------------------
     */
     
    __inline_call _const_it head (
        size_type _lpos
        ) const
    {   return _const_it(this->_lptr[_lpos], this) ;
    }
    
    __inline_call _write_it head (
        size_type _lpos
        )
    {   return _write_it(this->_lptr[_lpos], this) ;
    }
    
    __inline_call _const_it hend (
        size_type //_lpos
        ) const
    {   return _const_it(nullptr , this);
    }
    
    __inline_call _write_it hend (
        size_type //_lpos
        )
    {   return _write_it(nullptr , this);
    }
    
    __inline_call _const_it tend (
        size_type //_lpos
        ) const
    {   return _const_it(nullptr , this);
    }
    
    __inline_call _write_it tend (
        size_type //_lpos
        )
    {   return _write_it(nullptr , this);
    }
    
    
    /*
    --------------------------------------------------------
     * _POP: "_pop" data from container.
    --------------------------------------------------------
     */
     
    __inline_call void_type _pop (
        _write_it _prev,
        _write_it _item,
        size_type _lpos
        )
    {
        item_type *_iptr=_item.item() ;
        item_type *_pptr=_prev.item() ;
    
        this->_size -= +1 ;
    
        _pop_item(this->_lptr [_lpos] ,
            _iptr, _pptr) ;
            
        free_item( _iptr) ;
    }
    
    /*
    --------------------------------------------------------
     * PUSH: "push" data onto container.
    --------------------------------------------------------
     */
    
    __inline_call _write_it push (  // copy
        data_type const&_data,
        size_type       _lpos
        )
    {
        this->_size += +1 ;
        if (_lpos >= this->_lptr.count())
        this->_lptr.set_count( _lpos + 1, 
        containers::loose_alloc, nullptr) ;
            
        item_type*_item = nullptr ;
        make_item(_data , _item);
        push_item( this->
            _lptr[_lpos], _item);
            
        return _write_it (_item, this) ;
    }
    
    __inline_call _write_it push (  // move
        data_type &&    _data,
        size_type       _lpos
        )
    {
        this->_size += +1 ;
        if (_lpos >= this->_lptr.count())
        this->_lptr.set_count( _lpos + 1, 
        containers::loose_alloc, nullptr) ;
            
        item_type*_item = nullptr ;
        make_item(_data , _item);
        push_item( this->
            _lptr[_lpos], _item);
            
        return _write_it (_item, this) ;
    }
    
    } ;

#   undef  __cont


    }

#   endif   //__ARRAY_LIST__   
    
    
    
