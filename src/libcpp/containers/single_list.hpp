
/*
------------------------------------------------------------
 * a singly-linked list object.
------------------------------------------------------------
 *
 * SINGLE-LIST is a "singly-linked" list type, where 
 * list items comprise a singly-oriented chain of 
 * pointers. This variant maintains both "head" and 
 * "tail" objects, but only supports singly-oriented 
 * manipulation. Insertion/deletion is O(1), acccess is 
 * O(N). List count is stored explicitly, and is O(1) 
 * as a result. 
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
------------------------------------------------------------
 *
 * Last updated: 03 May, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __SINGLE_LIST__
#   define __SINGLE_LIST__

#   include "singleiter.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename D
             > 
    class  single_item
    { 
/*------------------------------- singly-linked list item */
    public  :
    
    typedef D                           data_type ; 
    typedef 
    __cont::single_item<data_type>      self_type ;
    
    public  :

    self_type          *_next ;
    
    data_type           _data ;

    public  :
    
/*------------------------------- c'tor/d'tor/assign etc. */
    __inline_call single_item (
        self_type *_nsrc
        ) : _next( _nsrc) 
    { /* void construct */
    }
    __inline_call single_item (
        self_type *_nsrc ,
        data_type const&_dsrc
        ) : _next( _nsrc), 
            _data(
        __copy(data_type,_dsrc))
    { /* copy construct */
    }
    __inline_call single_item (
        self_type *_nsrc ,
        data_type&&_dsrc
        ) : _next( _nsrc), 
            _data(
        __move(data_type,_dsrc))
    { /* move construct */
    }

    __inline_call~single_item()         = default ;

    __inline_call single_item (
        self_type const&_src
        )                               = default ;
    __inline_call single_item (
        self_type &&    _src
        )                               = default ;

    __inline_call 
        self_type & operator= (
        self_type const&_src
        )                               = default ;
    __inline_call 
        self_type & operator= (
        self_type &&    _src
        )                               = default ;

/*------------------------------ access to "next" pointer */
    __inline_call self_type*& next (
        ) 
    {   return ( this->_next ) ; 
    }

    } ;


    template <
    typename D ,
    typename A = allocators::basic_alloc
             >
    class single_list : public 
            allocators::_item_alloc <
                single_item < D>, A >
    {
/*------ double-ended linked-list for singly-linked items */
    public  :
    
    typedef D                           data_type ;
    typedef A                           allocator ;

    typedef containers::single_list <
                data_type , 
                allocator >             self_type ;

    typedef containers::single_item <
                data_type >             item_type ;
    typedef allocators::_item_alloc <
                item_type , 
                allocator >             obj_alloc ;

    typedef containers::const_single_iterator <
                self_type >             _const_it ;
    typedef containers::write_single_iterator <
                self_type >             _write_it ;
    
    typedef typename 
            allocator::size_type        size_type ;
    typedef typename 
            allocator::diff_type        diff_type ;

    private :

    item_type* _hptr ;          // head/tail pointers
    item_type* _tptr ;
    
    size_type  _nobj ;

    private :
    
/*------------------------------ helper - construct range */
    template <
    typename      iter_type
             >
    __normal_call void_type copy_list (
        iter_type _head,
        iter_type _tend,
    __cont::base_iterator_kind
        )
    { /* push the sequence onto list */
        for ( ; _head != _tend; ++_head) 
            push_tail(*_head) ;
    }
    
    __inline_call void_type copy_list (
        size_type _size ,
        data_type const& _dsrc,
    __cont::null_iterator_kind
        )
    {   copy_data(_size, _dsrc);
    }
    
    __normal_call void_type copy_data (
        size_type _size,
        data_type const& _dsrc
        )
    { /* push _size copies onto list */
        for ( ; _size-- != 0; ) 
            push_head(_dsrc);
    }

    public  :

/*--------------------------- default c'tor - do nothing! */
    __inline_call single_list (
        allocator const&_asrc = allocator()
        ) : obj_alloc(  _asrc )
    { /* construct pointers and allocator */
        this->_nobj = +0 ;
        this->_hptr = nullptr ;
        this->_tptr = nullptr ;
    }
/*--------------------------- default c'tor - initialisor */
    __inline_call single_list (
        size_type _size ,
        data_type const&_dsrc = data_type(),
        allocator const&_asrc = allocator()
        ) : obj_alloc(  _asrc )
    {
        this->_nobj = +0 ;
        this->_hptr = nullptr ;
        this->_tptr = nullptr ;
        copy_data(_size,_dsrc);
    }
/*--------------------------- default c'tor - initialisor */
    template <
    typename  iter_type
             >
    __inline_call single_list (
        iter_type _head,
        iter_type _tail,
        allocator const&_asrc = allocator()
        ) : obj_alloc(  _asrc )
    {
        this->_nobj = +0 ;
        this->_hptr = nullptr ;
        this->_tptr = nullptr ;
        copy_list(_head,_tail,
            __cont::iter_kind(_head)) ;
    }
    
/*-------------------------------------------- _def d'tor */
    __inline_call~single_list (
        ) { clear() ; }

/*-------------------------------------------- copy c'tor */
    __inline_call single_list (
        self_type const& _src
        ) : obj_alloc(  
        __copy(obj_alloc,_src))
    {
        this->_nobj = +0 ;
        this->_hptr = nullptr ;
        this->_tptr = nullptr ;
    /*------------------------------- copy data from _src */
        _const_it _head, _tail;
        _head = _src.head();
        _tail = _src.tend();
        copy_list(_head,_tail,
            __cont::iter_kind(_head)) ;
    }
    
/*-------------------------------------------- move c'tor */
    __inline_call single_list (
        self_type && _src
        ) : obj_alloc(  
        __move(obj_alloc,_src))
    {
    /*------------------------------- move data from _src */
        this->_nobj = _src._nobj ;
        this->_hptr = _src._hptr ;
        this->_tptr = _src._tptr ;
        _src. _nobj =  +0 ;
        _src. _hptr = nullptr ;
        _src. _tptr = nullptr ;
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
        static_cast<obj_alloc&>(*this),
        static_cast<obj_alloc&>(_copy)) ;

        swap(this->_hptr, 
            _copy. _hptr) ;
        swap(this->_tptr, 
            _copy. _tptr) ;
        swap(this->_nobj, 
            _copy. _nobj) ;
        }
        
        return ( *this )  ;
    }
    
/*-------------------------------------------- move a-op. */
    __inline_call self_type &operator = (
        self_type && _src
        )
    {
        if (this != &_src)
        {
        self_type _copy (
            __move(self_type, _src)) ;

        using std::swap ;
        swap (
        static_cast<obj_alloc&>(*this),
        static_cast<obj_alloc&>(_copy)) ;

        swap(this->_hptr, 
            _copy. _hptr) ;
        swap(this->_tptr, 
            _copy. _tptr) ;
        swap(this->_nobj, 
            _copy. _nobj) ;
        }
        
        return ( *this )  ;
    }
    
/*-------------------------------- return container count */
    __inline_call size_type count (
        ) const { return ( this->_nobj ) ; }
        
/*-------------------------------- true if sequence empty */
    __inline_call bool_type empty (
        ) const
    { /* return empty status */
        return ( nullptr ==this->_hptr ) ;
    }
    
/*-------------------------------- return container alloc */
    __inline_call allocator get_alloc(
        ) const
    {   return static_cast<allocator>( *this ) ;
    }
  
/*------------------------------ "const" access iterators */
    __inline_call _const_it head (
        ) const
    {/*---- return iterator for list head */
        self_type *_obj = (self_type *)this;
        return _const_it(this->_hptr, _obj);
    }
    
    __inline_call _const_it tail (
        ) const
    {/*---- return iterator for list tail */
        self_type *_obj = (self_type *)this;
        return _const_it(this->_tptr, _obj);
    }
    
    __inline_call _const_it hend (
        ) const
    {/* return iterator "past" head, via null terminator! */
        self_type *_obj = (self_type *)this;
        return _const_it(nullptr, _obj) ;
    }
    
    __inline_call _const_it tend (
        ) const
    {/* return iterator "past" tail, via null terminator! */
        self_type *_obj = (self_type *)this;
        return _const_it(nullptr, _obj) ;
    }

/*------------------------------ "write" access iterators */
    __inline_call _write_it head (
        )
    {/*---- return iterator for list head */
        self_type *_obj = (self_type *)this;
        return _write_it(this->_hptr, _obj);
    }
    
    __inline_call _write_it tail (
        )
    {/*---- return iterator for list tail */
        self_type *_obj = (self_type *)this;
        return _write_it(this->_tptr, _obj);
    }
    
    __inline_call _write_it hend (
        )
    {/* return iterator "past" head, via null terminator! */
        self_type *_obj = (self_type *)this;
        return _write_it(nullptr, _obj) ;
    }
    
    __inline_call _write_it tend (
        )
    {/* return iterator "past" tail, via null terminator! */
        self_type *_obj = (self_type *)this;
        return _write_it(nullptr, _obj) ;
    }
    
/*-------------------------------------- "free" container */
    __normal_call void_type clear (
        ) 
    {/* _pop all items in traversal from head */
        _write_it _head = head();
        _write_it _tend = tend();
        for ( ; _head != _tend; )
        {
            _write_it _prev(_head++) ;
            self_type:: 
             _destruct(_prev.item()) ;
            self_type::
            deallocate(_prev.item(),1) ;
        }
        this->_nobj = +0 ;
        this->_hptr = nullptr ;
        this->_tptr = nullptr ;
    }

/*---------------------------- push data (_def construct) */
    __inline_call _write_it push_head (
        )
    {/* allocate and construct a new raw item */
        item_type *_this_item = 
        self_type::allocate(1); 
        self_type::construct(_this_item, 
            this->_hptr);

        if (this->_hptr!= nullptr) 
        {/* push onto list head only */
            this->_hptr = _this_item ;
        }
        else
        {/* push onto list head/tail */
            this->_hptr = _this_item ;
            this->_tptr = _this_item ;
        }
        this->_nobj += +1 ;

        return 
        _write_it(_this_item,(self_type*)this) ;
    }
    
/*---------------------------- push data (copy construct) */
    __inline_call _write_it push_head (
        data_type const&_data
        )
    {/* allocate and construct a new raw item */
        item_type *_this_item = 
        self_type::allocate(1); 
        self_type::construct(_this_item, 
            this->_hptr,//link! 
                __copy(data_type,_data)) ;

        if (this->_hptr!= nullptr) 
        {/* push onto list head only */
            this->_hptr = _this_item ;
        }
        else
        {/* push onto list head/tail */
            this->_hptr = _this_item ;
            this->_tptr = _this_item ;
        }
        this->_nobj += +1 ;

        return 
        _write_it(_this_item,(self_type*)this) ;
    }
    
/*---------------------------- push data (move construct) */
    __inline_call _write_it push_head (
        data_type &&_data
        )
    {/* allocate and construct a new raw item */
        item_type *_this_item = 
        self_type::allocate(1); 
        self_type::construct(_this_item, 
            this->_hptr,//link! 
                __move(data_type,_data)) ;

        if (this->_hptr!= nullptr) 
        {/* push onto list head only */
            this->_hptr = _this_item ;
        }
        else
        {/* push onto list head/tail */
            this->_hptr = _this_item ;
            this->_tptr = _this_item ;
        }
        this->_nobj += +1 ;

        return 
        _write_it(_this_item,(self_type*)this) ;
    }

/*---------------------------- push data (_def construct) */
    __inline_call _write_it push_tail (
        )
    {/* allocate and construct a new raw item */
        item_type *_this_item = 
        self_type::allocate(1); 
        self_type::construct(_this_item , 
       (item_type *) nullptr) ;

        if (this->_tptr != nullptr) 
        {/* push onto list tail only */
            this->_tptr->next() 
                        = _this_item ;
            this->_tptr = _this_item ;
        }
        else
        {/* push onto list tail/head */
            this->_hptr = _this_item ;
            this->_tptr = _this_item ;
        }
        this->_nobj += +1 ;

        return 
        _write_it(_this_item,(self_type*)this) ;
    }
    
/*---------------------------- push data (copy construct) */
    __inline_call _write_it push_tail (
        data_type const&_data
        )
    {/* allocate and construct a new raw item */
        item_type *_this_item = 
        self_type::allocate(1); 
        self_type::construct(_this_item, 
            (item_type *)nullptr, 
                __copy(data_type,_data)) ;

        if (this->_tptr != nullptr) 
        {/* push onto list tail only */
            this->_tptr->next() 
                        = _this_item ;
            this->_tptr = _this_item ;
        }
        else
        {/* push onto list tail/head */
            this->_hptr = _this_item ;
            this->_tptr = _this_item ;
        }
        this->_nobj += +1 ;

        return 
        _write_it(_this_item,(self_type*)this) ;
    }
    
/*---------------------------- push data (move construct) */
    __inline_call _write_it push_tail (
        data_type &&_data
        )
    {/* allocate and construct a new raw item */
        item_type *_this_item = 
        self_type::allocate(1); 
        self_type::construct(_this_item, 
            (item_type *)nullptr, 
                __move(data_type,_data)) ;

        if (this->_tptr != nullptr) 
        {/* push onto list tail only */
            this->_tptr->next() 
                        = _this_item ;
            this->_tptr = _this_item ;
        }
        else
        {/* push onto list tail/head */
            this->_hptr = _this_item ;
            this->_tptr = _this_item ;
        }
        this->_nobj += +1 ;

        return 
        _write_it(_this_item,(self_type*)this) ;
    }

/*---------------------------- push data (copy obj range) */
    template <
    typename      iter_type
             >
    __normal_call void_type push_head (
        iter_type _hend,
        iter_type _tail
        )
    { /* push a full sequence at tail */
        for(; _tail != _hend; --_tail ) 
            push_head(*_tail) ;
    }
    
    template <
    typename      iter_type
             >
    __normal_call void_type push_tail (
        iter_type _head,
        iter_type _tend
        )
    { /* push a full sequence at tail */
        for(; _head != _tend; ++_head ) 
            push_tail(*_head) ;
    }

/*--------------------------------------------- _pop data */
    __inline_call void_type _pop_head (
        )
    {
        __assert( this->_hptr != nullptr &&
        "single_list._pop_head: null ptr");
    /* shuffle next item onto list _head */
        item_type *_head_item = 
        this->_hptr ;
        this->_hptr = _head_item->next();

        if (this->_hptr == nullptr)
        { /* list became empty */
            this->_tptr =  nullptr;
        }
        this->_count -= 1 ;

    /* _destruct and deallocate old head */
        self_type:: _destruct(_head_item) ;     
        self_type::deallocate(_head_item,1) ;
    }
    
/*--------------------------------------------- _pop data */
    __inline_call void_type _pop_head (
        data_type &_data
        )
    {
        __assert( this->_hptr != nullptr &&
        "single_list._pop_head: null ptr");
    /*-- move data before _pop */
        _data = __move (
            data_type, this->_hptr->_data);
    /*-- delegate to _pop impl.*/
        self_type::_pop_head () ;
    }

/*-------------------------- insert data (_def construct) */
    __normal_call _write_it insert_next (
        _write_it  _prev_iter
        )
    {
        item_type *_prev_item, *_next_item ;
        if ((_prev_item = 
             _prev_iter. item()) == nullptr)
        { /*----- push onto list head */
            return push_head() ;
        }
        else
        if ((_next_item = 
             _prev_item->next()) == nullptr)
        { /*----- push onto list tail */
            return push_tail() ;
        }
        else
        { /* insert in middle of list */
            item_type *_this_item = 
            self_type::allocate(1); 
            self_type::construct(
                _this_item, _next_item);
            
            _prev_item->next() = _this_item;
            this->_count += +1 ;

            return _write_it(
              _this_item, (self_type*)this);
        }
    }
    
/*-------------------------- insert data (copy construct) */
    __normal_call _write_it insert_next (
        _write_it  _prev_iter,
        data_type const&_data
        )
    {
        item_type *_prev_item, *_next_item ;
        if ((_prev_item = 
             _prev_iter. item()) == nullptr)
        { /*----- push onto list head */
            return push_head(__copy(data_type,_data)) ;
        }
        else
        if ((_next_item = 
             _prev_item->next()) == nullptr)
        { /*----- push onto list tail */
            return push_tail(__copy(data_type,_data)) ;
        }
        else
        { /* insert in middle of list */
            item_type *_this_item = 
            self_type::allocate(1); 
            self_type::construct(_this_item, 
                _next_item, //link!
                    __copy(data_type,_data)) ;
            
            _prev_item->next() = _this_item;
            this->_count += +1 ;

            return _write_it(
              _this_item, (self_type*)this);
        }
    }
    
/*-------------------------- insert data (move construct) */
    __normal_call _write_it insert_next (
        _write_it _prev_iter,
        data_type &&_data
        )
    {
        item_type *_prev_item, *_next_item ;
        if ((_prev_item = 
             _prev_iter. item()) == nullptr)
        { /*----- push onto list head */
            return push_head(__move(data_type,_data)) ;
        }
        else
        if ((_next_item = 
             _prev_item->next()) == nullptr)
        { /*----- push onto list tail */
            return push_tail(__move(data_type,_data)) ;
        }
        else
        { /* insert in middle of list */
            item_type *_this_item = 
            self_type::allocate(1); 
            self_type::construct(_this_item, 
                _next_item, //link!
                    __move(data_type,_data)) ;
            
            _prev_item->next() = _this_item;            
            this->_count += +1 ;

            return _write_it(
              _this_item, (self_type*)this);
        }
    }

/*-------------------------------------------- erase item */
    __normal_call void_type erase (
        _write_it _prev_iter,
        _write_it _this_iter
        )
    {
        __assert( _this_iter.item() != nullptr &&
            "single_list.erase: null iterator!");
    /* _pop item, re-link, _destruct//deallocate */
        item_type 
       *_this_item = _this_iter. item() ,
       *_prev_item = _prev_iter. item() ,
       *_prev_next = _prev_item->next() ,
       *_next_item = _this_item->next() ;

        if (_prev_item == nullptr)
        {/* re-link item neighbours at list head */
            __assert( 
                this->_hptr == _this_item &&
            "single_list.erase: _bad iterator!") ;                      
            this->_hptr = _next_item ;
            if (_next_item == nullptr) 
            this->_tptr = _next_item ;
        }
        else        
        {/* re-link item neighbours in list mid  */
            __assert( 
                _prev_next  == _this_item &&
            "single_list.erase: _bad iterator!") ;          
            _prev_item->next() 
                        = _next_item ;                      
            if (_next_item == nullptr) 
            this->_tptr = _prev_item ;
        }
        this->_nobj -= +1 ;

        self_type:: _destruct(_this_item) ;         
        self_type::deallocate(_this_item,1) ;
    }

    } ;

#   undef  __cont


    }

#   endif   //__SINGLE_LIST__



