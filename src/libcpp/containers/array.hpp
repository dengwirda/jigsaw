
/* 
------------------------------------------------------------
 * a contiguous dynamically allocated buffer.
------------------------------------------------------------
 *
 * ARRAY is a dynamically allocated, contiguous array 
 * object, where the linear array storage is guaranteed 
 * to consist of a single block, as returned by malloc, 
 * etc. Provides O(1) random access, amortarised 
 * O(log(n)) incremental insertion, etc. This is a 
 * striped-down implementation -- offering only 
 * "efficient" operations (no insert, delete, etc).
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

#   ifndef __ARRAY__
#   define __ARRAY__

#   include "array_iter.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename D , 
    typename A =
            allocators::basic_alloc
             >
    class array : public 
            allocators::_item_alloc < D, A >
    { 
/*------ a dynamically allocated, contiguous array object */
    public  :
    
    typedef D                               data_type ;
    typedef A                               allocator ;

    typedef typename 
            allocator::size_type            size_type ;
    typedef typename 
            allocator::diff_type            diff_type ;

    typedef __cont::array       <
                data_type, 
                allocator       >           self_type ;
        
    typedef __cont::write_array_iterator <
                self_type       >           _write_it ;
    typedef __cont::const_array_iterator <
                self_type       >           _const_it ;
   
    typedef allocators::_item_alloc <
                data_type, 
                allocator       >           obj_alloc ;     

    private :

        enum {_hptr, _tptr, _lptr} ;

    containers::
    fixed_array<__write_ptr(data_type),+3>  _ptrs ;

    private :

/*------------------------------ helper - construct range */

    __normal_call void_type ctor_iter (
        _write_it _head,
        _write_it _tail
        )
    { /* _def construct sequence */
        for ( ; _head != _tail ; ++_head)
        {
            self_type::construct(&*_head) ;
        }
    }   
    __normal_call void_type ctor_iter (
        _write_it _head,
        _write_it _tail,
        data_type const& _data
        )
    { /* copy construct sequence */
        for ( ; _head != _tail ; ++_head) 
        {
            self_type::construct(&*_head, 
                                   _data) ;
        }
    }
    
/*------------------------------ helper - _destruct range */

    __normal_call void_type dtor_iter (
        _write_it _head,
        _write_it _tail
        )
    {
        for ( ; _head != _tail ; ++_head)
        {
            self_type::_destruct(&*_head) ;
        }
    }
    
/*------------------------------ helper - inc. allocation */

    __inline_call void_type inc_alloc (
        size_type _new_count
        )
    { /* Grow allocation as geometric series */
        if (alloc() < _new_count)
        {
        if (alloc() <= +64)
        set_alloc(std::max((alloc()*2)/1,_new_count));
        else
        /* from x^n = 1 + ... + x^(n-2)  */
        set_alloc(std::max((alloc()*3)/2,_new_count));
        }
    }
    
/*------------------------------ helper - copy allocation */
    
    __inline_call void_type copy_iter (
        size_type _size ,
        data_type const&_dsrc,
        __cont::null_iterator_kind
        )
    { /* copy _data onto object */
        set_count(_size, __cont::tight_alloc, _dsrc) ;
    }
    
    template <
    typename      iter_type
             >
    __normal_call void_type copy_iter (
        iter_type _head,
        iter_type _tail,
    __cont::base_iterator_kind
        )
    { /* copy range onto object */
        for(;_head!=_tail;++_head) push_tail(*_head) ;
    }
    
    public  :
    
/*--------------------------- default c'tor - do nothing! */
    __inline_call array (
        allocator const&_asrc = allocator()
        ) : obj_alloc(  _asrc )
    {
        this->_ptrs[_hptr] = nullptr;
        this->_ptrs[_tptr] = nullptr;
        this->_ptrs[_lptr] = nullptr;
    }
    
/*--------------------------- default c'tor - initialisor */
    __inline_call array (
        size_type _size,
        data_type const&_dsrc = data_type(),
        allocator const&_asrc = allocator()
        ) : obj_alloc(  _asrc )
    {
        this->_ptrs[_hptr] = nullptr;
        this->_ptrs[_tptr] = nullptr;
        this->_ptrs[_lptr] = nullptr;
        copy_iter(_size,_dsrc,
        __cont::null_iterator_kind()) ;
    }
    
/*--------------------------- default c'tor - initialisor */
    template <
    typename  iter_type
             >
    __inline_call array (
        iter_type _head,
        iter_type _tail,
        allocator const&_asrc = allocator()
        ) : obj_alloc(  _asrc )
    {
        this->_ptrs[_hptr] = nullptr;
        this->_ptrs[_tptr] = nullptr;
        this->_ptrs[_lptr] = nullptr;
        copy_iter(_head,_tail,
            __cont::iter_kind(_head)) ;
    }
    
/*-------------------------------------------- _def d'tor */
    __inline_call~array (
        )
    {/* _destruct all items while allocator remains valid */
        set_count(+0, __cont::tight_alloc); 
    }
    
/*-------------------------------------------- copy c'tor */
    __inline_call array (
        self_type const& _src
        ) : obj_alloc(
        __copy(obj_alloc,_src))
    {
        this->_ptrs[_hptr] = nullptr;
        this->_ptrs[_tptr] = nullptr;
        this->_ptrs[_lptr] = nullptr;
    /*------------------------------------- alloc storage */
        set_alloc(_src.count());
    /*------------------------------------- copy sequence */
        _const_it _head, _tail ; 
        _head = _src.head();
        _tail = _src.tend();
        copy_iter(_head,_tail,
            __cont::iter_kind(_head)) ;
    }
    
/*-------------------------------------------- move c'tor */
    __inline_call array (
        self_type &&_src
        ) : obj_alloc(
        __move(obj_alloc,_src))
    {
        this->_ptrs[_hptr] = 
        _src. _ptrs[_hptr] ;
        this->_ptrs[_tptr] = 
        _src. _ptrs[_tptr] ;
        this->_ptrs[_lptr] = 
        _src. _ptrs[_lptr] ; 
        _src. _ptrs.fill(nullptr);
    }
    
/*-------------------------------------------- copy a-op. */
    __inline_call self_type& operator = (
        self_type const&_src
        )
    {
        if (this != &_src)
        {
        self_type _copy (
            __copy(self_type , _src)) ;

        using std::swap ;
        swap (
        static_cast<obj_alloc&>(*this),
        static_cast<obj_alloc&>(_copy)
            ) ;

        swap(this->_ptrs[_hptr], 
            _copy. _ptrs[_hptr]) ;
        swap(this->_ptrs[_tptr], 
            _copy. _ptrs[_tptr]) ;
        swap(this->_ptrs[_lptr], 
            _copy. _ptrs[_lptr]) ;
        }
        
        return ( *this ) ;
    }
    
/*-------------------------------------------- move a-op. */
    __inline_call self_type &operator = (
        self_type && _src
        )
    {
        if (this != &_src)
        {
        self_type _copy (
            __move(self_type , _src)) ;

        using std::swap ;
        swap (
        static_cast<obj_alloc&>(*this),
        static_cast<obj_alloc&>(_copy)
            ) ;

        swap(this->_ptrs[_hptr], 
            _copy. _ptrs[_hptr]) ;
        swap(this->_ptrs[_tptr], 
            _copy. _ptrs[_tptr]) ;
        swap(this->_ptrs[_lptr], 
            _copy. _ptrs[_lptr]) ;
        }
        
        return ( *this ) ;
    }
    
/*-------------------------------------- push array alloc */
    __normal_call void_type set_alloc (
        size_type _new_alloc
        )
    {
        size_type _cur_count = count();
        size_type _cur_alloc = alloc();
        
        if (_new_alloc  < _cur_count )
        {   
    /*---------------------------------------- dec. alloc */
            set_count(_new_alloc,__cont::loose_alloc) ;
            
            this->_ptrs[_hptr] = 
            self_type::reallocate(this->_ptrs[_hptr], 
                _cur_alloc, _new_alloc, _new_alloc );
            this->_ptrs[_tptr] = 
            this->_ptrs[_hptr] + _new_alloc ;
            this->_ptrs[_lptr] = 
            this->_ptrs[_hptr] + _new_alloc ;
        }
        else
        if (_new_alloc != _cur_alloc )
        {
    /*---------------------------------------- inc. alloc */
            this->_ptrs[_hptr] = 
            self_type::reallocate(this->_ptrs[_hptr], 
                _cur_alloc, _new_alloc, _cur_count );
            this->_ptrs[_tptr] = 
            this->_ptrs[_hptr] + _cur_count ;
            this->_ptrs[_lptr] = 
            this->_ptrs[_hptr] + _new_alloc ;
        }
    }
    
/*-------------------------------------- trim array alloc */
    __inline_call void_type  fix_alloc (
        ) { set_alloc(count()) ; }

/*---------------------- increment count (_def construct) */
    __normal_call void_type  inc_count (
        size_type _inc_count ,
        __cont::alloc_types _new_alloc = 
        __cont::loose_alloc
        )
    {
    /*--------- expand underlying buffer, based on policy */
        size_type _cur_count = count() ;
        switch (_new_alloc )
        {
        case __cont::loose_alloc : 
        { inc_alloc(_cur_count + _inc_count); break; }
        case __cont::tight_alloc : 
        { set_alloc(_cur_count + _inc_count); break; }
        }
    /*--------- update buffer pointers and ctor new items */
        _write_it _head = tend() ; 
        this->_ptrs[_tptr] += _inc_count ;
        ctor_iter(_head , tend());
    }

/*---------------------- increment count (copy construct) */
    __normal_call void_type  inc_count (
        size_type _inc_count ,
        __cont::alloc_types _new_alloc ,
        data_type const& _data
        )
    {
    /*--------- expand underlying buffer, based on policy */
        size_type _cur_count = count() ;
        switch (_new_alloc )
        {
        case __cont::loose_alloc : 
        { inc_alloc(_cur_count + _inc_count); break; }
        case __cont::tight_alloc : 
        { set_alloc(_cur_count + _inc_count); break; }
        }
    /*--------- update buffer pointers and ctor new items */
        _write_it _head = tend();
        this->_ptrs[_tptr] += _inc_count;
        ctor_iter(_head , tend(), _data);
    }
    
/*--------------------------------------- decrement count */
    __normal_call void_type  dec_count (
        size_type _dec_count ,
        __cont::alloc_types _new_alloc = 
        __cont::loose_alloc
        )
    {
    /*---- update buffer pointers and _destruct old items */
        size_type _cur_count = count() ;
        _write_it _head = tend() - _dec_count;
        dtor_iter(_head , tend());
        this->_ptrs[_tptr] -= _dec_count;
    /*---- shrink underlying buffer, based on size policy */
        switch (_new_alloc )
        {
        case __cont::loose_alloc : 
        { /* do nothing - preserve buffer */  break; }
        case __cont::tight_alloc : 
        { set_alloc(_cur_count - _dec_count); break; }
        }
    }
    
/*---------------------------- set count (_def construct) */
    __inline_call void_type  set_count (
        size_type _new_count ,
        __cont::alloc_types _new_alloc = 
        __cont::loose_alloc
        )
    {
    /*----- manipulate the ctor'd range within the buffer */
        size_type _cur_count = count() ;
        if(_new_count > _cur_count)
        {
            inc_count(
           _new_count - _cur_count, _new_alloc) ;
        }
        else
        if(_new_count <=_cur_count)
        {
            dec_count(
           _cur_count - _new_count, _new_alloc) ;
        }
    }
    
/*---------------------------- set count (copy construct) */
    __inline_call void_type  set_count (
        size_type _new_count ,
        __cont::alloc_types _new_alloc ,
        data_type const&_data
        )
    {
    /*----- manipulate the ctor'd range within the buffer */
        size_type _cur_count = count() ;
        if(_new_count > _cur_count)
        {
            inc_count(
           _new_count - _cur_count, _new_alloc, _data);
        }
        else
        if(_new_count <=_cur_count)
        {
            dec_count(
           _cur_count - _new_count, _new_alloc) ;
        }
    }
    
/*-------------------------------- return container alloc */
    __inline_call allocator get_alloc(
        ) const
    {   return static_cast<allocator>(  *this ) ;
    }
    
/*-------------------------------- return container count */
    __inline_call size_type count (
        ) const
    { /* take pointer difference */
        return  ( size_type )
       (this->_ptrs[_tptr] - this->_ptrs[_hptr]);
    }
    
/*-------------------------------- return container alloc */
    __inline_call size_type alloc (
        ) const
    { /* take pointer difference */
        return  ( size_type )
       (this->_ptrs[_lptr] - this->_ptrs[_hptr]);
    }
    
/*-------------------------------- true if sequence empty */
    __inline_call bool_type empty (
        ) const
    { /* take pointer difference */
        return 
        this->_ptrs[_hptr]== this->_ptrs[_tptr] ;
    }
    
/*-------------------------------- const access iterators */
    __inline_call _const_it head (
        ) const
    {/*---- return iterator for list head */
        data_type *_ptr = 
            this->_ptrs[_hptr] != nullptr ?
            this->_ptrs[_hptr]  : nullptr ;
        self_type *_obj = (self_type *)this ;
        return _const_it(_ptr, _obj) ;
    }
    
    __inline_call _const_it tail (
        ) const
    {/*---- return iterator for list head */
        data_type *_ptr = 
            this->_ptrs[_tptr] != nullptr ?
            this->_ptrs[_tptr]-1: nullptr ;
        self_type *_obj = (self_type *)this ;
        return _const_it(_ptr, _obj) ;
    }
    
    __inline_call _const_it hend (
        ) const
    {/*------ return iterator "past" head */
        data_type *_ptr = 
            this->_ptrs[_hptr] != nullptr ?
            this->_ptrs[_hptr]-1: nullptr ;
        self_type *_obj = (self_type *)this ;
        return _const_it(_ptr, _obj) ;
    }
    
    __inline_call _const_it tend (
        ) const
    {/*------ return iterator "past" tail */
        data_type *_ptr = 
            this->_ptrs[_tptr] != nullptr ?
            this->_ptrs[_tptr]  : nullptr ;
        self_type *_obj = (self_type *)this ;
        return _const_it(_ptr, _obj) ;
    }   
    
/*-------------------------------- write access iterators */
    __inline_call _write_it head (
        )
    { /* return iterator for list head */
        data_type *_ptr = 
            this->_ptrs[_hptr] != nullptr ?
            this->_ptrs[_hptr]  : nullptr ;
        self_type *_obj = (self_type *)this ;
        return _write_it(_ptr, _obj) ;
    }
    
    __inline_call _write_it tail (
        )
    { /* return iterator for list head */
        data_type *_ptr = 
            this->_ptrs[_tptr] != nullptr ?
            this->_ptrs[_tptr]-1: nullptr ;
        self_type *_obj = (self_type *)this ;
        return _write_it(_ptr, _obj) ;
    }
    
    __inline_call _write_it hend (
        )
    { /* return iterator "past" head */
        data_type *_ptr = 
            this->_ptrs[_hptr] != nullptr ?
            this->_ptrs[_hptr]-1: nullptr ;
        self_type *_obj = (self_type *)this ;
        return _write_it(_ptr, _obj) ;
    }
    
    __inline_call _write_it tend (
        )
    { /* return iterator "past" tail */
        data_type *_ptr = 
            this->_ptrs[_tptr] != nullptr ?
            this->_ptrs[_tptr]  : nullptr ;
        self_type *_obj = (self_type *)this ;
        return _write_it(_ptr, _obj) ;
    }
        
/*-------------------------------------- "free" container */
    __inline_call void_type clear (
        __cont::alloc_types _this_alloc = 
        __cont::loose_alloc
        )
    {
        if (_this_alloc == __cont::loose_alloc)
            set_count(0, _this_alloc) ;
        else
            set_count(0, _this_alloc) ;
    }
    
/*---------------------------- push data (_def construct) */
    __inline_call size_type push_tail (
        )
    {/* default construct obj and increment count */
        size_type result = count();
        if (result != alloc())
        {
            self_type::construct(
            this->_ptrs[_tptr]++) ;
        }
        else
        {
            inc_alloc(result + 1) ;
            self_type::construct(
            this->_ptrs[_tptr]++) ;
        }
        return result;
    }
    
/*---------------------------- push data (copy construct) */
    __inline_call size_type push_tail (
        data_type const&_data
        )
    {/* copy data onto object and increment count */
        size_type result = count();
        if (result != alloc())
        {
            self_type::construct(
            this->_ptrs[_tptr]++,
                __copy(data_type,_data)) ;
        }
        else
        {
            data_type _temp(
                __copy(data_type,_data)) ;
                 
            inc_alloc(result + 1) ;
            self_type::construct(
            this->_ptrs[_tptr]++,
                __move(data_type,_temp)) ;
        }
        return result;
    }
    
/*---------------------------- push data (move construct) */
    __inline_call size_type push_tail (
        data_type &&_data
        )
    {/* move data onto object and increment count */
        size_type result = count();
        if (result != alloc())
        {
            self_type::construct(
            this->_ptrs[_tptr]++,
                __move(data_type,_data)) ;
        }
        else
        {
            data_type _temp(
                __move(data_type,_data)) ; 
            
            inc_alloc(result + 1) ;
            self_type::construct(
            this->_ptrs[_tptr]++,
                __move(data_type,_temp)) ;
        }
        return result;
    }
    
/*---------------------------- push data (copy obj range) */
    template <
    typename      iter_type
             >
    __normal_call void_type push_tail (
        iter_type _head,
        iter_type _tend
        )
    {/* push full sequence at tail and inc. count */
        for(; _head != _tend; ++_head) 
        {
        self_type::push_tail ( *_head);
        }
    }
    
/*--------------------------------------------- _pop data */
    __inline_call void_type _pop_tail (
        )
    { /* _destruct tail item and dec count */
        __assert ( count() != +0 
            && "array._pop_tail: empty") ;
        
        self_type::
        _destruct(--this->_ptrs[_tptr] ) ;
    }
    
/*--------------------------------------------- _pop data */
    __inline_call void_type _pop_tail (
        data_type &_data
        )
    { /* _destruct tail item and dec count */
        __assert ( count() != +0 
            && "array._pop_tail: empty") ;
        
        _data    = 
        std::move(*(this->_ptrs[_tptr]-1)) ;
        self_type::
        _destruct(--this->_ptrs[_tptr] ) ;
    }
    
/*------------------------------------ operator[] (write) */
    __inline_call data_type      &operator[] (
        size_type _pos
        )
    {
        __assert ( _pos >=     +0 && 
                   _pos < count() &&
            "array[]: out of range!" ) ;
        
        return *(this->_ptrs[_hptr] + _pos) ;
    }
    
/*------------------------------------ operator[] (const) */
    __inline_call data_type const&operator[] (
        size_type _pos
        ) const
    {
        __assert ( _pos >=     +0 && 
                   _pos < count() &&
            "array[]: out of range!" ) ;
        
        return *(this->_ptrs[_hptr] + _pos) ;
    }

    } ;

#   undef  __cont


    }

#   endif//__ARRAY__



