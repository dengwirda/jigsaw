
/*
------------------------------------------------------------
 * "fixed-array" -- contiguous static buffer.
------------------------------------------------------------
 *
 * FIXED-ARRAY is a simple wrapper over a "c"-style
 * statically allocated array. Array size is fixed. It
 * provides O(1) random access, and that's about it!
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
------------------------------------------------------------
 *
 * Last updated: 03 May, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __FIXED_ARRAY__
#   define __FIXED_ARRAY__

#   include "array_iter.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename    D,
    std::size_t L
             >
    class fixed_array
    {
/*------------ simple wrapper for "c-style" static arrays */
    public  :
    typedef D                           data_type ;

    typedef std::   size_t              size_type ;
    typedef std::ptrdiff_t              diff_type ;

    size_type static const _size = L > +0 ?
                                   L : +1 ;

    typedef __cont::fixed_array <
                data_type, L    >       self_type ;

    typedef __cont::const_array_iterator <
                self_type       >       _const_it ;
    typedef __cont::write_array_iterator <
                self_type       >       _write_it ;

    private  :

    data_type    _data[ _size]  ;

    private :

/*-------------------------------- helper - copy sequence */
    __normal_call void_type copy_data (
        size_type _SIZE,
        data_type const& _dsrc
        )
    {
        for (auto _item  =  head() ;
                  _SIZE != +0 ;
                --_SIZE, ++_item )
           *_item = __copy(data_type, _dsrc) ;
    }

    template <
    typename      iter_type
             >
    __normal_call void_type copy_data (
        iter_type _head,
        iter_type _tend
        )
    {
        for (auto _item  =  head() ;
                  _head != _tend;
                ++_item, ++_head )
           *_item = __copy(data_type,*_head) ;
    }

/*-------------------------------- helper - copy sequence */
    template <
    typename      iter_type
             >
    __normal_call void_type copy_iter (
        iter_type _head,
        iter_type _tend,
    __cont::base_iterator_kind
        ) { copy_data(_head, _tend) ; }

    __inline_call void_type copy_iter (
        size_type _SIZE,
        data_type const& _dsrc ,
    __cont::null_iterator_kind
        ) { copy_data(_SIZE, _dsrc) ; }

    public  :

/*--------------------------- default c'tor - do nothing! */
    __inline_call fixed_array (
        ) {}

/*--------------------------- default c'tor - initialisor */
    __inline_call fixed_array (
        size_type _SIZE,
        data_type const&_dsrc = data_type ()
        )
    {   copy_data(_SIZE, _dsrc) ;
    }

/*--------------------------- default c'tor - initialisor */
    template <
        typename  iter_type
             >
    __inline_call fixed_array (
        iter_type _head,
        iter_type _tail
        )
    {
        copy_iter(_head,_tail ,
            __cont::iter_kind(_head)) ;
    }

/*--------------------------- default d'tor/copy/move ops */
    __inline_call~fixed_array () = default ;

    __inline_call fixed_array (
        self_type const& _src
            )                    = default ;
    __inline_call fixed_array (
        self_type &&     _src
            )                    = default ;

    __inline_call
        self_type& operator = (
        self_type const& _src
            )                    = default ;
    __inline_call
        self_type& operator = (
        self_type  &&    _src
            )                    = default ;

/*-------------------------------------------- fill array */
    __inline_call void_type fill (
        data_type const& _dsrc = data_type()
        )
    {   copy_data(_size, _dsrc) ;
    }

/*-------------------------------------------- operator[] */
    __inline_call data_type      &operator[] (
        size_type _pos
        )
    { /*------ index into underlying array */
        __assert ( _pos >=   +0 &&
                   _pos < _size &&
        "fixed_array[]: out of range!") ;

        return this->_data[_pos] ;
    }

    __inline_call data_type const&operator[] (
        size_type _pos
        ) const
    { /*------ index into underlying array */
        __assert ( _pos >=   +0 &&
                   _pos < _size &&
        "fixed_array[]: out of range!") ;

        return this->_data[_pos] ;
    }

/*--------------------------------------- container count */
    __inline_call size_type count (
        ) const { return ( _size ); }

/*--------------------------------------- container alloc */
    __inline_call size_type alloc (
        ) const { return ( _size ); }

/*--------------------------------------- write iterators */
    __inline_call _write_it head (
        )
    { /* return iterator to head of array */
        size_type  _off = 0;
        self_type *_obj = (self_type *)this ;
        return _write_it(this->_data + _off, _obj);
    }

    __inline_call _write_it tail (
        )
    { /* return iterator to tail of array */
        size_type  _off = _size - 1;
        self_type *_obj = (self_type *)this ;
        return _write_it(this->_data + _off, _obj);
    }

    __inline_call _write_it hend (
        )
    { /* return iterator at head - 1 */
        size_type _off = 1;
        self_type *_obj = (self_type *)this ;
        return _write_it(this->_data + _off, _obj);
    }

    __inline_call _write_it tend (
        )
    { /* return iterator at tail + 1 */
        size_type  _off = _size - 0;
        self_type *_obj = (self_type *)this ;
        return _write_it(this->_data + _off, _obj);
    }

/*--------------------------------------- const iterators */
    __inline_call _const_it head (
        ) const
    { /* return iterator to head of array */
        size_type  _off = 0;
        self_type *_obj = (self_type *)this ;
        return _const_it(this->_data + _off, _obj);
    }

    __inline_call _const_it tail (
        ) const
    { /* return iterator to tail of array */
        size_type  _off = _size - 1;
        self_type *_obj = (self_type *)this ;
        return _const_it(this->_data + _off, _obj);
    }

    __inline_call _const_it hend (
        ) const
    { /* return iterator at head - 1 */
        size_type _off = 1;
        self_type *_obj = (self_type *)this ;
        return _const_it(this->_data - _off, _obj);
    }

    __inline_call _const_it tend (
        ) const
    { /* return iterator at tail + 1 */
        size_type  _off = _size - 0;
        self_type *_obj = (self_type *)this ;
        return _const_it(this->_data + _off, _obj);
    }

    } ;

#   undef  __cont


    }

#   endif   //__FIXED_ARRAY__



