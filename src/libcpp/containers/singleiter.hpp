
/*
------------------------------------------------------------
 * "single" iterator -- for singly-linked sequences.
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

#   ifndef __SINGLE_ITER__
#   define __SINGLE_ITER__

#   include "iter_base.hpp"

    namespace containers {

    template <
    typename L,
    typename I
             >
    class single_iterator_base_ : 
        public containers::single_iterator_base
    { 
/*---------------------- iterator for singly-linked lists */
    public  :
    
    typedef L                           list_type ;
    typedef I                           iter_type ;
    
    typedef typename 
            list_type::data_type        data_type ;
    typedef typename 
            list_type::item_type        item_type ; 
    typedef typename 
            list_type::size_type        size_type ;
    typedef typename 
            list_type::diff_type        diff_type ;
    
    typedef single_iterator_base_ <
            list_type , 
            iter_type             >     self_type ;

    public  :

#   ifdef _DEBUG
    item_type* _ptr;    // list item
    list_type* _obj;    // list obj.
#   else
    item_type* _ptr;    // list item
#   endif

    public  :
    
    __inline_call single_iterator_base_ (
        item_type *_psrc = nullptr,
        list_type *_lsrc = nullptr
#   ifdef _DEBUG
        ) : _ptr(_psrc) , 
            _obj(_lsrc) {}
#   else
        ) : _ptr(_psrc) { __unreferenced(_lsrc) ; }
#   endif

    __inline_call~single_iterator_base_ (
        )                               = default ;
    __inline_call single_iterator_base_ (
        self_type const&_src
        )                               = default ;
    __inline_call single_iterator_base_ (
        self_type &&    _src
        )                               = default ;

    __inline_call 
        self_type & operator = (
        self_type const&_src
        )                               = default ;
    __inline_call 
        self_type & operator = (
        self_type &&    _src
        )                               = default ;

/*--------------------------------- translation operators */

    __inline_call iter_type &operator++ (
        )
    {
#   ifdef _DEBUG
        __assert( this->_ptr != nullptr &&
        "single_iterator: null link!" ) ;
#   endif
        this->_ptr = this->_ptr->next() ;
        return *(iter_type*)this ;
    }

    __inline_call 
        iter_type operator++ ( int )
    {
        iter_type _tmp(*(iter_type*)this);
        ++*this;
        return _tmp;
    }

    __inline_call iter_type &operator+= (
        size_type _off
        )
    {
        for ( ; _off > 0; --_off)++*this ;
        return *(iter_type*)this ;
    }

    __inline_call iter_type  operator + (
        size_type _off
        ) const
    {
        iter_type _tmp(*(iter_type*)this);
        return _tmp += _off;
    }

/*--------------------------------- comparative operators */

    __inline_call bool_type operator == (
        self_type const&_src
        ) const
    {
        return this->_ptr == _src._ptr ;
    }
    __inline_call bool_type operator != (
        self_type const&_src
        ) const
    {
        return this->_ptr != _src._ptr ;
    }

    } ;

    /*
    --------------------------------------------------------
     * "const" iterators for singly-linked sequences.
    --------------------------------------------------------
     */

    template <
    typename L
             > 
    class const_single_iterator: public 
            single_iterator_base_<L, 
            const_single_iterator<L> >
    {
/*---------------------- iterator for singly-linked lists */
    public  :
    
    typedef L                           list_type ;

    typedef const_single_iterator<
            list_type            >      self_type ;
    typedef single_iterator_base_< 
            list_type , 
            self_type            >      base_type ;

    typedef typename 
            base_type::data_type        data_type ;
    typedef typename 
            base_type::item_type        item_type ; 
    typedef typename 
            base_type::size_type        size_type ;
    typedef typename 
            base_type::diff_type        diff_type ;

    typedef data_type const&            const_ref ;
    typedef data_type const*            const_ptr ;

    public  :
    
    __inline_call const_single_iterator (
        item_type *_psrc = nullptr,
        list_type *_lsrc = nullptr
        ) : base_type(_psrc, _lsrc) {}

    __inline_call~const_single_iterator (
        )                              =  default ;
    __inline_call const_single_iterator (
        self_type const&_src
        )                              =  default ;
    __inline_call const_single_iterator (
        self_type &&    _src
        )                              =  default ;

    __inline_call 
        self_type & operator = (
        self_type const&_src
        )                              =  default ;
    __inline_call 
        self_type & operator = (
        self_type &&    _src
        )                              =  default ;

/*-------------------------------- "const" access to data */

    __inline_call const_ref operator * (
        ) const 
    {
#       ifdef _DEBUG
        __assert( this->_ptr != nullptr &&
        "single_iterator: null pointer!");
#       endif
        return this->_ptr->_data;
    }

    __inline_call const_ptr operator-> (
        ) const { return &**this ;   }

    __inline_call item_type const*item ( 
        ) const { return this->_ptr; }
    
    } ;

    /*
    --------------------------------------------------------
     * "write" iterators for singly-linked sequences.
    --------------------------------------------------------
     */

    template <
    typename L
             > 
    class write_single_iterator: public 
            single_iterator_base_<L, 
            write_single_iterator<L> >
    {
/*---------------------- iterator for singly-linked lists */
    public  :
    
    typedef L                           list_type ;

    typedef write_single_iterator<
            list_type            >      self_type ;
    typedef single_iterator_base_< 
            list_type , 
            self_type            >      base_type ;

    typedef typename 
            base_type::data_type        data_type ;
    typedef typename 
            base_type::item_type        item_type ; 
    typedef typename 
            base_type::size_type        size_type ;
    typedef typename 
            base_type::diff_type        diff_type ;

    typedef data_type &                 write_ref ;
    typedef data_type *                 write_ptr ;

    public  :
    
    __inline_call write_single_iterator (
        item_type *_psrc = nullptr,
        list_type *_lsrc = nullptr
        ) : base_type(_psrc, _lsrc) {}

    __inline_call~write_single_iterator (
        )                              =  default ;
    __inline_call write_single_iterator (
        self_type const&_src
        )                              =  default ;
    __inline_call write_single_iterator (
        self_type &&    _src
        )                              =  default ;

    __inline_call 
        self_type & operator = (
        self_type const&_src
        )                              =  default ;
    __inline_call 
        self_type & operator = (
        self_type &&    _src
        )                              =  default ;

/*-------------------------------- "write" access to data */

    __inline_call write_ref operator * (
        ) const 
    {
#       ifdef _DEBUG
        __assert( this->_ptr != nullptr &&
        "single_iterator: null pointer!");
#       endif
        return this->_ptr->_data;
    }

    __inline_call write_ptr operator-> (
        ) const { return &**this ;   }

    __inline_call item_type*item ( 
        ) const { return this->_ptr; }
    
    } ;


    }

#   endif   //__SINGLE_ITER__



