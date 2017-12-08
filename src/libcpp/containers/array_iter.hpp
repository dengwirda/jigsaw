
/* 
------------------------------------------------------------
 * "array" iterator -- for contiguous sequences.
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

#   ifndef __ARRAY_ITER__
#   define __ARRAY_ITER__

#   include "iter_base.hpp"

    namespace containers {

    template <
    typename C ,
    typename I 
             >
    class array_iterator_base:
        public containers::random_iterator_base
    { 
/* random access iterator for contiguous array containers */
    public	:
    
	typedef C                           container ;
	typedef I                           iter_type ;

	typedef typename 
            container::data_type        data_type ;
	typedef typename 
            container::diff_type        diff_type ;
	typedef typename 
            container::size_type        size_type ;
	
    typedef array_iterator_base <
            container , 
            iter_type           >       self_type ;

    public	:

#   ifdef _DEBUG
    data_type* _ptr;    // array item
	container* _obj;    // array obj.
#   else
    data_type* _ptr;    // array item
#   endif

    public	:
    
	__inline_call array_iterator_base (
		data_type *_psrc = nullptr,
		container *_osrc = nullptr
#   ifdef _DEBUG
		) : _ptr(_psrc), 
			_obj(_osrc) {}
    #else
		) : _ptr(_psrc) { __unreferenced(_osrc) ; }
#   endif

    __inline_call~array_iterator_base() = default ;

    __inline_call array_iterator_base (
        self_type const&_src
        )                               = default ;
    __inline_call array_iterator_base (
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

	__inline_call iter_type& operator++ (
		)
	{ /*----------------- pre-increment */
		++this->_ptr;
		return *(iter_type*)this;
	}
	__inline_call iter_type& operator-- (
		)
	{ /*----------------- pre-decrement */
		--this->_ptr;
		return *(iter_type*)this;
	}

	__inline_call iter_type& operator+= (
		size_type _off
		)
	{ /*--------------------- increment */
		this->_ptr += _off ;
		return *(iter_type*)this;
	}
	__inline_call iter_type& operator-= (
		size_type _off
		)
	{ /*--------------------- decrement */
		this->_ptr -= _off ;
		return *(iter_type*)this;
	}

	__inline_call 
        iter_type operator++ ( int )
	{ /*---------------- post-increment */
		iter_type _tmp(*(iter_type*)this);
		++*this;
		return _tmp;
	}
	__inline_call 
        iter_type operator-- ( int )
	{ /*---------------- post-decrement */
		iter_type _tmp(*(iter_type*)this);
		--*this;
		return _tmp;
	}

	__inline_call iter_type operator + (
		size_type _off
		)
	{ /*-------------------- add offset */
		iter_type _tmp(*(iter_type*)this);
		return _tmp += _off;
	}
	__inline_call iter_type operator - (
		size_type _off
		)
	{ /*-------------------- sub offset */
		iter_type _tmp(*(iter_type*)this);
		return _tmp -= _off;
	}

	__inline_call diff_type operator - (
		self_type const& _src
		) const
	{ /*----------- iterator difference */ 
		return this->_ptr - _src._ptr; 
	}

/*--------------------------------- comparative operators */

	__inline_call bool_type operator== (
		self_type const& _src
		) const
	{ /*--------------------- equal to */
		return this->_ptr == _src._ptr ;
	}
	__inline_call bool_type operator!= (
		self_type const& _src
		) const
	{ /*----------------- not equal to */
		return this->_ptr != _src._ptr ;
	}

	__inline_call bool_type operator < (	
		self_type const& _src
		) const
	{ /*-------------------- less than */
		return this->_ptr <  _src._ptr ;
	}
	__inline_call bool_type operator > (
		self_type const& _src
		) const
	{ /*-------------------- more than */
		return this->_ptr >  _src._ptr ;
	}

	__inline_call bool_type operator<= (
		self_type const& _src
		) const
	{ /*-------- less than or equal to */
		return this->_ptr <= _src._ptr ;
	}
	__inline_call bool_type operator>= (
		self_type const& _src
		) const
	{ /*-------- more than or equal to */
		return this->_ptr >= _src._ptr ;
	}

    } ;

    /*
    --------------------------------------------------------
     * "const" iterators for contiguous sequences.
    --------------------------------------------------------
     */

    template <
    typename C
             >
    class const_array_iterator: public  
            array_iterator_base <C, 
            const_array_iterator<C> >
    { 
/* random access iterator for contiguous array containers */
    public	:
    
	typedef C                           container ;

	typedef const_array_iterator<
            container           >       self_type ;
	typedef array_iterator_base < 
            container , 
            self_type           >       base_type ;

	typedef typename 
            base_type::data_type        data_type ;
	typedef typename 
            base_type::diff_type        diff_type ;
	typedef typename 
            base_type::size_type        size_type ;

    public	:
    
	__inline_call const_array_iterator (
		data_type *_psrc = nullptr,
		container *_osrc = nullptr
		) : base_type(_psrc, _osrc) {}

	__inline_call~const_array_iterator()= default ;

    __inline_call const_array_iterator (
        self_type const&_src
        )                              =  default ;
    __inline_call const_array_iterator (
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

    __inline_call data_type const* operator-> (
		) const
	{   return &**this ;
	}

	__inline_call data_type const& operator[] (
		size_type _pos
		) const
	{   return *(*this + _pos);
	}

	__inline_call data_type const& operator * (
		) const
	{ /* return reference to underlying data */ 
#       ifdef _DEBUG
		__assert ( 
            this->_ptr != nullptr   && 
	    "const_array_iterator: null pointer!" ) ;
		__assert ( 
                this->_ptr >= 
            this->_obj->head()._ptr &&
		        this->_ptr <= 
            this->_obj->tail()._ptr &&
        "const_array_iterator: out of range!" ) ;
#       endif
		return ( *this->_ptr ); 
	}

    } ;

    /*
    ------------------------------------------------------------
     * "write" iterators for contiguous sequences.
    ------------------------------------------------------------
     */

    template <
    typename C
             >
    class write_array_iterator: public  
            array_iterator_base <C, 
            write_array_iterator<C> >
    { 
/* random access iterator for contiguous array containers */
    public	:
    
	typedef C                           container ;
	
	typedef write_array_iterator<
            container           >       self_type ;
	typedef array_iterator_base < 
            container , 
            self_type           >       base_type ;

	typedef typename 
            base_type::data_type        data_type ;
	typedef typename 
            base_type::diff_type        diff_type ;
	typedef typename 
            base_type::size_type        size_type ;

    public	:
	
	__inline_call write_array_iterator (
		data_type *_psrc = nullptr,
		container *_osrc = nullptr
		) : base_type(_psrc, _osrc) {}

	__inline_call~write_array_iterator()= default ;

    __inline_call write_array_iterator (
        self_type const&_src
        )                              =  default ;
    __inline_call write_array_iterator (
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

    __inline_call data_type* operator-> (
		) const
	{   return &**this ;
	}

	__inline_call data_type& operator[] (
		size_type _pos
		) const
	{   return *(*this + _pos) ;
	}

	__inline_call data_type& operator * (
		) const
	{ /* return reference to underlying data */ 
#       ifdef _DEBUG
		__assert ( 
            this->_ptr != nullptr   && 
	    "write_array_iterator: null pointer!" ) ;
		__assert ( 
                this->_ptr >= 
            this->_obj->head()._ptr &&
		        this->_ptr <= 
            this->_obj->tail()._ptr &&
        "write_array_iterator: out of range!" ) ;
#       endif
		return ( *this->_ptr ) ; 
	}

    } ;


    }

#   endif   //__ARRAY_ITER__



