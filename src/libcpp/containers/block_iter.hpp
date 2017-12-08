
/* 
------------------------------------------------------------
 * "block" iterator -- for non-contiguous sequences.
------------------------------------------------------------
 *
 * BLOCK-ARRAY is a block-wise array implementation, 
 * where the linear array storage is broken across 
 * non-contiguous "blocks" in memory. It is otherwise 
 * equivalent to a standard array interface. The blocked 
 * data-structure can be helpful for very large 
 * dynamically-resizable arrays, since memory is only 
 * re-alloc'd on the block level. Peformance asymptotics 
 * are equivalent to the standard arry object, with 
 * larger constant factors. 
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

#   ifndef __BLOCK_ITER__
#   define __BLOCK_ITER__

#   include "iter_base.hpp"

    namespace containers {

    template <
    typename C ,
    typename I 
             >
    class block_iterator_base:
        public containers::random_iterator_base
    { 
/* random access iterator for block-wise array containers */
    public	:
	
	typedef C                           container ;
	typedef I                           iter_type ;

	typedef typename 
            container::data_type        data_type ;
	typedef typename 
            container::size_type        size_type ;
	typedef typename 
            container::diff_type        diff_type ;
    typedef typename 
            container::root_type        root_type ;
    typedef typename 
            container::leaf_type        leaf_type ;
	
    typedef block_iterator_base <
            container , 
            iter_type           >       self_type ;

    public	:

#   ifdef _DEBUG
    root_type *_ptr;    // block ptr.
    size_type  _off;    // offset
	container *_obj;    // block obj.
#   else
    root_type *_ptr;    // block ptr.
    size_type  _off;    // offset
#   endif

    public	:
    
	__inline_call block_iterator_base (
		container *_bsrc = nullptr,
		root_type *_psrc = nullptr, 
		size_type  _osrc = +0
#   ifdef _DEBUG
		) : _ptr(_psrc) , 
			_off(_osrc) , 
			_obj(_bsrc) {}
#   else
		) : _ptr(_psrc) , 
			_off(_osrc) { __unreferenced(_bsrc) ; }
#   endif

	__inline_call~block_iterator_base() = default ;

    __inline_call block_iterator_base (
        self_type const&_src
        )                               = default ;
    __inline_call block_iterator_base (
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
		++this->_off;
		return *(iter_type*)this;
	}
	__inline_call iter_type& operator-- (
		)
	{ /*----------------- pre-decrement */
		--this->_off;
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

	__inline_call iter_type& operator+= (
		size_type _off
		)
	{ /*--------------------- increment */
		this->_off += _off;
		return *(iter_type*)this;
	}
	__inline_call iter_type& operator-= (
		size_type _off
		)
	{ /*--------------------- decrement */
		this->_off -= _off;
		return *(iter_type*)this;
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
		self_type const&_src
		) const
	{ /*----------- iterator difference */
		return this->_off - _src._off ; 
	}

/*--------------------------------- comparative operators */

	__inline_call bool_type operator== (
		self_type const&_src
		) const
	{ /*--------------------- equal to */
		return this->_off == _src._off ;
	}
	__inline_call bool_type operator!= (
		self_type const&_src
		) const
	{ /*----------------- not equal to */
		return this->_off != _src._off ;
	}
	__inline_call bool_type operator < (
		self_type const&_src
		) const
	{ /*-------------------- less than */
		return this->_off <  _src._off ;
	}
	__inline_call bool_type operator > (
		self_type const&_src
		) const
	{ /*-------------------- more than */
		return this->_off >  _src._off ;
	}
	__inline_call bool_type operator<= (
		self_type const&_src
		) const
	{ /*-------- less than or equal to */
		return this->_off <= _src._off ;
	}
	__inline_call bool_type operator>= (
		self_type const&_src
		) const
	{ /*-------- more than or equal to */
		return this->_off >= _src._off ;
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
    class const_block_iterator: public  
            block_iterator_base <C, 
            const_block_iterator<C> >
    { 
/* random access iterator for block-wise array containers */
    public	:
    
	typedef C                           container ;

	typedef const_block_iterator<
            container           >       self_type ;
	typedef block_iterator_base < 
            container , 
            self_type           >       base_type ;

    typedef typename 
            base_type::data_type        data_type ;
	typedef typename 
            base_type::size_type        size_type ;
	typedef typename 
            base_type::diff_type        diff_type ;
    typedef typename 
            base_type::root_type        root_type ;
    typedef typename 
            base_type::leaf_type        leaf_type ;

    public	:
    
    __inline_call const_block_iterator (
		container *_bsrc = nullptr,
		root_type *_psrc = nullptr, 
		size_type  _osrc = 0
		) : base_type(_bsrc, _psrc, _osrc) {}

	__inline_call~const_block_iterator()= default ;

    __inline_call const_block_iterator (
        self_type const&_src
        )                              =  default ;
    __inline_call const_block_iterator (
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
	{   return &**this;
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
            this->_ptr != nullptr &&
	    "const_block_iterator: null pointer!" ) ;
		__assert ( 
            this->_off >= 0 && 
            this->_off < this->_obj->count() &&
        "const_block_iterator: out of range!" ) ;
#       endif
		return (*this->_ptr)
            [this->_off / container::_size]
            [this->_off % container::_size] ;
	}

    } ;

    /*
    --------------------------------------------------------
     * "write" iterators for contiguous sequences.
    --------------------------------------------------------
     */

    template <
    typename C
             >
    class write_block_iterator: public  
            block_iterator_base <C, 
            write_block_iterator<C> >
    { 
/* random access iterator for block-wise array containers */
    public	:
    
	typedef C                           container ;

	typedef write_block_iterator<
            container           >       self_type ;
	typedef block_iterator_base < 
            container , 
            self_type           >       base_type ;

    typedef typename 
            base_type::data_type        data_type ;
	typedef typename 
            base_type::size_type        size_type ;
	typedef typename 
            base_type::diff_type        diff_type ;
    typedef typename 
            base_type::root_type        root_type ;
    typedef typename 
            base_type::leaf_type        leaf_type ;

    public	:
    
    __inline_call write_block_iterator (
		container *_bsrc = nullptr,
		root_type *_psrc = nullptr, 
		size_type  _osrc = 0
		) : base_type(_bsrc, _psrc, _osrc) {}

	__inline_call~write_block_iterator()= default ;

    __inline_call write_block_iterator (
        self_type const&_src
        )                              =  default ;
    __inline_call write_block_iterator (
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
	{   return &**this;
	}

	__inline_call data_type& operator[] (
		size_type _pos
		) const
	{   return *(*this + _pos);
	}

	__inline_call data_type& operator * (
		) const
	{ /* return reference to underlying data */ 
#       ifdef _DEBUG
		__assert ( 
            this->_ptr != nullptr &&
	    "write_block_iterator: null pointer!" ) ;
		__assert ( 
            this->_off >= 0 && 
            this->_off < this->_obj->count() &&
        "write_block_iterator: out of range!" ) ;
#       endif
		return (*this->_ptr)
            [this->_off / container::_size]
            [this->_off % container::_size] ;
	}

    } ;


    }

#   endif   //__BLOCK_ITER__



