
    /* 
    --------------------------------------------------------
     * "pool'd" allocator -- "node"-based containers...
    --------------------------------------------------------
     *
     * _POOL-ALLOC provides "pool'd" allocation for "small"
     * objects, suitable for use with "node"-based 
     * containers, including lists, hash-tables, etc. 
     * Nothing too tricky(!), just a linked list of memory 
     * slabs, a "free" list of deallocated objects, and a 
     * fall-back onto the base allocator for alloc.'s of a 
     * larger size.
     *
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
     * Last updated: 02 May, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ALLOC_POOL___
#   define __ALLOC_POOL___

#   include <algorithm>
#   include <cstring>

#   include "alloc_base.hpp"

    namespace allocators {

    template <
    typename    A = basic_alloc ,
    std::size_t S = 64 * 1024
             >
    class  _pool_alloc : public A
	{ 
/*--- local memory pool: delegate to base alloc otherwise */
	public	:
	
	typedef A				            base_type ;

    std::size_t static constexpr pool_size  = S ;

	typedef typename 
            base_type::size_type	    size_type ;
	typedef typename 
            base_type::diff_type	    diff_type ;

	typedef _pool_alloc <
		     base_type  , 
             pool_size  >               self_type ;

	private	:
	
/*-------------------------------- local allocation lists */
	char_type*  _alloc ;        // pointer to block
	char_type*	_shift ;        // pointer to block 
	char_type*	_block ;
	char_type*	_cache ;        // cached item list
	
/*-------------------------------- local allocation sizes */
	size_type   _item_size ;    // size of alloc.
	size_type   _slab_size ;    // size of buffer

    /*
    --------------------------------------------------------
     * local helper: access encoded list pointers.
    --------------------------------------------------------
     */
 
	#define __nextblock(__block)    \
		*((char **)(__block  +      \
		this->_slab_size*this->_item_size))
	
	#define __nextcache(__cache)    \
			* (char **)(__cache)
	
	#define __slabbytes			    \
	    this->_slab_size *		    \
		this->_item_size + sizeof (char *)
		
    /*
    --------------------------------------------------------
     * local helper: allocate and push new buffer.
    --------------------------------------------------------
     */
 
    __inline_call void_type push_block ( // push new buffer
		)
	{
	/*------------------------- allocate new buffer item  */
		char_type *_this_slab = 
	        (char_type *)base_type
	            ::allocate(__slabbytes);
	            
	/*------------------------- push new buffer onto list */
		__nextblock(_this_slab)= this->_block;
		
		this->_alloc = _this_slab + 
		    this->_slab_size*this->_item_size;
		this->_block = _this_slab ;
		this->_shift = _this_slab ;
	}
	
    /*
    --------------------------------------------------------
     * local helper: _pop top buffer and de-alloc.
    --------------------------------------------------------
     */
     
	__inline_call void_type _pop_block ( // _pop top buffer
		)
	{
	/*------------------------- _pop top buffer from list */
		char_type *_tail_slab = this->_block ;
		this->_block =
		    __nextblock (_tail_slab) ;
		    
	/*------------------------- free top buffer data/item */
		base_type::deallocate(
		    _tail_slab, __slabbytes) ;
	}

	public	:
	
    /*
    --------------------------------------------------------
     * _POOL_ALLOC: construct pool'd allocator from base.
    --------------------------------------------------------
     */
     
	__inline_call _pool_alloc ( // default c'tor
		size_type const _size = +1 ,
		base_type const&_asrc = base_type()
		) : base_type(  _asrc), 
	/*--------------- construct with empty alloc. buffers */
			_alloc(nullptr) , 
			_shift(nullptr) ,
			_block(nullptr) , 
			_cache(nullptr)
	{ 
	/*--------------- construct to force alloc on request */
		this->_item_size = 
			std::max(_size, (size_type)sizeof(char*)) ;
		this->_slab_size = 
			(size_type)pool_size /  this->_item_size;
		this->_slab_size = 
			std::max((size_type)+1, this->_slab_size) ;
	}
	
	__inline_call~_pool_alloc ( // default d'tor
		)
	{ 
	/*--------------- _free items while base alloc. valid */
		clear() ;
	}
	
    /*
    --------------------------------------------------------
     * CLEAR: de-alloc. and free the allocator.
    --------------------------------------------------------
     */
     
	__normal_call void_type clear (
		)
	{
	/*-------------- _destruct and deallocate all buffers */
		for ( ; this->_block != nullptr; ) _pop_block() ;
		
	/*-------------- set in "empty" default c'tor'd state */
		this->_block = nullptr ;
		this->_cache = nullptr ;
		this->_alloc = nullptr ;
		this->_shift = nullptr ;
	}
	
    /*
    --------------------------------------------------------
     * return count/stat.'s for alloc.
    --------------------------------------------------------
     */
 
    __normal_call size_type nfree (
        ) const
    {
    /*------------------------------- traverse cache list */
        size_type  _size = +0 ;
        char_type* _cptr = this->_cache ;
        for ( ; _cptr != nullptr; )
        {
            _cptr  = __nextcache(_cptr) ;
            _size +=  +1 ;
        }
        return ( _size ) ;
    }

    __inline_call size_type nslab (
        ) const
    {
    /*------------------------------- traverse block list */
        size_type  _size = +0 ;
        char_type* _bptr = this->_block ;
        for ( ; _bptr != nullptr; )
        {
            _bptr  = __nextblock(_bptr) ;
            _size +=  +1 ;
        }
        return ( _size ) ;
    }

    /*
    --------------------------------------------------------
     * ALLOCATE: do buffer allocation.
    --------------------------------------------------------
     */
 
	__inline_call char_type*   allocate (
		size_type _new_count
		)
	{
	/*--------------------- catch any "empty" allocations */
		if (_new_count <= +0 ) return ( nullptr ) ;
	
	/*--------------------- deal with genuine allocations */
		char_type *_nptr;
		if (_new_count <= this->_item_size)
		{
			if((_nptr = this->_cache)!= nullptr )								
			{ 
			/*--------------- _pop cached item from stack */
				this->_cache = __nextcache(_nptr);
			}
			else
			{ 
			/*---------------------- alloc new item block */
				if (this->_shift >= this->_alloc) 
					push_block() ;
				_nptr = this->_shift ;
				
			/*---------------------- inc. offset in block */ 
				this->_shift += this->_item_size ;
			}
		}
		else
		{   
    /*------------------------ delegate to base allocator */
			_nptr = base_type
			    ::allocate(_new_count);
		}
		
		return ( _nptr );
	}
	
    /*
    --------------------------------------------------------
     * REALLOCATE: re-allocate buffer.
    --------------------------------------------------------
     */
 
	__inline_call char_type* reallocate (
        char_type*_addr ,
		size_type _old_count,
		size_type _new_count
		)
	{
	/*------------------ catch any first-time allocations */
		if (_addr == nullptr) 
		    return self_type::allocate(_new_count);
	
	/*------------------ deal with genuine reallocations  */
		char_type *_nptr;
		if (_old_count <= this->_item_size)    
		{
    /*------------------ current alloc is from local pool */
		if (_new_count <= this->_item_size)
		{ 
		/*--- keep current, block is already large enough */
			_nptr = _addr ;
		}
		else
		{ 
		/*--- deallocate back to pool, allocate from base */
			_nptr = 
			base_type::allocate( _new_count);
			
			size_type _nnum = std::min(
			    _new_count, _old_count) ;
			std::memcpy(_nptr, _addr, _nnum);
			
			self_type::
                deallocate(_addr,_old_count);
		}
		}
		else                                    
		{
    /*------------------ current alloc is from base alloc */
		if (_new_count > this->_item_size)
		{ 
		/*--- re-allocate within base, too large for pool */
			_nptr = base_type::reallocate(
			    _addr, _old_count, _new_count) ;
		}
		else
		{ 
		/*--- deallocate back to base, allocate from pool */
			_nptr = 
			self_type::allocate( _new_count);
			
			size_type _nnum = std::min(
			    _new_count, _old_count) ;
			std::memcpy(_nptr, _addr, _nnum);
			
			base_type::
                deallocate(_addr,_old_count);
		}
		}
		
		return ( _nptr ) ;
	}
	
    /*
    --------------------------------------------------------
     * DEALLOCATE: de-allocate buffer.
    --------------------------------------------------------
     */
     
    __inline_call void_type  deallocate (
        char_type*_addr ,
		size_type _old_count
		)
	{
		if (_old_count <= this->_item_size)
		{ 
			if (_addr != nullptr)
			{
    /*------------------------------ cache un-used buffer */
				__nextcache(_addr) = 
				this->_cache ;
				this->_cache = _addr ;
			}
		}
		else
		{
    /*------------------------------ deallocate from base */
			base_type::
            deallocate(_addr, _old_count) ;
		}
	}

    #undef  __nextblock
	#undef  __nextcache
    #undef  __slabbytes

	} ;
    
    
    }    //} allocators

#   endif//__ALLOC_POOL__



