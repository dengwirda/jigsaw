
    /* 
    --------------------------------------------------------
     * "wrap" an allocator -- "stateful" allocators...
    --------------------------------------------------------
     *
     * _WRAP-ALLOC is a wrapper for so-called "stateful" 
     * allocators, allowing, for example, _POOL-ALLOC to be 
     * used within containers without duplication. It just 
     * maintains a pointer-to-allocator internally, and 
     * marshals calls through to the base.
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

#   ifndef __ALLOC_WRAP___
#   define __ALLOC_WRAP___

    namespace allocators {


    template <
    typename A
             > 
    class  _wrap_alloc
	{ 
/*-------- pointer to base allocator: stateful allocators */
	public	:
	typedef A				                allocator ;

    typedef _wrap_alloc<allocator>	        self_type ;

	typedef typename 
            allocator::size_type	        size_type ;
	typedef typename 
            allocator::diff_type	        diff_type ;

	private :

    allocator     *_aptr ;

	public	:
	
    /*
    --------------------------------------------------------
     * _WRAP_ALLOC: construct alloc. wrapper.
    --------------------------------------------------------
     */
 
	__inline_call _wrap_alloc (
		allocator *_asrc = nullptr
		) : _aptr( _asrc ) {}

    __inline_call~_wrap_alloc () = default ;
    
    __inline_call _wrap_alloc (
        self_type const& _src  
            )                    = default ;
    __inline_call _wrap_alloc (
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
                
    /*
    --------------------------------------------------------
     * GET/PUT_ALLOCATOR: access to base obj.
    --------------------------------------------------------
     */
 
    __inline_call 
        allocator      & put_allocator (
        ) 
    {   
        return *this->_aptr ;
    }
    
    __inline_call 
        allocator const& get_allocator (
        )const
    {   
        return *this->_aptr ;
    }
    
    /*
    --------------------------------------------------------
     * ALLOCATE: do buffer allocation.
    --------------------------------------------------------
     */
     
	__inline_call char_type*   allocate (
		size_type _new_count
		) const
	{
		__assert( this->_aptr != nullptr &&
	    "_wrap_alloc: null allocator!" ) ;

		return this->
            _aptr->allocate( _new_count) ;
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
		) const
	{
		__assert( this->_aptr != nullptr &&
	    "_wrap_alloc: null allocator!" ) ;

		return this->_aptr->reallocate(_addr, 
		        _old_count, 
                _new_count) ;
	}
	
    /*
    --------------------------------------------------------
     * DEALLOCATE: de-allocate buffer.
    --------------------------------------------------------
     */
	
	__inline_call void_type  deallocate (
        char_type*_addr ,
		size_type _old_count
		) const
	{
		__assert( this->_aptr != nullptr &&
	    "_wrap_alloc: null allocator!" ) ;

	    this->_aptr->
            deallocate(_addr,_old_count) ;
	}
	
	} ;
    
    
    }

#   endif   //__ALLOC_WRAP__



