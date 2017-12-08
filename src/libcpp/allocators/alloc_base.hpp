
    /* 
    --------------------------------------------------------
     * basic allocator -- wrap "malloc" and friends.
    --------------------------------------------------------
     *
     * BASIC-ALLOC is a simple wrapper for "c"-style alloc.
     * via malloc and friends, intended to expose "realloc"
     * for c++ containers. This class deals only with bytes 
     * see _ITEM-ALLOC for the treatment of non-POD objects
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

#   ifndef __ALLOC_BASE___
#   define __ALLOC_BASE___

#   include <cstddef>
#   include <cstdlib>
#   include <new>

    namespace allocators {

    class   basic_alloc
	{
/*-------- raw byte allocator: wrap malloc, realloc, free */
	public  :
	
	typedef std::ptrdiff_t		        diff_type;
	typedef std::size_t 		        size_type;

    typedef basic_alloc                 self_type;

	public	:
    
    /*
    --------------------------------------------------------
     * ALLOCATE: do buffer allocation.
    --------------------------------------------------------
     */
          
	__static_call
	__inline_call char_type*   allocate (
		size_type _new_count
		)
	{ 
    /*---------------------- delegate to malloc & friends */
		char_type* _addr = 
	        (char_type*)std::malloc(_new_count) ;

    /*------------------------------ throw on alloc fails */
		if ( _new_count != +0 &&
             nullptr == _addr )
        throw std::bad_alloc ();

    /*------------------------------ return memory buffer */
		return ( _addr ) ;
	}
	
    /*
    --------------------------------------------------------
     * REALLOCATE: re-allocate buffer.
    --------------------------------------------------------
     */
     
	__static_call
	__inline_call char_type* reallocate (
        char_type*_addr ,
		size_type _old_count,
		size_type _new_count
		)
	{ 
        __unreferenced(_old_count) ;

    /*---------------------- delegate to malloc & friends */
		char_type* _nptr = 
	   (char_type*)std::realloc(_addr,_new_count) ;

    /*------------------------------ throw on alloc fails */
		if ( _new_count != +0 &&
             nullptr == _nptr )
        throw std::bad_alloc ();
    
    /*------------------------------ return memory buffer */
		return ( _nptr ) ;
	}
	
    /*
    --------------------------------------------------------
     * DEALLOCATE: de-allocate buffer.
    --------------------------------------------------------
     */
 
	__static_call
	__inline_call void_type  deallocate (
        char_type*_addr ,
		size_type _old_count
		)
	{ 
        __unreferenced(_old_count);

    /*---------------------- delegate to malloc & friends */
	    if (nullptr != _addr) 
            std::free( _addr) ;
	}

    } ;


    }

#   endif   //__ALLOC_BASE__



