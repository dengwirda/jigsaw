
    /* 
    --------------------------------------------------------
     * "item" allocator -- for general obj. types.
    --------------------------------------------------------
     *
     * _ITEM-ALLOC provides a "realloc"-able c++ allocator, 
     * disambiguating between POD and non-POD object types.
     * Only POD types benefit from calls through to realloc 
     * non-trivial objects induce alloc/copy/free patterns.
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

#   ifndef __ALLOC_ITEM___
#   define __ALLOC_ITEM___

#   include <type_traits>

    namespace allocators {
    
    
    template <
    typename D , 
    typename A
             > 
    class obj_alloc : public A
    { 
/*---------------- typed allocator wrapper for obj. types */
    public  :
    
    typedef D                           data_type ;
    typedef A                           allocator ;

    typedef typename 
            allocator::size_type        size_type ;
    typedef typename 
            allocator::diff_type        diff_type ;

    typedef obj_alloc <
            data_type ,
            allocator >                 self_type ;

    public  :
        
    /*
    --------------------------------------------------------
     * OBJ_ALLOC: construct allocator for ITEM obj.
    --------------------------------------------------------
     */
     
    __inline_call obj_alloc (
        allocator const&_asrc = allocator()
        ) : allocator(  _asrc )
    { /*------ construct base from _asrc */
    }

    __inline_call~obj_alloc ()   = default;
    
    __inline_call obj_alloc (
        self_type const& _src  
            )                    = default;
    __inline_call obj_alloc (
        self_type &&     _src  
            )                    = default;

    __inline_call 
        self_type& operator = (
        self_type const& _src  
            )                    = default;
    __inline_call 
        self_type& operator = (
        self_type  &&    _src  
            )                    = default;

    /*
    --------------------------------------------------------
     * ALLOCATE: do obj. buffer allocation.
    --------------------------------------------------------
     */
 
    __inline_call data_type*   allocate (
        size_type _new_count
        )
    {   return (data_type*)allocator::allocate (
            _new_count * sizeof(data_type));
    }
    
    /*
    --------------------------------------------------------
     * REALLOCATE: re-allocate obj. buffer.
    --------------------------------------------------------
     */
     
    __inline_call data_type* reallocate (
        data_type*_addr,
        size_type _old_alloc ,
        size_type _new_alloc ,
        size_type _old_count
        )
    {
    /*------------------------- allocate a new raw buffer */
        data_type *_rptr = 
            (data_type *)allocator:: allocate (
                _new_alloc * sizeof(data_type))  ;
    
    /*---- copy objects to new buffer and free old buffer */
        if (_addr != nullptr)
        {
            __write_ptr(data_type) _head = _addr ; 
            __write_ptr(data_type) _tail ; 
            __write_ptr(data_type) _item = _rptr ;
            
        /*--- number of ctor'd items to be moved to _rptr */
            size_type _new_count = std::min (
                _old_count, _new_alloc);
                
        /*-------------------- move items between buffers */
            _tail = _addr + _new_count;
            for ( ; _head != _tail; ++_head, 
                                    ++_item) 
            { 
            /*----- move ctor'd items from _addr to _rptr */
                construct(_item, std::move(*_head)) ;
                _destruct(_head) ;
            }
            _tail = _addr + _old_count;
            for ( ; _head != _tail; ++_head) 
            { 
            /*---- _destruct all remaining items in _addr */
                _destruct(_head) ;
            }
            
            allocator:: deallocate((char_type*)_addr, 
                _old_alloc * sizeof(data_type)) ;
        }
        
        return ( _rptr ) ;
    }
    
    /*
    --------------------------------------------------------
     * DEALLOCATE: de-allocate obj. buffer.
    --------------------------------------------------------
     */
     
    __inline_call void_type  deallocate (
        data_type*_addr,
        size_type _old_alloc
        )
    { 
    /*----------------------------- free non-null buffers */
        if (_addr != nullptr)
        allocator::deallocate (
            (char_type*)_addr , 
                _old_alloc * sizeof(data_type)) ;
    }
    
    
    /*
    --------------------------------------------------------
     * _DESTRUCT: destory allocated object.
    --------------------------------------------------------
     */
     
    __static_call
    __inline_call void_type _destruct (
        data_type* _xptr
        ) 
    { 
        __unreferenced(_xptr);
        
        _xptr->~data_type () ; 
    }
    
    /*
    --------------------------------------------------------
     * CONSTRUCT: regular object construction.
    --------------------------------------------------------
     */
     
    __static_call
    __inline_call void_type construct ( // _def c'tor obj.
        data_type* _xptr
        ) 
    {   new(_xptr) data_type ; 
    }
    
    __static_call
    __inline_call void_type construct ( // copy c'tor obj.
        data_type* _xptr,
        data_type const&_xval
        ) 
    {   new( _xptr)data_type (
            __copy(data_type, _xval)) ;
    }
    
    __static_call
    __inline_call void_type construct ( // move c'tor obj.
        data_type* _xptr,
        data_type&&_xval
        ) 
    {   new( _xptr)data_type (
            __move(data_type, _xval)) ;
    }
    
    template <
    typename...   data_part
             >
    __static_call
    __inline_call void_type construct (
        data_type  * _xptr,
        data_part const&... _part
        )
    {   new(_xptr)data_type(_part...) ;
    }
    
    } ;
    
    
    template <
    typename D , 
    typename A 
             > 
    class pod_alloc : public A
    { 
/*---------------------- buffer allocator for "pod" types */
    public  :
    
    typedef D                           data_type ;
    typedef A                           allocator ;

    typedef typename 
            allocator::size_type        size_type ;
    typedef typename 
            allocator::diff_type        diff_type ;

    typedef pod_alloc <
            data_type ,
            allocator >                 self_type ;

    public  :
        
    /*
    --------------------------------------------------------
     * POD_ALLOC: construct allocator for POD. obj.
    --------------------------------------------------------
     */
     
    __inline_call pod_alloc (
        allocator const&_asrc = allocator()
        ) : allocator(  _asrc )
    { /*------ construct base from _asrc */
    }

    __inline_call~pod_alloc ()   = default;
    
    __inline_call pod_alloc (
        self_type const& _src  
            )                    = default;
    __inline_call pod_alloc (
        self_type &&     _src  
            )                    = default;

    __inline_call 
        self_type& operator = (
        self_type const& _src  
            )                    = default;
    __inline_call 
        self_type& operator = (
        self_type  &&    _src  
            )                    = default;

    /*
    --------------------------------------------------------
     * ALLOCATE: do obj. buffer allocation.
    --------------------------------------------------------
     */
     
    __inline_call data_type*   allocate (
        size_type _new_alloc
        )
    {   return (data_type *)allocator::  allocate (
            _new_alloc * sizeof(data_type)) ;
    }
    
    /*
    --------------------------------------------------------
     * REALLOCATE: re-allocate obj. buffer.
    --------------------------------------------------------
     */
     
    __inline_call data_type* reallocate (
        data_type*_addr,
        size_type _old_alloc,
        size_type _new_alloc,
        size_type 
        )
    {   return (data_type *)allocator::reallocate (
               (char_type *)_addr ,
            _old_alloc * sizeof(data_type),
            _new_alloc * sizeof(data_type)) ; 
    }
    
    /*
    --------------------------------------------------------
     * DEALLOCATE: de-allocate obj. buffer.
    --------------------------------------------------------
     */
 
    __inline_call void_type  deallocate (
        data_type*_addr,
        size_type _old_alloc
        )
    {   if (_addr != nullptr) 
        allocator:: deallocate((char_type*) _addr , 
            _old_alloc * sizeof(data_type)) ;
    }
    
    /*
    --------------------------------------------------------
     * _DESTRUCT: destory allocated object.
    --------------------------------------------------------
     */
     
    __static_call
    __inline_call void_type _destruct (
        data_type*
        ) { /* do nothing */ }
    
    /*
    --------------------------------------------------------
     * CONSTRUCT: regular object construction.
    --------------------------------------------------------
     */
     
    __static_call
    __inline_call void_type construct ( // _def c'tor obj.
    __write_ptr  (data_type)
        ) {/* do nothing */}
        
    __static_call
    __inline_call void_type construct ( // copy c'tor obj.
    __write_ptr  (data_type) _xptr,
        data_type const&_xval
        ) { *_xptr = _xval ; }
        
    __static_call
    __inline_call void_type construct ( // move c'tor obj.
    __write_ptr  (data_type) _xptr,
        data_type  &&_xval
        ) { *_xptr = _xval ; }

    template <
    typename...   data_part
             >
    __static_call
    __inline_call void_type construct (
        data_type  * _xptr,
        data_part const&... _part
        )
    {   new(_xptr)data_type(_part...) ;
    }
    
    } ;


    template <
    typename D , 
    typename A , 
    bool     F 
             >
    class  loc_alloc : public  obj_alloc <D, A>
    { 
/*----- specialisation to derive from _POD if flag = TRUE */
    public  :
    
    typedef D                   data_type ;
    typedef A                   allocator ;

    typedef obj_alloc <
            data_type , 
            allocator >         base_type ;

    typedef loc_alloc <
            data_type ,
            allocator , F >     self_type ;

    public  :

    /*
    --------------------------------------------------------
     * LOC_ALLOC: construct obj. from _ASRC.
    --------------------------------------------------------
     */
     
    __inline_call loc_alloc (
        allocator const&_asrc = allocator()
        ) : base_type(  _asrc )
    { /*------ construct base from _asrc */
    }

    __inline_call~loc_alloc ()   = default;
    
    __inline_call loc_alloc (
        self_type const& _src  
            )                    = default;
    __inline_call loc_alloc (
        self_type &&     _src  
            )                    = default;

    __inline_call 
        self_type& operator = (
        self_type const& _src  
            )                    = default;
    __inline_call 
        self_type& operator = (
        self_type  &&    _src  
            )                    = default;

    } ;
    

    template <
    typename D , 
    typename A
             >
    class loc_alloc<D, A, true> : 
                    public pod_alloc<D, A>
    {
/*----- specialisation to derive from _POD if flag = TRUE */
    public  :
    
    typedef D                   data_type ;
    typedef A                   allocator ;

    typedef pod_alloc <
            data_type , 
            allocator >         base_type ;

    typedef loc_alloc <
            data_type ,
            allocator , 
            true      >         self_type ;

    public  :

    /*
    --------------------------------------------------------
     * LOC_ALLOC: construct obj. from _ASRC.
    --------------------------------------------------------
     */

    __inline_call loc_alloc (
        allocator const&_asrc = allocator()
        ) : base_type(  _asrc )
    { /*------ construct base from _asrc */
    }

    __inline_call~loc_alloc ()   = default;
    
    __inline_call loc_alloc (
        self_type const& _src  
            )                    = default;
    __inline_call loc_alloc (
        self_type &&     _src  
            )                    = default;

    __inline_call 
        self_type& operator = (
        self_type const& _src  
            )                    = default;
    __inline_call 
        self_type& operator = (
        self_type  &&    _src  
            )                    = default;

    } ;
    
    
    template <
    typename D , 
    typename A = basic_alloc
             >
    class _item_alloc : public  loc_alloc <
                                D , A , 
                    std::is_pod<D>::value >
    { 
/*-------------- type-aware item alloc - public interface */
    public  :
    
    typedef D                   data_type ;
    typedef A                   allocator ;

    typedef loc_alloc   <
            data_type   , 
            allocator   , 
    std::is_pod<D>::value>      base_type ;

    typedef _item_alloc <
            data_type   ,
            allocator   >       self_type ;

    public  :

    /*
    --------------------------------------------------------
     * _ITEM_ALLOC: construct obj. from _ASRC.
    --------------------------------------------------------
     */

    __inline_call _item_alloc (
        allocator const&_asrc = allocator()
        ) : base_type ( _asrc)
    { /*------ construct base from _asrc */
    }

    __inline_call~_item_alloc () = default;
    
    __inline_call _item_alloc (
        self_type const& _src  
            )                    = default;
    __inline_call _item_alloc (
        self_type &&     _src  
            )                    = default;

    __inline_call 
        self_type& operator = (
        self_type const& _src  
            )                    = default;
    __inline_call 
        self_type& operator = (
        self_type  &&    _src  
            )                    = default;

    } ;

    
    }

#   endif   //__ALLOC_ITEM__



