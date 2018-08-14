
/*
------------------------------------------------------------
 * a dynamically-allocated, "block-wise" array object.
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

#   ifndef __BLOCK_ARRAY__
#   define __BLOCK_ARRAY__

#   include "block_iter.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename D , 
    typename A = allocators::basic_alloc 
             >
    class block_array
    {
/*------ a dynamically allocated, block-wise array object */
    public  :
    typedef D                           data_type ;
    typedef A                           allocator ;

    typedef __cont::block_array <
                data_type , 
                allocator >             self_type ;

    typedef typename 
            allocator::size_type        size_type ;
    typedef typename 
            allocator::diff_type        diff_type ;

    typedef __cont::const_block_iterator <
                self_type >             _const_it ;
    typedef __cont::write_block_iterator <
                self_type >             _write_it ;
    
    typedef __cont::array <
                data_type , 
                allocator >             leaf_type ;
    typedef __cont::array <
                leaf_type ,
                allocator >             root_type ;

    size_type static const _sizt = sizeof(data_type) ;
    size_type static const _sizb = (64*1024)/_sizt;    
    size_type static const _size = (_sizb>4)?_sizb:4 ;

    private :

    root_type           _block ;
    size_type           _count ;
    size_type           _alloc ;

    /*
    --------------------------------------------------------
     * INC-ALLOC: helper - broaden storage.
    --------------------------------------------------------
     */

    __normal_call void_type inc_alloc (
        size_type _new_count
        )
    {
        if (_new_count <= alloc()) return ;
    /*-------------------------- alloc. storage in blocks */
        size_type _new_alloc = alloc();
        if (_new_count <= _size)
        { /* grow allocation as multiples of two */
            _new_alloc = 
            std::max(_new_alloc * 2, _new_count);
        }
        else
        { /* round to the nearest block boundary */
            _new_alloc = 
          ((_new_count / _size) + 1) * _size    ;
        }
        set_alloc(_new_alloc) ;
    }
    
    /*
    --------------------------------------------------------
     * COPY-ITER: helper - construct range.
    --------------------------------------------------------
     */
    
    template <
    typename      iter_type 
             >
    __normal_call void_type copy_iter (
        iter_type _head,
        iter_type _tail,
    __cont::base_iterator_kind
        )
    { /* copy range onto object */
        for(; _head != _tail; ++_head)
        {
            push_tail(*_head) ;
        }
    }

    __inline_call void_type copy_iter (
        size_type _size,
        data_type const& _dsrc,
    __cont::null_iterator_kind
        )
    { /* copy _data onto object */
        set_count(_size, 
            __cont::tight_alloc, _dsrc) ;
    }

    public  :
    
    /*
    --------------------------------------------------------
     * Container c'tor's/d'tor's/assignment op's, etc
    --------------------------------------------------------
     */

    __inline_call block_array (
        allocator const&_asrc = allocator()
        ) : _block( _asrc)
    {
    /*----------------------- default c'tor - do nothing! */
        this->_count = +0;
        this->_alloc = +0;
    }

    __inline_call block_array (
        size_type _size,
        data_type const&_dsrc = data_type(),
        allocator const&_asrc = allocator()
        ) : _block( _asrc)
    {
    /*----------------------- default c'tor - initialisor */
        this->_count = +0;
        this->_alloc = +0;
        copy_iter(_size,_dsrc,
        __cont::null_iterator_kind ()) ;
    }

    template <
    typename  iter_type
             >
    __inline_call block_array (
        iter_type _head,
        iter_type _tail,
        allocator const&_asrc = allocator()
        ) : _block( _asrc)
    {
    /*----------------------- default c'tor - initialisor */
        this->_count = +0;
        this->_alloc = +0;
        copy_iter(_head,_tail,
            __cont::iter_kind( _head)) ;
    }

    __inline_call~block_array (
        )
    /*---------------------------------------- _def d'tor */
    { set_count(0,__cont::tight_alloc) ; 
    }

/*-------------------------------------------- copy c'tor */
    __inline_call block_array (
        self_type const&_src
        )                       = default ;
/*-------------------------------------------- move c'tor */
    __inline_call block_array (
        self_type &&    _src
        ) : _block(std::move(_src._block) )
    {
        this->_alloc = 
        _src. _alloc ;    
        _src. _alloc = +0u ;
        
        this->_count = 
        _src. _count ;
        _src. _count = +0u ;
    }

/*-------------------------------------------- copy a-op. */
    __inline_call 
        self_type & operator = (
        self_type const&_src
        )                       = default ;
/*-------------------------------------------- move a-op. */
    __inline_call 
        self_type & operator = (
        self_type &&    _src
        )
    {
        if (this != &_src)
        {
        using std::swap ;
        self_type _copy (
                __move(self_type , _src)) ;

        swap(this->_block, 
            _copy. _block) ;
        swap(this->_count, 
            _copy. _count) ;
        swap(this->_alloc, 
            _copy. _alloc) ;
        }
        
        return ( *this ) ;
    }

    /*
    --------------------------------------------------------
     * SET-ALLOC: push new array alloc.
    --------------------------------------------------------
     */

    __normal_call void_type set_alloc (
        size_type _new_alloc
        )
    {
        size_type _cur_count = count();
        size_type _cur_alloc = alloc();
        
        size_type _cur_block = +0;
        size_type _new_block = +0;
        size_type _end_block = +0;
        size_type _end_shift = +0;

        size_type _siz_block = _size ;
        size_type _nil_block = +0    ;

        if (_new_alloc  < _cur_count )
        {
        /*--------------- _destruct un-needed range first */
            set_count(_new_alloc, 
                __cont::loose_alloc) ;
        }

        if (_cur_alloc!=+0)
        {
        /*-------------- calc. block-wise storage offsets */
            _cur_block = 
           (_cur_alloc - 1) / _size + 0 ;
        }
        if (_new_alloc!=+0)
        {
            _new_block = 
           (_new_alloc - 1) / _size + 0 ;
            _end_block = 
           (_new_alloc - 1) / _size + 1 ;
            _end_shift = 
           (_new_alloc - 1) % _size + 1 ;
        }

        if (_new_alloc  < _cur_alloc )
        { 
        /*--------------- _pop old blocks and realloc end */
            size_type _pos_block ;
            for(_pos_block=_cur_block; 
                _pos_block>_new_block; 
              --_pos_block )
            {
            this->_block [_pos_block].
                set_alloc(_nil_block);
            }           
            this->_block [_new_block].
                set_alloc(_end_shift);

            this->_block.set_count (
                _end_block, __cont::loose_alloc) ;

            this->_alloc =_new_alloc ;
        }
        else
        if (_new_alloc  > _cur_alloc )
        { 
        /*--------------- push new blocks and realloc end */
            this->_block.set_count (
                _end_block, __cont::loose_alloc, 
            leaf_type(this->_block.get_alloc())) ;   

            size_type _pos_block ;
            for(_pos_block=_cur_block; 
                _pos_block<_new_block; 
              ++_pos_block )
            {
            this->_block [_pos_block].
                set_alloc(_siz_block);
            }
            this->_block [_new_block].
                set_alloc(_end_shift);

            this->_alloc =_new_alloc ;
        }
    }
    
    /*
    --------------------------------------------------------
     * FIX-ALLOC: trim current alloc to count.
    --------------------------------------------------------
     */

    __inline_call void_type  fix_alloc (
        ) { set_alloc(count()) ; }

    /*
    --------------------------------------------------------
     * INC-COUNT: increment container count.
    --------------------------------------------------------
     */

    __normal_call void_type  inc_count (   // _def-construct
        size_type _inc_count ,
        __cont::alloc_types _new_alloc = 
        __cont::loose_alloc
        )
    {
        size_type _cur_count = count() ;
        size_type _new_count = count() + 
                  _inc_count ;
        
        size_type _cur_block = +0;
        size_type _new_block = +0;
        size_type _end_shift = +0;

        size_type _siz_block = _size ;
    
    /*------------------ calc. block-wise storage offsets */
        if (_cur_count!=+0)
        {
            _cur_block = 
           (_cur_count - 1) / _size + 0;
        }
        if (_new_count!=+0)
        {
            _new_block = 
           (_new_count - 1) / _size + 0;
            _end_shift = 
           (_new_count - 1) % _size + 1;
        }

    /*--------- expand underlying buffer, based on policy */
        switch (_new_alloc )
        {
        case __cont::loose_alloc : 
        { inc_alloc(_cur_count + _inc_count); break ; }
        case __cont::tight_alloc : 
        { set_alloc(_cur_count + _inc_count); break ; }
        }
    /*--------- update block structure and ctor new items */
        for(size_type _pos_block = _cur_block; 
                      _pos_block < _new_block; 
                    ++_pos_block )
        {
            this->_block [_pos_block].set_count(
                _siz_block, __cont::loose_alloc) ;
        }
        if (_new_count!=+0)
        {
            this->_block [_new_block].set_count(
                _end_shift, __cont::loose_alloc) ;
        }

        this->_count =_new_count ;
    }

    /*
    --------------------------------------------------------
     * INC-COUNT: increment container count.
    --------------------------------------------------------
     */

    __normal_call void_type  inc_count (   // copy-construct
        size_type _inc_count ,
        __cont::alloc_types _new_alloc ,
        data_type const&_data
        )
    {
        size_type _cur_count = count() ;
        size_type _new_count = count() + 
                  _inc_count ;
        
        size_type _cur_block = +0;
        size_type _new_block = +0;
        size_type _end_shift = +0;

        size_type _siz_block = _size ;
        
    /*------------------ calc. block-wise storage offsets */
        if (_cur_count!=+0)
        {
            _cur_block = 
           (_cur_count - 1) / _size + 0;
        }
        if (_new_count!=+0)
        {
            _new_block = 
           (_new_count - 1) / _size + 0;
            _end_shift = 
           (_new_count - 1) % _size + 1;
        }

    /*--------- expand underlying buffer, based on policy */
        switch (_new_alloc )
        {
        case __cont::loose_alloc : 
        { inc_alloc(_cur_count + _inc_count); break ; }
        case __cont::tight_alloc : 
        { set_alloc(_cur_count + _inc_count); break ; }
        }
    /*--------- update block structure and ctor new items */
        for(size_type _pos_block = _cur_block ;
                      _pos_block < _new_block ;
                    ++_pos_block )
        {
            this->_block [_pos_block].set_count(
                _siz_block, __cont::loose_alloc, _data) ;
        }
        if (_new_count!=+0)
        {
            this->_block [_new_block].set_count(
                _end_shift, __cont::loose_alloc, _data) ;
        }

        this->_count =_new_count ;
    }

    /*
    --------------------------------------------------------
     * DEC-COUNT: decrement container count.
    --------------------------------------------------------
     */

    __normal_call void_type  dec_count (
        size_type _dec_count,
        __cont::alloc_types _new_alloc = 
        __cont::loose_alloc
        )
    {
        size_type _cur_count = count() ;
        size_type _new_count = count() - 
                  _dec_count ;
        
        size_type _nil_block = +0;
        
        size_type _cur_block = +0;
        size_type _new_block = +0;
        size_type _end_shift = +0;
        
    /*------------------ calc. block-wise storage offsets */
        if (_cur_count!=+0)
        {
            _cur_block = 
           (_cur_count - 1) / _size + 0 ;
        }
        if (_new_count!=+0)
        {
            _new_block = 
           (_new_count - 1) / _size + 0 ;
            _end_shift = 
           (_new_count - 1) % _size + 1 ;
        }

    /*--------- update block structure and ctor new items */
        for(size_type _pos_block = _cur_block; 
                      _pos_block > _new_block; 
                    --_pos_block )
        {
            this->_block [_pos_block].set_count(
                _nil_block, __cont::loose_alloc) ;
        }
        if (_cur_count!=+0)
        {
            this->_block [_new_block].set_count(
                _end_shift, __cont::loose_alloc) ;
        }

        this->_count =_new_count ;

    /*--------- shrink underlying buffer, based on policy */
        switch (_new_alloc )
        {
        case __cont::loose_alloc : 
        { /* do nothing - preserve buffer */  break ; }
        case __cont::tight_alloc : 
        { set_alloc(_cur_count - _dec_count); break ; }
        }
    }

    /*
    --------------------------------------------------------
     * SET-COUNT: set new container count.
    --------------------------------------------------------
     */

    __inline_call void_type  set_count (   // _def-construct
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
    
    __inline_call void_type  set_count (   // copy-construct
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

    /*
    --------------------------------------------------------
     * GET-STATS: return container statistics.
    --------------------------------------------------------
     */

    __inline_call allocator get_alloc ( // return alloc obj.
        ) const
    { return static_cast<allocator>(  *this ) ;
    }
    
    __inline_call size_type count (// return container count
        ) const { return this->_count ; }
        
    __inline_call size_type alloc (// return container alloc
        ) const { return this->_alloc ; }
        
    __inline_call bool_type empty (// true if sequence empty
        ) const { return this->_count == +0 ; }
        
    /*
    --------------------------------------------------------
     * (const.-access) container iterators.
    --------------------------------------------------------
     */
     
    __inline_call _const_it head (
        ) const
    {/*--- return iterator for list head */
        size_type  _off = +0;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _const_it(_obj, _ptr, _off);
    }
    
    __inline_call _const_it tail (
        ) const
    {/*--- return iterator for list tail */
        size_type  _off =  this->_count -1;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _const_it(_obj, _ptr, _off);
    }
    
    __inline_call _const_it hend (
        ) const
    {/*----- return iterator "past" head */
        size_type  _off = -1;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _const_it(_obj, _ptr, _off);
    }
    
    __inline_call _const_it tend (
        ) const
    {/*----- return iterator "past" tail */
        size_type  _off =  this->_count -0;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _const_it(_obj, _ptr, _off);
    }

    /*
    --------------------------------------------------------
     * (write.-access) container iterators.
    --------------------------------------------------------
     */

    __inline_call _write_it head (
        )
    {/*--- return iterator for list head */
        size_type  _off = +0;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _write_it(_obj, _ptr, _off);
    }
    
    __inline_call _write_it tail (
        )
    {/*--- return iterator for list tail */
        size_type  _off =  this->_count -1;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _write_it(_obj, _ptr, _off);
    }
    
    __inline_call _write_it hend (
        )
    {/*----- return iterator "past" head */
        size_type  _off = -1;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _write_it(_obj, _ptr, _off);
    }
    
    __inline_call _write_it tend (
        )
    {/*----- return iterator "past" tail */
        size_type  _off =  this->_count -0;
        root_type *_ptr =(root_type*)&this->_block;
        self_type *_obj =(self_type*) this;
        return _write_it(_obj, _ptr, _off);
    }

    /*
    --------------------------------------------------------
     * CLEAR: "free" container storage.
    --------------------------------------------------------
     */

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
    
    /*
    --------------------------------------------------------
     * PUSH-TAIL: append data to tail .
    --------------------------------------------------------
     */
    
    __inline_call size_type push_tail ( // _def-construct
        )
    {/* _def construct object and increment count */
        size_type result = count();
        size_type offset = count()/_size ;  

        if (result == alloc()) 
        inc_alloc(result + 1);

        this-> _block[offset].push_tail();
        this-> _count += + 1 ;

        return result ;
    }
    
    __inline_call size_type push_tail ( // copy-construct
        data_type const&_data
        )
    {/* copy data onto object and increment count */
        size_type result = count();
        size_type offset = count()/_size ;

        if (result == alloc()) 
        inc_alloc(result + 1);

        this-> _block[offset].
        push_tail(__copy(data_type,_data)) ;
        this-> _count += + 1 ;

        return result ;
    }
    
    __inline_call size_type push_tail ( // move-construct
        data_type &&  _data
        )
    {/* copy data onto object and increment count */
        size_type result = count();
        size_type offset = count()/_size ;  

        if (result == alloc()) 
        inc_alloc(result + 1);

        this-> _block[offset].
        push_tail(__move(data_type,_data)) ;
        this-> _count += + 1 ;

        return result ;
    }
    
    template <
    typename      iter_type
             >
    __normal_call void_type push_tail ( // copy obj-range
        iter_type _head,
        iter_type _tend
        )
    {/* push full sequence at tail and inc. count */
        for(; _head != _tend; ++_head) 
        {
        self_type::push_tail ( *_head);
        }
    }
    
    /*
    --------------------------------------------------------
     * _POP-TAIL: erase data about tail.
    --------------------------------------------------------
     */
    
    __inline_call void_type _pop_tail (
        )
    { /* _destruct tail item and dec count */
        dec_count(+1, __cont::loose_alloc) ;        
    }
    
    __inline_call void_type _pop_tail (
        data_type &_data
        )
    { /* _destruct tail item and dec count */
        _data  = __move(data_type,*tail()) ;
        dec_count(+1, __cont::loose_alloc) ;
    }

    /*
    --------------------------------------------------------
     * OPERATOR[]: index into container.
    --------------------------------------------------------
     */

    __inline_call data_type      &operator[] (  // write
        size_type _pos
        )
    { /*------------ subscript operator */
        __assert ( _pos >=   +0   && 
                   _pos < count() &&
        "::block_array[]: out of range!" );
        
        return this->_block[_pos / _size]
                           [_pos % _size] ;
    }
    
    __inline_call data_type const&operator[] (  // const
        size_type _pos
        ) const
    { /*------------ subscript operator */
        __assert ( _pos >=   +0   && 
                   _pos < count() &&
        "::block_array[]: out of range!" );
        
        return this->_block[_pos / _size]
                           [_pos % _size] ;
    }

    } ;

#   undef  __cont


    }

#endif  //__BLOCK_ARRAY__



