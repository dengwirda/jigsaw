
/*
------------------------------------------------------------
 * an "n-ary"-heap based priority queue.
------------------------------------------------------------
 *
 * PRIORITY-IDX utilises an "n-ary" heap data-structure,
 * designed to maintain sorted precedence in mutable
 * collections. This "n-ary" heap simply generalises
 * a standard binary heap (i.e. a 2-heap) to trees with
 * "n" children per level.
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
 * Last updated: 15 Jun., 2022
 *
 * Copyright 2013-2022
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __PRIORITY_IDX__
#   define __PRIORITY_IDX__

#   include "array.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename T ,
    typename I
             >
    class _idx_pair
    {
/*----------- local data-pair for an indexed "n"-ary heap */
    public  :
    typedef T       data_type ;
    typedef I       iptr_type ;

    public  :
    data_type       _data ; // heap data
    iptr_type       _kidx ; // keys idx.
    } ;

#   define  D _idx_pair<T, K>

    template <
    typename T ,
    typename K ,
    typename P = std::less<T> ,
    typename A = allocators::basic_alloc ,
    typename C = containers::array <D,A>
             >
    class priorityidx
    {
/*----------- a dynamic priority queue as an "n"-ary heap */
    public  :

    #undef  D

    typedef T                           data_type ;
    typedef K                           iptr_type ;
    typedef C                           container ;
    typedef P                           pred_type ;
    typedef A                           allocator ;

    typedef containers::priorityidx <
            data_type ,
            iptr_type ,
            pred_type ,
            allocator ,
            container   >               self_type ;

    typedef typename
            container::data_type        pair_type ;
    typedef typename
            container::size_type        size_type ;

    typedef typename
            container::_write_it        _write_it ;
    typedef typename
            container::_const_it        _const_it ;

    typedef containers::array <
            iptr_type ,
            allocator         >         iptr_list ;

    size_type static constexpr _nfan = +8 ; // fan out

    public  :

    container           _heap ;

    iptr_list           _keys ;

    pred_type           _pred ;

    public  :

/*------------- helper - set "null" item marker as maxint */
    __static_call
    __inline_call iptr_type const null (
        ) 
    {   return std::numeric_limits<iptr_type>::max() ;
    }

    private :

/*------------- helper - push "hole" into sorted position */
    __normal_call _write_it push_upper (
        _write_it _head,
        _write_it _ipos,
        data_type const&_data
        )
    {
        for (; _ipos != _head ; )
        {
        /*----------------------- find position of parent */
            _write_it _ppos =
            _head + (_ipos-_head-1) / _nfan;
        /*----------------------- swap parent with "hole" */
            if (this->_pred(_data,_ppos->_data))
            {
            /*------------------ update heap-keys mapping */
                _ipos->_data =
                    std::move(_ppos->_data);
                _ipos->_kidx =
                    std::move(_ppos->_kidx);
            /*------------------ update keys-heap mapping */
                 this->_keys [_ipos->_kidx] =
                  (iptr_type)(_ipos- _head) ;
            /*------------------ traverse to upper levels */
                _ipos =_ppos ;
            }
            else break ;
        }

        return (_ipos) ;
    }

/*------------- helper - push "hole" into sorted position */
    __normal_call _write_it push_lower (
        _write_it _head,
        _write_it _tail,
        _write_it _ipos,
        data_type const&_data
        )
    {
        if (_head == _tail) return _ipos ;
    /*-- find position of root of last partial (sub)-tree */
        _write_it _cpos ;
        _write_it _imin ;
        _write_it _iend =  _head;
        if (_tail - _head > 0)
        {
            _iend =
            _head + (_tail-_head-1) / _nfan ;
        }
    /*-- push "hole" into sorted position on lower levels */
        for ( ; _ipos < _iend; )
        {
        /*---------------------------- pos of right child */
            _cpos =  _head +
            _nfan * (_ipos - _head) + _nfan ;
        /*---------------------------- find minimum child */
            _imin = _cpos;
            switch (_nfan)
            {
            case 8 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 7 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 6 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 5 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 4 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 3 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 2 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            }
        /*---------------------- update heap-keys mapping */
            _ipos->_data =
                std::move(_imin->_data);
            _ipos->_kidx =
                std::move(_imin->_kidx);
        /*---------------------- update keys-heap mapping */
             this->_keys [_ipos->_kidx] =
              (iptr_type)(_ipos- _head) ;
        /*---------------------- move onto next child pos */
            _ipos =_imin ;
        }

    /*-- deal with special case - last partial sub - tree */
        if (_ipos == _iend)
        {
            size_type _inum =
            _tail-_head-_nfan*(_ipos-_head);
            if (_inum > +0)
            {
        /*---------------------------- pos of right child */
            _cpos =  _head +
            _nfan * (_ipos -_head) + _inum ;
        /*---------------------------- find minimum child */
            _imin = _cpos;
            switch (_inum)
            {
            case 8 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 7 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 6 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 5 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 4 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 3 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            case 2 :
            if (this->_pred((--_cpos)->_data ,
                               _imin ->_data))
                _imin=_cpos;  // falls through
            }
        /*---------------------- update heap-keys mapping */
            _ipos->_data =
                std::move(_imin->_data);
            _ipos->_kidx =
                std::move(_imin->_kidx);
        /*---------------------- update keys-heap mapping */
             this->_keys [_ipos->_kidx] =
              (iptr_type)(_ipos- _head) ;
        /*---------------------- move onto next child pos */
            _ipos =_imin ;
            }
        }

    /*--- push "hole" into final position sorted position */
        return push_upper(_head, _ipos, _data) ;
    }

    public  :

/*--------------------------- default c'tor - do nothing! */
    __inline_call priorityidx (
        allocator const&_asrc = allocator(),
        pred_type const&_psrc = pred_type()
        ) : _heap(_asrc),
            _keys(_asrc),
            _pred(_psrc)   {}

/*--------------------------- default c'tor - initialisor */
    __inline_call priorityidx (
        pred_type const&_psrc
        ) : _heap(allocator()),
            _keys(allocator()),
            _pred(_psrc)   {}

    __inline_call~priorityidx() = default ;

    __inline_call priorityidx (
        self_type const&_src
        )                       = default ;
    __inline_call priorityidx (
        self_type &&    _src
        )                       = default ;

    __inline_call
        self_type & operator= (
        self_type const&_src
        )                       = default ;
    __inline_call
        self_type & operator= (
        self_type &&    _src
        )                       = default ;

/*------------------------------------- peek at heap root */

    __inline_call data_type const& root (
        ) const
    {   return this->_heap[  +0 ]._data ;
    }

/*------------------------------------- peek at heap item */
    __inline_call data_type const& peek (
        size_type  _hpos   = +0
        ) const
    {   return this->_heap[_hpos]._data ;
    }

/*------------------------------------- peek at item idx. */
    __inline_call iptr_type const  keys (
        size_type  _kpos   = +0
        ) const
    {   return this->_keys[_kpos] ;
    }

/*-------------------------------------------- empty test */
    __inline_call bool_type empty (
        ) const
    { return ( this->_heap. empty() ) ;
    }

/*-------------------------------------------- heap count */
    __inline_call size_type count (
        ) const
    { return ( this->_heap. count() ) ;
    }

/*-------------------------------------------- heap alloc */
    __inline_call size_type alloc (
        ) const
    { return ( this->_heap. alloc() ) ;
    }

/*-------------------------------------------- _set alloc */
    __inline_call size_type set_alloc (
        size_type _asiz
        )
    { this->_heap.set_alloc( _asiz) ;
    }

/*------------------------------------- (const) iterators */
    __inline_call _const_it head (
        ) const
    {   return this->_heap. head() ;
    }
    __inline_call _const_it tail (
        ) const
    {   return this->_heap. tail() ;
    }
    __inline_call _const_it hend (
        ) const
    {   return this->_heap. hend() ;
    }
    __inline_call _const_it tend (
        ) const
    {   return this->_heap. tend() ;
    }

    /*
    --------------------------------------------------------
     * PUSH: push data onto heap
    --------------------------------------------------------
     */

    __inline_call void_type push ( // copy
        iptr_type _kpos,
        data_type const&_data
        )
    {/*---------------- make room for new item in mapping */
        iptr_type _size = 
       (iptr_type) this->_keys .count() ;
        if (_size <= _kpos)
            this->_keys.set_count( 
                _kpos + 1 , 
        containers::loose_alloc, null());

        if (this->_keys[_kpos] != null())
            throw std::invalid_argument(
        "priorityidx.push: key in use") ;

    /*----------------- push _data onto tail of container */
        size_type _tpos =
            this->_heap.push_tail()  ;
    /*----------------- sort corresponding "hole" to pos. */
        _write_it _ipos = push_upper (
            this->_heap.head() ,
            this->_heap.head() + _tpos ,
        __copy(data_type, _data)) ;

    /*------------------------- copy new data into "hole" */
        _ipos->_kidx =    _kpos   ;
        _ipos->_data =
        __copy(data_type, _data)  ;

    /*------------------------- map _item to _data "hole" */
        this->_keys[_kpos] = 
            (iptr_type)(
            _ipos - this->_heap.head()) ;
    }

    __inline_call void_type push ( // move
        iptr_type _kpos,
        data_type && _data
        )
    {/*---------------- make room for new item in mapping */
        iptr_type _size = 
       (iptr_type) this->_keys .count() ;
        if (_size <= _kpos)
            this->_keys.set_count( 
                _kpos + 1 , 
        containers::loose_alloc, null());

        if (this->_keys[_kpos] != null())
            throw std::invalid_argument(
        "priorityidx.push: key in use") ;

    /*----------------- push _data onto tail of container */
        size_type _tpos =
            this->_heap.push_tail()  ;
    /*----------------- sort corresponding "hole" to pos. */
        _write_it _ipos = push_upper (
            this->_heap.head() ,
            this->_heap.head() + _tpos ,
        __copy(data_type, _data)) ;

    /*------------------------- move new data into "hole" */
        _ipos->_kidx =    _kpos   ;
        _ipos->_data =
        __move(data_type, _data)  ;

    /*------------------------- map _item to _data "hole" */
        this->_keys[_kpos] =
            (iptr_type)(
            _ipos - this->_heap.head()) ;
    }

    /*
    --------------------------------------------------------
     * _POP: _pop data from heap
    --------------------------------------------------------
     */

    __inline_call void_type _pop_root (
        )
    {
    /*---------------------------- _pop root, cache temp. */
        data_type  _temp ;
        _pop( +0,  _temp);
    }
    __inline_call void_type _pop_root (
        iptr_type &_kpos ,
        data_type &_data
        )
    {
    /*---------------------------- _pop item, return data */
        size_type _hpos =   +0 ;
        _kpos =    this->_heap[_hpos]._kidx ;
        _data =
        std::move (this->_heap[_hpos]._data);
    /*----------------------- push "hole" to lower levels */
        if (this->_heap.head() + _hpos !=
            this->_heap.tail() )
        {
    /*----------- sort "hole" at root to updated position */
        _write_it _ipos =  push_lower (
            this->_heap.head() ,
            this->_heap.tail() - 1    ,
            this->_heap.head() + _hpos,
            this->_heap.tail()-> _data) ;

    /*----------- copy current tail into updated position */
        _ipos->_data = std::move(
            this->_heap.tail()-> _data) ;
        _ipos->_kidx = std::move(
            this->_heap.tail()-> _kidx) ;

    /*---------------------- update mapping for tail item */
        this->_keys[_ipos->_kidx] =
            (iptr_type)(
            _ipos  - this->_heap.head());
        }
        this->_heap._pop_tail();

        this->_keys[_kpos] = null() ;
    }

    __inline_call void_type _pop (
        iptr_type  _kpos
        )
    {
    /*---------------------------- _pop data, cache temp. */
         data_type _temp ;
        _pop(_kpos,_temp);
    }
    __normal_call void_type _pop (
        iptr_type  _kpos ,
        data_type &_data
        )
    {
        if (this->_keys[_kpos] == null())
            return ;
    /*---------------------------- _pop item, return data */
        size_type _hpos = this->_keys[_kpos];
        _data =
        std::move (this->_heap[_hpos]._data);
    /*----------------------- push "hole" to lower levels */
        if (this->_heap.head() + _hpos !=
            this->_heap.tail() )
        {
    /*----------- sort "hole" at root to updated position */
        _write_it _ipos =  push_lower (
            this->_heap.head() ,
            this->_heap.tail() - 1    ,
            this->_heap.head() + _hpos,
            this->_heap.tail()-> _data) ;

    /*----------- copy current tail into updated position */
        _ipos->_data = std::move(
            this->_heap.tail()-> _data) ;
        _ipos->_kidx = std::move(
            this->_heap.tail()-> _kidx) ;

    /*---------------------- update mapping for tail item */
        this->_keys[_ipos->_kidx] =
            (iptr_type)(
            _ipos  - this->_heap.head());
        }
        this->_heap._pop_tail();

        this->_keys[_kpos] = null() ;
    }

    /*
    --------------------------------------------------------
     * UPDATE: update data in heap
    --------------------------------------------------------
     */

    __inline_call void_type update ( // copy
        iptr_type _kpos,
        data_type const&_data
        )
    {/*------------------ move "hole" to updated position */
        size_type _hpos =
            this->_keys[_kpos];
        _write_it _ipos =
            this->_heap.tend();

    /*------------------------ push "hole" to lower level */
        _ipos = push_lower (
            this-> _heap.head(),
            this-> _heap.tail(),
            this-> _heap.head()+_hpos ,
        __copy(data_type,_data)) ;

    /*------------------------ copy this data into "hole" */
        _ipos->_kidx =   _kpos   ;
        _ipos->_data =
        __copy(data_type,_data)  ;

    /*------------------------ copy position into mapping */
        this->_keys[_kpos] =
            (iptr_type)
            _ipos - this->_heap.head()  ;
    }

    __inline_call void_type update ( // move
        iptr_type _kpos,
        data_type &&_data
        )
    {/*------------------ move "hole" to updated position */
        size_type _hpos =
            this->_keys[_kpos];
        _write_it _ipos =
            this->_heap.tend();

    /*------------------------ push "hole" to lower level */
        _ipos = push_lower (
            this-> _heap.head(),
            this-> _heap.tail(),
            this-> _heap.head()+_hpos ,
        __copy(data_type,_data)) ;

    /*------------------------ copy this data into "hole" */
        _ipos->_kidx =   _kpos   ;
        _ipos->_data =
        __move(data_type,_data)  ;

    /*------------------------ copy position into mapping */
        this->_keys[_kpos] =
            (iptr_type)(
            _ipos - this->_heap.head()) ;
    }

    /*
    --------------------------------------------------------
     * REDUCE: update data in heap ("lower" priority)
    --------------------------------------------------------
     */

    __inline_call void_type reduce ( // copy
        iptr_type _kpos,
        data_type const&_data
        )
    {/*------------------ move "hole" to updated position */
        size_type _hpos =
            this->_keys[_kpos];
        _write_it _ipos = push_upper (
            this->_heap.head(),
            this->_heap.head()+_hpos ,
        __copy(data_type,_data)) ;

    /*------------------------ copy this data into "hole" */
        _ipos->_kidx =   _kpos   ;
        _ipos->_data =
        __copy(data_type,_data)  ;

    /*------------------------ copy position into mapping */
        this->_keys[_kpos] =
            (iptr_type)(
            _ipos - this->_heap.head()) ;
    }

    __inline_call void_type reduce ( // move
        iptr_type _kpos,
        data_type &&_data
        )
    {/*------------------ move "hole" to updated position */
        size_type _hpos =
            this->_keys[_kpos];
        _write_it _ipos = push_upper (
            this->_heap.head(),
            this->_heap.head()+_hpos ,
        __copy(data_type,_data)) ;

    /*------------------------ copy this data into "hole" */
        _ipos->_kidx =   _kpos   ;
        _ipos->_data =
        __move(data_type,_data)  ;

    /*------------------------ copy position into mapping */
        this->_keys[_kpos] =
            (iptr_type)(
            _ipos - this->_heap.head()) ;
    }

    /*
    --------------------------------------------------------
     * test heap validity (debug only!)
    --------------------------------------------------------
     */

    __normal_call bool_type test_heap (
        )
    {
#       ifdef _DEBUG
    /*------------------ test relationships for all nodes */
        size_type _ipos = +1 ;
        size_type _iend =
        this->_heap.count () ;
        for ( ; _ipos < _iend; ++_ipos)
        {
    /*------------------ heap is invalid if lower < upper */
            size_type _ppos =
                (_ipos - 1)/_nfan ;

            if (this->_pred(
                this->_heap[_ipos]._data,
                this->_heap[_ppos]._data
                )   )
                return false ;
        }
#       endif

    /*------------------ must be valid if we got this far */
        return true ;
    }

    } ;

#   undef  __cont


    }

#   endif   //__PRIORITY_IDX__



