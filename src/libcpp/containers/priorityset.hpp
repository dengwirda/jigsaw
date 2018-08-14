
/*
------------------------------------------------------------
 * an "n-ary"-heap based priority queue.
------------------------------------------------------------
 *
 * PRIORITY-MAP utilises an "n-ary" heap data-structure, 
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

#   ifndef __PRIORITYSET__
#   define __PRIORITYSET__

#   include "array.hpp"

    namespace containers {

#   define __cont containers
 
    template <
    typename D ,
    typename P = std::less<D> ,
    typename A = allocators::basic_alloc ,
    typename C = containers::array <D,A>
             >
    class priorityset
    { 
/*----------- a dynamic priority queue as an "n"-ary heap */
    public  :
    
    typedef D                           data_type ;
    typedef C                           container ;
    typedef P                           pred_type ;
    typedef A                           allocator ;
    
    typedef typename 
            container::size_type        size_type ;

    typedef priorityset <
            data_type , 
            pred_type ,
            allocator , 
            container   >               self_type ;
    
    typedef typename 
            container::_write_it        _write_it ;
    typedef typename 
            container::_const_it        _const_it ;

    size_type static const _nfan = +4;  // fan out

    public  :

    container              _heap ;
    
    pred_type              _pred ;

    private :
        
    /*
    --------------------------------------------------------
     * helper: push "hole" into sorted position
    --------------------------------------------------------
     */
     
    __normal_call _write_it push_upper (
        _write_it _head,
        _write_it _ipos,
        data_type const&_data
        )
    {
        for ( ; _ipos != _head ; )
        {
        /*----------------------- find position of parent */
            _write_it _ppos = 
            _head+(_ipos-_head-1) / _nfan ;

        /*----------------------- swap parent with "hole" */
            if (this->_pred(_data, *_ppos))
            {
               *_ipos = std::move( *_ppos);
                _ipos = _ppos ;
            }
            else break;
        }
        
    /*---------------------------- return sorted position */
        return (_ipos);
    }
    
    /*
    --------------------------------------------------------
     * helper: push "hole" into sorted position
    --------------------------------------------------------
     */
     
    __normal_call _write_it push_lower (
        _write_it _head,
        _write_it _tail,
        _write_it _ipos,
        data_type const&_data
        )
    {
        if (_head == _tail) return _ipos ;
    /*-- find position of "root" within trailing sub-tree */
        _write_it _cpos ;
        _write_it _imin ;
        _write_it _iend = _head;
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
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 7 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 6 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 5 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 4 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 3 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 2 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            }
        /*------------------------------ swap with "hole" */
           *_ipos = std::move(*_imin) ;
            _ipos = _imin ;
        }
    /*-- deal with special cases - sort trailing sub-tree */
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
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 7 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 6 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 5 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 4 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 3 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            case 2 :
            if (this->_pred(*--_cpos,*_imin))
                _imin=_cpos;
            }
        /*------------------------------ swap with "hole" */
            *_ipos = std::move( *_imin) ;
             _ipos = _imin ;
            }
        }
    /*------------------ push "hole "into sorted position */
        return push_upper(_head, _ipos, _data) ;
    }

    public  :

/*--------------------------- default c'tor - do nothing! */
    __inline_call priorityset (
        allocator const&_asrc = allocator(),
        pred_type const&_psrc = pred_type()
        ) : _heap(_asrc),
            _pred(_psrc)   {}
/*--------------------------- default c'tor - initialisor */
    __inline_call priorityset (
        pred_type const&_psrc
        ) : _heap(allocator()),
            _pred(_psrc)   {}

    __inline_call~priorityset() = default ;

    __inline_call priorityset (
        self_type const&_src
        )                       = default ;
    __inline_call priorityset (
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
    {   return this->_heap[   +0] ; 
    }
    
/*------------------------------------- peek at heap item */
    __inline_call data_type const& peek (
        size_type  _hpos
        ) const
    {   return this->_heap[_hpos] ;
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
    __inline_call void_type set_alloc (
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
        data_type const&_data
        )
    {
    /*---------------- push data onto "tail" of container */
        size_type _tpos = 
        this->_heap.push_tail() ;
    /*---------------- sort "hole" for data into position */
        _write_it _ipos = push_upper (
            this->_heap.head() , 
            this->_heap.head() + _tpos , 
        __copy(data_type,_data));
    /*---------------- copy new data into sorted position */
        *_ipos  = 
        __copy(data_type,_data) ;
    }
    
    __inline_call void_type push ( // move
        data_type &&_data
        )
    {
    /*---------------- push data onto "tail" of container */
        size_type _tpos = 
        this->_heap.push_tail() ;
    /*---------------- sort "hole" for data into position */
        _write_it _ipos = push_upper (
            this->_heap.head() , 
            this->_heap.head() + _tpos , 
        __copy(data_type,_data));
    /*---------------- copy new data into sorted position */
        *_ipos  = 
        __move(data_type,_data) ;
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
        data_type &_data
        )
    {
    /*---------------------------- _pop root, return data */
        _pop( +0,  _data);
    }
    
    __inline_call void_type _pop (
        size_type  _hpos
        )
    {
    /*---------------------------- _pop data, cache temp. */
         data_type _temp ;
        _pop(_hpos,_temp);
    }
    __normal_call void_type _pop (
        size_type  _hpos,
        data_type &_data
        )
    {
    /*---------------------------- _pop data, return data */
        _data = std::move (
            *(this->_heap.head() + _hpos));
    /*--------- push "hole" at HPOS into updated position */
        if (this->_heap.count() > 1)
        {
        /*---------- find new position for "hole" at HPOS */
            _write_it _ipos = push_lower (
                this->_heap.head() ,  
                this->_heap.tail() -    1 , 
                this->_heap.head() + _hpos, 
               *this->_heap.tail());
        /*---------- copy tail data into updated position */
           *_ipos = 
            std::move(*this->_heap.tail());
        } 
        this->_heap._pop_tail() ;
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
        size_type _iend = this->_heap.count() ;
        for ( ; _ipos < _iend; ++_ipos)
        {
            size_type _ppos = (_ipos-1)/_nfan ;
        /*-------------- heap is invalid if lower < upper */
            if (this->_pred(this->_heap[_ipos], 
                            this->_heap[_ppos])
                        )
                return false ;
        }
#       endif

    /*------------------ must be valid if we got this far */
        return true ;
    }
    
    } ;

#   undef  __cont


    }

#   endif//__PRIORITYSET__



