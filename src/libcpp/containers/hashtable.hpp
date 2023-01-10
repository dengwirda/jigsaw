
/*
------------------------------------------------------------
 * a chained hash-table.
------------------------------------------------------------
 *
 * HASH-TABLE is a dynamically-sized, chained hash-table
 * implementation, essentially a linear array of singly-
 * linked hash buckets.
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
 * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant
 * this code in any way whatsoever.  This code is
 * provided "as-is" to be used at your own risk.
 *
 * THE CONTRIBUTORS include:
 * (a) The University of Sydney
 * (b) The Massachusetts Institute of Technology
 * (c) Columbia University
 * (d) The National Aeronautics & Space Administration
 * (e) Los Alamos National Laboratory
 *
------------------------------------------------------------
 *
 * Last updated: 28 Apr., 2020
 *
 * Copyright 2013-2020
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __HASH_TABLE__
#   define __HASH_TABLE__

#   include "arraylist.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename D ,
    typename H ,
    typename P ,
    typename A = allocators::basic_alloc
             >
    class hash_table  :
        public containers::array_list<D,A>
    {
/*------------------------- a dynamic, chained hash-table */
    public  :

    typedef D                           data_type ;
    typedef H                           hash_type ;
    typedef P                           pred_type ;
    typedef A                           allocator ;

    typedef containers::array_list  <
                    data_type   ,
                    allocator   >       base_type ;

    typedef typename
            base_type::size_type        size_type ;
    typedef typename
            base_type::diff_type        diff_type ;
    typedef typename
            base_type::item_type        item_type ;
    typedef typename
            base_type::lptr_list        lptr_list ;

    typedef typename
            base_type::_write_it        _write_it ;
    typedef typename
            base_type::_const_it        _const_it ;

    typedef containers::hash_table  <
                    data_type   ,
                    hash_type   ,
                    pred_type   ,
                    allocator   >       self_type ;

    size_type static constexpr _mini_count = +8 ;

    public  :

    hash_type const     _hash;
    pred_type const     _pred;

    double              _load;

    public  :

/*------------------------ update table count and re-hash */
    __inline_call size_type hash_mask (
        ) const
    {   return this->_lptr.count() -  1 ;
    }

    __static_call
    __inline_call bool_type   is_pwr2 (
        size_type _xval
        )
    {   return (_xval & (_xval-1)) == 0 ;
    }

    __static_call
    __inline_call size_type next_pwr2 (
        size_type _xval
        )
    {
        if (_xval<= 1) return 1 ;
        int _pwr2 = 2;
        _xval--;
        while (_xval >>= 1) _pwr2 <<= 1 ;
        return _pwr2 ;
    }

    __static_call
    __inline_call size_type prev_pwr2 (
        size_type _xval
        )
    {
        int _pwr2 = 1;
        while (_xval >>= 1) _pwr2 <<= 1 ;
        return _pwr2 ;
    }

    __normal_call void_type set_slots (
        size_type _slots,
        __cont::alloc_types _alloc =
        __cont::loose_alloc
        )
    {
    /*------------------------------- round to next 2^pwr */
        _slots = next_pwr2 (_slots);

    /*------------------------------- inc/dec. table size */
        if (_slots <
                this->_mini_count )
            this->_lptr.set_count (
            this->_mini_count ,
                       _alloc ,  nullptr) ;
        else
            this->_lptr.set_count (
               _slots, _alloc ,  nullptr) ;

    /*------------------------------- re-hash all objects */
        redo_hash() ;
    }

    __normal_call void_type grow_hash (
        )
    {
    /*------------------------------- increase table size */
        if (this->_lptr.count() <
                this->_mini_count )
            this->_lptr.set_count (
            this->_mini_count ,
        containers::tight_alloc, nullptr) ;
        else
            this->_lptr.set_count (
            this->_lptr.count()*+2,
        containers::tight_alloc, nullptr) ;

    /*------------------------------- re-hash all objects */
        redo_hash() ;
    }

    __normal_call void_type redo_hash (
        )
    {
        if (this->_size == +0) return ;

        __assert (
            is_pwr2(this->_lptr.count() ) &&
        "hashtable: count must be 2^n!" ) ;

        item_type*_head = nullptr ;
        item_type*_tail = nullptr ;
        item_type*_next = nullptr ;
    /*--------------- splice together a list of all items */
        typename lptr_list::_write_it
        _iter = this->_lptr.head(),
        _tend = this->_lptr.tend();
        for( ; _iter != _tend; ++_iter)
        {
        /*--------------------- find first non-empty list */
            if (*_iter != nullptr)
            {
                _head = *_iter;
                _tail = *_iter;
               *_iter = nullptr;
                _iter++; break;
            }
        }
        for( ; _iter != _tend; ++_iter)
        {
        /*--------------- splice remaining lists together */
            if (*_iter != nullptr)
            {
                _tail =  *_iter;
                for(; nullptr != _tail->_next; )
                {
                    _tail=_tail->_next ;
                }
            /*------ re-attach at _tail and jump to _head */
                _tail->
                _next =  _head;
                _head = *_iter;
               *_iter = nullptr ;
            }
        }
    /*------------------ re-hash all items onto new table */
        for ( ; _head != nullptr; _head = _next)
        {
            _next = _head->_next;
        /*--------------------------- calc new hash value */
            size_type _hpos =
                this->_hash(_head->_data)
                    & this-> hash_mask () ;
        /*--------------------------- push onto new table */
            base_type::push_item(
                this->_lptr[_hpos],_head) ;
        }
    }

    public  :

/*--------------------------- default c'tor - do nothing! */
    __inline_call  hash_table (
        hash_type const&_hsrc = hash_type(),
        pred_type const&_psrc = pred_type(),
        double    const&_lsrc = double(+.8),
        allocator const&_asrc = allocator(),
        size_type const&_nsrc = size_t(+ 0)
    /*----------------------------- c'tor alloc from obj. */
        ) : base_type(  _asrc),
    /*----------------------------- c'tor other from obj. */
            _hash(_hsrc) ,
            _pred(_psrc) ,
            _load(_lsrc)
    {
        set_slots(_nsrc,
                containers::tight_alloc) ;
    }

/*--------------------------- default d'tor/copy/move ops */
    __inline_call~hash_table ()  = default ;

    __inline_call hash_table (
        self_type const& _src
            )                    = default ;
    __inline_call hash_table (
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

/*------------------------------- calc. table load-factor */
    __inline_call double    load_fact (
        size_type _iinc = +0
        )
    {   return   (double)(this->_size +_iinc) /
                 (double) this->_lptr.count() ;
    }

/*----------------------------- push data onto hash table */
    __inline_call _write_it push ( // copy construct
        data_type const&_data
        )
    {
    /*------------- re-size table if load factor exceeded */
        if (load_fact(1)> this->_load)
            grow_hash() ;
    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask () ;

    /*------------------------------- push data onto list */
        return base_type::push( _data,_hpos ) ;
    }

/*----------------------------- push data onto hash table */
    __inline_call _write_it push ( // move construct
        data_type &&    _data
        )
    {
    /*------------- re-size table if load factor exceeded */
        if (load_fact(1)> this->_load)
            grow_hash() ;
    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask () ;

    /*------------------------------- push data onto list */
        return base_type::push( _data,_hpos ) ;
    }

/*-------------------------- scan table for exact matches */
    __normal_call bool_type find (
        data_type const&_data,
        item_type     *&_same
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask ();

    /*------------------------------- scan list from head */
        _same = this->_lptr[_hpos];
    /*------------------------------- check exact matches */
        for( ; _same != nullptr;
                    _same = _same->_next)
        {
            if (this->_pred(_same->_data,_data))
            {
                return true ;
            }
        }

    /*---------------------------------- no matches found */
        return false ;
    }

    __normal_call bool_type find (
        data_type const&_data,
        _write_it      &_same
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask ();

    /*------------------------------- scan list from head */
        item_type*_spos = this->_lptr[_hpos] ;
    /*------------------------------- check exact matches */
        for( ; _spos != nullptr;
                    _spos = _spos->_next)
        {
            if (this->_pred(_spos->_data,_data))
            {
                _same = _write_it (_spos,
                        (base_type*)this) ;

                return true ;
            }
        }

    /*---------------------------------- no matches found */
        return false ;
    }

/*-------------------------- scan table for exact matches */
    template <
    typename      list_type
             >
    __normal_call bool_type find (
        data_type const&_data,
        list_type      &_list
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask ();

    /*------------------------------- scan list from head */
        item_type*_same = this->_lptr[_hpos] ;
    /*------------------------------- check exact matches */
        for( ; _same != nullptr;
                    _same = _same->_next)
        {
            if (this->_pred(_same->_data,_data))
            {
                _list.push_tail(_same);
            }
        }

    /*---------------------------------- no matches found */
        return ( !_list.empty() );
    }

/*--------------------- _pop any exact matches from table */
    __normal_call bool_type _pop (
        data_type const&_data,
        data_type      &_same
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask ();

    /*------------------------------- scan list from head */
        item_type*_prev = nullptr;
        item_type*_item = this->_lptr[_hpos] ;
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if (this->_pred(_item->_data,_data))
            {
            /*-------------------- steal data before _pop */
                _same = std::move(_item->_data);
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev,this),
                        _write_it(_item,this),_hpos) ;
            /*------------------------- exact match found */
                return true ;
            }
        }

    /*--------------------------------- couldnt find data */
        return false ;
    }

/*--------------------- _pop any exact matches from table */
    __normal_call bool_type _pop (
        data_type const&_data
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_data )
                        & this-> hash_mask () ;

    /*------------------------------- scan list from head */
        item_type*_prev = nullptr;
        item_type*_item = this->_lptr[_hpos] ;
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if (this->_pred(_item->_data,_data))
            {
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev,this),
                        _write_it(_item,this),_hpos) ;
            /*------------------------- exact match found */
                return true ;
            }
        }

    /*--------------------------------- couldnt find data */
        return false ;
    }

/*------------------- _pop any pointer matches from table */
    __normal_call bool_type _pop (
        item_type *_iptr
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        size_type _hpos = this->_hash(_iptr->
                                      _data )
                        & this-> hash_mask ();

    /*------------------------------- scan list from head */
        item_type*_prev = nullptr;
        item_type*_item = this->_lptr[_hpos] ;
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if(_iptr == _item)
            {
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev,this),
                        _write_it(_item,this),_hpos) ;
            /*------------------------- exact match found */
                return true ;
            }
        }

    /*--------------------------------- couldnt find data */
        return false ;
    }

/*---------------------- scan table and report statistics */
    __normal_call void_type _get_info (
        size_type &_min_count,
        size_type &_max_count,
        double    &_nil_ratio,
        double    &_bad_ratio,
        double    &_ave_count
        )
    {
        _min_count =
         std::numeric_limits<size_type>::max();
        _max_count =
         std::numeric_limits<size_type>::min();

        size_type _bad_limit =
            (size_type) (+3 * this->_load);

        _nil_ratio = (double)+.0;
        _bad_ratio = (double)+.0;
        _ave_count = (double)+.0;
    /*-------------------------- iter. all lists in table */
        typename lptr_list::_write_it
        _iter = this->_lptr.head(),
        _tend = this->_lptr.tend();
        for( ; _iter != _tend; ++_iter)
        {
        /*---------------------- count items in each list */
            item_type*_next =*_iter;
            size_type _lsiz = +0;
            for( ; _next != nullptr;
                        _next = _next->_next)
            {
                _lsiz += +1  ;
            }
        /*----------------------------- update statistics */
            if (_min_count > _lsiz)
                _min_count = _lsiz;
            if (_max_count < _lsiz)
                _max_count = _lsiz;

            if (_lsiz == (size_type) +0)
                _nil_ratio += +1. ;

            if (_bad_limit < _lsiz)
                _bad_ratio += +1. ;

            _ave_count += (double)_lsiz;
        }

        _nil_ratio /= this->_lptr.count() ;
        _bad_ratio /= this->_lptr.count() ;
        _ave_count /= this->_lptr.count() ;
    }

    } ;

#   undef  __cont


    }

#   endif   //__HASH_TABLE__



