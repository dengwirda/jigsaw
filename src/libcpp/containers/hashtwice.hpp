
/*
------------------------------------------------------------
 * a "2-choice" chained hash-table.
------------------------------------------------------------
 *
 * HASH-TABLE is a dynamically-sized, chained hashtable
 * implementation, essentially a linear array of singly-
 * linked hash buckets. Here, the so-called "2-choice"
 * variant is implemented, in which each operation uses
 * a pair of hash functions. Obj. are inserted into the
 * slot with the lower collision count.
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

#   ifndef __HASH_TWICE__
#   define __HASH_TWICE__

#   include "arraylist.hpp"

    namespace containers {

#   define __cont containers

    template <
    typename D ,
    typename H ,
    typename P ,
    typename A = allocators::basic_alloc
             >
    class hash_twice :
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

    typedef containers::array   <
                    size_type   ,
                    allocator   >       size_list ;

    typedef containers::hash_twice  <
                    data_type   ,
                    hash_type   ,
                    pred_type   ,
                    allocator   >       self_type ;

    size_type static constexpr _mini_count = +8 ;

    public  :

    hash_type const     _hfun;
    pred_type const     _pred;

    size_list           _lsiz;

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
        this->_lsiz.set_count(
        this->_lptr.count(), containers::tight_alloc, 0);

        typename size_list::_write_it
        _ITER = this->_lsiz.head(),
        _TEND = this->_lsiz.tend();
        for( ; _ITER != _TEND; ++_ITER)
        {
            *_ITER = (size_type)+0;
        }

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
            uint32_t _hash[2] ;
             this->_hfun(
            _head->_data, _hash[0], _hash[1]) ;

            _hash[0] &= this->hash_mask();
            _hash[1] &= this->hash_mask();

        /*--------------------------- push onto new table */
            if (this->_lsiz[_hash[0]] <=
                this->_lsiz[_hash[1]] )
            {
                base_type::push_item(
                this->_lptr[_hash[0]], _head) ;

                this->_lsiz[_hash[0]] += 1 ;
            }
            else
            {
                base_type::push_item(
                this->_lptr[_hash[1]], _head) ;

                this->_lsiz[_hash[1]] += 1 ;
            }
        }
    }

    public  :

/*--------------------------- default c'tor - do nothing! */
    __inline_call  hash_twice (
        hash_type const&_hsrc = hash_type(),
        pred_type const&_psrc = pred_type(),
        double    const&_lsrc = double(+.8),
        allocator const&_asrc = allocator(),
        size_type const&_nsrc = size_t(+ 0)
    /*----------------------------- c'tor alloc from obj. */
        ) : base_type(  _asrc),
    /*----------------------------- c'tor other from obj. */
        _hfun(_hsrc), _pred(_psrc),
        _lsiz(_asrc), _load(_lsrc)
    {
        set_slots(_nsrc,
                containers::tight_alloc) ;
    }

/*--------------------------- default d'tor/copy/move ops */
    __inline_call~hash_twice ()  = default ;

    __inline_call hash_twice (
        self_type const& _src
            )                    = default ;
    __inline_call hash_twice (
        self_type &&     _src
            )                    = default ;

    __inline_call self_type& operator = (
        self_type const& _src
            )                    = default ;
    __inline_call self_type& operator = (
        self_type  &&    _src
            )                    = default ;

/*------------------------------- calc. table load-factor */
    __inline_call double    load_fact (
        size_type _iinc = +0
        )
    {   return   (double) (this->_size +_iinc) /
                 (double)  this->_lptr.count() ;
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
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- push data onto list */
        if (this->_lsiz[_hash[0]] <= +1 ||
            this->_lsiz[_hash[0]] <=
            this->_lsiz[_hash[1]] )
        {
        this->_lsiz[_hash[0]] += +1 ;
        return base_type::push(_data,_hash[0]) ;
        }
        else
        {
        this->_lsiz[_hash[1]] += +1 ;
        return base_type::push(_data,_hash[1]) ;
        }
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
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- push data onto list */
        if (this->_lsiz[_hash[0]] <= +1 ||
            this->_lsiz[_hash[0]] <=
            this->_lsiz[_hash[1]] )
        {
        this->_lsiz[_hash[0]] += +1 ;
        return base_type::push(_data,_hash[0]) ;
        }
        else
        {
        this->_lsiz[_hash[1]] += +1 ;
        return base_type::push(_data,_hash[1]) ;
        }
    }

/*-------------------------- scan table for exact matches */
    __normal_call bool_type find (
        data_type const&_data,
        item_type     *&_same
        )
    {
        if (this->_lptr.empty()) return false;

    /*------------------------------- evaluate hash value */
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- scan list from head */
        _same = this->_lptr[_hash[0]];
    /*------------------------------- check exact matches */
        for( ; _same != nullptr;
                    _same = _same->_next)
        {
            if (this->_pred(_same->_data,_data))
            {
                return true ;
            }
        }

    /*------------------------------- scan list from head */
        _same = this->_lptr[_hash[1]];
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
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- scan list from head */
        item_type *_spos ;
        _spos = this->_lptr[_hash[0]];
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

    /*------------------------------- scan list from head */
        _spos = this->_lptr[_hash[1]];
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
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- scan list from head */
        item_type *_same;
        _same = this->_lptr[_hash[0]] ;
    /*------------------------------- check exact matches */
        for( ; _same != nullptr;
                    _same = _same->_next)
        {
            if (this->_pred(_same->_data,_data))
            {
                _list.push_tail(_same);
            }
        }

    /*------------------------------- scan list from head */
        _same = this->_lptr[_hash[1]] ;
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
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- scan list from head */
        item_type *_prev , *_item;
        _prev = nullptr;
        _item = this->_lptr[_hash[0]];
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if (this->_pred(_item->_data,_data))
            {
                this->_lsiz[_hash[0]] -= 1 ;
            /*-------------------- steal data before _pop */
                _same = std::move(_item->_data);
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev, this) ,
                    _write_it(_item, this) ,
                    _hash[0]) ;
            /*------------------------- exact match found */
                return true ;
            }
        }

    /*------------------------------- scan list from head */
        _prev = nullptr;
        _item = this->_lptr[_hash[1]];
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if (this->_pred(_item->_data,_data))
            {
                this->_lsiz[_hash[1]] -= 1 ;
            /*-------------------- steal data before _pop */
                _same = std::move(_item->_data);
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev, this) ,
                    _write_it(_item, this) ,
                    _hash[1]) ;
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
        uint32_t _hash[2] ;
        this->_hfun(_data, _hash[0], _hash[1]) ;

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- scan list from head */
        item_type *_prev , *_item;
        _prev = nullptr;
        _item = this->_lptr[_hash[0]];
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if (this->_pred(_item->_data,_data))
            {
                this->_lsiz[_hash[0]] -= 1 ;
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev, this) ,
                    _write_it(_item, this) ,
                    _hash[0]) ;
            /*------------------------- exact match found */
                return true ;
            }
        }

    /*------------------------------- scan list from head */
        _prev = nullptr;
        _item = this->_lptr[_hash[1]];
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if (this->_pred(_item->_data,_data))
            {
                this->_lsiz[_hash[1]] -= 1 ;
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev, this) ,
                    _write_it(_item, this) ,
                    _hash[1]) ;
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
        uint32_t _hash[2] ;
        this->_hfun(_iptr->_data, _hash[0], _hash[1]);

        _hash[0] &= this->hash_mask();
        _hash[1] &= this->hash_mask();

    /*------------------------------- scan list from head */
        item_type *_prev , *_item;
        _prev = nullptr;
        _item = this->_lptr[_hash[0]];
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if(_iptr == _item)
            {
                this->_lsiz[_hash[0]] -= 1 ;
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev, this) ,
                    _write_it(_item, this) ,
                    _hash[0]) ;
            /*------------------------- exact match found */
                return true ;
            }
        }

    /*------------------------------- scan list from head */
        _prev = nullptr;
        _item = this->_lptr[_hash[1]];
    /*------------------------------- check exact matches */
        for( ; _item != nullptr;
               _prev  = _item, _item = _item->_next)
        {
            if(_iptr == _item)
            {
                this->_lsiz[_hash[1]] -= 1 ;
            /*------------------- _pop and _destruct item */
                base_type::_pop(
                    _write_it(_prev, this) ,
                    _write_it(_item, this) ,
                    _hash[1]) ;
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
            size_type _lnum = +0;
            for( ; _next != nullptr;
                        _next = _next->_next)
            {
                _lnum += +1  ;
            }
        /*----------------------------- update statistics */
            if (_min_count > _lnum)
                _min_count = _lnum;
            if (_max_count < _lnum)
                _max_count = _lnum;

            if (_lnum == (size_type) +0)
                _nil_ratio += +1. ;

            if (_bad_limit < _lnum)
                _bad_ratio += +1. ;

            _ave_count += (double)_lnum;
        }

        _nil_ratio /= this->_lptr.count() ;
        _bad_ratio /= this->_lptr.count() ;
        _ave_count /= this->_lptr.count() ;
    }

    } ;

#   undef  __cont


    }

#   endif   //__HASH_TWICE__



