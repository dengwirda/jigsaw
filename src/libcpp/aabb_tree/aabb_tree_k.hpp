
/*
------------------------------------------------------------
 * AABB-TREE-K: AABB-tree construction in R^k.
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
 * Last updated: 20 Apr., 2023
 *
 * Copyright 2013-2023
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __AABB_TREE_K__
#   define __AABB_TREE_K__

    namespace geom_tree {

    template <
    typename I ,  // item-to-be-stored typedef
    size_t   K ,  // number of geometric dimensions
    typename N = aabb_node_base_k,  // inner node typedef
    typename A = allocators::basic_alloc
             >
    class aabb_tree
    {
/*---------------------- a static d-dimensional AABB-tree */
    public  :

    typedef I                           item_type ;
    typedef N                           node_user ;
    typedef A                           allocator ;

    typedef typename
            item_type::real_type        real_type ;
    typedef typename
            item_type::iptr_type        iptr_type ;

    iptr_type static constexpr         _dims =  K ;

    iptr_type static constexpr middle_split  = +0 ;
    iptr_type static constexpr median_split  = +1 ;
    iptr_type static constexpr window_split  = +2 ;

    typedef geom_tree::aabb_tree  <
                item_type ,
                    _dims ,
                node_user ,
                allocator >             tree_type ;

    typedef containers::single_item <
                item_type >             item_data ;

    class node_type ;
    class node_data ;

    class node_type : public node_user
        {
    /*------------------------------- aabb-tree node type */
        public  :
        real_type           _pmin[ K] ;
        real_type           _pmax[ K] ;
    /*------------------------------- linkedlist of items */
        item_data          *_hptr;
        iptr_type           _size = 0 ;
    /*------------------------------- tree-based pointers */
        node_type          *_pptr;
        node_data          *_lptr;
    /*----------------------------- tree-pointer handlers */
        __inline_call item_data* items (
            ) const
        {   return  this->_hptr ;
        }
        __inline_call node_type* upper (
            ) const
        {   return  this->_pptr ;
        }
        __inline_call node_type* lower (
            iptr_type _ipos
            ) const
        {
            if (this->_lptr != nullptr)
        /*------------------------------ pointer to child */
            return &this->_lptr
                        ->_node[_ipos];
            else
        /*------------------------------ must be a "leaf" */
            return nullptr;
        }
        } ;
    class node_data
        {
    /*----------------- keep node-pair together in memory */
        public  :
        node_type           _node[ 2] ;
        } ;

/*----------- two-layer pool'd allocator -- items + nodes */
    typedef allocators::_pool_alloc<
                allocator>              base_pool ;
    typedef allocators::_wrap_alloc<
                base_pool>              wrap_pool ;

    typedef allocators::_item_alloc<
                node_data,
                wrap_pool>              node_pool ;
    typedef allocators::_item_alloc<
                item_data,
                wrap_pool>              item_pool ;

    typedef containers::array <
                node_type* ,
                allocator>              work_list ;

    public  :
/*------------------------------------- tree root pointer */
    node_type                   _rdat ;
    node_type                  *_root ;
    iptr_type                   _size ;
    iptr_type                   _imax ;

    work_list                   _work ;

/*------------------------------------- pool'd allocators */
    base_pool              _node_base ;
    base_pool              _item_base ;

    node_pool              _node_pool ;
    item_pool              _item_pool ;

/*------------------------------------- tree shape params */
    iptr_type                   _npop ;
    real_type                   _long ;
    real_type                   _vtol ;

    private :
/*---------------------- helper - delete a block of nodes */
    __inline_call void_type free_node (
      __write_ptr(node_data)&_data
        )
    {
        this->
       _node_pool._destruct (_data) ;
        this->
       _node_pool.deallocate(_data,1) ;
       _data = nullptr ;
    }

/*---------------------- helper - create a block of nodes */
    __inline_call void_type make_node (
      __write_ptr(node_data)&_data
        )
    {
       _data =
        this->_node_pool.allocate (1) ;
        this->
       _node_pool.construct (_data) ;
    }

/*------------------------------- helper - delete an item */
    __inline_call void_type free_item (
      __write_ptr(item_data) _idat
        )
    {
        this->
       _item_pool._destruct (_idat) ;
        this->
       _item_pool.deallocate(_idat,1) ;
       _idat = nullptr ;
    }

/*------------------------------- helper - create an item */
    __inline_call void_type make_item (
        item_type const&_item,
      __write_ptr(item_data)&_idat
        )
    {
       _idat =
        this->_item_pool.allocate (1) ;
        this->_item_pool.construct(
       _idat, nullptr, _item);
    }

/*--------------------- helper - push an item onto a list */
    __static_call
    __inline_call void_type push_item (
      __write_ptr(item_data)&_head,
      __write_ptr(item_data) _item
        )
    {
        _item->_next = _head ;
        _head        = _item ;
    }

/*--------------------- helper - _pop an item from a list */
    __static_call
    __inline_call void_type _pop_item (
      __write_ptr(item_data)&_head,
      __write_ptr(item_data) _item,
      __write_ptr(item_data) _prev
        )
    {
        if (_prev != nullptr )
        {
    /*---------------------------- re-link in list middle */
        _prev->
        _next = _item->_next ;
        }
        else
        {
    /*---------------------------- re-link upon list head */
        _head = _item->_next ;
        }
    }

    public  :
/*----------------------------------------- default c'tor */
    __inline_call  aabb_tree  (
        allocator const&_asrc = allocator()
        ) : _root(nullptr) ,
            _size(+0) , _imax(+0),
            _work(  _asrc) ,
    /*-------------------------------------- "base" pools */
           _node_base (
        sizeof(node_data), _asrc),
           _item_base (
        sizeof(item_data), _asrc),
    /*-------------------------------------- "obj." pools */
           _node_pool (
        wrap_pool(&_node_base))  ,
           _item_pool (
        wrap_pool(&_item_base)) {}

/*----------------------------------------- default d'tor */

    // pool'd alloc. takes care of dealloc...

    __inline_call void_type clear (
        )
    {
        this->_root = nullptr ;
        this->_size = +0;
        this->_imax = +0;
        this->_node_base.clear() ;
        this->_item_base.clear() ;
    }

/*------------- return container count / empty statistics */
    __inline_call bool_type empty (
        ) const { return
              nullptr == this->_root ; }

    __inline_call iptr_type count (
        ) const { return this->_size ; }

/*------------- helper - find min-enclosing aabb for node */
    __static_call
    __normal_call void_type init_aabb_node (
        node_type *_nptr,
        node_type *_pptr=nullptr
        )
    {
    /*------------------------ init nodal aabb at -+ inf. */
        for(auto _idim = _dims ; _idim-- != +0; )
        {
            if (_pptr == nullptr)
            {
                _nptr->_pmin[_idim] =
        +std::numeric_limits<real_type>::infinity() ;
                _nptr->_pmax[_idim] =
        -std::numeric_limits<real_type>::infinity() ;
            }
            else
            {
                _nptr->_pmin[_idim] =
                _pptr->_pmax[_idim] ;
                _nptr->_pmax[_idim] =
                _pptr->_pmin[_idim] ;
            }
        }
    /*------------------ take min/max of item-wise aabb's */
        for(item_data*_iptr  = _nptr->_hptr;
                      _iptr !=  nullptr;
                      _iptr  = _iptr->_next)
        {
            push_aabb_node(_nptr, _iptr) ;
        }
    }

    __static_call
    __inline_call void_type push_aabb_node (
        node_type *_nptr, // node to update
        item_data *_iptr  // item to update with
        )
    {
        for (auto _idim = _dims; _idim-- != +0; )
        {
        /*------------------------ find max. on each axis */
            if (_nptr->_pmax[_idim] <
                    _iptr->_data.pmax(_idim))
            {
                _nptr->_pmax[_idim] =
                    _iptr->_data.pmax(_idim);
            }
        /*------------------------ find min. on each axis */
            if (_nptr->_pmin[_idim] >
                    _iptr->_data.pmin(_idim))
            {
                _nptr->_pmin[_idim] =
                    _iptr->_data.pmin(_idim);
            }
        }
    }

/*------------- helper - calc. best "split" axis for node */
    __normal_call void_type find_best_axis (
        node_type *_node,
        iptr_type &_bdim,
        real_type &_blen
        )
    {
        class rect_data
            {
    /*------------------ helper: sort aabb's by dimension */
            public  :
                real_type       _alen ;
                iptr_type       _axis ;
            } ;

        class rect_less
            {
    /*------------------ helper: sort aabb's by dimension */
            public  :
            __inline_call
                bool_type operator () (
                rect_data const&_xdat ,
                rect_data const&_ydat
                ) const
            {   return ( _xdat. _alen <
                         _ydat. _alen ) ;
            }
            } ;

    /*------------------------------ form aabb dimensions */
        rect_data _rect [_dims];
        iptr_type _best = +0 ;
        for(auto _idim = _dims ; _idim-- != +0; )
        {
    /*------------------------------ find long bbox. axis */
            real_type _alen;
            _alen = _node->_pmax[_idim] -
                    _node->_pmin[_idim] ;

            _rect[_idim]._alen = _alen ;
            _rect[_idim]._axis = _idim ;
        }

    /*------------------------------ sort aabb dimensions */
        algorithms::ssort (
           &_rect[0], &_rect[_dims], rect_less()) ;

    /*------------------------------ shortcuts: item=node */
        if (!item_type::_is_a_rect)
        {
            auto _last = _dims - 1;

            _bdim = _rect[_last]._axis ;
            _blen = this->_long *
                    _rect[_last]._alen ;

            return ;              // exit, no "long" node
        }

    /*------------------------------ scan aabb dimensions */
        for(auto _idim = _dims ; _idim-- != +0; )
        {
            iptr_type _axis ;
            _axis = _rect[_idim]._axis ;
            real_type _llen ;
            _llen = _rect[_idim]._alen ;

            _llen *= this->_long ;

            iptr_type _pnum = +0 ;
            iptr_type _cnum = +0 ;

        /*------------------ partition "long"/inner items */
            for (auto _iptr  = _node->_hptr;
                      _iptr !=  nullptr ;
                      _iptr  = _iptr->_next)
            {
                if (_iptr->
                    _data.plen(_axis) >= _llen)
                {
                    _pnum  += +1 ;
                }
                else
                {
                    _cnum  += +1 ;
                }
            }

        /*------------------ keep longest, non-empty dim. */
            if (_best < _cnum )
            {
                _best = _cnum ;
                _bdim = _axis ;
                _blen = _llen ;

            if (_pnum == +0) break ;
            }
        }

    }

/*--- top-down assembly of aabb-tree via recursive splits */
    template <
    typename      iter_type
             >
    __normal_call void_type load (
        iter_type _head ,  // head of items array to load
        iter_type _tend ,  // last of items array to load
        iptr_type _NPOP = + 32 ,  // leaf population lim.
        iptr_type _RULE = window_split ,  // split scheme
        real_type _LONG = +.75 ,  // "too-long" item lim.
        real_type _VTOL = +.50    // leaf vol. split lim.
        )
    {
        if (_tend <= _head) return ;

        this->_root = &this->_rdat ;
        this->_root->_size = +0 ;
        this->_size = +1 ;

    /*------------------------------ set node fill params */
        this->_npop = _NPOP ;
        this->_long = _LONG ;
        this->_vtol = _VTOL ;

    /*------------------------------ set null "tree" ptrs */
        this->_root->_pptr = nullptr;
        this->_root->_lptr = nullptr;
        this->_root->_hptr = nullptr;

    /*------------------------------ push items onto root */
        init_aabb_node (this->_root);
        
        {        
        item_data *_hptr   = nullptr;
        item_data *_idat   = nullptr;
        for(; _head != _tend; ++_head)
        {
    /*------------------------------ alloc. and push item */
            make_item(*_head, _idat);
            push_item( _hptr, _idat);

            this->_root->_size += +1;

            this->_imax = std::max(
            this->_imax,_head->ipos ()) ;
            
            push_aabb_node(this->_root, _hptr ) ;
        }
        this->_root->_hptr  = _hptr ;
        }

    /*-------------------------- a list of un-split nodes */
        this->_work.clear() ;
        this->_work.
            push_tail(this->_root);

    /*--- refine tree until all nodes satisfy constraints */
        containers::array<real_type> _sort;
        for( ; !this->_work.empty() ; )
        {
            node_type *_pnod = nullptr;
          __write_ptr(node_data) _ndat = nullptr;
            node_type *_rnod = nullptr;
            node_type *_lnod = nullptr;

            real_type  _blen ;
            iptr_type  _bdim = -1 ;
            iptr_type  _pnum = +0 ;
            iptr_type  _cnum = +0 ;
            iptr_type  _lnum = +0 ;
            iptr_type  _rnum = +0 ;

          __write_ptr(item_data) _next = nullptr;
          __write_ptr(item_data) _pptr = nullptr;
          __write_ptr(item_data) _cptr = nullptr;
          __write_ptr(item_data) _lptr = nullptr;
          __write_ptr(item_data) _rptr = nullptr;
          __write_ptr(item_data) _iptr = nullptr;

        /*-------------------------- _pop node from stack */
            this->_work._pop_tail (_pnod);

        /*---------- find best "split" axis for this node */
            find_best_axis (_pnod, _bdim ,
                            _blen) ;

            if (_bdim == -1)  continue;

        /*---------- compute split-pos. via user-def RULE */
            real_type _spos = (real_type)+0. ;
            iptr_type _nnum = _pnod->_size ;
            item_data*_hptr = _pnod->_hptr ;
            if (_RULE == middle_split)
            {
        /*---------- split pos. - middle of items in node */
            double _SPOS = +0. ;
            for(_iptr  = _hptr ; _iptr != nullptr ;
                _iptr  = _next )
            {
                _next  = _iptr-> _next ;
        /*------------------------------- partition items */
                if (_iptr->_data
                    .plen (_bdim) > _blen)
                {
                     push_item(_pptr, _iptr);
                    _pnum += +1 ;
                }
                else
                {
                     push_item(_cptr, _iptr);
                    _cnum += +1 ;

                    _SPOS +=
                    _iptr->_data.pmid(_bdim);
                }
            }

            _spos = (real_type)(_SPOS/_cnum);
            }
            else
            if (_RULE == median_split)
            {
        /*---------- split pos. - median of items in node */

        //  (1+-a)-approx. with probability p
        //  require: 1/a^2 + log(1/p) samples
            size_t _nsel =
           (size_t)(+500 * std::log2(_nnum + 1) + 1) ;

            size_t _filt = _nnum / _nsel + 1;
            size_t _njmp = _filt ;
            size_t _inum = +0;
            _sort.clear();
            for(_iptr  = _hptr ; _iptr != nullptr ;
                _iptr  = _next ,
                _inum += +1  )
            {
                _next  = _iptr-> _next ;
        /*------------------------------- partition items */
                if (_iptr->_data
                    .plen (_bdim) > _blen)
                {
                     push_item(_pptr, _iptr);
                    _pnum += +1 ;
                }
                else
                {
                     push_item(_cptr, _iptr);
                    _cnum += +1 ;

                    if(_inum >= _njmp)
                    {
                    _sort.push_tail(
                    _iptr->_data.pmid(_bdim)) ;

                    _njmp +=    _filt;
                    }
                }
            }

            if(!_sort.empty())
            {
            algorithms::qsort(
                _sort.head(), _sort.tend(),
                    std::less<real_type>()) ;

            _spos = _sort[_sort.count() / 2];
            }
            }
            else
            if (_RULE == window_split)
            {
        /*---------- split pos. - window of items in node */

        //  (1+-a)-approx. with probability p
        //  require: 1/a^2 + log(1/p) samples
            size_t _nsel =
           (size_t)(+500 * std::log2(_nnum + 1) + 1) ;

            size_t _filt = _nnum / _nsel + 1;
            size_t _njmp = _filt ;
            size_t _inum = +0;
            double _SPOS = +0. ;
            _sort.clear();
            for(_iptr  = _hptr ; _iptr != nullptr ;
                _iptr  = _next ,
                _inum += +1  )
            {
                _next  = _iptr-> _next ;
        /*------------------------------- partition items */
                if (_iptr->_data
                    .plen (_bdim) > _blen)
                {
                     push_item(_pptr, _iptr);
                    _pnum += +1 ;
                }
                else
                {
                     push_item(_cptr, _iptr);
                    _cnum += +1 ;

                    _SPOS +=
                    _iptr->_data.pmid(_bdim);

                    if(_inum >= _njmp)
                    {
                    _sort.push_tail(
                    _iptr->_data.pmid(_bdim)) ;

                    _njmp +=    _filt;
                    }
                }
            }

            if(!_sort.empty())     // middle-median blend
            {
            algorithms::qsort(
                _sort.head(), _sort.tend(),
                    std::less<real_type>()) ;

            _spos = _sort[_sort.count() / 2];
            _SPOS/= _cnum ;

            _spos = (real_type)(+1./3.*_spos) +
                    (real_type)(+2./3.*_SPOS) ;
            }
            }

            if (_cnum == +0) continue ;

        /*-------------------- alloc. new child node data */
            make_node(_ndat);
            _lnod = &_ndat->_node[ 0] ;
            _rnod = &_ndat->_node[ 1] ;

            this->_size  += +2 ;

            _pnod->_hptr  = _pptr ;
            _pnod->_lptr  = _ndat ;

            _lnod->_hptr  = nullptr ;
            _lnod->_lptr  = nullptr ;
            _lnod->_pptr  = _pnod ;

            init_aabb_node (_lnod , _pnod);

            _rnod->_hptr  = nullptr ;
            _rnod->_lptr  = nullptr ;
            _rnod->_pptr  = _pnod ;

            init_aabb_node (_rnod , _pnod);

        /*---------- partition list - split on left|right */
            for(_iptr  = _cptr ; _iptr != nullptr ;
                _iptr  = _next )
            {
                _next  = _iptr-> _next ;
        /*------------------------------- partition items */
                if (_iptr->_data
                    .pmid (_bdim) > _spos)
                {
                    push_aabb_node( _rnod, _iptr) ;
                    push_item(_rptr,_iptr);

                   _rnum += +1 ;
                }
                else
                {
                    push_aabb_node( _lnod, _iptr) ;
                    push_item(_lptr,_iptr);

                   _lnum += +1 ;
                }
            }

            _rnod->_hptr = _rptr ;
            _rnod->_size = _rnum ;

            _lnod->_hptr = _lptr ;
            _lnod->_size = _lnum ;

        /*------------------ push new children onto stack */
            if (_cnum < this->_npop)
            {
                double _volp, _voll, _volr ;

                _volp = +1. ;
                _voll = +1. ;
                _volr = +1. ;

                for (auto _idim = _dims; _idim-- != +0; )
                {
        /*--------------- parent // child (hyper) volumes */
                    _volp *=
                    _pnod->_pmax[_idim] -
                    _pnod->_pmin[_idim] ;

                    _voll *=
                    _lnod->_pmax[_idim] -
                    _lnod->_pmin[_idim] ;

                    _volr *=
                    _rnod->_pmax[_idim] -
                    _rnod->_pmin[_idim] ;
                }

                if (_voll + _volr <= _volp * this->_vtol)
                {
        /*--------------- push children due to vol. ratio */
                if (_lnum < _cnum && _rnum < _cnum)
                {
                    _work.push_tail (_rnod);
                    _work.push_tail (_lnod);
                }
                }

            }
            else
            {
        /*--------------- push children due to items cap. */
                if (_lnum < _cnum && _rnum < _cnum)
                {
                    _work.push_tail (_rnod);
                    _work.push_tail (_lnod);
                }
            }
        }

    }

/*-------- _pop items from tree while preserving geometry */
    __normal_call void_type _pop (
        )
    {
        if (this->_root == nullptr) return ;

        this->_work.clear() ;
        this->_work.push_tail(this->_root) ;

    /*---------------------------- nullify lists per node */
        for ( ; !this->_work.empty() ; )
        {
            node_type *_node = nullptr ;
            this->_work._pop_tail(_node) ;

            if (_node->items () != nullptr)
            {
                _node->_hptr = nullptr ;
            }

            if (_node->lower(0) != nullptr)
            {
                {
                    this->_work.push_tail (
                        _node->lower( 0)) ;

                    this->_work.push_tail (
                        _node->lower( 1)) ;
                }
            }
        }
    /*---------------------------- clear item linked list */
        this->_item_base.clear() ;
    }

/*-------- form a biased, spatially-local insertion order */
    template <
        typename  iptr_list
             >
    __normal_call void_type brio (
        iptr_list &_iset  // return brio indices for item
        )
    {
        containers::array<item_data*>_next ;

        if (this->_root == nullptr) return ;

        this->_work.clear() ;
        this->_work.push_tail(this->_root) ;

    /*---------------------------- init. leading item ptr */
        for ( ; !this->_work.empty() ; )
        {
            node_type *_node = nullptr ;
            this->_work._pop_tail(_node) ;

            if (_node->items () != nullptr)
            {
                _next.push_tail(
                    _node->items () ) ;
            }

            if (_node->lower(0) != nullptr)
            {
                uint32_t  _rsiz =
                    sizeof(real_type) * +K ;
                uint32_t  _usiz =
                    sizeof(uint32_t ) * +1 ;

                uint32_t  _hash ;
                _hash = hash::hashword (
               (uint32_t*)_node->_pmin,
                    _rsiz / _usiz, +13) ;

                if (_hash % 2 == +0 )
                {
                    this->_work.push_tail (
                        _node->lower( 0)) ;

                    this->_work.push_tail (
                        _node->lower( 1)) ;
                }
                else
                {
                    this->_work.push_tail (
                        _node->lower( 1)) ;

                    this->_work.push_tail (
                        _node->lower( 0)) ;
                }
            }
        }

    /*---------------------------- build order from lists */
        bool_type _push  = true ;
        for (auto _ilim  = +2 ;
                _push; _ilim *= 4 )
        {
            _push = false;
        for (auto _iter  = _next.head() ;
                  _iter != _next.tend() ;
                ++_iter  )
        {
            iptr_type _inum  = +0 ;
            for (auto _ipos  = *_iter ;
                      _ipos != nullptr;
                _ipos  = _ipos->_next )
            {
                _push  =  true ;

                _iset.push_tail(
                    _ipos->_data.ipos());

               *_iter  = _ipos->_next ;

                if (++_inum>_ilim) break;
            }
        }
        }
    }

/*-------- get items-per-nodes, as a sparse list-of-lists */
    template <
        typename  lptr_list,
        typename  iptr_list
             >
    __normal_call void_type dcmp (
        lptr_list &_lptr,  // ith list: _lp[i]:_lp[i+1]-1
        iptr_list &_iset   // unrolled list of item index
        )
    {
        _lptr.push_tail(_iset.count()) ;

        if (this->_root == nullptr) return;

        this->_work.clear() ;
        this->_work.push_tail(this->_root);

        for ( ; !this->_work.empty() ; )
        {
            node_type *_node = nullptr ;
            this->_work._pop_tail(_node) ;

    /*---------------------------- push items to CCS list */
            if (_node->items () != nullptr)
            {
                for (auto
                    _item  = _node->items() ;
                          _item != nullptr;
                    _item  = _item->_next )
                {
                    _iset.push_tail(
                        _item->_data.ipos());
                }
                _lptr.
                 push_tail(_iset.count()) ;
            }

    /*---------------------------- descend to child nodes */
            if (_node->lower(0) != nullptr)
            {
                uint32_t  _rsiz =
                    sizeof(real_type) * +K ;
                uint32_t  _usiz =
                    sizeof(uint32_t ) * +1 ;

                uint32_t  _hash ;
                _hash = hash::hashword (
               (uint32_t*)_node->_pmin,
                    _rsiz / _usiz, +13) ;

                if (_hash % 2 == +0 )
                {
                    this->_work.push_tail (
                        _node->lower( 0)) ;

                    this->_work.push_tail (
                        _node->lower( 1)) ;
                }
                else
                {
                    this->_work.push_tail (
                        _node->lower( 1)) ;

                    this->_work.push_tail (
                        _node->lower( 0)) ;
                }
            }
        }
    }

/*-------- helper: calc. "distance" between pos. and aabb */
    __static_call
    __normal_call real_type calc_rect_dist (
      __const_ptr(real_type) _ppos,
      __const_ptr(real_type) _bmin,
      __const_ptr(real_type) _bmax
        )
    {
        real_type _dist = (real_type)0.;
        for (auto _idim = _dims; _idim-- != +0; )
        {
        if (_ppos[_idim] < _bmin[_idim])
        {
            real_type _dloc =
            _bmin[_idim] - _ppos[_idim];
            _dist +=  _dloc ;
        }
        else
        if (_ppos[_idim] > _bmax[_idim])
        {
            real_type _dloc =
            _ppos[_idim] - _bmax[_idim];
            _dist +=  _dloc ;
        }
        }

        return ( _dist * _dist ) ;
    }

/*-------- search collection via recursive aabb traversal */
    template <
        typename  tree_pred ,// tree intersections
        typename  item_pred  // item intersections
             >
    __normal_call bool_type find (
        tree_pred &_pred,
        item_pred &_fout
        )
    {
        if (this->_root == nullptr) return false ;

        containers::array<
                node_type* ,
                allocator>                 _wset ;

    /*----------------- maintain stack of unvisited nodes */
        _wset.set_alloc (+64)  ;
        if (_pred(this->_root->_pmin,
                  this->_root->_pmax)
                 )
        _wset.push_tail (this->_root)  ;

    /*----------------- traverse while _pred remains true */
        bool_type _find =  false ;
        for ( ; !_wset.empty() ; )
        {
            node_type *_node = nullptr ;
            _wset._pop_tail(_node) ;

            if (_node->_hptr!= nullptr )
            {
                _find = true ;
        /*-------------------- leaf: push onto output */
                _fout( _node->_hptr) ;
            }

            if (_node->lower(0) != nullptr)
            {
        /*-------------------- traverse into children */
            if (_pred(
                _node->lower(0)->_pmin ,
                _node->lower(0)->_pmax )
                     )
                _wset.push_tail (
                      _node->lower(0)) ;

            if (_pred(
                _node->lower(1)->_pmin ,
                _node->lower(1)->_pmax )
                     )
                _wset.push_tail (
                      _node->lower(1)) ;
            }
        }

        return ( _find ) ;
    }

/*-------- check for nearsest in collection via traversal */
    template <
        typename  projector
             >
    __normal_call bool_type near (
      __const_ptr(real_type) _ppos ,  // point to index
        projector &_proj    // closest distance projector
        )
    {
        class node_dist
            {
    /*----------------------------- node/dist type for PQ */
            public  :
                real_type       _dsqr ;
                node_type      *_node ;
            } ;

        class node_pred
            {
    /*----------------------------- node/dist less for PQ */
            public  :
            __inline_call
                bool_type operator () (
                node_dist const&_adat ,
                node_dist const&_bdat
                ) const
            {   return ( _adat. _dsqr <
                         _bdat. _dsqr ) ;
            }
            } ;

        if (this->_root == nullptr) return false ;

    /*----------------- maintain stack of unvisited nodes */
        containers::priorityset<
            node_dist ,
            node_pred    > _nnpq ;

        _nnpq.set_alloc(+64) ;

        node_dist _ndat ;
        _ndat._node =  _root ;
        _ndat._dsqr = calc_rect_dist(_ppos,
           &_root->_pmin[ 0] ,
           &_root->_pmax[ 0]);
        _nnpq.push    (_ndat);

    /*----------------- traverse tree while len. reducing */
        bool_type _find = false ;

        real_type _dsqr =
            +std::numeric_limits
                <real_type>::infinity() ;

        for ( ; !_nnpq.empty()  ; )
        {
        /*------------------------ test next closest node */
            _nnpq._pop_root (_ndat) ;

            if (_ndat._dsqr<=_dsqr)
            {
        /*------------------------ descend if maybe close */

            if (_ndat.
                _node-> _hptr   != nullptr)
            {
        /*------------------------ leaf: update item-dist */
                _find =  true ;

                _dsqr = _proj(
                        _ndat._node ,
                        _ndat._node->_hptr) ;
            }

            if (_ndat.
                _node->lower(0) != nullptr)
            {
        /*------------------------ traverse into children */
                node_type*_inod =
                    _ndat._node->lower(0) ;

                node_type*_jnod =
                    _ndat._node->lower(1) ;

                _ndat._node =  _inod ;
                _ndat._dsqr =
                    calc_rect_dist(_ppos,
                       &_inod->_pmin[ 0],
                           &_inod->_pmax[ 0]) ;

                _nnpq.push(_ndat)  ;

                _ndat._node =  _jnod ;
                _ndat._dsqr =
                    calc_rect_dist(_ppos,
                       &_jnod->_pmin[ 0],
                           &_jnod->_pmax[ 0]) ;

                _nnpq.push(_ndat)  ;
            }

            }
        }

    /*---------------------------- must have found a node */
        return ( _find )  ;
    }

    } ;


    }

#   endif   //__AABB_TREE_K__



