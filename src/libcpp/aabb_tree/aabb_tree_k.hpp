
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
 * Disclaimer:  Neither I nor: Columbia University, The
 * Massachusetts Institute of Technology, The 
 * University of Sydney, nor The National Aeronautics
 * and Space Administration warrant this code in any 
 * way whatsoever.  This code is provided "as-is" to be 
 * used at your own risk.
 *
------------------------------------------------------------
 *
 * Last updated: 03 October, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __AABB_TREE_K__
#   define __AABB_TREE_K__

    namespace geom_tree {
    
    template <
    typename I ,
    size_t   K ,
    typename N = aabb_node_base_k,
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
    
    iptr_type static constexpr         _dims =  K * +1 ;

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
    /*------------------------------- tree-based pointers */
        node_type          *_pptr;
        node_data          *_lptr;
    /*------------------------------- marker: search algs */
        iptr_type           _mark;
        iptr_type           _size;
    /*----------------------------- tree-pointer handlers */
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

    typedef     node_type*              node_ptrt ;

    typedef containers::array      <      
                node_ptrt,
                allocator>              work_list ;

    public  :
/*------------------------------------- tree root pointer */
    node_type                   _rdat ;
    node_type                  *_root ;
    iptr_type                   _size ;

    work_list                   _work ;

/*------------------------------------- pool'd allocators */
    base_pool              _node_base ;
    base_pool              _item_base ;
    
    node_pool              _node_pool ;
    item_pool              _item_pool ;

/*------------------------------------- tree shape params */
    iptr_type                   _imax ;
    real_type                   _long ;
    real_type                   _vtol ;
        
    private :
/*---------------------- helper - delete a block of nodes */
    __inline_call void_type free_node (
        node_data*&_data
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
        node_data*&_data
        )
    {
       _data = 
        this->_node_pool.allocate (1) ;
        this->
       _node_pool.construct (_data) ;
    }

/*------------------------------- helper - delete an item */
    __inline_call void_type free_item (
        item_data *_idat
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
        item_data     *&_idat
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
        item_data*&_head,
        item_data *_item
        )
    {
        _item->_next = _head ;
        _head        = _item ;
    }
    
/*--------------------- helper - _pop an item from a list */
    __static_call
    __inline_call void_type _pop_item (
        item_data*&_head,
        item_data *_item,
        item_data *_prev
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
        for(auto _idim = _dims ; _idim-- != +0; )
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
        iptr_type _best = 0;
        for(auto _idim = _dims ; _idim-- != +0; )
        {
            real_type _alen;
            _alen = _node->_pmax[_idim]
                  - _node->_pmin[_idim];
        
            _rect[_idim]._alen = _alen ;
            _rect[_idim]._axis = _idim ; 
        }
    
    /*------------------------------ sort aabb dimensions */
        algorithms::ssort (
            &_rect[0], &_rect[_dims], rect_less()) ;
   
    /*------------------------------ scan aabb dimensions */
        for(auto _idim = _dims ; _idim-- != +0; )
        {  
            iptr_type _axis ;
            _axis = _rect[_idim]._axis ;
            real_type _long ; 
            _long = _rect[_idim]._alen ;
            
            _long *= this->_long ;
            
            iptr_type _pnum = +0 ;
            iptr_type _cnum = +0 ;
        
        /*------------------ partition "long"/inner items */
            for(item_data*_iptr  = _node->_hptr; 
                          _iptr !=  nullptr ; 
                          _iptr  = _iptr->_next)
            {
                if (_iptr->
                    _data.plen(_axis) >= +_long)
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
                _blen = _long ;
                
            if (_pnum == +0) break ;
            }                   
        }
    }
    
/*--- top-down assembly of aabb-tree via recursive splits */
    template <
    typename      iter_type
             >
    __normal_call void_type load (
        iter_type _head ,
        iter_type _tend ,
        iptr_type _imax = + 32 ,
        real_type _long = +.75 ,
        real_type _vtol = +.50
        )
    {
        this->_root = &this->_rdat ;

        this->_size = +1 ;

    /*------------------------------ set node fill params */
        this->_imax = _imax ;
        this->_long = _long ;
        this->_vtol = _vtol ;

    /*------------------------------ set null "tree" ptrs */
        this->_root->_pptr = nullptr;
        this->_root->_lptr = nullptr;
        this->_root->_hptr = nullptr;

    /*------------------------------ push items onto root */
        item_data *_hptr   = nullptr;
        item_data *_idat   = nullptr;
        {
        iptr_type  _rsiz = +0 ;
        for(; _head != _tend; ++_head, ++_rsiz)
        {
    /*------------------------------ alloc. and push item */
            make_item(*_head, _idat);
            push_item( _hptr, _idat);
        }
        this->_root->_hptr  = _hptr ;
        this->_root->_size  = _rsiz ;
        }
        
        init_aabb_node (this->_root);
        
    /*-------------------------- a list of un-split nodes */
        this->_work.set_count( +0);
        this->_work.
            push_tail(this->_root);
        
    /*--- refine tree until all nodes satisfy constraints */
        for( ; !this->_work.empty() ; )
        {
            node_type *_pnod = nullptr;
            node_data *_ndat = nullptr;
            node_type *_rnod = nullptr;
            node_type *_lnod = nullptr;
            
            real_type  _blen ;      
            iptr_type  _bdim = -1 ;
            iptr_type  _pnum = +0 ;
            iptr_type  _cnum = +0 ;
            iptr_type  _lnum = +0 ;
            iptr_type  _rnum = +0 ;
            
            item_data *_next = nullptr;
            item_data *_pptr = nullptr;
            item_data *_cptr = nullptr;
            item_data *_lptr = nullptr;
            item_data *_rptr = nullptr;
            
        /*-------------------------- _pop node from stack */
            this->_work._pop_tail(_pnod) ;

        /*---------- find best "split" axis for this node */
            find_best_axis (_pnod, _bdim , 
                            _blen) ;
            
            if (_bdim == -1) continue ;
                   
        /*---------- partition list - remove "long" items */
            for(item_data *_iptr  = _pnod->
                                    _hptr ; 
                           _iptr != nullptr ; 
                           _iptr  = _next )
            {
                _next = _iptr->_next ;
        /*------------------------------- partition items */
                if (_iptr->_data
                    .plen (_bdim) > _blen)
                {
                    push_item(_pptr,_iptr);
                   _pnum += +1 ;
                }
                else
                {
                    push_item(_cptr,_iptr); 
                   _cnum += +1 ;
                }
            }
        
        /*---------- split pos. - mean of non-long aabb's */
            real_type _spos = (real_type)+0.;
            
            for(item_data *_iptr  = _cptr ; 
                           _iptr != nullptr ; 
                           _iptr  = _next )
            {
                _next = _iptr->_next ;
         
                _spos +=_iptr->
                        _data.pmid (_bdim);
            }
            
            _spos  = _spos / _cnum ;
  
        /*---------- partition list - split on left|right */
            for(item_data *_iptr  = _cptr ; 
                           _iptr != nullptr ; 
                           _iptr  = _next )
            {
                _next = _iptr->_next ;
        /*------------------------------- partition items */
                if (_iptr->_data
                    .pmid (_bdim) > _spos)
                {
                    push_item(_rptr,_iptr); 
                   _rnum += +1 ;
                }
                else
                {
                    push_item(_lptr,_iptr); 
                   _lnum += +1 ;
                }
            }
        
        /*-------------------- alloc. new child node data */
            make_node(_ndat);
            _lnod = &_ndat->_node[ 0] ;
            _rnod = &_ndat->_node[ 1] ;
            
            this->_size  += +2 ;
            
            _pnod->_hptr  = _pptr ;
            _pnod->_size  = _pnum ;
            _pnod->_lptr  = _ndat ;
            
            _lnod->_hptr  = _lptr ;
            _lnod->_size  = _lnum ;
            _lnod->_lptr  = nullptr ;
            _lnod->_pptr  = _pnod ;
            
            _rnod->_hptr  = _rptr ;
            _rnod->_size  = _rnum ;
            _rnod->_lptr  = nullptr ;
            _rnod->_pptr  = _pnod ;
        
            init_aabb_node (_lnod, _pnod);
            init_aabb_node (_rnod, _pnod);
            
        /*------------------ push new children onto stack */
            if (_cnum < this->_imax)
            {       
                real_type _volp, _voll, _volr ;
        
                _volp = (real_type) +1. ;
                _voll = (real_type) +1. ;
                _volr = (real_type) +1. ;
        
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
        
                if (_voll + _volr <= this->_vtol * _volp)
                {
        /*--------------- push children due to vol. ratio */
                    this->_work
                       .push_tail(_rnod);
                    this->_work
                       .push_tail(_lnod);
                }
                
            }
            else
            {
        /*--------------- push children due to items cap. */
                    this->_work
                       .push_tail(_rnod);
                    this->_work
                       .push_tail(_lnod);         
            }
        }
    
    }
    
/*------- helper: calc. "distance" between point and AABB */
    __normal_call real_type calc_rect_dist (
        real_type *_ppos,
        real_type *_bmin,
        real_type *_bmax
        )
    {
        real_type _dist = (real_type)0.;
        for (auto _idim = _dims; _idim-- != +0; )
        {
        if (_ppos[_idim] < _bmin[_idim])
        {
            real_type _dloc = 
            _bmin[_idim] - _ppos[_idim];
        
            _dist = 
                std::max (_dist, _dloc);
        }
        else
        if (_ppos[_idim] > _bmax[_idim])
        {
            real_type _dloc = 
            _ppos[_idim] - _bmax[_idim];
        
            _dist = 
                std::max (_dist, _dloc);
        }
        }
        
        return ( _dist ) ;
    }
    
/*-------- search collection via recursive aabb traversal */
    template <
    typename      tree_pred ,// tree intersections
    typename      item_pred  // item intersections
             >
    __normal_call bool_type find (
        tree_pred &_pred,
        item_pred &_fout
        )
    {
        if (this->_root == nullptr) return false ;

    /*----------------- maintain stack of unvisited nodes */
        this->_work.set_count( +0);
        this->_work.
            push_tail(this->_root);

    /*------------ traverse tree while _pred remains true */
        bool_type _find =  false ;
        for ( ; !this->_work.empty() ; )
        {
            node_type *_node = nullptr ;
            this->_work.
                _pop_tail(_node) ;

        /*------------- eval. node intersection predicate */
            if (_pred(_node->_pmin, _node->_pmax))
            {
                if (_node->_hptr    != nullptr)
                {
            /*-------------------- leaf: push onto output */
                    _fout( _node->_hptr) ;

                    _find = true ;
                }

                if (_node->lower(0) != nullptr)
                {
            /*-------------------- traverse into children */
                this->_work.push_tail (
                        _node->lower(0)) ;
                this->_work.push_tail (
                        _node->lower(1)) ;
                }
            }
        }
        
        return ( _find ) ;
    }
    
/*-------- check for nearsest in collection via traversal */
    __normal_call bool_type near (
        real_type *_ppos ,
        item_data*&_near
        )
    {
        class node_dist
            {
    /*----------------------------- node/dist type for PQ */
            public  :
                real_type       _dist ;
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
            {   return ( _adat. _dist <
                         _bdat. _dist ) ;
            }
            } ;
    
        bool_type _find = false;
    
        if (this->_root 
                == nullptr) return _find;
  
        real_type _dist = 
            +std::numeric_limits
                <real_type>::infinity() ;
        
        real_type _dloc = 
            -std::numeric_limits
                <real_type>::infinity() ;

    /*----------------- maintain stack of unvisited nodes */
        containers::priorityset<
            node_dist ,
            node_pred    > _nnpq ;
    
        node_dist _ndat ;
        _ndat._node =  _root ;
        _ndat._dist = 
        calc_rect_dist(_ppos ,
           &_root->_pmin[0],
               &_root->_pmax[0]) ;
        _nnpq.push    (_ndat);
        
    /*----------------- traverse tree while len. reducing */
        for ( ; !_nnpq.empty()  ; )
        {
        /*------------------------ test next closest node */
            _nnpq._pop_root (_ndat) ;

            if (_ndat._dist<=_dist)
            {
        /*------------------------ descend if maybe close */
            
            if (_ndat.
                _node->_hptr  != nullptr )
            {
        /*------------------------ leaf: update item-dist */
                for (item_data  *_iptr = 
                    _ndat._node->_hptr ; 
                        _iptr != nullptr ; 
                    _iptr = _iptr->_next )
                {
                    _dloc = 
                    calc_rect_dist(_ppos ,
                   &_iptr->_data.pmin(0) ,
                   &_iptr->_data.pmax(0)
                        ) ;
                  
                    if (_dloc < _dist)
                    {
                /*---------------- keep track of min-len. */
                        _dist = _dloc;
                        _near = _iptr;
                    }
                }

                _find =  true ;
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
                _ndat._dist = 
                calc_rect_dist(_ppos ,
                   &_inod->_pmin[ 0],
                       &_inod->_pmax[ 0]) ;
                _nnpq.push(_ndat)  ;
                
                _ndat._node =  _jnod ;
                _ndat._dist = 
                calc_rect_dist(_ppos ,
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



