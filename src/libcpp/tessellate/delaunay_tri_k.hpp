
    /* 
    --------------------------------------------------------
     * a "k"-dimensional delaunay triangulation object.
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
     * Last updated: 02 June, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __DELAUNAY_TRI_K__
#   define __DELAUNAY_TRI_K__

    namespace mesh {
        
    /*
    --------------------------------------------------------
     * DELAUNAY-TRI-K: a k-dim. delaunay tri. object.
    --------------------------------------------------------
     * NODE-TYPE - node-type class, from DEL-TRI-NODE-K.
     * TRIA-TYPE - tria-type class, from DEL-TRI-TRIA-K.
     * TRIA-PRED - predicate class, from DEL-TRI-PRED-K.
     * ALLOCATOR - memory allocator.
    --------------------------------------------------------
     */

    template <
    typename N ,
    typename T ,
    typename P ,
    typename A = allocators::basic_alloc 
             >
    class delaunay_tri_k
    { 
/*------------------------- delaunay triangulation in R^k */
    public  :
    
    typedef N                               node_type ;
    typedef T                               tria_type ;
    typedef A                               allocator ;
    typedef P                               tria_pred ;

    typedef typename 
            node_type::real_type            real_type ;
    typedef typename 
            tria_type::iptr_type            iptr_type ;

    typedef char_type                       fptr_type ;

/*------------------------- alternative to "static const" */
    
    iptr_type static constexpr null_mark =  +0 ;
    iptr_type static constexpr tria_dims =
                               tria_pred::_dims;
    
    __static_call
    __inline_call iptr_type null_flag (
        ) 
    {   return std::numeric_limits<iptr_type>::max()- 2 ;
    }

    typedef mesh::delaunay_tri_k    <
            node_type, 
            tria_type, 
            tria_pred, 
            allocator   >                   self_type ;

    typedef typename 
            tria_pred::face_ptrs            face_ptrs ;
    typedef typename 
            tria_pred::face_hash            face_hash ;
    typedef typename 
            tria_pred::face_same            face_same ;

    iptr_type static constexpr pool_byte_size = 96*1024 ;

    typedef allocators::_pool_alloc <
            allocators::basic_alloc ,
                pool_byte_size      >       pool_base ;

    typedef allocators::_wrap_alloc <
                pool_base           >       pool_wrap ;

    typedef containers::hash_table  <
                face_ptrs ,
                face_hash ,
                face_same ,
                pool_wrap           >       hash_list ;

    typedef containers::block_array <
                node_type , 
                allocator           >       node_list ;

    typedef containers::block_array <
                tria_type , 
                allocator           >       tria_list ;

    typedef containers::array       <
                iptr_type , 
                allocator           >       iptr_list ;

    public  :
    
    pool_base                     _fpol ;

    hash_list                     _fset ;

    node_list                     _nset ;
    tria_list                     _tset ;
    
    iptr_list                     _ftri ;
    iptr_list                     _fnod ;

    iptr_list                     _work ;

    public  :

#   define __mark(__ti) tria((__ti))->mark()
           
    /*
    --------------------------------------------------------
     * NODE/TRIA "get" and "set" access.
    --------------------------------------------------------
     */
     
    __inline_call node_type      * node (
        iptr_type  _ipos
        )
    {   return &this->_nset[_ipos] ;
    }
    
    __inline_call tria_type      * tria (
        iptr_type _ipos
        )
    {   return &this->_tset[_ipos] ;
    }
    
    __inline_call node_type const* node (
        iptr_type  _ipos
        ) const
    {   return &this->_nset[_ipos] ;
    }
    
    __inline_call tria_type const* tria (
        iptr_type _ipos
        ) const
    {   return &this->_tset[_ipos] ;
    }
       
    /*
    --------------------------------------------------------
     * PUT-NODE: "delete" node from active set.
    --------------------------------------------------------
     */
 
    __inline_call void_type _put_node (
        iptr_type _ipos
        )
    {
        this->_fnod.push_tail(_ipos);
        this->_nset[_ipos].mark()=-1;
    }
    
    /*
    --------------------------------------------------------
     * PUT-TRIA: "delete" tria from active set.
    --------------------------------------------------------
     */
 
    __inline_call void_type _put_tria (
        iptr_type _ipos
        )
    {
        this->_ftri.push_tail(_ipos);
        this->_tset[_ipos].mark()=-1;
    }

    /*
    --------------------------------------------------------
     * GET-NODE: "create" node into active set.
    --------------------------------------------------------
     */
     
    __inline_call iptr_type _get_node (
        )
    {
        iptr_type _ipos = -1;
        if (this->_fnod.count() != +0 )
        {
    /*------------------------ recycle from free list */
            this->
           _fnod._pop_tail(_ipos) ;
        }
        else
        {
    /*------------------------ alloc. from underlying */
            _ipos = (iptr_type)
                this->_nset.count() ;
            this->_nset.push_tail() ;
        }

        node(_ipos)->mark() = +0 ;

        return ( _ipos ) ;
    }
    
    /*
    --------------------------------------------------------
     * GET-TRIA: "create" tria into active set.
    --------------------------------------------------------
     */
     
    __inline_call iptr_type _get_tria (
        )
    {
        iptr_type _ipos = -1;
        if (this->_ftri.count() != +0 )
        {
    /*------------------------ recycle from free list */
            this->
           _ftri._pop_tail(_ipos) ;
        }
        else
        {
    /*------------------------ alloc. from underlying */
            _ipos = (iptr_type)
                this->_tset.count() ;
            this->_tset.push_tail() ;
        }

        tria(_ipos)->mark() = +0 ;

        return ( _ipos ) ;
    }
    
    /*
    --------------------------------------------------------
     * FIND-PAIR: return adjacent tria/face indexing.
    --------------------------------------------------------
     */
     
    __inline_call void_type find_pair (
        iptr_type  _ipos,
        iptr_type &_jpos,
        iptr_type  _ifac,
        iptr_type &_jfac,
        iptr_type &_flag
        ) const
    {
    /*---------------- unwind face adjacency pointers */
        _jfac = __unflip(
            tria(_ipos)->fpos(_ifac)) ;
        
        _jpos = __unflip(
            tria(_ipos)->next(_ifac)) ;
        
        _flag = __isflip(
            tria(_ipos)->next(_ifac)) 
              ? -1 : +1 ;
    }  
    
    /*
    --------------------------------------------------------
     * PUSH-PAIR: encode adjacent tria/face indexing.
    --------------------------------------------------------
     */
    
    __inline_call void_type push_pair (
        iptr_type  _itri,
        iptr_type  _jtri,
        iptr_type  _ifac,
        iptr_type  _jfac,
        iptr_type  _flag
        )
    {
    /*---------------- encode face adjacency pointers */
        if (_itri != this->null_flag())
        {
            if (_flag >= +0)
        tria(_itri)->
            next(_ifac) = _jtri ;
            else
        tria(_itri)->
            next(_ifac)=__doflip(_jtri) ;
        
        tria(_itri)->
        fpos(_ifac) = (char_type)_jfac;
        }
        
        if (_jtri != this->null_flag())
        {
            if (_flag >= +0)
        tria(_jtri)->
            next(_jfac) = _itri ;
            else
        tria(_jtri)->
            next(_jfac)=__doflip(_itri) ;
        
        tria(_jtri)->
        fpos(_jfac) = (char_type)_ifac;
        }   
    } 
 
    /*
    --------------------------------------------------------
     * WALK-TRIA-LIST: breadth-first-search about tria.
    --------------------------------------------------------
     */
     
#   include "delaunay_walk_tria.inc" 
        
        
    /*
    --------------------------------------------------------
     * WALK-MESH-NODE: enclosing tria. via walking.
    --------------------------------------------------------
     */
 
#   include "delaunay_walk_node.inc"


    /*
    --------------------------------------------------------
     * STAR-TRIA-VOID: re-tria. star-shaped cavity.
    --------------------------------------------------------
     */
     
#   include "delaunay_star_void.inc"

     
    public  :
    
    /*
    --------------------------------------------------------
     * CONSTRUCT: create a del.-tri. obj. from alloc.
    --------------------------------------------------------
     */
    
    __inline_call  delaunay_tri_k   (
        allocator const& _asrc = allocator ()
        ) : _fpol(sizeof(
      typename hash_list::item_type)) ,

            _fset(face_hash(), 
                  face_same(), 
            +.9, (pool_wrap(&_fpol))) ,

            _nset( _asrc),
            _tset( _asrc),
            _ftri( _asrc),
            _fnod( _asrc), 
            _work( _asrc)
    {   tria_pred::exactinit() ;    // init. predicates
    }
    
    /*
    --------------------------------------------------------
     * PUSH-ROOT: push the bounding "super" triangle.
    --------------------------------------------------------
     */
    
    __normal_call void_type push_root (
        real_type *_pmin,
        real_type *_pmax
        )
    {
    /*---------------------------- de-alloc. existing */
        this->_nset.clear() ;
        this->_tset.clear() ;

    /*---------------------------- scale initial tria */
        real_type static constexpr _scal = 
            (real_type)+tria_pred::_dims ;

        real_type _pmid[tria_pred::_dims];
        real_type _pdel = (real_type)+.0 ;
        for (auto _idim = tria_pred::_dims + 0 ; 
                  _idim-- != + 0 ; )
        {
            real_type _xdel=_pmax[_idim] - 
                            _pmin[_idim] ;   
            
            _pdel = 
                std::max(_pdel,_xdel) ;
        
            _pmid[_idim] = 
            _pmax[_idim] *(real_type)+.5 +
            _pmin[_idim] *(real_type)+.5 ;
        }

    /*---------------------------- push tria indexing */
        iptr_type _itri = _get_tria() ;
        for (auto _inod = tria_pred::_dims + 1 ; 
                  _inod-- != + 0 ; )
        {
            tria(_itri)->node(_inod) =_inod;
            
            tria(_itri)->fpos(_inod) = -1;
            
            tria(_itri)->next(_inod) =
                __doflip(this->null_flag());
        }

    /*------------------------- push node coordinates */
        for (auto _inod = tria_pred::_dims + 1 ; 
                  _inod-- != + 0 ; )
        {
        iptr_type _jnod = _get_node() ;
        node(_jnod)->next() = _itri;
        for (auto _idim = tria_pred::_dims + 0 ; 
                  _idim-- != + 0 ; )
        {
            if (_idim != _jnod - 1 )
                node(_jnod)->pval(_idim) = 
                    _pmid[_idim] - 
               (real_type).5*_pdel ;
            else
                node(_jnod)->pval(_idim) = 
                    _pmid[_idim] + 
               (real_type)1.*_pdel*_scal ;
        }
        }
    }
    
    /*
    --------------------------------------------------------
     * FIND-NODE: find enclosing triangle for a point.
    --------------------------------------------------------
     */
     
    __inline_call bool_type find_node (
        real_type *_ppos,
        iptr_type &_tpos,
        iptr_type  _hint  = null_flag()
        )
    {   
        return walk_mesh_node (
                   _ppos, _tpos, _hint) ;
    }

    /*
    --------------------------------------------------------
     * CIRC-LIST: find "conflict" triangles for point.
    --------------------------------------------------------
     */
    
    template <
    typename      list_type
             >
    __normal_call bool_type circ_list (
        real_type *_ppos,
        iptr_type &_elem,
        list_type &_list,
        iptr_type  _hint  = null_flag()
        )
    {
        this->_work.clear();
    /*--------------------------- find enclosing tria */
        if (walk_mesh_node(_ppos, _elem, _hint))
        {
    /*------------------------------ bfs about cavity */
        typename tria_pred:: 
        template circ_pred<
             self_type >_pred( _ppos) ;
             
        walk_tria_list (_elem, +1, 
                        _pred, _work) ;
                        
    /*------------------------------ push index lists */
        _list.push_tail(_work.head(),
                        _work.tend()) ;
                        
        return (  true ) ;
        }
        
        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * PUSH-NODE: push a new node into a triangulation.
    --------------------------------------------------------
     */
        
    __inline_call bool_type push_node (
        real_type *_ppos,
        iptr_type &_node,
        iptr_type  _hint  = null_flag()
        )
    {
        bool_type _pass   = push_node (
                _ppos, _node, _hint,
                (iptr_list*)nullptr,
                (iptr_list*)nullptr,
                (iptr_list*)nullptr)  ;
        
        return (  _pass )  ;
    }

    template <
    typename      list_type
             >
    __normal_call bool_type push_node (
        real_type *_ppos,
        iptr_type &_node,
        iptr_type  _hint  = null_flag() ,
        list_type *_tnew  = nullptr   ,
        list_type *_told  = nullptr   ,
        list_type *_circ  = nullptr
        )
    {
        this->_work.clear ();

    /*--------------------------- _find enclosing element */
        iptr_type _elem = -1;
        if (walk_mesh_node(_ppos, _elem, _hint))
        {

    /*--------------------------- push new node onto list */
        _node = _get_node() ;
        for (auto _idim = tria_pred::_dims + 0 ; 
                  _idim-- != +0 ; )
        {
        node(_node)->pval(_idim) = _ppos[_idim];
        }   

    /*--------------------------- grab enclosing indexing */
        iptr_type _tnod[ +1 + tria_pred::_dims];
        for (auto _inod = tria_pred::_dims + 1 ; 
                  _inod-- != +0 ; )
        {
        _tnod[_inod] = tria(_elem)->node(_inod);
        }

    /*--------------------------- test for node duplicate */
        real_type _dist  = 
            std::numeric_limits
                <real_type>::infinity() ;
                
        for (auto _inod = tria_pred::_dims + 1 ; 
                  _inod-- != +0 ; )
        {
        _dist = std::min(_dist,
            tria_pred::lensqr_kd(_ppos,
                &node(_tnod[_inod])->pval(+0)));
        }

        if (_dist == (real_type)0.)
        {
            _put_node(_node) ;
        /*----------------------- bail-out on duplicates! */
            return ( false ) ;
        }

    /*-------------------- retriangulate enclosing cavity */
        typename tria_pred::
        template circ_pred<
            self_type> _pred( _ppos) ;
        
        if (_circ == nullptr)
        walk_tria_list(_elem, +1   , 
                       _pred, _work) ;
        else
        _work.push_tail(_circ->head(), 
                        _circ->tend()
                       ) ;
            
        star_tria_void(_work, _node, 
                   +1, _tnew, _told) ;

    /*-------------------- delaunay topology is recovered */
        return (  true ) ;
        }
        
        return ( false ) ;
    }
    
    /*
    --------------------------------------------------------
     * ROLL-BACK: "roll-back" an update.
    --------------------------------------------------------
     */
     
    template <
    typename      list_type
             >
    __normal_call void_type roll_back (
        list_type &_tnew,
        list_type &_told
        )
    {
        iptr_type _npos ;
        iptr_type _fpos ;
        for (auto _tpos  = _tnew.head();
                  _tpos != _tnew.tend();
                ++_tpos )
        {
        for (_npos = tria_pred::_dims+1; 
             _npos-- != +0 ; )
        {
            node(tria(*_tpos )
               ->node( _npos))
               ->next() = null_flag();
        }
        }
        
        for (auto _tpos  = _told.head();
                  _tpos != _told.tend();
                ++_tpos )
        {
        for (_npos = tria_pred::_dims+1; 
             _npos-- != +0 ; )
        {
            node(tria(*_tpos )
               ->node( _npos))
               ->next() = *_tpos ;
        }

        for (_fpos = tria_pred::_dims+1; 
             _fpos-- != +0 ; )
        {
            iptr_type  _tadj;
            iptr_type  _fadj;
            iptr_type  _flag;
            find_pair(*_tpos, _tadj , 
                _fpos, _fadj, _flag )  ;

            push_pair(*_tpos, _tadj , 
                       _fpos, _fadj , 
                       _flag) ;
        }
        }

        for (auto _tpos  = _tnew.head();
                  _tpos != _tnew.tend();
                ++_tpos )
        {
        for (_npos = tria_pred::_dims+1; 
             _npos-- != +0 ; )
        {
            if (node(tria(*_tpos )
                   ->node( _npos)) 
                   ->mark()>= +0 )
            if (node(tria(*_tpos )
                   ->node( _npos))
                   ->next() == null_flag())
            {
           _put_node(tria(*_tpos )
                   ->node( _npos)) ;
            }
        }
        }
        
        for (auto _tpos  = _tnew.head();
                  _tpos != _tnew.tend();
                ++_tpos )
        {
            _put_tria(*_tpos) ;
        }
    }
    
    /*
    --------------------------------------------------------
     * WALK-NODE: do breadth-first-search about node.
    --------------------------------------------------------
     */
    
    template <
    typename      tria_pred ,
    typename      list_type
             >
    __normal_call void_type walk_node (
        iptr_type  _node ,
        tria_pred  _pred ,
        list_type &_list
        )
    {
        if (node(_node)->
            next() != this->null_flag())
        {
        this->_work.clear()  ;
    /*----------------------------------- bfs about _node */
        walk_tria_list (
            node(_node)->next(), 
                 +1   , _pred, _work) ;

    /*----------------------------------- push index list */
        _list.push_tail(_work.head(),
                        _work.tend()) ;
        }
    }
    
    /*
    --------------------------------------------------------
     * WALK-TRIA: do breadth-first-search about tria.
    --------------------------------------------------------
     */
     
    template <
    typename      tria_pred,
    typename      list_type
             >
    __normal_call void_type walk_tria (
        iptr_type  _tria ,
        tria_pred  _pred ,
        list_type &_list
        )
    {
        this->_work.clear()  ;  
    /*----------------------------------- bfs about _tria */
        walk_tria_list (_tria,    +1, 
                        _pred, _work) ;
                        
    /*----------------------------------- push index list */
        _list.push_tail(_work.head(),
                        _work.tend()) ;
    }

 
#   undef  __mark   
    
    } ;
    
    }
    
#   endif  //__DELAUNAY_TRI_K__



