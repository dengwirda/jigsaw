
    /*
    --------------------------------------------------------
     * _ZIP-MESH-2: merge/split operations.
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
    --------------------------------------------------------
     *
     * Last updated: 11 Dec., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    /*
    --------------------------------------------------------
     * _ZIP-MESH: edge merge/split operations.
    --------------------------------------------------------
     */

    #include "iter_zips_2.inc"
    #include "iter_divs_2.inc"

    __static_call
    __normal_call void_type _zip_fcel (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        char_type  _kern ,
        real_list &_hval ,  // cached node spacing array
        iptr_list &_nset ,  // list of nodes active
        mark_list &_mark ,  // sets of node/cell flags
        iptr_type  _imrk ,  // outer iteration marker
        iter_opts &_opts ,  // user options
        real_type  _QLIM ,  // cell quality threshold
        real_type  _DLIM ,  // dual quality threshold
        iptr_type &_nzip ,  // number of successful zip's
        iptr_type &_ndiv ,  // number of successful div's
        iter_stat &_tcpu    // CPU timers/info
        )
    {
        class sort_data
            {
            public  :
    /*------------------------ tuple for edge re-ordering */
            iptr_type           _inod , _jnod ,
                                _knod ;
            float               _cost ;
            public  :
            __inline_call sort_data (
                iptr_type _isrc ,
                iptr_type _jsrc ,
                iptr_type _ksrc ,
                float     _csrc
                ) : _inod(_isrc), _jnod(_jsrc), 
                    _knod(_ksrc),
                    _cost(_csrc)   {}
            } ;
        class sort_less
            {
            public  :
    /*------------------------ less-than op. for cost-tup */
            __inline_call
                bool_type operator () (
                sort_data const&_idat ,
                sort_data const&_jdat
                ) const
            {  return  _idat._cost < _jdat._cost;
            }
            } ;

    #   define MARKNODE(_NODE) _mark._node[_NODE]

    #   define PUSHMARK              \
            init_mark( _mesh, _mark,        \
                std::max(_imrk - 1, +0));   \
            if (std::abs(           \
                _mark._node[_nnew])!= _imrk+1)  \
            {                       \
            if (_mark._node[_nnew] >= +0)   \
            {                       \
                _mark._node[_nnew] = +_imrk+1;  \
            }                       \
            else                    \
            {                       \
                _mark._node[_nnew] = -_imrk-1;  \
            }                       \
                _nset.push_tail(_nnew) ;    \
            }                       \

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*------------------------ add/pop nodes to fix topo. */
        typedef containers::
            array<sort_data> sort_list ;

        __unreferenced ( _DLIM );

        sort_list _sort ;
        conn_list _iadj , _jadj, _kadj ;
        conn_list _conn , _cset; 
        conn_list _iset , _jset, _kset ;
        real_list _qold , _qnew, _qtmp ;

        _nzip = +0 ; _ndiv = +0 ;

    //  assemble list of faces attached to "recent" nodes

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        _sort.set_alloc(
            _mesh.tri3().count()) ;

        for (auto _iter  = _mesh.tri3().head() ;
                  _iter != _mesh.tri3().tend() ;
                ++_iter  )
        {
             auto _inod =  _iter->node(0) ;
             auto _jnod =  _iter->node(1) ;
             auto _knod =  _iter->node(2) ;

             auto _iptr = _mesh.
             node().head()+_iter->node(0) ;
             auto _jptr = _mesh.
             node().head()+_iter->node(1) ;
             auto _kptr = _mesh.
             node().head()+_iter->node(2) ;

             auto _flag = _imrk - 2 ; // recent

            if (_iter->mark() >= +0 &&
               (    std::abs(
                _mark._node[_inod]) > _flag ||
                    std::abs(
                _mark._node[_jnod]) > _flag ||
                    std::abs(
                _mark._node[_knod]) > _flag ))
            {
                float _area  =
               (float)pred_type::tri3_mass (
                    & _iptr->pval(0) ,
                    & _jptr->pval(0) ,
                    & _kptr->pval(0) ) ;

                auto _ndeg = _iptr->ndeg() + 
                             _jptr->ndeg() +
                             _kptr->ndeg() ;

                _area *=  (float) _ndeg;

                _sort.push_tail(
                    sort_data(
                _inod, _jnod, _knod, _area)) ;
            }
        }

        if (_sort.empty()) return ;

        algorithms::qsort(      // sort edge list by lsqr
            _sort.head() ,
            _sort.tend() , sort_less());

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_zips += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    //  scan cells largest-to-smallest and try to div any
    //  unvisited cells

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _iter  = _sort.tail();
                  _iter != _sort.hend();
                --_iter  )
        {
    /*--------------------------- try to "div" local edge */
            iptr_type  _tadj, _tnod[3] ;
            _tnod[0] = _iter->_inod;
            _tnod[1] = _iter->_jnod;
            _tnod[2] = _iter->_knod;

            if (MARKNODE(_tnod[0]) < +0 &&
                MARKNODE(_tnod[1]) < +0 &&
                MARKNODE(_tnod[2]) < +0 ) continue ;

            if (MARKNODE(_tnod[0])>_imrk||
                MARKNODE(_tnod[1])>_imrk||
                MARKNODE(_tnod[2])>_imrk) continue ;

            if(!_mesh.find_tri3(
                          _tnod, _tadj) ) continue ;

            if (_opts.div_())
            {
    /*--------------------------- "div" for topo. + score */
                iptr_type  _nnew = -1;

                bool_type  _move = false ;
            //  _div_tri3( ...

                if (_move)
                {
                    PUSHMARK; _ndiv += +1;
                }
            }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_divs += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    //  scan cells smallest-to-largest and try to zip any
    //  unvisited cells

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _iter  = _sort.head();
                  _iter != _sort.tend();
                ++_iter  )
        {
    /*--------------------------- try to "zip" local edge */
            iptr_type  _tadj, _tnod[3] ;
            _tnod[0] = _iter->_inod;
            _tnod[1] = _iter->_jnod;
            _tnod[2] = _iter->_knod;

            if (MARKNODE(_tnod[0]) < +0 ||
                MARKNODE(_tnod[1]) < +0 ||
                MARKNODE(_tnod[2]) < +0 ) continue ;

            if (MARKNODE(_tnod[0])>_imrk||
                MARKNODE(_tnod[1])>_imrk||
                MARKNODE(_tnod[2])>_imrk) continue ;

            if(!_mesh.find_tri3(
                          _tnod, _tadj) ) continue ;

            if (_opts.zip_())
            {
    /*--------------------------- "zip" for topo. + score */
                iptr_type  _nnew = -1;

                bool_type  _move = false ;
                _zip_tri3( _geom, _mesh,
                    _hfun, _hval, _opts,
                    _mark. _node,
                    _imrk, _tadj,
                    _kern, _move, _nnew,
                    _iset, _jset, _kset,
                    _iadj, _jadj, _kadj, 
                    _cset, _qold, _qnew,
                    _qtmp, _QLIM, _tcpu) ;

                if (_move)
                {
                    PUSHMARK; _nzip += +1;
                }
            }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_zips += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    #   undef   PUSHMARK
    #   undef   MARKNODE
    }

    __static_call
    __normal_call void_type _zip_ecel (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        char_type  _kern ,
        real_list &_hval ,  // cached node spacing array
        iptr_list &_nset ,  // list of nodes active
        mark_list &_mark ,  // sets of node/cell flags
        iptr_type  _imrk ,  // outer iteration marker
        iter_opts &_opts ,  // user options
        real_type  _QLIM ,  // cell quality threshold
        real_type  _DLIM ,  // dual quality threshold
        iptr_type &_nzip ,  // number of successful zip's
        iptr_type &_ndiv ,  // number of successful div's
        iter_stat &_tcpu    // CPU timers/info
        )
    {
        class sort_data
            {
            public  :
    /*------------------------ tuple for edge re-ordering */
            iptr_type           _inod ;
            iptr_type           _jnod ;
            float               _cost ;
            public  :
            __inline_call sort_data (
                iptr_type _isrc ,
                iptr_type _jsrc ,
                float     _csrc
                ) : _inod(_isrc), _jnod(_jsrc),
                    _cost(_csrc)   {}
            } ;
        class sort_less
            {
            public  :
    /*------------------------ less-than op. for cost-tup */
            __inline_call
                bool_type operator () (
                sort_data const&_idat ,
                sort_data const&_jdat
                ) const
            {  return  _idat._cost < _jdat._cost;
            }
            } ;

    #   define MARKNODE(_NODE) _mark._node[_NODE]

    #   define PUSHMARK              \
            init_mark( _mesh, _mark,        \
                std::max(_imrk - 1, +0));   \
            if (std::abs(           \
                _mark._node[_nnew])!= _imrk+1)  \
            {                       \
            if (_mark._node[_nnew] >= +0)   \
            {                       \
                _mark._node[_nnew] = +_imrk+1;  \
            }                       \
            else                    \
            {                       \
                _mark._node[_nnew] = -_imrk-1;  \
            }                       \
                _nset.push_tail(_nnew) ;    \
            }                       \

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*------------------------ add/pop nodes to fix topo. */
        typedef containers::
            array<sort_data> sort_list ;

        __unreferenced ( _DLIM );

        sort_list _sort ;
        conn_list _aset , _bset, _cset ;
        conn_list _conn , _iset, _jset ;
        real_list _qold , _qnew, _qtmp ;

        _nzip = +0 ; _ndiv = +0 ;

    //  assemble list of edges attached to "recent" nodes

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        _sort.set_alloc(
            _mesh.edge().count()) ;

        for (auto _iter  = _mesh.edge().head() ;
                  _iter != _mesh.edge().tend() ;
                ++_iter  )
        {
             auto _inod =  _iter->node(0) ;
             auto _jnod =  _iter->node(1) ;

             auto _iptr = _mesh.
             node().head()+_iter->node(0) ;
             auto _jptr = _mesh.
             node().head()+_iter->node(1) ;

             auto _flag = _imrk - 2 ; // recent

            if (_iter->mark() >= +0 &&
               (    std::abs(
                _mark._node[_inod]) > _flag ||
                    std::abs(
                _mark._node[_jnod]) > _flag ))
            {
                float _lsqr  =
               (float)pred_type::length_sq (
                    & _iptr->pval(0) ,
                    & _jptr->pval(0) ) ;

                auto _ndeg = _iptr->ndeg() + 
                             _jptr->ndeg() ;

                _lsqr *=  (float) _ndeg;

                _sort.push_tail(
                 sort_data(_inod, _jnod, _lsqr)) ;
            }
        }

        if (_sort.empty()) return ;

        algorithms::qsort(      // sort edge list by lsqr
            _sort.head() ,
            _sort.tend() , sort_less());

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_zips += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    //  scan edges longest-to-shortest and try to div any
    //  unvisited edges

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _iter  = _sort.tail();
                  _iter != _sort.hend();
                --_iter  )
        {
    /*--------------------------- try to "div" local edge */
            iptr_type  _eadj, _enod[2] ;
            _enod[0] = _iter->_inod;
            _enod[1] = _iter->_jnod;

            if (MARKNODE(_enod[0]) < +0 &&
                MARKNODE(_enod[1]) < +0 ) continue ;

            if (MARKNODE(_enod[0])>_imrk||
                MARKNODE(_enod[1])>_imrk) continue ;

            if(!_mesh.find_edge(
                          _enod, _eadj) ) continue ;

            if (_opts.div_())
            {
    /*--------------------------- "div" for topo. + score */
                iptr_type  _nnew = -1;

                bool_type  _move = false ;
                _div_edge( _geom, _mesh,
                    _hfun, _hval, _opts,
                    _imrk, _eadj,
                    _kern, _move, _nnew,
                    _iset, _jset,
                    _aset, _bset,
                    _qold, _qnew,
                    _qtmp, _QLIM, _tcpu) ;

                if (_move)
                {
                    PUSHMARK; _ndiv += +1;
                }
            }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_divs += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    //  scan edges shortest-to-longest and try to zip any
    //  unvisited edges

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _iter  = _sort.head();
                  _iter != _sort.tend();
                ++_iter  )
        {
    /*--------------------------- try to "zip" local edge */
            iptr_type  _eadj, _enod[2] ;
            _enod[0] = _iter->_inod;
            _enod[1] = _iter->_jnod;

            if (MARKNODE(_enod[0]) < +0 ||
                MARKNODE(_enod[1]) < +0 ) continue ;

            if (MARKNODE(_enod[0])>_imrk||
                MARKNODE(_enod[1])>_imrk) continue ;

            if(!_mesh.find_edge(
                          _enod, _eadj) ) continue ;

            if (_opts.zip_())
            {
    /*--------------------------- "zip" for topo. + score */
                iptr_type  _nnew = -1;

                bool_type  _move = false ;
                _zip_edge( _geom, _mesh,
                    _hfun, _hval, _opts,
                    _mark. _node,
                    _imrk, _eadj,
                    _kern, _move, _nnew,
                    _iset, _jset,
                    _aset, _bset, _cset,
                    _qold, _qnew,
                    _qtmp, _QLIM, _tcpu) ;

                if (_move)
                {
                    PUSHMARK; _nzip += +1;
                }
            }
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._core_zips += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    #   undef   PUSHMARK
    #   undef   MARKNODE
    }

    __static_call
    __normal_call void_type _zip_mesh (
        geom_type &_geom ,  // geometry object
        mesh_type &_mesh ,  // mesh object
        hfun_type &_hfun ,  // spacing object
        char_type  _kern ,
        real_list &_hval ,  // cached node spacing array
        iptr_list &_nset ,  // list of nodes active
        mark_list &_mark ,  // sets of node/cell flags
        iptr_type  _imrk ,  // outer iteration marker
        iter_opts &_opts ,  // user options
        real_type  _QLIM ,  // cell quality threshold
        real_type  _DLIM ,  // dual quality threshold
        iptr_type &_nzip ,  // number of successful zip's
        iptr_type &_ndiv ,  // number of successful div's
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*--------------------- zip//div topo. cell hierarchy */
        _zip_fcel( _geom, _mesh, _hfun,
            _kern, _hval, 
            _nset, _mark, _imrk, 
            _opts, _QLIM, _DLIM,
            _nzip, _ndiv, _tcpu) ;

        _zip_ecel( _geom, _mesh, _hfun,
            _kern, _hval, 
            _nset, _mark, _imrk, 
            _opts, _QLIM, _DLIM,
            _nzip, _ndiv, _tcpu) ;

    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        for (auto _iter  = _mark._node.head() ;
                  _iter != _mark._node.tend() ;
                ++_iter  )
        {
    /*--------------------- undo local inc. on node flags */
            if (*_iter > + _imrk)
                *_iter = + _imrk;
            else
            if (*_iter < - _imrk)
                *_iter = - _imrk;
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._init_zips += _tcpu.time_span(_ttic, _ttoc) ;
    #   endif//__use_timers

    }



