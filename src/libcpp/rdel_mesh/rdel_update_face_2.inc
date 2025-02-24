
    /*
    --------------------------------------------------------
     * RDEL-UPDATE-2: update restricted delaunay in R^2.
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
     * Last updated: 11 Nov., 2024
     *
     * Copyright 2013-2024
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_2.hpp


    /*
    --------------------------------------------------------
     * _POP-EDGE: delete edge from restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type _pop_edge (
        mesh_type &_mesh ,
        iptr_type  _tpos
        )
    {
        if(_mesh._eset._size == +0) return ;

        for(auto _fpos =+3; _fpos-- != +0; )
        {
    /*--------------------------------- get edge indexing */
            iptr_type _tnod[ +3] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_tnod, _fpos, +2, +1);
            _tnod[0] = _mesh._tria.
                tria(_tpos)->node(_tnod[0]);
            _tnod[1] = _mesh._tria.
                tria(_tpos)->node(_tnod[1]);

            algorithms::isort (
                &_tnod[0], &_tnod[2],
                    std::less<iptr_type>());

            edge_data _edat;
            _edat._node[0] = _tnod[0];
            _edat._node[1] = _tnod[1];

    /*--------------------------------- remove if present */
            _mesh._pop_edge(_edat);
        }
    }

    /*
    --------------------------------------------------------
     * _POP-TRIA: delete tria from restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type _pop_tria (
        mesh_type &_mesh ,
        iptr_type  _tpos
        )
    {
        if (_mesh._tset._size == 0) return;

        {
    /*--------------------------------- get tria indexing */
            tria_data _tdat;
            _tdat._node[0] = _mesh.
                _tria.tria(_tpos)->node(0);
            _tdat._node[1] = _mesh.
                _tria.tria(_tpos)->node(1);
            _tdat._node[2] = _mesh.
                _tria.tria(_tpos)->node(2);

    /*--------------------------------- remove if present */
            _mesh._pop_tria(_tdat);
        }
    }

    /*
    --------------------------------------------------------
     * INIT-BALL: add new ball to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_ball (
        mesh_type &_mesh,
        geom_type &_geom,
        hfun_type &_hfun,
        iptr_type  _npos,
        ball_list &_bset,
        ball_list &_bscr,
        char_type  _kind,
        iptr_type  _pass,
        rdel_opts &_opts
        )
    {
        __unreferenced ( _geom ) ;
        __unreferenced ( _opts ) ;

        if (_mesh._tria.
             node(_npos)->feat()==hard_feat)
        {
    /*---------- push protecting ball for "hard" features */
            ball_data _ball ;
            _ball._node[0] = _npos;

            _ball._pass    = _pass;
            _ball._kind    = _kind;

            _ball._ball[0] = _mesh.
                _tria.node(_npos)->pval(0) ;
            _ball._ball[1] = _mesh.
                _tria.node(_npos)->pval(1) ;

            real_type _rbal = _hfun.eval (
               &_mesh._tria.
            node(_npos)->pval(0) ,
                _mesh._tria.
            node(_npos)->idxh()) ;

            _ball._ball[2] = _rbal * _rbal ;

            _bset.push_tail (_ball) ;
            _bscr.push_tail (_ball) ;
        }
    }

    /*
    --------------------------------------------------------
     * PUSH-EDGE: add new edge to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_edge (
        mesh_type &_mesh ,
        geom_type &_geom ,
        hfun_type &_hfun ,
        iptr_type  _tpos ,
        edat_list &_eset ,
        escr_list &_escr ,
        typename
    mesh_type::edge_list & _edge_test ,
        iptr_type &_nedg ,
        iptr_type &_ndup ,
        iptr_type  _pass ,
        rdel_opts &_opts
        )
    {
        if(!_geom.have_feat(1))return ;

    /*-------------------------------- correct node dims? */
        iptr_type _fdim =+0;
        for (auto _node =+3; _node-- != +0; )
        {
        if (_mesh._tria.node (
            _mesh._tria.tria (
            _tpos)->node(_node))->fdim() <= +1)
            _fdim += +1 ;
        }
    /*-------------------------------- quick break if not */
        if (_fdim  < +2 ) return ;

    /*-------------------------------- check "restricted" */
        for (auto _fpos =+3; _fpos-- != +0; )
        {
        /*---------------------------- extract face nodes */
            iptr_type _tnod[ +3] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_tnod, _fpos, +2, +1) ;
            _tnod[0] = _mesh._tria.
            tria(_tpos)->node(_tnod[0]);
            _tnod[1] = _mesh._tria.
            tria(_tpos)->node(_tnod[1]);

            if (_mesh._tria.node(
                _tnod[0])->fdim() < 2 &&
                _mesh._tria.node(
                _tnod[1])->fdim() < 2 )
            {
        /*--------------- face contains correct k^d nodes */
            algorithms::isort (
                &_tnod[0], &_tnod[2],
                    std::less<iptr_type>()) ;

            edge_data _fdat;
            _fdat._node[0] = _tnod[ 0] ;
            _fdat._node[1] = _tnod[ 1] ;

            edge_cost _cdat;
            _cdat._node[0] = _tnod[ 0] ;
            _cdat._node[1] = _tnod[ 1] ;

            bool_type _test = true;
            typename mesh_type::
                     edge_list::
                item_type *_mptr = nullptr  ;
        //  pragma omp critical(update_face_edg2)
            if(_edge_test.
                find( _fdat, _mptr) )
            {
        /*--------------------------- count bnd. repeats! */
                _ndup +=
                _mptr->_data._dups;

        /*--------------------------- don't test repeats! */
                _test = false;
            }
            
            _cdat._pass    = _pass;
            _fdat._pass    = _pass;

            _fdat._tadj    = _tpos;
            _fdat._eadj    =
                (char_type)  _fpos;
            _fdat._dups    = +0; // count num. dup's
                                 // only in hash-set

            if (_test)
            {
        /*--------------------------- call face predicate */
            char_type _topo[ 2] , _feat;
            real_type _fbal[ 3] ;
            real_type _sbal[ 3] ;
            mesh_pred::edge_cost (
                    _geom, _hfun, _mesh,
                    _fdat. _tadj,
                    _fdat. _eadj,
                    _opts, _cdat,
                    _fdat. _part,
                    _feat, _topo,
                    _fdat. _kind,
                    _fbal, _sbal) ;

        /*--------------------------- push edge onto mesh */
        //  pragma omp critical(update_face_edg2)
            {
            if (_fdat._kind == mesh::ring_item)
                _escr. push_tail(_cdat) ;

            if (_fdat._kind != mesh::null_item)
                _nedg += +1 ;

            if (_fdat._kind != mesh::null_item)
                _fdat._dups = +1 ;

            if (_fdat._kind != mesh::null_item)
                _eset. push_tail(_fdat) ;

            _edge_test.push (_fdat) ;
            }
            }

            }
        } // for (auto _fpos = +3; _fpos-- != +0; )
    }

    /*
    --------------------------------------------------------
     * PUSH-TRIA: add new tria to restricted-tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_tria (
        mesh_type &_mesh ,
        geom_type &_geom ,
        hfun_type &_hfun ,
        iptr_type  _tpos ,
        iptr_type &_sign ,
        tdat_list &_tset ,
        tscr_list &_tscr ,
        /* typename
    mesh_type::tria_list & _tria_test , */
        iptr_type &_ntri ,
        iptr_type  _pass ,
        rdel_opts &_opts
        )
    {
    /*-------------------------------- check "restricted" */
        {
            if(!_geom.have_feat(1)) return ;

            iptr_type  _tnod[ +3] ;
            _tnod[0] = _mesh.
            _tria.tria(_tpos)->node(0);
            _tnod[1] = _mesh.
            _tria.tria(_tpos)->node(1);
            _tnod[2] = _mesh.
            _tria.tria(_tpos)->node(2);

            if (_mesh._tria.node(
                _tnod[0])->fdim() < 3 &&
                _mesh._tria.node(
                _tnod[1])->fdim() < 3 &&
                _mesh._tria.node(
                _tnod[2])->fdim() < 3 )
            {
        /*--------------- face contains correct k^d nodes */

            tria_data _tdat;
            _tdat._node[0] = _tnod[ 0] ;
            _tdat._node[1] = _tnod[ 1] ;
            _tdat._node[2] = _tnod[ 2] ;

            _tdat._tadj    = _tpos;

            tria_cost _cdat;
            _cdat._node[0] = _tnod[ 0] ;
            _cdat._node[1] = _tnod[ 1] ;
            _cdat._node[2] = _tnod[ 2] ;

        //  _tria_test.push( _tdat) ; guarantee no repeat

        /*--------------------------- call tria predicate */
            _tdat._part =  _sign ;

            _cdat._pass =  _pass ;
            _tdat._pass =  _pass ;

            mesh_pred::tria_cost (
                    _geom, _hfun, _mesh,
                    _tdat. _tadj,
                    _opts, _cdat,
                    _tdat. _part,
                    _tdat. _kind);

            _sign = _tdat. _part ;

        /*--------------------------- push edge onto mesh */
        //  pragma omp critical(update_face_tri3)
            {
            if (_tdat._kind == mesh::ring_item)
                _tscr. push_tail(_cdat) ;

            if (_tdat._kind != mesh::null_item)
                _ntri += +1 ;

            if (_tdat._kind != mesh::null_item)
                _tset. push_tail(_tdat) ;
            }

            }
        }
    }

    /*
    --------------------------------------------------------
     * TRIA-CIRC: calc. circumball for tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type tria_circ ( // get "exact" tbal?
        mesh_type &_mesh ,
        iptr_type  _tpos
        )
    {
        iptr_type _tnod[3] = {
            _mesh.
        _tria.tria( _tpos)->node(0) ,
            _mesh.
        _tria.tria( _tpos)->node(1) ,
            _mesh.
        _tria.tria( _tpos)->node(2)
            } ;

        algorithms::isort(
            &_tnod[0], &_tnod[3],
                std::less<iptr_type>()) ;

    /*---------------------- calc. ball in floating-point */
        real_type _tbal[3] ;
        geometry::perp_ball_2d (
            _tbal ,
       &_mesh._tria.
            node(_tnod[0])->pval(0) ,
       &_mesh._tria.
            node(_tnod[1])->pval(0) ,
       &_mesh._tria.
            node(_tnod[2])->pval(0)
            ) ;

        _mesh._tria.tria(
        _tpos)->circ(0) =  _tbal[0] ;
        _mesh._tria.tria(
        _tpos)->circ(1) =  _tbal[1] ;
    }

    /*
    --------------------------------------------------------
     * PUSH-RDEL: push faces onto restricted tria.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type push_rdel (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        bool_type  _init ,      // true to init. circball
        iptr_list &_nnew ,
        iptr_list &_tnew ,
        escr_list &_escr ,      // lists of new rDT faces
        edat_list &_edat ,
        tscr_list &_tscr ,
        tdat_list &_tdat ,
        ball_list &_bscr ,
        ball_list &_bdat ,
        iptr_type  _sign ,
        iptr_type  _pass ,
        mode_type  _dim0 ,      // lo rDT dim. to examine
        mode_type  _dim1 ,      // hi rDT dim. to examine
        rdel_opts &_opts
        )
    {
    /*------------------------- no. "restricted" subfaces */
        iptr_type _nedg = +0 ;
        iptr_type _ntri = +0 ;
        iptr_type _ndup = +0 ;

    //  iptr_type _ncpu = +1 ;  //!! sequential for now...

    //  ifdef  __use_openmp
    //  omp_set_num_threads(_ncpu);
    //  endif//__use_openmp

    /*------------------------- calc. voronoi-dual points */
        //  pragma omp parallel default(shared)
        {
        //  pragma omp for schedule(static)
            for( auto _iter  = _tnew.head();
                      _iter != _tnew.tend();
                    ++_iter  )
            {
            if (_init) tria_circ(_mesh, *_iter) ;
            }
        }

    /*------------- push any new protecting balls created */
        if (_dim0 <= node_mode &&
            _dim1 >= node_mode )
        {
            for( auto _iter  = _nnew.head();
                      _iter != _nnew.tend();
                    ++_iter  )
            {
                char_type _kind = feat_ball;

                init_ball( _mesh, _geom, _hfun,
                   *_iter, _bdat, _bscr,
                           _kind, _pass, _opts) ;
            }
        }
    /*------------- push any new restricted edges created */
        if (_dim0 <= edge_mode &&
            _dim1 >= edge_mode )
        {
            _mesh._etwo.clear();
            _mesh._etwo.set_slots(_tnew.count() * +3) ;

        //  pragma omp parallel default(shared)
            {
        //  pragma omp for schedule(static)
            for( auto _iter  = _tnew.head();
                      _iter != _tnew.tend();
                    ++_iter  )
            {
                push_edge( _mesh, _geom, _hfun,
                   *_iter, _edat, _escr,
                           _mesh. _etwo,
                    _nedg, _ndup, _pass, _opts) ;
            }
            }
        }
    /*------------- push any new restricted trias created */
        if (_dim0 <= tria_mode &&
            _dim1 >= tria_mode )
        {
        //  _mesh._ttwo.clear();
        //  _mesh._ttwo.set_slots(_tnew.count() * +1) ;

            bool_type _safe = true ;
            if (_dim0 >  +0) _safe = false ;
        //  if (_nedg >= +1) _safe = false ;
            if (_ndup >= +1) _safe = false ;

        //  pragma omp parallel default(shared)
            {
        //  pragma omp for schedule(static)
            for( auto _iter  = _tnew.head();
                      _iter != _tnew.tend();
                    ++_iter )
            {
               _sign =   (!_safe) ? -1 : _sign;

                push_tria( _mesh, _geom, _hfun,
                   *_iter, _sign, _tdat, 
                    _tscr, _ntri, _pass, _opts) ;
            }
            }
        }
    }



