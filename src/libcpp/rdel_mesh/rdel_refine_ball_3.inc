
    /*
    --------------------------------------------------------
     * RDEL-REFINE-3: refine restricted subfaces in R^3.
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
     * Last updated: 30 April, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_3.hpp


    /*
    --------------------------------------------------------
     * _BAD-BALL: refine a "bad" ball.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind _bad_ball (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        ball_heap &_nbpq ,
        edat_list &_eold ,
        edat_list &_ecav ,
        escr_list &_escr ,
        fdat_list &_fold ,
        fdat_list &_fcav ,
        fscr_list &_fscr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_bcav ,
        ball_list &_bscr ,
        char_type &_tdim ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
        class node_pred
        {
    /*--------------------- find adj. set of tria-to-node */
        public  :

        typedef typename
            mesh_type::tria_type tria_type ;

        public  :
            real_type _dist =
                std::numeric_limits
                    <real_type>::infinity();

            iptr_type _npos ;

        public  :
        __inline_call node_pred (
            iptr_type _nsrc
            ) : _npos(_nsrc) {}
    /*--------------------- find adj. set of tria-to-node */
        __inline_call bool_type operator() (
            tria_type&_tria ,
            iptr_type _tpos ,
            iptr_type _fpos
            )
        {
        #   define CLOSEST(_pptr, _qptr)   \
            real_type _dsqr ;       \
            _dsqr = geometry::lensqr_3d (  \
               &_pptr->pval(0),     \
               &_qptr->pval(0)) ;   \
            if (this->_dist > _dsqr)    \
            {                       \
                this->_dist = _dsqr;    \
                return true ;       \
            }                       \

    /*--------------------- add if minimise dist. to node */
             auto _nptr =
            _tria .node( this->_npos ) ;

             auto _iptr = _tria.node(
            _tria .tria(_tpos)->node(0)) ;
             auto _jptr = _tria.node(
            _tria .tria(_tpos)->node(1)) ;
             auto _kptr = _tria.node(
            _tria .tria(_tpos)->node(2)) ;
             auto _lptr = _tria.node(
            _tria .tria(_tpos)->node(3)) ;

            __unreferenced(_fpos) ;

            if (_iptr != _nptr)
            {
                CLOSEST( _iptr, _nptr )
            }
            if (_jptr != _nptr)
            {
                CLOSEST( _jptr, _nptr )
            }
            if (_kptr != _nptr)
            {
                CLOSEST( _kptr, _nptr )
            }
            if (_lptr != _nptr)
            {
                CLOSEST( _lptr, _nptr )
            }

            return false;

        #   undef CLOSEST
        }
        } ;

    /*--------------------- refine next "protecting ball" */
        typename rdel_opts::node_kind
        _kind =  rdel_opts::null_kind ;

        __unreferenced( _hfun) ;
        __unreferenced( _mode) ;
        __unreferenced( _tdim) ;
        __unreferenced( _pass) ;
        __unreferenced( _args) ;

        _nnew.set_count(  +0 ) ;
        _nold.set_count(  +0 ) ;
        _tnew.set_count(  +0 ) ;
        _told.set_count(  +0 ) ;

        _eold.set_count(  +0 ) ;
        _escr.set_count(  +0 ) ;
        _ecav.set_count(  +0 ) ;
        _fold.set_count(  +0 ) ;
        _fscr.set_count(  +0 ) ;
        _fcav.set_count(  +0 ) ;
        _tscr.set_count(  +0 ) ;
        _tcav.set_count(  +0 ) ;

        _bscr.set_count(  +0 ) ;
        _bcav.set_count(  +0 ) ;

    /*--------------------- pop() leading element from PQ */
        typename mesh_type::
                 ball_list::
        item_type*_bptr = nullptr;
        ball_data _ball ;

        bool_type _find = false;
        for ( ; !_nbpq.empty() ; )
        {
    /*--------------------- cycles until valid ball found */
            _nbpq._pop_root(_ball) ;

            if (_mesh.find_ball(_ball,_bptr))
            {
            if (_bptr->_data._pass ==
                       _ball._pass )
            {
                _find = true; break;
            }
            }
        }
        if ( !_find ) return _kind ;

    /*--------------------- check dist. to existing nodes */
        real_type static constexpr _HTOL =
            (real_type) +3./8. * +3./8. ;

        real_type static constexpr _STOL =
            (real_type) +1./2. * +1./2. ;

        iptr_type _base = _ball._node[0];
        real_type _rmax = _ball._ball[3];

        iptr_list _tset ;
        _tset.set_alloc ( +32 );

        _mesh._tria.walk_node (
        _base, node_pred(_base), _tset) ;

        for (auto _tpos  = _tset.head() ;
                  _tpos != _tset.tend() ;
                ++_tpos  )
        {
            iptr_type _npos;
            for (_npos = +3; _npos-- != +0; )
            {
                iptr_type _node =
                    _mesh._tria.
                tria(*_tpos)->node(_npos) ;

                if (_node == _base) continue;

                iptr_type _feat =
                    _mesh.
                _tria.node(_node)->feat() ;

                real_type _dsqr =
                    geometry::lensqr_3d (
                   &_mesh._tria.
                     node(_base)->pval(0) ,
                   &_mesh._tria.
                     node(_node)->pval(0)
                    ) ;

                if (_feat == hard_feat)
                {
                    _rmax = std::min (
                    _rmax, _HTOL * _dsqr) ;
                }
                else
                {
                    _rmax = std::min (
                    _rmax, _STOL * _dsqr) ;
                }
            }
        }

    /*--------------------- otherwise form "collar" radii */
        _ball._ball[3] = std:: min (
        _ball._ball[3],_rmax) ;

        _ball._ball[3] = std::sqrt (
            _ball._ball [ 3]) ;

         auto _ndeg =
        _mesh._tria.node(_base)->topo() ;

    /*--------------------- find intersections with geom. */
        typename
        geom_type::ball_type _pdat, _hdat ;
        _pdat._pmid[0] = _ball._ball[0] ;
        _pdat._pmid[1] = _ball._ball[1] ;
        _pdat._pmid[2] = _ball._ball[2] ;

        _hdat._pmid[0] = _ball._ball[0] ;
        _hdat._pmid[1] = _ball._ball[1] ;
        _hdat._pmid[2] = _ball._ball[2] ;

        mesh::keep_all_3d <
            real_type, iptr_type> _pred ;
        mesh::keep_all_3d <
            real_type, iptr_type> _halo ;

        for(auto _iter = 64; _iter-- != +0 ; )
        {
            bool_type _okay = true ;

            _pred.clear() ;
            _halo.clear() ;

            _pdat._rrad =
           (real_type)+1.0 * _ball._ball[3] ;

            _hdat._rrad =
           (real_type)+1.5 * _ball._ball[3] ;

            _bptr->_data._ball[3] = // save to mesh "proper"
           (real_type)+1.0 * _ball._ball[3]
                           * _ball._ball[3] ;

            _okay = _okay &
                _geom.intersect (_pdat, _pred) ; // actual!

            if(!_okay)
            {
                _ball._ball[3] *=   // reduce rad. for geom.
                    (real_type) +.67 ;
                continue ;
            }
            else
            if (_pred._inum !=  _ndeg)
            {
                _ball._ball[3] *=   // reduce rad. for topo.
                    (real_type) +.67 ;
                continue ;
            }

            _okay = _okay &
                _geom.intersect (_hdat, _halo) ; // bracket

            if(!_okay)
            {
                _ball._ball[3] *=   // reduce rad. for geom.
                    (real_type) +.67 ;
                continue ;
            }
            else
            if (_halo._inum !=  _ndeg)
            {
                _ball._ball[3] *=  // reduce rad. for topo.
                    (real_type) +.67 ;
                continue ;
            }
            else break ;
        }

        return ( _kind ) ;
    }

    /*
    --------------------------------------------------------
     * MARK-FEAT: "feature" nodes in "collar".
    --------------------------------------------------------
     */

    template <
        typename  halo_type
             >
    __static_call
    __normal_call void_type mark_feat (
        mesh_type &_mesh ,
        iptr_type  _base ,
        halo_type &_halo ,
        iptr_list &_mark ,
        rdel_opts &_args
        )
    {
    /*------------ MARK[i] = 1 for nodes in "sharp-feat." */
        real_type _DtoR =
       (real_type)+3.141592653589793 / +180.;

        real_type _PHI1 =
       (real_type)+180. - _args.phi1() ;

        real_type _HARD =
            std::cos(_PHI1 * _DtoR) ;

        auto _BASE = _mesh._tria.node(_base);

        _mark.set_count(_halo._list.count());
        _mark.fill( +0 ) ;

        for (auto _ipos  = _halo._list.head() ;
                  _ipos != _halo._list.tend() ;
                ++_ipos  )
        for (auto _jpos  = _halo._list.head() ;
                  _jpos != _halo._list.tend() ;
                ++_jpos  )
        {
    /*--------------------- calc. angle from IPOS to JPOS */
            if (_ipos == _jpos) continue ;

             auto _inum  =
            _ipos - _halo._list.head();
             auto _jnum  =
            _jpos - _halo._list.head();

            real_type _ivec[3], _jvec[3] ;
            geometry::vector_3d (
               &_BASE->pval(0),
               &_ipos->pval(0), _ivec);

            geometry::vector_3d (
               &_BASE->pval(0),
               &_jpos->pval(0), _jvec);

            real_type _acos =
            geometry::cosine_3d(_ivec, _jvec) ;

            if (_acos >= _HARD)
            {
    /*--------------------- MARK pair if feat. is "sharp" */
                _mark[_inum] = +1 ;
                _mark[_jnum] = +1 ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * INIT-BALL: set-up ball "collar".
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_ball (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        typename
    mesh_type::edge_list &_epro ,
        typename
    mesh_type::face_list &_fpro ,
        iptr_type  _mode ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
    /*--------------------- push finalised "collar" nodes */
        __unreferenced(_hfun) ;
        __unreferenced(_mode) ;
        __unreferenced(_pass) ;
        __unreferenced(_fpro) ;

        iptr_list _mark  ;
        for (auto _bpos  =
            _mesh._bset._lptr.head () ;
                  _bpos !=
            _mesh._bset._lptr.tend () ;
                ++_bpos  )
        {
            if ( *_bpos == nullptr) continue ;

            for (auto _ball  = *_bpos ;
                _ball != nullptr;
                _ball  = _ball->_next )
            {
    /*--------------------- build "collar" for given ball */
            typename
                geom_type::ball_type _hdat ;
            _hdat._pmid[0] =
                _ball->_data._ball[0] ;
            _hdat._pmid[1] =
                _ball->_data._ball[1] ;
            _hdat._pmid[2] =
                _ball->_data._ball[2] ;

            _hdat._rrad =
            std::sqrt(_ball->_data._ball[3]) ;

            iptr_type _base =
                _ball->_data._node[0] ;

            mesh::keep_all_3d <
                real_type, iptr_type>  _halo ;

            _geom.intersect( _hdat, _halo) ;

            mark_feat(_mesh,
                _base,_halo, _mark, _args) ;

            iptr_type _inum  = +0 ;
            for (auto _iter  = _halo._list.head() ;
                      _iter != _halo._list.tend() ;
                    ++_iter, ++_inum)
            {
    /*--------------------- push nodes + constraint faces */
                if (_mark[_inum] <= +0) continue;

                iptr_type _hint =
                    _mesh._tria.node(_base)->next() ;

                real_type  _PPOS[4] ; // weight = 0
                _PPOS[0] = _iter->pval(0) ;
                _PPOS[1] = _iter->pval(1) ;
                _PPOS[2] = _iter->pval(2) ;
                _PPOS[3] = (real_type)+0. ;

                iptr_type _node = -1;
                if (_mesh._tria.push_node(
                        &_PPOS[ 0], _node, _hint) )
                {
                    _mesh._tria.node(
                        _node)->idxh() =
                     hfun_type::null_hint();

                    _mesh._tria.node(
                        _node)->fdim() = 1 ;

                    _mesh._tria.node(
                        _node)->feat() = 0 ;

                    _mesh._tria.node(
                        _node)->topo() = 2 ;

                    edge_data  _edat;
                    _edat._node[0] = _base ;
                    _edat._node[1] = _node ;

                    _epro.push(_edat)  ;
                }
            }

            }
        }
    }

    /*
    --------------------------------------------------------
     * INIT-DISC: set-up disc "collar".
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_disc (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        typename
    mesh_type::edge_list &_epro ,
        typename
    mesh_type::face_list &_fpro ,
        iptr_type  _mode ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
        real_type static constexpr
            _DTOL = (real_type) 7./4. ;

    /*--------------------- push finalised "collar" nodes */
        __unreferenced(_hfun) ;
        __unreferenced(_mode) ;
        __unreferenced(_pass) ;
        __unreferenced(_args) ;
        __unreferenced(_epro) ;

        for (auto _epos  =
            _mesh._eset._lptr.head () ;
                  _epos !=
            _mesh._eset._lptr.tend () ;
                ++_epos  )
        {
            if ( *_epos == nullptr) continue ;

            for (auto _edge  = *_epos ;
                _edge != nullptr;
                _edge  = _edge->_next )
            {
            if (_edge->_data._feat == hard_feat)
            {

            auto _inod = _edge->_data._node[0] ;
            auto _jnod = _edge->_data._node[1] ;

            auto _iptr =
                 _mesh._tria.node(_inod) ;
            auto _jptr =
                 _mesh._tria.node(_jnod) ;

            real_type _ebal[4];
            geometry::circ_ball_3d(_ebal ,
               &_iptr->pval(0),
               &_jptr->pval(0)) ;

            real_type _evec[3];
            geometry::vector_3d (
               &_iptr->pval(0),
               &_jptr->pval(0), _evec) ;


            typename
            geom_type::disc_type _hdat ;
            _hdat._pmid[0] = _ebal[0] ;
            _hdat._pmid[1] = _ebal[1] ;
            _hdat._pmid[2] = _ebal[2] ;

            _hdat._nvec[0] = _evec[0] ;
            _hdat._nvec[1] = _evec[1] ;
            _hdat._nvec[2] = _evec[2] ;

            _hdat._rrad =
                _DTOL * std::sqrt(_ebal[3]) ;

            mesh::keep_all_3d <
                real_type, iptr_type> _halo ;

            real_type _junk[3] = { +0. } ;

            _geom.intersect(_hdat, _junk, _halo) ;


            for (auto _iter  =
                      _halo._list.head() ;
                      _iter !=
                      _halo._list.tend() ;
                    ++_iter  )
            {
                iptr_type _hint =
                    _mesh._tria.node(_inod)->next() ;

                /*
                iptr_type _near = -1;
                _mesh._tria.find_node(
                   &_iter->pval(0), _near, _inod);

                if (_near != _inod &&
                    _near != _jnod )
                    continue ;
                */

                real_type  _PPOS[4] ; // weight = 0
                _PPOS[0] = _iter->pval(0) ;
                _PPOS[1] = _iter->pval(1) ;
                _PPOS[2] = _iter->pval(2) ;
                _PPOS[3] = (real_type)+0. ;

                iptr_type _node = -1;
                if (_mesh._tria.push_node(
                        &_PPOS[ 0], _node, _hint) )
                {
                    _mesh._tria.node(
                        _node)->idxh() =
                     hfun_type::null_hint();

                    _mesh._tria.node(
                        _node)->fdim() = 2 ;

                    _mesh._tria.node(
                        _node)->feat() = 0 ;

                    _mesh._tria.node(
                        _node)->topo() = 2 ;

                    face_data  _fdat;
                    _fdat._node[0] =
                    _edge ->_data._node [0];
                    _fdat._node[1] =
                    _edge ->_data._node [1];
                    _fdat._node[2] = _node ;

                    _fpro.push(_fdat)  ;
                }
            }

            }

            }
        }
    }



