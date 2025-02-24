
    /*
    --------------------------------------------------------
     * RDEL-REFINE-2: refine restricted subfaces in R^2.
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
     * Last updated: 15 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_2.hpp


    /*
    --------------------------------------------------------
     * ETOP-DISK: return TRUE if 1-manifold disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type etop_disk (
        mesh_type &_mesh,
        iptr_type  _npos,
        iptr_list &_tset,
        edat_list &_eset
        )
    {
        typename
            mesh_type::edge_list _edup (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(),
            .8, _mesh._eset.get_alloc()) ;

        _edup.set_slots (
            _tset.count() * +3 ,
                containers::tight_alloc) ;

        iptr_type _feat =
        _mesh._tria.node(_npos)->feat() ;
        uint_type _topo =
        _mesh._tria.node(_npos)->topo() ;

        _topo =
       (_feat != null_feat) ? _topo : 2 ;

        for (auto _tpos  = _tset.head() ;
                  _tpos != _tset.tend() ;
                ++_tpos  )
        {
    /*--------------------- assemble 1-dim. topo disk */
        for (auto _epos = +3; _epos-- != +0; )
        {
            iptr_type _enod [ +3];
            mesh_type::tria_type::tria_type::
            face_node(_enod, _epos, 2, 1);
            _enod[0] = _mesh._tria.
             tria(*_tpos)->node(_enod[0]);
            _enod[1] = _mesh._tria.
             tria(*_tpos)->node(_enod[1]);

            iptr_type _same = +0 ;
            if (_enod[0] == _npos)
                _same += +1;
            if (_enod[1] == _npos)
                _same += +1;
            if (_same != +1) continue ;

            algorithms::isort (
                &_enod[0], &_enod[2],
                    std::less<iptr_type>()) ;

            edge_data _edat;
            _edat._node[0] = _enod[0] ;
            _edat._node[1] = _enod[1] ;

            typename mesh_type::
                     edge_list::
                 item_type*_mptr=nullptr;
            typename mesh_type::
                     edge_list::
                 item_type*_eptr=nullptr;
            if(!_mesh.
                 find_edge(_edat, _mptr))
                continue ;
            if( _edup.find(_edat, _eptr))
                continue ;

           _edup.push (_edat ) ;

           _eset.push_tail(
                _mptr->_data ) ;
        }
        }

        return ( _eset.count() <= _topo ) ;
    }

    /*
    --------------------------------------------------------
     * TOPO-NODE: refinement point for bad disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind topo_node (
        geom_type &_geom,
        mesh_type &_mesh,
        edat_list &_eset,
        typename
    mesh_type::edge_list &_epro ,
        iptr_type  _npos,
        real_type *_ppos,
        real_type *_pmax,
        char_type &_tdim,
        iptr_type &_tadj,
        rdel_opts &_opts
        )
    {
        typename rdel_opts::node_kind
        _kind =  rdel_opts::null_kind ;

        real_type _rEPS = std::pow (
            std::numeric_limits
       <real_type>::epsilon (), +.75) ;

        iptr_type _enod[ +3] = {0} ;
        real_type _vnrm[ +3] ;
        _vnrm[0] = (real_type) +0. ;
        _vnrm[1] = (real_type) +0. ;
        _vnrm[2] = (real_type) +0. ;

        _pmax[2] = (real_type) +0. ;

    /*--------------------- find max. SDB in edge-set */
        for (auto _edge  = _eset.head() ;
                  _edge != _eset.tend() ;
                ++_edge  )
        {
            mesh_type::tria_type::
            tria_type::face_node(
                _enod, _edge->_eadj, 2, 1) ;
            _enod[0] =_mesh._tria.
                tria(_edge->_tadj)
                    ->node( _enod[ 0]) ;
            _enod[1] =_mesh._tria.
                tria(_edge->_tadj)
                    ->node( _enod[ 1]) ;

             auto _iptr =
            _mesh._tria.node(_enod[0]) ;
             auto _jptr =
            _mesh._tria.node(_enod[1]) ;

            algorithms::isort (
                &_enod[0], &_enod[2],
                    std::less<iptr_type>()) ;

            edge_data _edat;
            _edat._node[0] = _enod[ 0] ;
            _edat._node[1] = _enod[ 1] ;

            typename mesh_type::    // skip protected
                     edge_list::
                item_type *_same = nullptr  ;
            if(_epro.find(
                    _edat, _same)) continue ;

            char_type _feat;
            char_type _topo[ +2] ;
            iptr_type _part;
            real_type _fbal[ +3] ;
            real_type _sbal[ +3] ;
            real_type _enrm[ +2] ;
            if (!mesh_pred::edge_ball  (
                _geom, _mesh ,
                _edge->_tadj ,
                _edge->_eadj ,
                _fbal, _sbal ,
                _feat, _topo , _part ) )
            __assert( false &&
                "TOPO-NODE: edge-ball" ) ;

    /*--------------------- try to est. local normals */
            geometry::line_norm_2d (
               &_iptr->pval(0),
               &_jptr->pval(0), _enrm
                ) ;

            geometry::normalise_2d (_enrm) ;

            real_type _sign =
            geometry::dot_2d(_enrm, _vnrm) ;

            if (_sign >= (real_type) +0. )
            {
                _vnrm[0]+= _enrm[0];
                _vnrm[1]+= _enrm[1];
            }
            else
            {
                _vnrm[0]-= _enrm[0];
                _vnrm[1]-= _enrm[1];
            }

            if (_pmax[2] < _sbal[2])
            {
    /*--------------------- keep max. SDB in edge-set */
                _pmax[2] = _sbal[2];

                _kind    =
                 rdel_opts::circ_kind  ;

                _tdim    = +1;

                _tadj =  _edge->_tadj  ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
            }
        }

        if(__chkbit(
        _opts.rule(),rdel_opts::offT_kind) )
        {
    /*--------------------- calc. "skinny" off-centre */
            _vnrm[2] =
            geometry::normalise_2d (_vnrm) ;

            real_type _rrad =
                std::sqrt(_pmax[ 2]) ;

            if (_vnrm[2] >= _rEPS * _rrad)
            {
    /*--------------------- project on local "normal" */
            real_type _pvec[2] ;
            geometry::
            vector_2d(_ppos, _pmax, _pvec) ;

            real_type _sign =
            geometry::dot_2d(_pvec, _vnrm) ;

            if (_sign > (real_type)+0.)
            {
            _vnrm[0]*=
                (real_type)+2. * _rrad ;
            _vnrm[1]*=
                (real_type)+2. * _rrad ;
            }
            else
            {
            _vnrm[0]*=
                (real_type)-2. * _rrad ;
            _vnrm[1]*=
                (real_type)-2. * _rrad ;
            }

    /*--------------------- calc. local intersections */
            typename
            geom_type::line_type _ldat ;
            _ldat._ipos[0] = _ppos[0];
            _ldat._ipos[1] = _ppos[1];

            _ldat._jpos[0] = _ppos[0];
            _ldat._jpos[1] = _ppos[1];
            _ldat._jpos[0]+= _vnrm[0];
            _ldat._jpos[1]+= _vnrm[1];

            mesh::keep_all_2d <
                real_type ,
                iptr_type >  _halo   ;

            _geom.intersect( _ldat ,
                             _halo ) ;

            if (_halo._list.count() <=  +2)
            {
    /*--------------------- keep "best" intersections */
            auto _ioff = _halo._list.tend();

            real_type _dmin  =
                _rEPS * std::pow(_rrad, +2);

            real_type _best  =
                +std::numeric_limits
                    <real_type>::infinity();

            for (auto
                _iter  = _halo._list.head();
                _iter != _halo._list.tend();
              ++_iter  )
            {
                iptr_type _near = -1;
                _mesh._tria.find_node (
               &_iter->pval(0), _near , _npos) ;

                real_type _dsqr =
                    geometry::lensqr_2d (
                       &_iter->pval(0), _ppos) ;

                if (_near == _npos &&
                    _dsqr >= _dmin &&
                    _dsqr <= _best)
                {
                    _kind =
                     rdel_opts::offT_kind  ;

                    _best = _dsqr ;
                    _ioff = _iter ;
                }
            }

            if (_ioff != _halo._list.tend())
            {
                _pmax[0] = _ioff->pval (0) ;
                _pmax[1] = _ioff->pval (1) ;
            }
            }

            }
        }

        return  _kind ;
    }

    /*
    --------------------------------------------------------
     * _BAD-ETOP: refine a "bad" edge-disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind _bad_etop (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        typename
    mesh_type::edge_list &_epro ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        node_heap &_etpq ,
        iptr_list &_emrk ,
        edat_list &_etmp ,
        edat_list &_ecav ,
        escr_list &_escr ,
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
            iptr_type _npos;

        public  :
        __inline_call node_pred (
            iptr_type _nsrc
            ) : _npos(_nsrc) {}
    /*--------------------- find adj. set of tria-to-node */
        __inline_call bool_type operator() (
            tria_type&_tria,
            iptr_type _tpos,
            iptr_type _fpos
            )
        {
            iptr_type _tnod[3] = {
            _tria.tria(_tpos)->node(0) ,
            _tria.tria(_tpos)->node(1) ,
            _tria.tria(_tpos)->node(2)
                } ;

            __unreferenced(_fpos) ;

            if (_tnod[0] == this->_npos)
            return  true ;
            if (_tnod[1] == this->_npos)
            return  true ;
            if (_tnod[2] == this->_npos)
            return  true ;

            return false ;
        }
        } ;

        typename rdel_opts::node_kind
        _kind =  rdel_opts::null_kind ;

    /*--------------------- pop() leading element from PQ */
        _nnew.set_count(  + 0) ;
        _nold.set_count(  + 0) ;
        _tnew.set_count(  + 0) ;
        _told.set_count(  + 0) ;

        _etmp.set_count(  + 0) ;
        _escr.set_count(  + 0) ;
        _ecav.set_count(  + 0) ;
        _tscr.set_count(  + 0) ;
        _tcav.set_count(  + 0) ;

        _bscr.set_count(  + 0) ;
        _bcav.set_count(  + 0) ;

        containers::array<iptr_type> _tset;
        containers::array<edge_data> _eset;

        _tset.set_alloc(  +32) ;
        _eset.set_alloc(  +32) ;

    /*--------------------- cycles until a bad disk found */
        bool_type _find = false;
        iptr_type _npos = -1 ;
        for ( ; !_etpq.empty() ; )
        {
            node_data _ndat;
            _etpq._pop_root(_ndat) ;

            _npos = _ndat._node[0] ;

            _emrk[_npos] =   -1;

            _tset.set_count( +0) ;
            _eset.set_count( +0) ;

            _mesh._tria.walk_node(
            _npos, node_pred(_npos), _tset) ;

            if (!etop_disk(
                _mesh, _npos, _tset, _eset) )
            {
                _find = true; break ;
            }
        }
        if (!_find) return _kind ;

    /*------------------------- calc. best node in cavity */
        real_type _pmax [ +3] ;
        char_type _tmax = -1;
        iptr_type _hint = -1;

        auto _ppos =
            &_mesh .
        _tria.node (_npos)->pval(0) ;

        _kind = topo_node( _geom,
            _mesh , _eset, _epro,
            _npos , _ppos, _pmax,
            _tmax , _hint, _args) ;

    /*------------------------- push node via constraints */
        _kind = push_node( _geom,
            _hfun , _mesh, _mode,
            _tmax , _pmax, _kind,
            _epro , _nnew, _nold,
            _tnew , _told, _etmp,
            _ecav , _escr,
            _tcav , _tscr,
            _bcav , _bscr, _hint,
            _tdim , _pass, _args) ;
         
        return ( _kind ) ;
    }



