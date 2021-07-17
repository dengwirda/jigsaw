
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
     * Last updated: 15 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_3.hpp


    /*
    --------------------------------------------------------
     * FTOP-SCAN: return no. connected components.
    --------------------------------------------------------
     */

    __static_call
    __normal_call iptr_type ftop_scan (
        mesh_type &_mesh,
        fdat_list &_fset,
        typename
    mesh_type::edge_list  &_edup
        )
    {
        containers::array <
        typename mesh_type::
            edge_list::item_type*> _aset;

        iptr_list _seen, _fbfs;
        _seen.set_count(_fset.count(),
            containers::loose_alloc, 0) ;

        iptr_type _fcon  = +0 ;
        iptr_type _fnum  = +0 ;
        iptr_type _fpos  ;

        for (auto _face  = _fset.head() ;
                  _face != _fset.tend() ;
                ++_face, ++_fnum)
        {
            if (_seen [_fnum] == +0)
            {

            _seen[_fnum] = +1 ;

            _fbfs.push_tail(_fnum) ;

            _fcon       += +1 ;

            for ( ; !_fbfs.empty() ; )
            {

            _fbfs._pop_tail(_fpos) ;

            for(auto _epos = +3; _epos-- != +0; )
            {
                iptr_type _tadj =
                    _fset[_fpos]._tadj ;
                iptr_type _fadj =
                    _fset[_fpos]._fadj ;

                iptr_type _fnod [ +4];
                mesh_type::tria_type::tria_type::
                face_node(_fnod, _fadj, 3, 2);
                _fnod[0] = _mesh._tria.
                 tria( _tadj)->node(_fnod[0]);
                _fnod[1] = _mesh._tria.
                 tria( _tadj)->node(_fnod[1]);
                _fnod[2] = _mesh._tria.
                 tria( _tadj)->node(_fnod[2]);

                iptr_type _enod [ +3] ;
                mesh_type::tria_type::tria_type::
                face_node(_enod, _epos, 2, 1);
                _enod[0] = _fnod[_enod[ 0] ] ;
                _enod[1] = _fnod[_enod[ 1] ] ;

                algorithms::isort (
                    &_enod[0], &_enod[2],
                        std::less<iptr_type>()) ;

                edge_data _edat;
                _edat._node[0] = _enod[0] ;
                _edat._node[1] = _enod[1] ;

                _aset.set_count(+0,
                    containers::loose_alloc);

                if(!_edup.find(_edat, _aset))
                    continue ;

                for (auto _iter  = _aset.head() ;
                          _iter != _aset.tend() ;
                        ++_iter  )
                {
                    iptr_type _iadj =
                      (*_iter)->_data._tadj;

                    if (_seen[_iadj] == +0)
                    {
                        _seen[_iadj]  = +1 ;

                        _fbfs.
                         push_tail( _iadj );
                    }
                }
            }

            }

            }
        }

        return ( _fcon ) ;
    }

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
           +.8, _mesh._eset.get_alloc()) ;

        _edup.set_slots (
            _tset.count() * +6 ,
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
        for (auto _epos = +6; _epos-- != +0; )
        {
            iptr_type _enod [ +4];
            mesh_type::tria_type::tria_type::
            face_node(_enod, _epos, 3, 1);
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
            else
            {   _edup.push(_edat) ;
            }

           _eset.push_tail(_mptr->_data ) ;
        }
        }

        return ( _eset.count() <= _topo ) ;
    }

    /*
    --------------------------------------------------------
     * FTOP-DISK: return TRUE if 2-manifold disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type ftop_disk (
        mesh_type &_mesh,
        iptr_type  _npos,
        iptr_list &_tset,
        edat_list &_eset,
        fdat_list &_fset
        )
    {
        typename
            mesh_type::edge_list _edup (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(),
           +.8, _mesh._eset.get_alloc()) ;
        typename
            mesh_type::face_list _fdup (
        typename mesh_type::face_hash(),
        typename mesh_type::face_pred(),
           +.8, _mesh._fset.get_alloc()) ;

        containers::array <
        typename mesh_type::
            edge_list::item_type*> _aset ;

        _edup.set_slots (
            _tset.count() * +6 ,
                containers::tight_alloc) ;

        _fdup.set_slots (
            _tset.count() * +4 ,
                containers::tight_alloc) ;

        for (auto _tpos  = _tset.head() ;
                  _tpos != _tset.tend() ;
                ++_tpos  )
        {
    /*--------------------- assemble 2-dim. topo disk */
        for (auto _fpos = +4; _fpos-- != +0; )
        {
            iptr_type _fnod [ +4];
            mesh_type::tria_type::tria_type::
            face_node(_fnod, _fpos, 3, 2);
            _fnod[0] = _mesh._tria.
             tria(*_tpos)->node(_fnod[0]);
            _fnod[1] = _mesh._tria.
             tria(*_tpos)->node(_fnod[1]);
            _fnod[2] = _mesh._tria.
             tria(*_tpos)->node(_fnod[2]);

            iptr_type _same = +0 ;
            if (_fnod[0] == _npos)
                _same += +1;
            if (_fnod[1] == _npos)
                _same += +1;
            if (_fnod[2] == _npos)
                _same += +1;
            if (_same != +1) continue ;

            algorithms::isort (
                &_fnod[0], &_fnod[3],
                    std::less<iptr_type>()) ;

            face_data _fdat;
            _fdat._node[0] = _fnod[0] ;
            _fdat._node[1] = _fnod[1] ;
            _fdat._node[2] = _fnod[2] ;

            typename mesh_type::
                     face_list::
                 item_type*_mptr=nullptr;
            typename mesh_type::
                     face_list::
                 item_type*_fptr=nullptr;
            if(!_mesh.
                 find_face(_fdat, _mptr))
                continue ;
            if( _fdup.find(_fdat, _fptr))
                continue ;
            else
            {   _fdup.push(_fdat) ;
            }

           _fset.push_tail(_mptr->_data ) ;
        }
        }

        iptr_type _fpos  = +0;
        for (auto _face  = _fset.head() ;
                  _face != _fset.tend() ;
                ++_face, ++_fpos)
        {
    /*--------------------- assemble 1-dim. topo disk */
        for (auto _epos = +3; _epos-- != +0; )
        {
            iptr_type _tadj = _face->_tadj ;
            iptr_type _fadj = _face->_fadj ;

            iptr_type _fnod [ +4];
            mesh_type::tria_type::tria_type::
            face_node(_fnod, _fadj, 3, 2);
            _fnod[0] = _mesh._tria.
             tria( _tadj)->node(_fnod[0]);
            _fnod[1] = _mesh._tria.
             tria( _tadj)->node(_fnod[1]);
            _fnod[2] = _mesh._tria.
             tria( _tadj)->node(_fnod[2]);

            iptr_type _enod [ +3] ;
            mesh_type::tria_type::tria_type::
            face_node(_enod, _epos, 2, 1);
            _enod[0] = _fnod[_enod[ 0] ] ;
            _enod[1] = _fnod[_enod[ 1] ] ;

            iptr_type _same = +0;
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

            _edat._tadj =
                (iptr_type) _fpos ; //!! no ??
            _edat._eadj =
                (char_type) _epos ;

            typename mesh_type::
                     edge_list::
                 item_type*_eptr=nullptr;
            if(!_edup.find(_edat, _eptr))
            {
                _edup.push(_edat) ;
                _eset.
                 push_tail(_edat) ;
            }
            else
            {
                _edup.push(_edat) ;
            }
        }
        }

        for (auto _edge  = _eset.head() ;
                  _edge != _eset.tend() ;
                ++_edge  )
        {
            _aset.set_count ( +0,
                containers::loose_alloc);

            if(_edup.find(*_edge, _aset))
            {
        /*------------------- "internal" topo. degree */
            uint_type _topo = +2;

            typename mesh_type::
                     edge_list::
                 item_type*_eptr=nullptr;
            if(_mesh.
                find_edge(*_edge, _eptr))
            {
        /*------------------- "external" topo. degree */
                _topo =
                    _eptr->_data. _topo ;
            }

            if (_aset.count () != _topo )
                return false ;
            }
        }

    /*--------------------- count no. connected comp. */
        return ftop_scan(_mesh, _fset, _edup) == 1;
    }

    /*
    --------------------------------------------------------
     * TOPO-NODE: refinement point for bad disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind etop_node (
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

        iptr_type _enod[ +4] = {0} ;
        real_type _vdir[ +4] ;
        _vdir[0] = (real_type) +0. ;
        _vdir[1] = (real_type) +0. ;
        _vdir[2] = (real_type) +0. ;
        _vdir[3] = (real_type) +0. ;

        _pmax[3] = (real_type) +0. ;

    /*--------------------- find max. SDB in edge-set */
        for (auto _edge  = _eset.head() ;
                  _edge != _eset.tend() ;
                ++_edge  )
        {
            mesh_type::tria_type::
            tria_type::face_node(
                _enod, _edge->_eadj, 3, 1) ;
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
                &_enod[0], &_enod[1],
                    std::less<iptr_type>()) ;

            edge_data _edat;
            _edat._node[0] = _enod[ 0] ;
            _edat._node[1] = _enod[ 1] ;

            typename mesh_type::    // skip protected
                     edge_list::
                item_type *_same = nullptr  ;
            if(_epro.find(
                    _edat, _same)) continue ;

            char_type _feat, _topo ;
            iptr_type _part;
            real_type _fbal[ +4] ;
            real_type _sbal[ +4] ;
            real_type _edir[ +3] ;
            if (!mesh_pred::edge_ball  (
                _geom, _mesh ,
                _edge->_tadj ,
                _edge->_eadj ,
                _fbal, _sbal ,
                _feat, _topo ,  _part) )
            __assert( false &&
                "TOPO-NODE: edge-ball" ) ;

    /*--------------------- try to est. local normals */
            geometry::vector_3d (
               &_iptr->pval(0),
               &_jptr->pval(0), _edir
                ) ;

            geometry::normalise_3d (_edir) ;

            real_type _sign =
            geometry::dot_3d(_edir, _vdir) ;

            if (_sign >= (real_type) +0. )
            {
                _vdir[0]+= _edir[0];
                _vdir[1]+= _edir[1];
                _vdir[2]+= _edir[2];
            }
            else
            {
                _vdir[0]-= _edir[0];
                _vdir[1]-= _edir[1];
                _vdir[2]-= _edir[2];
            }

            if (_pmax[3] < _sbal[3])
            {
    /*--------------------- keep max. SDB in edge-set */
                _pmax[3] = _sbal[3];

                _kind    =
                 rdel_opts::circ_kind  ;

                _tdim    = +1;

                _tadj =  _edge->_tadj  ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
                _pmax[2] = _sbal[2];
            }
        }

        if(__chkbit(
        _opts.rule(),rdel_opts::offT_kind) )
        {
    /*--------------------- calc. "skinny" off-centre */
            _vdir[3] =
            geometry::normalise_3d (_vdir) ;

            real_type _rrad =
                std::sqrt(_pmax[ 3]);

            if (_vdir[3] >= _rEPS * _rrad)
            {
    /*--------------------- project on local "normal" */
            typename
            geom_type::flat_type _fdat ;
            _fdat._ppos[0] = _ppos[0];
            _fdat._ppos[1] = _ppos[1];
            _fdat._ppos[2] = _ppos[2];

            _fdat._nvec[0] = _vdir[0];
            _fdat._nvec[1] = _vdir[1];
            _fdat._nvec[2] = _vdir[2];

            real_type _RMIN[3] = {
            _ppos[0]-(real_type) +2.*_rrad ,
            _ppos[1]-(real_type) +2.*_rrad ,
            _ppos[2]-(real_type) +2.*_rrad ,
                } ;
            _fdat._bnds.push_tail(_RMIN) ;

            real_type _RMAX[3] = {
            _ppos[0]+(real_type) +2.*_rrad ,
            _ppos[1]+(real_type) +2.*_rrad ,
            _ppos[2]+(real_type) +2.*_rrad ,
                } ;
            _fdat._bnds.push_tail(_RMAX) ;

    /*--------------------- calc. local intersections */
            mesh::keep_all_3d <
                real_type ,
                iptr_type >  _halo   ;

            _geom.intersect( _fdat ,
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
                    geometry::lensqr_3d (
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
                _pmax[2] = _ioff->pval (2) ;
            }
            }

            }
        }

        return _kind ;
    }

    __static_call
    __normal_call
    typename rdel_opts::node_kind ftop_node (
        geom_type &_geom,
        mesh_type &_mesh,
        fdat_list &_fset,
        typename
    mesh_type::face_list &_fpro ,
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

        iptr_type _fnod[ +4] = {0} ;
        real_type _vnrm[ +4] ;
        _vnrm[0] = (real_type) +0. ;
        _vnrm[1] = (real_type) +0. ;
        _vnrm[2] = (real_type) +0. ;
        _vnrm[3] = (real_type) +0. ;

        _pmax[3] = (real_type) +0. ;

    /*--------------------- find max. SDB in face-set */
        for (auto _face  = _fset.head() ;
                  _face != _fset.tend() ;
                ++_face  )
        {
            mesh_type::tria_type::
            tria_type::face_node(
                _fnod, _face->_fadj, 3, 2) ;
            _fnod[0] =_mesh._tria.
                tria(_face->_tadj)
                    ->node( _fnod[ 0]) ;
            _fnod[1] =_mesh._tria.
                tria(_face->_tadj)
                    ->node( _fnod[ 1]) ;
            _fnod[2] =_mesh._tria.
                tria(_face->_tadj)
                    ->node( _fnod[ 2]) ;

             auto _iptr =
            _mesh._tria.node(_fnod[0]) ;
             auto _jptr =
            _mesh._tria.node(_fnod[1]) ;
             auto _kptr =
            _mesh._tria.node(_fnod[2]) ;

            algorithms::isort (
                &_fnod[0], &_fnod[2],
                    std::less<iptr_type>()) ;

            face_data _fdat;
            _fdat._node[0] = _fnod[ 0] ;
            _fdat._node[1] = _fnod[ 1] ;
            _fdat._node[2] = _fnod[ 2] ;

            typename mesh_type::    // skip protected
                     face_list::
                item_type *_same = nullptr  ;
            if(_fpro.find(
                    _fdat, _same)) continue ;

            char_type _feat, _topo ;
            iptr_type _part;
            real_type _fbal[ +4] ;
            real_type _sbal[ +4] ;
            real_type _fnrm[ +3] ;
            if (!mesh_pred::face_ball  (
                _geom, _mesh ,
                _face->_tadj ,
                _face->_fadj ,
                _fbal, _sbal ,
                _feat, _topo ,  _part) )
            __assert( false &&
                "TOPO-NODE: face-ball" ) ;

    /*--------------------- try to est. local normals */
            geometry::tria_norm_3d (
               &_iptr->pval(0),
               &_jptr->pval(0),
               &_kptr->pval(0), _fnrm) ;

            geometry::normalise_3d (_fnrm) ;

            real_type _sign =
            geometry::dot_3d(_fnrm, _vnrm) ;

            if (_sign >= (real_type) +0. )
            {
                _vnrm[0]+= _fnrm[0];
                _vnrm[1]+= _fnrm[1];
                _vnrm[2]+= _fnrm[2];
            }
            else
            {
                _vnrm[0]-= _fnrm[0];
                _vnrm[1]-= _fnrm[1];
                _vnrm[2]-= _fnrm[2];
            }

            if (_pmax[3] < _sbal[3])
            {
    /*--------------------- keep max. SDB in edge-set */
                _pmax[3] = _sbal[3];

                _kind    =
                 rdel_opts::circ_kind  ;

                _tdim    = +2;

                _tadj =  _face->_tadj  ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
                _pmax[2] = _sbal[2];
            }
        }

        if(__chkbit(
        _opts.rule(),rdel_opts::offT_kind) )
        {
    /*--------------------- calc. "skinny" off-centre */
            _vnrm[3] =
            geometry::normalise_3d (_vnrm) ;

            real_type _rrad =
                std::sqrt(_pmax[ 3]);

            if (_vnrm[3] >= _rEPS * _rrad)
            {
    /*--------------------- project on local "normal" */
            real_type _pvec[3] ;
            geometry::
            vector_3d(_ppos, _pmax, _pvec) ;

            real_type _sign =
            geometry::dot_3d(_pvec, _vnrm) ;

            if (_sign > (real_type)+0.)
            {
            _vnrm[0]*=
                (real_type)+2. * _rrad ;
            _vnrm[1]*=
                (real_type)+2. * _rrad ;
            _vnrm[2]*=
                (real_type)+2. * _rrad ;
            }
            else
            {
            _vnrm[0]*=
                (real_type)-2. * _rrad ;
            _vnrm[1]*=
                (real_type)-2. * _rrad ;
            _vnrm[2]*=
                (real_type)-2. * _rrad ;
            }

    /*--------------------- calc. local intersections */
            typename
            geom_type::line_type _ldat ;
            _ldat._ipos[0] = _ppos[0];
            _ldat._ipos[1] = _ppos[1];
            _ldat._ipos[2] = _ppos[2];

            _ldat._jpos[0] = _ppos[0];
            _ldat._jpos[1] = _ppos[1];
            _ldat._jpos[2] = _ppos[2];
            _ldat._jpos[0]+= _vnrm[0];
            _ldat._jpos[1]+= _vnrm[1];
            _ldat._jpos[2]+= _vnrm[2];

            mesh::keep_all_3d <
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
                    geometry::lensqr_3d (
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
                _pmax[2] = _ioff->pval (2) ;
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
        typename
    mesh_type::face_list &_fpro ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        node_heap &_etpq ,
        iptr_list &_emrk ,
        edat_list &_etmp ,
        edat_list &_ecav ,
        escr_list &_escr ,
        fdat_list &_ftmp ,
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
            iptr_type _tnod[4] = {
            _tria.tria(_tpos)->node(0) ,
            _tria.tria(_tpos)->node(1) ,
            _tria.tria(_tpos)->node(2) ,
            _tria.tria(_tpos)->node(3)
                } ;

            __unreferenced(_fpos) ;

            if (_tnod[0] == this->_npos)
            return  true ;
            if (_tnod[1] == this->_npos)
            return  true ;
            if (_tnod[2] == this->_npos)
            return  true ;
            if (_tnod[3] == this->_npos)
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
        _ftmp.set_count(  + 0) ;
        _fscr.set_count(  + 0) ;
        _fcav.set_count(  + 0) ;
        _tscr.set_count(  + 0) ;
        _tcav.set_count(  + 0) ;

        _bscr.set_count(  + 0) ;
        _bcav.set_count(  + 0) ;

        containers::array<iptr_type> _tset;
        containers::array<edge_data> _eset;
        containers::array<face_data> _fset;

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
        real_type _pmax [ +4] ;
        char_type _tmax = -1;
        iptr_type _hint = -1;

        auto _ppos =
            &_mesh .
        _tria.node (_npos)->pval(0) ;

        _kind = etop_node( _geom,
            _mesh , _eset, _epro,
            _npos , _ppos, _pmax,
            _tmax , _hint, _args) ;

    /*------------------------- push node via constraints */
        _kind = push_node( _geom,
            _hfun , _mesh, _mode,
            _tmax , _pmax, _kind,
            _epro , _fpro,
            _nnew , _nold,
            _tnew , _told,
            _etmp , _ecav, _escr,
            _ftmp , _fcav, _fscr,
            _tcav , _tscr,
            _bcav , _bscr, _hint,
            _tdim , _pass, _args) ;

        return ( _kind ) ;
    }

    /*
    --------------------------------------------------------
     * _BAD-FTOP: refine a "bad" face-disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind _bad_ftop (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        typename
    mesh_type::edge_list &_epro ,
        typename
    mesh_type::face_list &_fpro ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        node_heap &_ftpq ,
        iptr_list &_fmrk ,
        edat_list &_etmp ,
        edat_list &_ecav ,
        escr_list &_escr ,
        fdat_list &_ftmp ,
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
            iptr_type _tnod[4] = {
            _tria.tria(_tpos)->node(0) ,
            _tria.tria(_tpos)->node(1) ,
            _tria.tria(_tpos)->node(2) ,
            _tria.tria(_tpos)->node(3)
                } ;

            __unreferenced(_fpos) ;

            if (_tnod[0] == this->_npos)
            return  true ;
            if (_tnod[1] == this->_npos)
            return  true ;
            if (_tnod[2] == this->_npos)
            return  true ;
            if (_tnod[3] == this->_npos)
            return  true ;

            return false ;
        }
        } ;

        typename rdel_opts::node_kind
        _kind =  rdel_opts::null_kind ;

    /*--------------------- pop() leading element from PQ */
        _tnew.set_count(  + 0) ;
        _told.set_count(  + 0) ;

        _etmp.set_count(  + 0) ;
        _escr.set_count(  + 0) ;
        _ecav.set_count(  + 0) ;
        _ftmp.set_count(  + 0) ;
        _fscr.set_count(  + 0) ;
        _fcav.set_count(  + 0) ;
        _tscr.set_count(  + 0) ;
        _tcav.set_count(  + 0) ;

        _bscr.set_count(  + 0) ;
        _bcav.set_count(  + 0) ;

        containers::array<iptr_type> _tset;
        containers::array<edge_data> _eset;
        containers::array<face_data> _fset;

        _tset.set_alloc(  +32) ;
        _eset.set_alloc(  +32) ;
        _fset.set_alloc(  +32) ;

    /*--------------------- cycles until a bad disk found */
        bool_type _find = false;
        iptr_type _npos = -1 ;
        for ( ; !_ftpq.empty() ; )
        {
            node_data _ndat;
            _ftpq._pop_root(_ndat) ;

            _npos = _ndat._node[0] ;

            _fmrk[_npos] =   -1;

            _tset.set_count( +0) ;
            _eset.set_count( +0) ;
            _fset.set_count( +0) ;

            _mesh._tria.walk_node(
            _npos, node_pred(_npos), _tset) ;

            if (!ftop_disk(_mesh,
                _npos, _tset, _eset, _fset) )
            {
                _find = true; break ;
            }
        }
        if (!_find) return _kind ;

    /*------------------------- calc. best node in cavity */
        real_type _pmax [ +4] ;
        char_type _tmax = -1;
        iptr_type _hint = -1;

        auto _ppos =
            &_mesh .
        _tria.node (_npos)->pval(0) ;

        _kind = ftop_node( _geom,
            _mesh , _fset, _fpro,
            _npos , _ppos, _pmax,
            _tmax , _hint, _args) ;

    /*------------------------- push node via constraints */
        _kind = push_node( _geom,
            _hfun , _mesh, _mode,
            _tmax , _pmax, _kind,
            _epro , _fpro,
            _nnew , _nold,
            _tnew , _told,
            _etmp , _ecav, _escr,
            _ftmp , _fcav, _fscr,
            _tcav , _tscr,
            _bcav , _bscr, _hint,
            _tdim , _pass, _args) ;

        return ( _kind ) ;
    }



