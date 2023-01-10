
    /*
    --------------------------------------------------------
     * RDEL-PRED-DELFRONT-2: "frontal-DEL" kernel in R^2.
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
     * Last updated: 28 Dec., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_pred_delfront_2.hpp


    /*
    --------------------------------------------------------
     * EDGE-OFFH: "size"-optimal off-centre .
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind edge_offh (
        geom_type &_geom,
        hfun_type &_hfun,
        mesh_type &_mesh,
        iptr_type *_enod,
        real_type *_sbal,
        real_type *_ppos,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind
            _kind = rdel_opts::null_kind;

        typename hfun_type::hint_type
            _hint = _hfun  .null_hint() ,
            _HINT = _hfun  .null_hint() ;

        real_type static const
            _rtol = (real_type)+1.00E-003 ;

        real_type static const
            _dtol = (real_type)+8.33E-001 ;

        real_type static const
            _lgap = (real_type)+8./9. ;

        real_type static const
            _ugap = (real_type)+9./8. ;

        __unreferenced(_args) ;

        real_type  _hsiz [ 4] ;
        _hsiz[ 0] = _hfun.eval (
       &_mesh._tria.node(_enod[0])->pval(0) ,
        _mesh._tria.node(_enod[0])->idxh()) ;
        _hsiz[ 1] = _hfun.eval (
       &_mesh._tria.node(_enod[1])->pval(0) ,
        _mesh._tria.node(_enod[1])->idxh()) ;

        if (_hsiz[1] < _hsiz[0])
        {
            std::swap( _enod[0], _enod[1]);
            std::swap( _hsiz[0], _hsiz[1]);
        }

        real_type _dvec[3] ;
        real_type _base[2] = {
        _mesh._tria.node(_enod[0])->pval(0) ,
        _mesh._tria.node(_enod[0])->pval(1) ,
            } ;

        _hint = _HINT =
        _mesh._tria.node (_enod[0])->idxh() ;

        _dvec[ 0] =
        _sbal[ 0] - _base[ 0] ;
        _dvec[ 1] =
        _sbal[ 1] - _base[ 1] ;

        _dvec[ 2] =
        geometry:: normalise_2d (_dvec);

        _ppos[ 0] = _sbal[ 0] ;
        _ppos[ 1] = _sbal[ 1] ;

        typename
            geom_type::ball_type _ball ;
        _ball.
        _pmid[ 0] = _base[ 0] ;
        _ball.
        _pmid[ 1] = _base[ 1] ;

        real_type static const  _seps =
            (real_type)std::sqrt(
       +std::numeric_limits<real_type>::epsilon());

        uint32_t  _hash = hash_ball(_sbal) ;

        real_type _pert =
           ((real_type)(_hash % +4096)) / +4096 ;
        _pert  +=   + 1 ;
        _pert  *= _seps ;

        for (auto _iter = +0; _iter++ != +4 ; )
        {
            _hsiz[2] = _hfun.eval(_ppos, _hint) ;

            real_type _PMID[2] = {
           (real_type)+.5 * (_base[0]+_ppos[0]) ,
           (real_type)+.5 * (_base[1]+_ppos[1]) ,
                } ;

            _hsiz[3] = _hfun.eval(_PMID, _HINT) ;

            real_type _hval =
            _hsiz[0] * (real_type)1./4. +
            _hsiz[2] * (real_type)1./4. +
            _hsiz[3] * (real_type)2./4. ;

            real_type _dist = _hval ;

            if (_dist >= _dvec[2] * _lgap &&
                _dist <= _dvec[2] * _ugap )
            {                       // circumball cushion
                _dist  =
                    .5 *(_dist + _dvec[2] ) ;

                _kind  =  rdel_opts::offH_kind;
            }

            if (_dist >= _dvec[2] * _ugap )
            {                       // circumball limiter
                _ppos[0]=_sbal[0] ;
                _ppos[1]=_sbal[1] ;

                _kind  =  rdel_opts::circ_kind;
                 return  _kind ;
            }
            else                    // adv.-front limiter
            {
                _kind  =  rdel_opts::offH_kind;
            }

            _dist -= _pert * _dist ;

            iptr_type _dual = -1 ;
            cosine_intersect _pred(_base,_dvec) ;
            _ball._rrad  = _dist ;
            if (!_geom.intersect (
                 _ball, _pred) )
                return rdel_opts::null_kind;
            if (!_pred. _find  )
                return rdel_opts::null_kind;

            real_type _move =
                geometry::length_2d(
                    _ppos,&_pred._proj.pval(0)) ;

            _ppos[ 0] = _pred._proj.pval(0);
            _ppos[ 1] = _pred._proj.pval(1);

            _mesh._tria.find_node ( // test voro. limiter
                &_ppos[0], _dual, _enod[ 0]) ;

            if (_dual != _enod [ +0] )
            {
            real_type _dlen =
                geometry::length_2d(
                    _ppos, &_mesh._tria.
                       node(_dual)->pval(0)) ;

            if (_dlen <= _dtol*_hval )
            {
                return rdel_opts::null_kind;
            }
            }

            if (_move <= _rtol*_hval ) break ;
        }

        return ( _kind ) ;
    }

    /*
    --------------------------------------------------------
     * TRIA-OFFH: "size"-optimal off-centre .
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind tria_offh (
        geom_type &_geom,
        hfun_type &_hfun,
        mesh_type &_mesh,
        iptr_type *_fnod,
        real_type *_fbal,
        real_type *_tbal,
        real_type *_dvec,
        real_type *_ppos,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind
            _kind = rdel_opts::null_kind;

        typename hfun_type::hint_type
            _hint = _hfun  .null_hint() ,
            _HINT = _hfun  .null_hint() ;

        real_type static const
            _rtol = (real_type)+1.00E-003 ;

        real_type static const _alth =
            (real_type)std::sqrt(3.)/2. ;

        real_type static const
            _dtol = (real_type)+8.33E-001 ;

        real_type static const
            _lgap = (real_type)+8./9. ;

        real_type static const
            _ugap = (real_type)+9./8. ;

        __unreferenced(_geom) ;
        __unreferenced(_args) ;

        real_type  _hsiz [ 4] ;
        _hsiz[ 0] = _hfun.eval (
       &_mesh._tria.node(_fnod[0])->pval(0) ,
        _mesh._tria.node(_fnod[0])->idxh()) ;
        _hsiz[ 1] = _hfun.eval (
       &_mesh._tria.node(_fnod[1])->pval(0) ,
        _mesh._tria.node(_fnod[1])->idxh()) ;

        _hint = _HINT =
        _mesh._tria.node(_fnod[0])->idxh () ;

        real_type  _bias[2] = {
        _hsiz[ 1] / (_hsiz[0] + _hsiz[1]) ,
        _hsiz[ 0] / (_hsiz[0] + _hsiz[1]) ,
            } ;

        real_type  _base[2] ;
        _base[0] = _bias[0] *
        _mesh._tria.node(_fnod[0])->pval(0) +
                   _bias[1] *
        _mesh._tria.node(_fnod[1])->pval(0) ;
        _base[1] = _bias[0] *
        _mesh._tria.node(_fnod[0])->pval(1) +
                   _bias[1] *
        _mesh._tria.node(_fnod[1])->pval(1) ;

        _ppos[ 0] = _tbal[ 0] ;
        _ppos[ 1] = _tbal[ 1] ;

        real_type static const  _seps =
            (real_type)std::sqrt(
       +std::numeric_limits<real_type>::epsilon());

        uint32_t  _hash = hash_ball(_tbal) ;

        real_type _pert =
           ((real_type)(_hash % +4096)) / +4096 ;
        _pert  +=   + 1 ;
        _pert  *= _seps ;

        for (auto _iter = +0; _iter++ != +4 ; )
        {
            _hsiz[2] = _hfun.eval(_ppos, _hint) ;

            real_type _PMID[2] = {
           (real_type)+.5 * (_base[0]+_ppos[0]) ,
           (real_type)+.5 * (_base[1]+_ppos[1]) ,
                } ;

            _hsiz[3] = _hfun.eval(_PMID, _HINT) ;

            real_type _hmid =
            _hsiz[0] * (real_type)2./9. +
            _hsiz[1] * (real_type)2./9. +
            _hsiz[2] * (real_type)2./9. +
            _hsiz[3] * (real_type)3./9. ;

            _kind    = rdel_opts::offH_kind ;

            real_type _dist ;
            real_type _dsqr =
                _hmid*_hmid - _fbal[ 2] ;

            real_type _near =
            _fbal[2] * (real_type)1./3. ;

            if (_dsqr < _near)
            {                       // min.-space limiter
                _dist = std::numeric_limits
                    <real_type>::infinity() ;
            }
            else
            {
                _dist =
                (real_type)std::sqrt(_dsqr) ;
            }

            if (_dist >= _alth * _hmid)
            {                       // adv.-front limiter
                _dist  = _alth * _hmid ;
                _kind  =  rdel_opts::offH_kind;
            }

            if (_dist >= _dvec[3])  // off-centre limiter
            {
                _dist  =
                    2. * _dist * _dvec[3]
                      / (_dist + _dvec[3] ) ;

                _kind  =  rdel_opts::offC_kind;
            }

            if (_dist >= _dvec[2] * _lgap &&
                _dist <= _dvec[2] * _ugap )
            {                       // circumball cushion
                _dist  =
                    .5 *(_dist + _dvec[2] ) ;

                _kind  =  rdel_opts::offH_kind;
            }

            if (_dist >= _dvec[2] * _ugap )
            {                       // circumball limiter
                _ppos[0]=_tbal[0] ;
                _ppos[1]=_tbal[1] ;

                _kind  =  rdel_opts::circ_kind;
                 return  _kind ;
            }

            _dist -= _pert * _dist ;

            real_type _proj[2] ;
            _proj[ 0] =
            _base[ 0] + _dist * _dvec[0];
            _proj[ 1] =
            _base[ 1] + _dist * _dvec[1];

            iptr_type _dual = -1 ;
            real_type _move =
            geometry::length_2d(_ppos, _proj);

            _ppos[ 0] = _proj[ 0];
            _ppos[ 1] = _proj[ 1];

            _mesh._tria.find_node ( // test voro. limiter
                &_ppos[0], _dual, _fnod[ 0]) ;

            if (_dual != _fnod [ +0] &&
                _dual != _fnod [ +1] )
            {
            real_type _dlen =
                geometry::length_2d(
                    _ppos, &_mesh._tria.
                       node(_dual)->pval(0)) ;

            if (_dlen <= _dtol*_hmid )
            {
                return rdel_opts::null_kind;
            }
            }

            if (_move < _rtol * _hmid) break ;
        }

        return ( _kind )  ;
    }



