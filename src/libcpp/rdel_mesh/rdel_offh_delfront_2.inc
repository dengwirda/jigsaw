
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The 
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any 
     * way whatsoever.  This code is provided "as-is" to be 
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 04 October, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
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
        _kind  = rdel_opts::null_kind ;
        
        typename hfun_type::hint_type 
        _hint  =_hfun   .null_hint () ;

        real_type static const _rtol = 
            (real_type)+1.E-03;
            
        real_type static const _epsb = 
            (real_type)+9.E-01;

        __unreferenced(_args) ;

        if (_enod[1] < _enod[0])
        std::swap(_enod[0], _enod[1]) ;

        real_type  _nsiz [ 3] ;
        _nsiz[ 0] = _hfun.eval (
       &_mesh._tria.
         node(_enod[0])->pval(0) , 
        _mesh._tria.
         node(_enod[0])->idxh()) ;
        _nsiz[ 1] = _hfun.eval (
       &_mesh._tria.
         node(_enod[1])->pval(0) , 
        _mesh._tria.
         node(_enod[1])->idxh()) ;
        
        real_type _dvec[3] ;
        real_type _ipos[2] = {
        _mesh._tria.
         node(_enod[0])->pval(0) ,
        _mesh._tria.
         node(_enod[0])->pval(1)
            } ;
        
        _hint = _mesh._tria.
            node (_enod[0])->idxh() ;

        _dvec[ 0] = 
        _sbal[ 0] - _ipos[ 0] ;
        _dvec[ 1] = 
        _sbal[ 1] - _ipos[ 1] ;
  
        _dvec[ 2] = 
         geometry::length_2d(_dvec) ;
        _dvec[ 0]/= _dvec[ 2] ;
        _dvec[ 1]/= _dvec[ 2] ;

        _ppos[ 0] = _sbal[ 0] ;
        _ppos[ 1] = _sbal[ 1] ;
        
        real_type static const  _seps = 
            (real_type)std::sqrt(
       +std::numeric_limits<real_type>::epsilon());

        uint32_t  _hash = hash_ball(_sbal) ;

        real_type _pert =
           ((real_type)(_hash % +4096)) / +4096 ;
        _pert  *= _seps ;

        for(iptr_type _iter = +0; _iter++ != +8 ; )
        {
            _nsiz[ 2] = _hfun.eval (_ppos, _hint) ;

            real_type _hval = 
            _nsiz[0] * (real_type)+1./2. +
            _nsiz[2] * (real_type)+1./2. ;

            real_type _dist = _hval ;
            if (_dist >= _epsb * _dvec[2])
            {                       // circumball limiter
                _ppos[0]=_sbal[0];
                _ppos[1]=_sbal[1];
                
                _kind = 
                 rdel_opts::circ_kind ;
                 
                 return  _kind ;
            }
            else                    // adv.-front limiter
            {
                _kind = 
                 rdel_opts::offh_kind ;
            }
            
            _dist -= _pert * _dist ;
            
            cosine_intersect _pred(_ipos,_dvec);
            if (!_geom.intersect (
                 _ipos, 
                 _dist, _pred) )
                return  rdel_opts::null_kind ;
            if (!_pred. _find  )
                return  rdel_opts::null_kind ;
            
            real_type _move = 
                geometry::length_2d(
                    _ppos,&_pred._proj.pval(0));

            _ppos[ 0] = _pred._proj.pval(0);
            _ppos[ 1] = _pred._proj.pval(1);

            if (_move < _rtol * _hval) break ;
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
        _kind  = rdel_opts::null_kind ;
        
        typename hfun_type::hint_type 
        _hint  =_hfun   .null_hint () ;
    
        real_type static const _rtol = 
            (real_type)+1.E-03;

        real_type static const _epsb = 
            (real_type)+9.E-01;
            
        real_type static const _alth = 
            (real_type)std::sqrt(3.)/2. ;

        __unreferenced(_geom) ;
        __unreferenced(_args) ;

        real_type  _nsiz [ 3] ;
        _nsiz[ 0] = _hfun.eval (
       &_mesh._tria.
         node(_fnod[0])->pval(0) , 
        _mesh._tria.
         node(_fnod[0])->idxh()) ;
        _nsiz[ 1] = _hfun.eval (
       &_mesh._tria.
         node(_fnod[1])->pval(0) , 
        _mesh._tria.
         node(_fnod[1])->idxh()) ;
        
        _hint = 
        _mesh._tria.node(_fnod[0])->idxh() ;

        _ppos[ 0] = _tbal[ 0] ;
        _ppos[ 1] = _tbal[ 1] ;
        
        real_type static const  _seps = 
            (real_type)std::sqrt(
       +std::numeric_limits<real_type>::epsilon());

        uint32_t  _hash = hash_ball(_tbal) ;

        real_type _pert =
           ((real_type)(_hash % +4096)) / +4096 ;
        _pert  *= _seps ;
 
        for(iptr_type _iter = +0; _iter++ != +8 ; )
        {
            _kind    = rdel_opts::offh_kind ;

            _nsiz[2] = _hfun.eval(_ppos, _hint) ;
       
            real_type _hmid = 
            _nsiz[0] * (real_type)1./4. +
            _nsiz[1] * (real_type)1./4. +
            _nsiz[2] * (real_type)2./4. ;
            
            real_type _dist ;
            real_type _dsqr = 
                _hmid*_hmid - _fbal[ 2] ;

            real_type _near = 
            _fbal[2] * (real_type)1./3. ;

            if (_dsqr < _near)
            {                       // min.-space limiter
                _dist = 
           +std::numeric_limits<real_type>::infinity();
            }
            else
            {
                _dist = (real_type) std::sqrt( _dsqr) ;
            }
            
            if (_dist > _alth * _hmid)  
            {                       // adv.-front limiter
                _dist = _alth * _hmid;
                _kind = 
                 rdel_opts::offh_kind;
            }
            if (_dist >= _dvec[3])  // off-centre limiter
            {
                _dist =  _dvec[3];
                _kind = 
                 rdel_opts::offc_kind;
            }
            if (_dist >= _dvec[2] * _epsb) 
            {                       // circumball limiter
                _ppos[0]=_tbal[0];
                _ppos[1]=_tbal[1];
                
                _kind = 
                 rdel_opts::circ_kind;
                 
                 return  _kind;
            }
            
            _dist -= _pert * _dist ;

            real_type _proj[2] ;
            _proj[ 0] = 
            _fbal[ 0] + _dist*_dvec[ 0];
            _proj[ 1] = 
            _fbal[ 1] + _dist*_dvec[ 1];

            real_type _move = 
            geometry::length_2d(_ppos, _proj);

            _ppos[ 0] = _proj[ 0];
            _ppos[ 1] = _proj[ 1];

            if (_move < _rtol * _hmid) break ;
        }

        return ( _kind )  ;
    }
    
    
    
