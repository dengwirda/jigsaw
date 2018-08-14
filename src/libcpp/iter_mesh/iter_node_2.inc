
    /*
    --------------------------------------------------------
     * ITER-NODE-2: optim. schemes to reposition nodes.
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
     * Last updated: 18 April, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
    
    // from iter_mesh_k.hpp
    

    /*
    --------------------------------------------------------
     * _PVT-MOVE: point-wise voro. tessellation update.
    --------------------------------------------------------
     */
   
    // Shift node to weighted mean of adj. dual vertex
    // positions... Like ODT, but sans |t_i| terms?
    
    // xnew = xold + SUM(1.0/h_i * c_i) / SUM(1.0/h_i)
   
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _pvt_move_2 (
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iptr_list &_tset ,
        node_iter  _node ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        real_type _move[_dims] = {
            (real_type) +0.0 } ;

        _ladj = 
            (real_type) +0.0 ;
    
        __unreferenced (_pred) ; // for MSVC...

        real_type _wval = 
       +std::numeric_limits<real_type>::epsilon() ;
        real_type _wsum = 
       +std::numeric_limits<real_type>::epsilon() ;
    
        iptr_type _tnum = +0 ;
    
        for (auto _tria  = _tset.head() ,
                  _tend  = _tset.tend() ;
                  _tria != _tend ; 
                ++_tria, ++_tnum )
        {     
             auto _tptr = 
            _mesh._set3.head()+*_tria ;
            
            iptr_type _tnod[3] ;
            _tnod[ 0] = _tptr->node( 0) ;
            _tnod[ 1] = _tptr->node( 1) ;
            _tnod[ 2] = _tptr->node( 2) ;
            
             auto _inod = 
            _mesh._set1.head()+_tnod[0] ;
             auto _jnod = 
            _mesh._set1.head()+_tnod[1] ;
             auto _knod = 
            _mesh._set1.head()+_tnod[2] ;
                 
            real_type _ball [_dims + 1] ;
            _pred.perp_ball (_ball , 
                &_inod->pval(0),
                    &_jnod->pval(0),
                        &_knod->pval(0), true) ;
        
            if (_hval[_tnod[0]] < (real_type)+0.)
            {
                _hval[_tnod[0]] = _hfun.eval (
                   &_inod->pval(0) , 
                    _inod->hidx()) ;
            }            
            
            if (_hval[_tnod[1]] < (real_type)+0.)
            {
                _hval[_tnod[1]] = _hfun.eval (
                   &_jnod->pval(0) , 
                    _jnod->hidx()) ;
            }
            
            if (_hval[_tnod[2]] < (real_type)+0.)
            {
                _hval[_tnod[2]] = _hfun.eval (
                   &_knod->pval(0) , 
                    _knod->hidx()) ;
            }
            
            real_type _tsqr = std::max(
                (real_type) +0., _ball[_dims]);
            
            real_type _hbar = _hval[_tnod[0]] +
                              _hval[_tnod[1]] +
                              _hval[_tnod[2]] ;
            _hbar /= (real_type) +3.0 ; 
         
            _wval  = (real_type) +1.0 / _hbar ;
            
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _move[_idim] 
                   += _wval * _ball [_idim] ;
            }
            
            _wsum  += _wval ;
            _ladj  += _tsqr ;
        }
        
        if (_tnum > +0)
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] = 
                    _move[_idim] / _wsum
                        - _node->pval(_idim) ;
            }
            
            _ladj = std::sqrt(_ladj / _tnum) ;
        }
        
    }
    
    /*
    --------------------------------------------------------
     * _ODT-MOVE: optimal delaunay tessellation update.
    --------------------------------------------------------
     */
   
    // Shift node to weighted mean of adj. dual vertex
    // positions...
    
    // xnew = xold + SUM(|t_i|_r * c_i) / SUM(|t_i|_r)
    //
    // where |t_i|_r = INT|i rho(x) dA
    //
    //               ~ |t_i| / (h_i)^2
   
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _odt_move_2 (
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iptr_list &_tset ,
        node_iter  _node ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        real_type _move[_dims] = {
            (real_type) +0.0 } ;

        _ladj = 
            (real_type) +0.0 ;
    
        __unreferenced (_pred) ; // for MSVC...

        real_type _wsum = 
       +std::numeric_limits<real_type>::epsilon() ;
    
        iptr_type _tnum = +0 ;
    
        for (auto _tria  = _tset.head() ,
                  _tend  = _tset.tend() ;
                  _tria != _tend ; 
                ++_tria, ++_tnum )
        {     
             auto _tptr = 
            _mesh._set3.head()+*_tria ;
            
            iptr_type _tnod[3] ;
            _tnod[ 0] = _tptr->node( 0) ;
            _tnod[ 1] = _tptr->node( 1) ;
            _tnod[ 2] = _tptr->node( 2) ;
            
             auto _inod = 
            _mesh._set1.head()+_tnod[0] ;
             auto _jnod = 
            _mesh._set1.head()+_tnod[1] ;
             auto _knod = 
            _mesh._set1.head()+_tnod[2] ;
            
            real_type _ball [_dims + 1] ;
            _pred.perp_ball (_ball , 
                &_inod->pval(0),
                    &_jnod->pval(0),
                        &_knod->pval(0), true) ;
                 
            real_type _tmag = 
            _pred.mass_tria (
                &_inod->pval(0),
                    &_jnod->pval(0),
                        &_knod->pval(0) ) ;
            
            if (_hval[_tnod[0]] < (real_type)+0.)
            {
                _hval[_tnod[0]] = _hfun.eval (
                   &_inod->pval(0) , 
                    _inod->hidx()) ;
            }            
            
            if (_hval[_tnod[1]] < (real_type)+0.)
            {
                _hval[_tnod[1]] = _hfun.eval (
                   &_jnod->pval(0) , 
                    _jnod->hidx()) ;
            }
            
            if (_hval[_tnod[2]] < (real_type)+0.)
            {
                _hval[_tnod[2]] = _hfun.eval (
                   &_knod->pval(0) , 
                    _knod->hidx()) ;
            }
            
            real_type _tsqr = std::max(
                (real_type) +0., _ball[_dims]);
            
            real_type _irho = (real_type)+1. 
                            / _hval[_tnod[0]] ;
            real_type _jrho = (real_type)+1. 
                            / _hval[_tnod[1]] ;
            real_type _krho = (real_type)+1. 
                            / _hval[_tnod[2]] ;
   
            _irho = _irho * _irho  ;
            _jrho = _jrho * _jrho  ;
            _krho = _krho * _krho  ;
           
            real_type _rbar = 
               (_irho+_jrho+_krho) / (real_type)+3. ;
            
            real_type _wval = _tmag* _rbar  ;
            
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _move[_idim] 
                   += _wval * _ball [_idim] ;
            }
            
            _wsum  += _wval ;
            _ladj  += _tsqr ;
        }
        
        if (_tnum > +0)
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] = 
                    _move[_idim] / _wsum
                        - _node->pval(_idim) ;
            }
            
            _ladj = std::sqrt(_ladj / _tnum) ;
        }
        
    }
    
    /*
    --------------------------------------------------------
     * GRAD-MOVE: "local-ascent" node movement vector. 
    --------------------------------------------------------
     */
    
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type grad_move_2 (
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        iptr_list &_tset ,
        node_iter  _node ,
        real_list &_cost ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        real_type static const _ZERO = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.80) ;
            
        real_type static const _HINC = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.50) ;
                
        real_type static const _RMIN = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.75) ;
        
        real_type static const _RMAX = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.25) ;
    
        real_type _qmin =
            +std::numeric_limits
                <real_type>::infinity() ;
        
        __unreferenced(_hfun);
        __unreferenced(_pred); // for MSVC...

        iptr_type _tnum = +0 ;
        
        _ladj = (real_type)0.;
        
        real_type _bmin[_dims] = {
           +std::numeric_limits
                <real_type>::infinity()
            } ;
        real_type _bmax[_dims] = {
           -std::numeric_limits
                <real_type>::infinity()
            } ;
        
        real_type _dqdx[_dims] = {
            (real_type) +0.0 } ;
        
        real_type _save  [_dims] ;
        for (auto _idim = _dims; _idim-- != +0; )
        {
            _save[_idim] = 
                _node->pval(_idim) ;
        }
            
        for (auto _tria  = _tset.head(),
                  _tend  = _tset.tend();
                  _tria != _tend;
                ++_tria, ++_tnum)
        {     
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
         
             auto _inod = _mesh.
            _set1 .head()+_tptr->node(0) ;
             auto _jnod = _mesh.
            _set1 .head()+_tptr->node(1) ;
             auto _knod = _mesh.
            _set1 .head()+_tptr->node(2) ;
         
            _qmin = std::min(
                _qmin, _cost[_tnum]) ;
         
            real_type _pmid[_dims] = {
                (real_type) +0.0 } ;
         
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _pmid[_idim] += 
                    _inod->pval(_idim) ;
                _pmid[_idim] += 
                    _jnod->pval(_idim) ;
                _pmid[_idim] += 
                    _knod->pval(_idim) ;
            }
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _pmid[_idim] 
                    /= (real_type) +3. ;
            }
            
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _bmin[_idim] = std::min(
                    _bmin[_idim], 
                    _pmid[_idim])  ;
                
                _bmax[_idim] = std::max(
                    _bmax[_idim], 
                    _pmid[_idim])  ;
            }
            
            real_type _llen = 
                std::sqrt(_pred.length_sq (
                    _pmid, &_node->pval(0))) ;
            
            _ladj +=  _llen ;
        }

        real_type _qbar , _qlow ;
        _qbar = (real_type) +1. ;
        _qlow = (real_type) +0. ;
        
        iptr_type _lnum   = +0  ;
        iptr_type _hnum   = +1  ;
        
        real_type _qlim = _qmin + 
            (real_type) +1.0E-004 ;

        _tnum = (iptr_type) +0  ;
        for (auto _tria  = _tset.head() ,
                  _tend  = _tset.tend() ;
                  _tria != _tend;
                ++_tria, ++_tnum)
        {     
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
        
            real_type _qtri = _cost[_tnum];
            
            if (_qtri <=  _qlim)
            {
                for (auto _idim = _dims ; 
                          _idim-- != 0; )
                {      
                real_type _binc = _bmax[_idim] -
                                  _bmin[_idim] ;
           
                real_type _hdel = _HINC*_binc;
                
                real_type _hsum = (real_type)0.;
                
                real_type _sdel = (real_type)0.;
                real_type _sabs = (real_type)0.;
                real_type _sbar = (real_type)0.;
                
                for (auto _iter = +0; _iter++ != +8; )
                {
                    _node->pval(_idim) = 
                       _save[_idim] + _hdel;
                        
                    _hsum  = (real_type)0. ;
                    _hsum += _node->
                        pval(_idim) - _save[_idim] ;

                    real_type _scr1 = 
                        _pred .cost_tria (
                   &_mesh._set1[
                    _tptr->node(0)].pval(0),
                   &_mesh._set1[
                    _tptr->node(1)].pval(0),
                   &_mesh._set1[
                    _tptr->node(2)].pval(0)
                        ) ;
                    
                    _node->pval(_idim) = 
                       _save[_idim] - _hdel;

                    _hsum -= _node->
                        pval(_idim) - _save[_idim] ;

                    real_type _scr0 = 
                        _pred .cost_tria (
                   &_mesh._set1[
                    _tptr->node(0)].pval(0),
                   &_mesh._set1[
                    _tptr->node(1)].pval(0),
                   &_mesh._set1[
                    _tptr->node(2)].pval(0)
                        ) ;
                
                    _sdel = _scr1 - _scr0 ;
                    _sabs = 
                        std::abs(_sdel);
                    
                    _sbar = std::max(
                        std::abs(_scr1),
                            std::abs(_scr0)) ;
                    
                    _node->pval(_idim) 
                        = _save[_idim] ;
 
                    if (_sabs > (real_type)0.)
                    {                   
                    if (_sabs > _RMAX * _sbar)
                    {
                        _hdel *= 
                            (real_type) +.10 ;
                    }    
                    else
                    if (_sabs < _RMIN * _sbar)
                    {
                        _hdel *= 
                            (real_type) +10. ;
                    }
                    else { break ; }
                    }
                    else { break ; }
                }
            
                _dqdx[_idim]+= _sdel / _hsum ;
            
                }
            
                _qlow += _qtri ; _lnum += +1 ;
            
            }
            else
            {
                _qbar += _qtri ; _hnum += +1 ; 
            }
        }
 
        if (_tnum > +0)  
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _dqdx[_idim]  /= _lnum ;
            }
            
            _qlow /=  _lnum ;
            _ladj /=  _tnum ;
        }

        if (_hnum > +0)
        {
            _qbar /=  _hnum ;
        }
        
        real_type _scal = 
            std::sqrt(_pred.length_sq(_dqdx)) ;
        
        if (_scal <= _ZERO * _ladj)
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] 
                    = (real_type) +0. ;
            }
        }
        else
        {    
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] = 
                    _dqdx[_idim] / _scal ;
            }
        
            _scal = (_qbar - _qlow) /
            _pred .length_sq(_line, _dqdx) ;
            
            _scal = std::min(_scal, _ladj) ;
            
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] *= _scal ;
            }
        }
    
    }
    
    
    
