
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
     * Last updated: 27 October, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
    
    // from iter_mesh_k.hpp
    

    /*
    --------------------------------------------------------
     * _ODT-MOVE: optimal delaunay tessellation update.
    --------------------------------------------------------
     */
   
    // Shift node to weighted mean of adj. dual vertex
    // positions...
    
    // xnew <=== SUM( |t_i|_r * c_i ) / SUM( |t_i|_r )
    //
    // with |t_i|_r = INT|i rho(x) dA
    //
    //              ~ |t_i| / (h_i)^3

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

            algorithms::isort(
                &_tnod[0], &_tnod[3], 
                    std::less<iptr_type>()) ;
            
             auto _inod = 
            _mesh._set1.head()+_tnod[0] ;
             auto _jnod = 
            _mesh._set1.head()+_tnod[1] ;
             auto _knod = 
            _mesh._set1.head()+_tnod[2] ;
            
            real_type _ball [_dims + 1] ;
            real_type _kpos [_dims + 1] ;
            real_type _jpos [_dims + 1] ;
            real_type _ipos [_dims + 1] ;

            _ipos[_dims] = 
                (real_type)+.75 * 
                    _inod->pval (_dims) ;
            _jpos[_dims] = 
                (real_type)+.75 * 
                    _jnod->pval (_dims) ;
            _kpos[_dims] = 
                (real_type)+.75 * 
                    _knod->pval (_dims) ;

            for (auto _idim = _dims; _idim-- != +0; )
            {
                _ipos[_idim] = 
                    _inod->pval (_idim) ;
                _jpos[_idim] = 
                    _jnod->pval (_idim) ;
                _kpos[_idim] = 
                    _knod->pval (_idim) ;
            }

            _pred.perp_ball (_ball , 
                &_ipos [+0] ,
                &_jpos [+0] ,
                &_kpos [+0] , true);
                 
            real_type _tmag =
            _pred.mass_tria (
                &_ipos [+0] ,
                &_jpos [+0] ,
                &_kpos [+0] ) ;
                
            _tmag = std::abs(_tmag);

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

            real_type _rbar =
               (_irho+_jrho+_krho) / (real_type)+3. ;
            
            real_type _wval = 
                _tmag * std::pow(_rbar, +3) ;
            
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
     * _CVT-MOVE: centroidal voronoi tessellation move.
    --------------------------------------------------------
     */
   
    // Shift vertex to weighted mean of adj. dual cell
    
    // xnew <=== INT|i rho(x) * x dA / INT|i rho(x) dA
   
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type _cvt_move_2 (
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

            algorithms::isort(
                &_tnod[0], &_tnod[3], 
                    std::less<iptr_type>()) ;
            
             auto _inod = 
            _mesh._set1.head()+_tnod[0] ;
             auto _jnod = 
            _mesh._set1.head()+_tnod[1] ;
             auto _knod = 
            _mesh._set1.head()+_tnod[2] ;
            
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
            
            real_type _irho = (real_type)+1. 
                            / _hval[_tnod[0]] ;
            real_type _jrho = (real_type)+1. 
                            / _hval[_tnod[1]] ;
            real_type _krho = (real_type)+1. 
                            / _hval[_tnod[2]] ;
       
            node_iter _anod, _bnod, _cnod;
            real_type _arho, _brho, _crho;

            if (_inod == _node)
            {
                _anod = _inod; _arho = _irho;
                _bnod = _jnod; _brho = _jrho;
                _cnod = _knod; _crho = _krho;
            }
            else
            if (_jnod == _node)
            {
                _anod = _jnod; _arho = _jrho;
                _bnod = _knod; _brho = _krho;
                _cnod = _inod; _crho = _irho;
            }
            else
            // (_knod == _node)
            {
                _anod = _knod; _arho = _krho;
                _bnod = _inod; _brho = _irho;
                _cnod = _jnod; _crho = _jrho;
            }

            real_type _cpos [_dims + 1] ;
            real_type _bpos [_dims + 1] ;
            real_type _apos [_dims + 1] ;

            _apos[_dims] = 
                (real_type)+.75 * 
                    _anod->pval (_dims) ;
            _bpos[_dims] = 
                (real_type)+.75 * 
                    _bnod->pval (_dims) ;
            _cpos[_dims] = 
                (real_type)+.75 * 
                    _cnod->pval (_dims) ;

            for (auto _idim = _dims; _idim-- != +0; )
            {
                _apos[_idim] = 
                    _anod->pval (_idim) ;
                _bpos[_idim] = 
                    _bnod->pval (_idim) ;
                _cpos[_idim] = 
                    _cnod->pval (_idim) ;
            }

            real_type _0bal [_dims + 1] ;
            _pred.perp_ball (_0bal , 
                &_apos [+0] ,
                &_bpos [+0] ,
                &_cpos [+0] , true);

            real_type _1bal [_dims + 1] ;
            real_type _2bal [_dims + 1] ;
            _pred.perp_ball (_1bal , 
                &_apos [+0] ,
                &_bpos [+0] , true);

            _pred.perp_ball (_2bal , 
                &_apos [+0] ,
                &_cpos [+0] , true);

            real_type _1mag =
            std::abs(_pred.mass_tria (
                &_apos [+0] ,
                &_1bal [+0] ,
                &_0bal [+0] ) ) ;

            real_type _2mag =
            std::abs(_pred.mass_tria (
                &_apos [+0] ,
                &_0bal [+0] ,
                &_2bal [+0] ) ) ;

            real_type _0rho = 
               (_arho+_brho+_crho) / (real_type)+3. ;

            real_type _1rho = 
               (_arho+_brho) / (real_type)+2. ;
            real_type _2rho = 
               (_arho+_crho) / (real_type)+2. ;

            _0rho = std::pow(_0rho, +3) ;
            _1rho = std::pow(_1rho, +3) ;
            _2rho = std::pow(_2rho, +3) ;
            _arho = std::pow(_arho, +3) ;

            real_type _tsqr = std::max(
                (real_type)+0., _0bal[_dims]) ;


            for (auto _idim = _dims; _idim-- != +0; )
            {
                _move[_idim] +=
            // 1st sub-tria in voro. cell 
                _1mag*_0rho * _0bal[_idim] +
                _1mag*_1rho * _1bal[_idim] + 
                _1mag*_arho * _apos[_idim] +

            // 2nd sub-tria in voro. cell
                _2mag*_0rho * _0bal[_idim] +
                _2mag*_2rho * _2bal[_idim] + 
                _2mag*_arho * _apos[_idim] ;
            }
            
            _wsum += _1mag*(_0rho+_1rho+_arho) ;
            _wsum += _2mag*(_0rho+_2rho+_arho) ;

            _ladj += _tsqr ;
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
     * DQDX-MOVE: "local-ascent" node movement vector. 
    --------------------------------------------------------
     */
    
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdx_move_2 (
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
        real_type static const _HINC = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.50) ;
                
        real_type static const _RMIN = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.75) ;
        
        real_type static const _RMAX = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.25) ;
    
        __unreferenced(_hfun);
        __unreferenced(_pred); // for MSVC...
    
    /*------------------ calc. local characteristic scale */
        real_type _qmin =
            +std::numeric_limits
                <real_type>::infinity();
        
        real_type _bmin[_dims] = {
           +std::numeric_limits
                <real_type>::infinity()
            } ;
        real_type _bmax[_dims] = {
           -std::numeric_limits
                <real_type>::infinity()
            } ;
        
        iptr_type _tnum = +0 ;
        
        _ladj = (real_type)0.;
        
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
         
             iptr_type _tnod[3] ;
            _tnod[ 0] = _tptr->node( 0) ;
            _tnod[ 1] = _tptr->node( 1) ;
            _tnod[ 2] = _tptr->node( 2) ;

            algorithms::isort(
                &_tnod[0], &_tnod[3], 
                    std::less<iptr_type>()) ;
            
             auto _inod = 
            _mesh._set1.head()+_tnod[0] ;
             auto _jnod = 
            _mesh._set1.head()+_tnod[1] ;
             auto _knod = 
            _mesh._set1.head()+_tnod[2] ;
         
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
            
            real_type _lsqr = 
                _pred.length_sq (
                    _pmid, &_node->pval(0))  ;
            
            _ladj +=  _lsqr ;
        }

        _ladj = std::sqrt (_ladj/_tnum) ;

    /*------------------ adj. gradients w.r.t. W: dQ / dw */
        real_type _qbar , _qlow ;
        _qbar = (real_type) +1. ;
        _qlow = (real_type) +0. ;
        
        iptr_type _lnum =   +0  ;
        iptr_type _hnum =   +1  ;
        
        real_type _qlim = _qmin + 
            (real_type) +1.0E-002 ;

        _tnum = (iptr_type) +0  ;
        for (auto _tria  = _tset.head() ,
                  _tend  = _tset.tend() ;
                  _tria != _tend;
                ++_tria, ++_tnum)
        {     
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
        
            real_type _qtri = _cost[_tnum];
            
        /*-------------- only do gradients for 'poor' set */
            if (_qtri <=  _qlim)
            {
                for (auto _idim = _dims ; 
                          _idim-- != 0; )
                {
            /*---------- local iter. on finite-diff. step */      
                real_type _binc = _bmax[_idim] -
                                  _bmin[_idim] ;
           
                real_type _hdel = _HINC*_binc;
                
                real_type _hsum = (real_type)0.;
                
                real_type _sdel = (real_type)0.;
                real_type _sabs = (real_type)0.;
                real_type _sbar = (real_type)0.;
                
                for (auto _iter = +0; _iter++ != +2; )
                {
                
                /*------ centred finite-diff. for dQ / dw */
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
                
                    _sbar = _cost[_tnum];

                    _sdel = _scr1-_scr0 ;
                    
                    _sbar = std::abs(_sbar) ;
                    _sabs = std::abs(_sdel) ;
                    
                    _node->pval(_idim) 
                        = _save[_idim] ;
 
                /*------ try to adjust step on rel. diff. */
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
            
            /*---------- finalise gradient and accumulate */
                _dqdx[_idim]+= _sdel / _hsum ;
            
                }
            
                _qlow += _qtri ; _lnum += +1 ;         
            }
            else
            {
            /*---------- accumulate metrics in 'good' set */
                _qbar += _qtri ; _hnum += +1 ;
            }
        }
 
        if (_lnum > +0)  
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _dqdx[_idim]  /= _lnum ;
            }
            
            _qlow /=  _lnum ;
            _qbar /=  _hnum ;
        }
        
    /*------------------ 1st ord. Taylor-series step est. */
        real_type _scal = 
            std::sqrt(_pred.length_sq(_dqdx)) ;
        
        if (_scal*_ladj <= _HINC)
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
            _pred .innerprod(_line, _dqdx) ;
            
            _scal*= (real_type)+.50 ;

            _scal = std::min(_scal, _ladj) ;
            
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] *= _scal ;
            }
        }
    
    }
    
    
    
