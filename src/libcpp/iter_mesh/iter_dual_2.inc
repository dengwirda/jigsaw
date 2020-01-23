
    /*
    --------------------------------------------------------
     * ITER-DUAL-2: optim. schemes to reposition dual.
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
     * DQDW-MOVE: "local-ascent" weight movement vector. 
    --------------------------------------------------------
     */
    
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type dqdw_move_2 (
        geom_type &_geom ,
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        iptr_list &_tset ,
        node_iter  _node ,
        real_list &_cost ,
        real_type &_step ,
        real_type &_wadj
        )
    {
        real_type static const _WINC = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.50) ;
                
        real_type static const _RMIN = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.75) ;
        
        real_type static const _RMAX = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(), +.25) ;
    
        __unreferenced(_geom);
        __unreferenced(_hfun);
  
    /*------------------ calc. local characteristic scale */
        real_type _qmin =
            +std::numeric_limits
                <real_type>::infinity();

        iptr_type _tnum = +0 ;
        
        _step = (real_type)0.;
        _wadj = (real_type)0.;
        
        for (auto _tria  = _tset.head(),
                  _tend  = _tset.tend();
                  _tria != _tend;
                ++_tria, ++_tnum)
        {     
             auto _tptr  = 
            _mesh._set3.head()+*_tria;
         
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
                
            real_type _ball[_dims+1] ;
            _pred.circ_ball(_ball,
                &_inod->pval(0),
                &_jnod->pval(0),
                &_knod->pval(0) ) ;
                
            _wadj += _ball[_dims] ; // ball-rad.^2
        }
    
        _wadj /= _tset.count () ;
    
    /*------------------ adj. gradients w.r.t. W: dQ / dw */
        real_type _qbar, _qlow, _dqdw;
        _qbar = (real_type) +1. ;
        _qlow = (real_type) +0. ;
        _dqdw = (real_type) +0. ;
        
        iptr_type _lnum =   +0  ;
        iptr_type _hnum =   +1  ;
        
        real_type _save = 
            _node->pval  (_dims);
        
        real_type _qlim = _qmin + 
            (real_type) +1.0E-002 ;
        
        _tnum = (iptr_type) +0  ;
        for (auto _tria  = _tset.head(),
                  _tend  = _tset.tend();
                  _tria != _tend;
                ++_tria, ++_tnum)
        {     
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
        
            real_type _qtri = _cost[_tnum] ;
            real_type _DQDW ;
        
        /*-------------- only do gradients for 'poor' set */
            if (_qtri <=  _qlim)
            {
                real_type _hsum = (real_type)0.;
            
                real_type _wdel = _WINC*_wadj;
                
                real_type _sdel = (real_type)0.;
                real_type _sabs = (real_type)0.;
                real_type _sbar = (real_type)0.;
        
            /*---------- local iter. on finite-diff. step */        
                for (auto _iter = +0; _iter++ != +2; )
                {
                
                /*------ centred finite-diff. for dQ / dw */
                    _node->pval(_dims) = 
                        _save + _wdel;
                        
                    _hsum  = (real_type)+0.;
                    _hsum += _node->
                        pval(_dims) - _save;
            
                    real_type _scr1 = 
                        _pred.cost_dual (
                       &_mesh._set1[
                        _tptr->node(0)].pval(0),
                       &_mesh._set1[
                        _tptr->node(1)].pval(0),
                       &_mesh._set1[
                        _tptr->node(2)].pval(0)) ;                
                    
                    _node->pval(_dims) = 
                        _save - _wdel;

                    _hsum -= _node->
                        pval(_dims) - _save;

                    real_type _scr0 = 
                        _pred.cost_dual (
                       &_mesh._set1[
                        _tptr->node(0)].pval(0),
                       &_mesh._set1[
                        _tptr->node(1)].pval(0),
                       &_mesh._set1[
                        _tptr->node(2)].pval(0)) ;
                    
                    _sbar = _cost[_tnum];

                    _sdel = _scr1-_scr0 ;
                    
                    _sbar = std::abs(_sbar) ;
                    _sabs = std::abs(_sdel) ;
                    
                    _node->pval(_dims) =
                                _save;
 
                /*------ try to adjust step on rel. diff. */
                    if (_sabs > (real_type)0.)
                    {                   
                    if (_sabs > _RMAX * _sbar)
                    {
                        _wdel *= 
                            (real_type) +.10 ;
                    }  
                    else  
                    if (_sabs < _RMIN * _sbar)
                    {
                        _wdel *= 
                            (real_type) +10. ;
                    }
                    else { break ; }
                    }
                    else { break ; }
                }
   
                _node->pval(_dims) = _save ;
   
            /*---------- finalise gradient and accumulate */
                _DQDW = _sdel / _hsum ;
   
                _dqdw += _DQDW;
                 
                _qlow += _qtri; _lnum += 1 ;
            }
            else
            {
            /*---------- accumulate metrics in 'good' set */
                _qbar += _qtri; _hnum += 1 ;
            }
        }
    
        if (_tnum > +0)  
        {       
            _dqdw /=  _lnum ;
            _qlow /=  _lnum ;
            _qbar /=  _hnum ;
        }
        
    /*------------------ 1st ord. Taylor-series step est. */
        real_type _scal = std::abs (_dqdw) ;
        
        if (_scal*_wadj <= _WINC)
        {
            _step = (real_type) +0. ;
        }
        else
        {    
            _step = (_qbar-_qlow) / _dqdw  ;

            _step*= (real_type)+.50 ;
        }

        if (_step > (real_type) +0. )
        {
            _step = std::min(_step,+_wadj) ;
        }
        else
        if (_step < (real_type) +0. )
        {
            _step = std::max(_step,-_wadj) ;
        }
   
    }
    
    
    
