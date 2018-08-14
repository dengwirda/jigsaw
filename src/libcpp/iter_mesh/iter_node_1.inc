
    /*
    --------------------------------------------------------
     * ITER-NODE-1: optim. schemes to reposition nodes.
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
     * Last updated: 12 August, 2018
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
     * SPRG-MOVE: "spring"-based node movement vector. 
    --------------------------------------------------------
     */
   
    template <
        typename  node_iter
             >
    __static_call
    __normal_call void_type sprg_move_1 (
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval , 
        iptr_list &_eset ,
        node_iter  _node ,
        real_type *_line ,
        real_type &_ladj
        )
    {
        real_type _move[_dims] = {
            (real_type) +0.0 } ;

        _ladj = 
            (real_type) +0.0 ;
    
        real_type _wval = 
       +std::numeric_limits<real_type>::epsilon();
        real_type _wsum = 
       +std::numeric_limits<real_type>::epsilon();
   
        iptr_type _enum = (iptr_type)+0  ;

        real_type _smin = (real_type)-1. ;
        real_type _smax = (real_type)+1. ;
   
        for (auto _edge  = _eset.head(),
                  _eend  = _eset.tend();
                  _edge != _eend; 
                ++_edge, ++_enum )
        {     
             auto _eptr = 
            _mesh._set2.head()+*_edge ;
        
            iptr_type  _enod[2] ;
            _enod[0] = _eptr->node(0) ==
                       _node->node(0) ?
                       _eptr->node(0) :
                       _eptr->node(1) ;
                       
            _enod[1] = _eptr->node(0) ==
                       _node->node(0) ?
                       _eptr->node(1) :
                       _eptr->node(0) ;
        
             auto _inod = 
            _mesh._set1.head()+_enod[0] ;
             auto _jnod = 
            _mesh._set1.head()+_enod[1] ;
        
            real_type _ipos  [_dims] ;
            real_type _jpos  [_dims] ;
            real_type _evec  [_dims] ;
            real_type _proj  [_dims] ;
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _ipos[_idim] =_inod->
                               pval(_idim) ;
                _jpos[_idim] =_jnod->
                               pval(_idim) ;
                
                _evec[_idim] =_ipos[_idim] -
                              _jpos[_idim] ;
            }
  
            real_type _elen = 
            std::sqrt(_pred.length_sq(_evec));
  
            if (_hval[_enod[0]] 
                    < (real_type)+0.)
            {
                _hval[_enod[0]] = 
            _hfun.eval(_ipos, _inod->hidx()) ;
            }            
            
            if (_hval[_enod[1]] 
                    < (real_type)+0.)
            {
                _hval[_enod[1]] = 
            _hfun.eval(_jpos, _jnod->hidx()) ;
            }
            
            _ladj +=   _elen ;

            real_type  _hbar = 
           (real_type) +.5 * _hval[_enod[0]] +
           (real_type) +.5 * _hval[_enod[1]] ;
  
            real_type  _scal = 
               (_hbar - _elen) / _elen ;
              
            _scal = std::max(_scal, 
                             _smin) ;
            _scal = std::min(_scal, 
                             _smax) ;
            
          //_wval = std::abs(_scal) ;
          //_wval = (real_type) +1. ;
            _wval = 
            std::sqrt(std::abs(_scal)) ;
            
            _wsum   += _wval ;
            
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _proj[_idim] = 
                _ipos[_idim] + 
                    _scal * _evec[_idim] ;
            
                _move[_idim]+= 
                    _wval * _proj[_idim] ;
            }
        }
            
        if (_enum > +0)
        {
            for (auto _idim = _dims; _idim-- != +0; )
            {
                _line[_idim] = 
                    _move[_idim] / _wsum
                        - _node->pval(_idim) ;   
            }
            
            _ladj /= _enum ;
        }

    }
    
    
    
