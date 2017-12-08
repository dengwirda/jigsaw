
    /*
    --------------------------------------------------------
     * ITER-ZIPS-2: optim. schemes to merge nodes.
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
     * Last updated: 17 November, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
    
    // from iter_mesh_2.hpp
    

    /*
    --------------------------------------------------------
     * _ZIP-EDGE: try edge merge to improve adj. cost-fun.
    --------------------------------------------------------
     */
    
    __static_call 
    __normal_call void_type _zip_edge (
        geom_type &_geom ,
        mesh_type &_mesh , 
        size_type &_hfun ,
        pred_type &_pred ,
        iptr_type  _edge ,
        bool_type &_okay ,
        iptr_list &_iset ,
        iptr_list &_jset ,
        iptr_list &_aset ,
        iptr_list &_bset ,
        real_list &_csrc ,
        real_list &_cdst ,
        real_type  _lmax 
            = (real_type) +0.50E+00 
            *  std::sqrt( +3.00E+00),
        real_type  _good 
            = (real_type) +9.00E-01 ,
        real_type  _qinc 
            = (real_type) +1.00E-08
        )
    {
         auto _eptr = 
        _mesh._set2.head() + _edge;
        
        typename mesh_type::
            edge_type _edat(*_eptr) ;
        
        _okay = false ;
        
         auto _iptr = _mesh.
        _set1.head()+ _eptr->node(0);        
         auto _jptr = _mesh.
        _set1.head()+ _eptr->node(1);
        
        _iset.set_count(+0) ;
        _jset.set_count(+0) ;
        _aset.set_count(+0) ;
        _bset.set_count(+0) ;
        
        _csrc.set_count(+0) ;
        _cdst.set_count(+0) ;
        
        real_type _ipos[_dims] ;
        real_type _jpos[_dims] ;
        iptr_type _idim = +0;
        for (_idim = _dims+0; _idim-- != 0; )
        {
            _ipos[_idim] =
                _iptr->pval(_idim) ;
            _jpos[_idim] =
                _jptr->pval(_idim) ;
        }

        real_type  _isiz = 
        _hfun.eval(_ipos, _iptr->hidx ()) ;
        real_type  _jsiz = 
        _hfun.eval(_jpos, _jptr->hidx ()) ;
        
        real_type  _elen = std::
        sqrt(_pred.length_sq(_ipos, _jpos)) ;

        real_type  _hbar = 
       (real_type) (1./2.) * _isiz + 
       (real_type) (1./2.) * _jsiz ;

        if (_elen >= _hbar * _lmax)
            return  ;
        
        _mesh.node_tri3(_eptr->node (0), 
             _iset) ;
        _mesh.node_tri3(_eptr->node (1), 
             _jset) ;
                
        if (_iset.count() <= +1 ||
            _jset.count() <= +1 )
            return  ;
        if (_iset.count() >= +6 &&
            _jset.count() >= +6 )
            return  ;
        if (_iset.count() +
            _jset.count() > +12 )
            return  ;
        
        for (auto _tria  = _iset.head();
                  _tria != _iset.tend();
                ++_tria  )
        {
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
        
            iptr_type _same =  +0 ;
            if (_tptr->node(0) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(0) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_tptr->node(1) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(1) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_tptr->node(2) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(2) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_same >= +2) continue ;
            
            _aset.push_tail(*_tria) ;
        }
   
        for (auto _tria  = _jset.head();
                  _tria != _jset.tend();
                ++_tria  )
        {
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
        
            iptr_type _same =  +0 ;
            if (_tptr->node(0) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(0) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_tptr->node(1) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(1) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_tptr->node(2) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(2) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_same >= +2) continue ;
            
            _bset.push_tail(*_tria) ;
        }
    
        loop_tscr(_mesh, _pred , 
                  _iset, 
                  _csrc) ;
                   
        loop_tscr(_mesh, _pred , 
                  _bset, 
                  _csrc) ;
        
        real_type _pbal[_dims+1] = {   
            (real_type) +0.000 } ;
        
        for (auto _tria  = _iset.head();
                  _tria != _iset.tend();
                ++_tria  )
        {
            real_type _ball[_dims+1] ;
                _pred. mini_ball(_ball ,
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(0)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(1)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(2)].pval(0)
                    ) ;
         
            for(_idim = _dims; 
                    _idim-- != +0 ; )
            {
                _pbal [_idim] +=
                    _ball [_idim] ;
            }
        }
        
        for (auto _tria  = _jset.head();
                  _tria != _jset.tend();
                ++_tria  )
        {
            real_type _ball[_dims+1] ;
                _pred. mini_ball(_ball ,
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(0)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(1)].pval(0),
               &_mesh._set1[
                _mesh._set3[
               *_tria].node(2)].pval(0)
                    ) ;
         
            for(_idim = _dims; 
                    _idim-- != +0 ; )
            {
                _pbal [_idim] +=
                    _ball [_idim] ;
            }
        }
 
        for(_idim = _dims+0; _idim-- != +0; )
        {
        _pbal[_idim] /= (_iset.count() + 
                         _jset.count() )  ;
        }
        
        real_type _isav[_dims+1] ;
        real_type _jsav[_dims+1] ;        
        for(_idim = _dims+1; _idim-- != +0; )
        {
            _isav[_idim] = 
                _iptr->pval(_idim) ;
                
            _jsav[_idim] = 
                _jptr->pval(_idim) ;
        }
        
        real_type _pmid[_dims+1] = {   
            (real_type) +0.000 } ;
            
        for(_idim = _dims+0; _idim-- != +0; )
        {
            _pmid[_idim] = 
        (real_type)+.5*_iptr->pval(_idim) +
        (real_type)+.5*_jptr->pval(_idim) ;
        }
    
        iptr_type static 
            constexpr _INUM = (iptr_type)+16 ;
        
        real_type _pnew[_dims+1] = {   
            (real_type) +0.000 } ;
        real_type _circ[_dims+1] = {   
            (real_type) +0.000 } ;
        
        _pnew[_dims+0]  = 
        (real_type) +.5 * _iptr->pval (_dims)+
        (real_type) +.5 * _jptr->pval (_dims);
               
        for (auto _iloc = +0; _iloc != _INUM ; 
                ++_iloc )
        {   
            real_type _bias = 
           (real_type)_iloc / (_INUM - 1) ;
            
            real_type _apos = 
                (real_type)+0.0+_bias;    
            real_type _bpos = 
                (real_type)+1.0-_bias;
        
            for(_idim = _dims+0; _idim-- != 0; )
            {
                _pnew[_idim] = 
                _apos*_pmid[_idim] +
                _bpos*_pbal[_idim] ;
            }
        
            _pred. proj_node(
                _geom, _circ, _pnew) ;
       
            for(_idim = _dims+1; _idim-- != 0; )
            {
                _iptr->pval(_idim) 
                    = _pnew[_idim] ;
                        
                _jptr->pval(_idim) 
                    = _pnew[_idim] ;
            }
            
            _cdst.set_count(0) ;
         
            loop_tscr(_mesh, _pred , 
                      _aset, 
                      _cdst) ;
                      
            loop_tscr(_mesh, _pred , 
                      _bset, 
                      _cdst) ;
       
            move_okay(_cdst, _csrc , 
                      _okay, _good ,
                      _qinc) ;
             
            if (_okay) break ;         
        }
                  
        if(!_okay) 
        {
            for(_idim = _dims+1; _idim-- != 0; )
            {
                _iptr->pval(_idim) 
                    = _isav[_idim] ;
                        
                _jptr->pval(_idim) 
                    = _jsav[_idim] ;
            }
            return  ;
        }
            
        for (auto _tria  = _bset.head();
                  _tria != _bset.tend();
                ++_tria  )
        {              
             auto _tptr  = 
            _mesh._set3.head()+*_tria ;
        
            typename mesh_type
                ::tri3_type _tdat ;
            _tdat.node(0) = 
                _tptr->node(0) ;
            _tdat.node(1) = 
                _tptr->node(1) ;
            _tdat.node(2) = 
                _tptr->node(2) ;
       
            if (_tdat. node(0) == 
                _edat. node(1) ) 
                _tdat. node(0) = 
                _edat. node(0) ;
            
            if (_tdat. node(1) == 
                _edat. node(1) ) 
                _tdat. node(1) = 
                _edat. node(0) ;
            
            if (_tdat. node(2) == 
                _edat. node(1) ) 
                _tdat. node(2) = 
                _edat. node(0) ;                
                
            _mesh.push_tri3(_tdat) ;
        }
        
        for (auto _tria  = _jset.head();
                  _tria != _jset.tend();
                ++_tria  )
        {              
            _mesh.
                _pop_tri3(* _tria) ;
        }
          
    }
    
    
    
