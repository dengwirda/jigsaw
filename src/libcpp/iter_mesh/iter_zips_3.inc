
    /*
    --------------------------------------------------------
     * ITER-ZIPS-3: optim. schemes to merge nodes.
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
    
    // from iter_mesh_3.hpp
    

    /*
    --------------------------------------------------------
     * _ZIP-FACE: try face merge to improve adj. cost-fun.
    --------------------------------------------------------
     */

    /*
    __static_call 
    __normal_call void_type _zip_face_3 (
     */


    /*
    --------------------------------------------------------
     * _ZIP-EDGE: try edge merge to improve adj. cost-fun.
    --------------------------------------------------------
     */
    
    __static_call 
    __normal_call void_type _zip_edge_3 (
        geom_type &_geom ,
        mesh_type &_mesh , 
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iter_opts &_opts ,
        iptr_type  _edge ,
        bool_type &_okay ,
        iptr_list &_iset ,
        iptr_list &_jset ,
        iptr_list &_aset ,
        iptr_list &_bset ,
        iptr_list &_cset ,
        real_list &_tsrc ,
        real_list &_tdst ,
        real_list &_ttmp ,
        real_list &_dsrc ,
        real_list &_ddst ,
        real_type  _lmax 
            = (real_type) +0.80E+00 ,
        real_type  _good 
            = (real_type) +9.25E-01 ,
        real_type  _qinc 
            = (real_type) +1.00E-08
        )
    {
        __unreferenced(_pred) ;         // for MSVC...
    
    /*--------------------------------- get edge indexing */
         auto _eptr = 
        _mesh._set2.head() + _edge;
        
        typename mesh_type::
            edge_type _edat(*_eptr) ;
        
        iptr_type _inod, _jnod ; 
        _inod = _eptr->node(0) ;
        _jnod = _eptr->node(1) ;
        
         auto _iptr = _mesh.
        _set1.head()+ _eptr->node(0);        
         auto _jptr = _mesh.
        _set1.head()+ _eptr->node(1);
        
        _okay = false ;
        
        _iset.set_count(+0) ;
        _jset.set_count(+0) ;
        _aset.set_count(+0) ;
        _bset.set_count(+0) ;
        _cset.set_count(+0) ;
        
        _tsrc.set_count(+0) ;
        _tdst.set_count(+0) ;
        _ttmp.set_count(+0) ;
        _dsrc.set_count(+0) ;
        _ddst.set_count(+0) ;
    
    /*--------------------------------- exit if FEAT node */    
        if (_iptr->feat() != mesh::null_feat)
            return ;
        if (_jptr->feat() != mesh::null_feat)
            return ;
        
    /*--------------------------------- get edge h-sizing */
        real_type _ipos[_dims + 1] ;
        real_type _jpos[_dims + 1] ;
        iptr_type _idim = +0;
        for (_idim = _dims+1; _idim-- != 0; )
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
        
        real_type  _lsqr =
            _pred.length_sq(_ipos, _jpos) ;

        real_type  _hbar = 
            std::min(_isiz , _jsiz);

    /*--------------------------------- exit if too large */
        if (_lsqr >= _hbar * _lmax *
                     _hbar * _lmax )
            return  ;
        
    /*--------------------------------- get adjacent cell */
        _mesh.node_tri4(_eptr->node (0), 
             _iset) ;
        _mesh.node_tri4(_eptr->node (1), 
             _jset) ;
                
    /*--------------------------------- exit if poor topo */
        /*
        if (_iset.count() <= +1 ||
            _jset.count() <= +1 )
            return  ;
        if (_iset.count() >= +6 &&
            _jset.count() >= +6 )
            return  ;
        if (_iset.count() +
            _jset.count() > +12 )
            return  ;
        */
        
    /*--------------------------------- get disjoint face */
        for (auto _tria  = _iset.head();
                  _tria != _iset.tend();
                ++_tria  )
        {
             auto _tptr  = 
            _mesh._set4.head()+*_tria ;
        
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
            if (_tptr->node(3) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(3) ==
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
            _mesh._set4.head()+*_tria ;
        
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
            if (_tptr->node(3) ==
                _eptr->node(0) )
                _same += +1 ;
            if (_tptr->node(3) ==
                _eptr->node(1) )
                _same += +1 ;
            if (_same >= +2) continue ;
            
            _bset.push_tail(*_tria) ;
        }
    
    /*--------------------------------- get adjacent cost */
        real_type _minA =
        loop_tscr(_mesh, _pred , 
                  _iset, 
                  _tsrc) ;
        real_type _minB = 
        loop_tscr(_mesh, _pred , 
                  _bset, 
                  _tsrc) ;
                  
        real_type _TMIN = 
        std::min( _minA, _minB);
        
        if (_TMIN>_good) return;
        
    /*--------------------------------- get adjacent ball */
        real_type _pbal[_dims+1] = {   
            (real_type) +0.000 } ;
        
        for (auto _tria  = _iset.head();
                  _tria != _iset.tend();
                ++_tria  )
        {
            real_type _ball[_dims+1] ;
                _pred. perp_ball(_ball ,
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(0)].pval(0),
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(1)].pval(0),
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(2)].pval(0),
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(3)].pval(0),
                    true ) ;
         
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
                _pred. perp_ball(_ball ,
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(0)].pval(0),
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(1)].pval(0),
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(2)].pval(0),
               &_mesh._set1[
                _mesh._set4[
               *_tria].node(3)].pval(0),
                    true ) ;
         
            for(_idim = _dims; 
                    _idim-- != +0 ; )
            {
                _pbal [_idim] +=
                    _ball [_idim] ;
            }
        }
 
        for(_idim = _dims+0; _idim-- != 0; )
        {
        _pbal[_idim] /= (_iset.count() + 
                         _jset.count() ) ;
        }
        
        _pbal[_dims+0] = (real_type) +0. ;
    
    /*--------------------------------- try to merge edge */
        typename mesh_type
               ::node_type _ndat ;
        typename mesh_type
               ::tri4_type _tdat ;
        
        _ndat.hidx() = 
            size_type::null_hint() ;
            
        _ndat.fdim() =   +3;
        _ndat.feat() = 
            mesh::null_feat;
               
        for(_idim = _dims+1; _idim-- != 0; )
        {
            _ndat.pval(_idim) = 
                    _pbal[_idim] ;
        }
  
    /*--------------------------------- push node at PBAL */      
        iptr_type _inew =
            _mesh.push_node(_ndat) ;
            
         auto _nptr  = 
        _mesh._set1.head() + _inew ;
    
        _hval.set_count(std::max (
            _inew + 1, 
        (iptr_type)_hval.count())) ;
        
        _hval[_inew] = (real_type)-1. ;
    
    /*--------------------------------- push a new cavity */
        for (auto _tria  = _aset.head();
                  _tria != _aset.tend();
                ++_tria  )
        {              
             auto _tptr  = 
            _mesh._set4.head() +*_tria ;
        
            _tdat.node(0) = 
                _tptr->node(0) ;
            _tdat.node(1) = 
                _tptr->node(1) ;
            _tdat.node(2) = 
                _tptr->node(2) ;
            _tdat.node(3) = 
                _tptr->node(3) ;
       
            if (_tdat. node(0) == 
                _edat. node(0) ) 
                _tdat. node(0) = _inew ;
            
            if (_tdat. node(1) == 
                _edat. node(0) ) 
                _tdat. node(1) = _inew ;
            
            if (_tdat. node(2) == 
                _edat. node(0) ) 
                _tdat. node(2) = _inew ;
                
            if (_tdat. node(3) == 
                _edat. node(0) ) 
                _tdat. node(3) = _inew ;           
              
            _cset.push_tail(  
            _mesh.push_tri4(_tdat)) ;
        }
    
        for (auto _tria  = _bset.head();
                  _tria != _bset.tend();
                ++_tria  )
        {              
             auto _tptr  = 
            _mesh._set4.head() +*_tria ;
        
            _tdat.node(0) = 
                _tptr->node(0) ;
            _tdat.node(1) = 
                _tptr->node(1) ;
            _tdat.node(2) = 
                _tptr->node(2) ;
            _tdat.node(3) = 
                _tptr->node(3) ;
       
            if (_tdat. node(0) == 
                _edat. node(1) ) 
                _tdat. node(0) = _inew ;
            
            if (_tdat. node(1) == 
                _edat. node(1) ) 
                _tdat. node(1) = _inew ;
            
            if (_tdat. node(2) == 
                _edat. node(1) ) 
                _tdat. node(2) = _inew ;
                
            if (_tdat. node(3) == 
                _edat. node(1) ) 
                _tdat. node(3) = _inew ;                
                
            _cset.push_tail(
            _mesh.push_tri4(_tdat)) ;
        }
    
    /*--------------------------------- optim. node coord */
        iptr_type static 
            constexpr _INUM = (iptr_type) +8 ;
               
        for (auto _iloc = +0; _iloc != _INUM ; 
                ++_iloc )
        {
            bool_type _move = false ;
        
            _ttmp.set_count(0) ;
            _dsrc.set_count(0) ;
        
            real_type  _minC = 
            loop_tscr( _mesh, _pred , 
                       _cset, 
                       _ttmp ) ;
                       
            real_type  _minD = 
            loop_dscr( _mesh, _pred , 
                       _cset, 
                       _dsrc ) ;
            
            move_node( _geom, _mesh ,
                _hfun, _pred, _hval , 
                _opts, _nptr, +1    , 
                _move, _cset, 
                _ttmp, _tdst,
                _dsrc, _ddst, 
                _minC, _good,
                _minD, _good ) ;
               
            if (_move) continue;
            
            move_node( _geom, _mesh ,
                _hfun, _pred, _hval , 
                _opts, _nptr, +2    , 
                _move, _cset, 
                _ttmp, _tdst,
                _dsrc, _ddst, 
                _minC, _good,
                _minD, _good ) ;
               
            if (_move) continue;
            
            break ;
        }
    
    /*--------------------------------- compare cost scr. */
        _tdst.set_count(0) ;
    
        loop_tscr(_mesh, _pred, _cset, 
                  _tdst) ;
    
        move_okay(_tdst, _tsrc, _okay, 
                  _good,
                  _qinc) ;
              
        if (_okay)
        {  
    /*--------------------------------- delete old cavity */   
        for (auto _tria  = _aset.head();
                  _tria != _aset.tend();
                ++_tria  )
        {              
            _mesh.
                _pop_tri4(* _tria) ;
        }
        for (auto _tria  = _jset.head();
                  _tria != _jset.tend();
                ++_tria  )
        {              
            _mesh.
                _pop_tri4(* _tria) ;
        }   
        
            _mesh.
                _pop_node(& _inod) ;
            _mesh.
                _pop_node(& _jnod) ;
        }
        else
        {
    /*--------------------------------- delete new cavity */
        for (auto _tria  = _cset.head();
                  _tria != _cset.tend();
                ++_tria  )
        {              
            _mesh.
                _pop_tri4(* _tria) ;
        }
        
            _mesh.
                _pop_node(& _inew) ;
        }
     
    }
    
    
    
