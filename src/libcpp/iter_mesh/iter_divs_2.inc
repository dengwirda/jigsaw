
    /*
    --------------------------------------------------------
     * ITER-DIVS-2: optim. schemes to split edges.
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
     * Last updated: 20 August, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
    
    // from iter_mesh_2.hpp
    

    /*
    --------------------------------------------------------
     * _DIV-EDGE: try edge split to improve adj. cost-fun.
    --------------------------------------------------------
     */
 
    __static_call
    __normal_call void_type _div_edge (
        geom_type &_geom ,
        mesh_type &_mesh , 
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iter_opts &_opts ,
        iptr_type  _enum ,
        bool_type &_okay ,
        iptr_list &_tset ,
        iptr_list &_tnew ,
        real_list &_tsrc ,
        real_list &_tdst ,
        real_list &_ttmp ,
        real_list &_dsrc ,
        real_list &_ddst ,
        real_type  _lmin
            = (real_type) +0.67E+00,
        real_type  _good 
            = (real_type) +9.25E-01,
        real_type  _qinc 
            = (real_type) +1.00E-08
        )
    {    
        __unreferenced(_pred) ;         // for MSVC...
        __unreferenced(_opts) ;

    /*--------------------------------- get edge indexing */
         auto _edge = 
        _mesh._set2.head() + _enum ;
        
        iptr_type _enod[2] ;
        _enod[ 0] = _edge->node(0) ;
        _enod[ 1] = _edge->node(1) ;
        
         auto _iptr = _mesh.
        _set1.head()+ _edge->node(0) ;        
         auto _jptr = _mesh.
        _set1.head()+ _edge->node(1) ;
        
        _tset.set_count(+0);
        _tnew.set_count(+0);
        
        _mesh.edge_tri3(_enum,_tset) ;
    
        _okay = false ;
    
        if (_tset.count()!=2) return ;
    
    /*--------------------------------- get edge h-sizing */
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
        
        real_type  _lsqr =
            _pred.length_sq(_ipos, _jpos) ;

        real_type  _hbar = 
            std::min(_isiz , _jsiz);

    /*--------------------------------- exit if too small */
        if (_lsqr <= _hbar * _lmin *
                     _hbar * _lmin )
            return  ;
         
    /*--------------------------------- get adjacent face */
         auto _itri = 
        _mesh._set3.head()+_tset[0];
         auto _jtri = 
        _mesh._set3.head()+_tset[1];
    
        _tsrc.set_count(+2) ;
        
        _tsrc[0] = _pred.cost_tria (
           &_mesh._set1[
            _itri->node(0)].pval(0),
           &_mesh._set1[
            _itri->node(1)].pval(0),
           &_mesh._set1[
            _itri->node(2)].pval(0)
                ) ;
     
        _tsrc[1] = _pred.cost_tria (
           &_mesh._set1[
            _jtri->node(0)].pval(0),
           &_mesh._set1[
            _jtri->node(1)].pval(0),
           &_mesh._set1[
            _jtri->node(2)].pval(0)
                ) ;
            
        if (_tsrc[1] < _tsrc[0])
        std::swap(_itri, _jtri ) ;
        
        real_type _TMIN = std::min (
            _tsrc[0] , _tsrc[1]) ;
        
    /*--------------------------------- exit if too good! */
        if (_TMIN > _good) return;
        
    /*--------------------------------- get adjacent ball */
        real_type _ibal[_dims+1] ;
        for(_idim = _dims+0; _idim-- != 0; )
        {
            _ibal[_idim] = 
            _mesh._set1[
        _itri->node(0)].pval(_idim) + 
            _mesh._set1[
        _itri->node(1)].pval(_idim) + 
            _mesh._set1[
        _itri->node(2)].pval(_idim) ;
            
            _ibal[_idim]/= 
                (real_type) +3.0 ;
        }

        real_type _jbal[_dims+1] ;
        for(_idim = _dims+0; _idim-- != 0; )
        {
            _jbal[_idim] = 
            _mesh._set1[
        _jtri->node(0)].pval(_idim) + 
            _mesh._set1[
        _jtri->node(1)].pval(_idim) + 
            _mesh._set1[
        _jtri->node(2)].pval(_idim) ;
            
            _jbal[_idim]/= 
                (real_type) +3.0 ;
        }
        
    /*--------------------------------- get adjacent apex */
        iptr_type _inod [3] ;
        iptr_type _jnod [3] ;
        for(auto _inum = 3; _inum-- != 0; )
        {
            mesh_type::tri3_type::
                face_node(
                    _inod, _inum, 2, 1) ;
            _inod[0] = 
            _itri->node(_inod[0]);
            _inod[1] = 
            _itri->node(_inod[1]);
            _inod[2] = 
            _itri->node(_inod[2]);
        
            if (_inod[2] != _enod[0])
            if (_inod[2] != _enod[1])
                break ;
        }
        for(auto _inum = 3; _inum-- != 0; )
        {
            mesh_type::tri3_type::
                face_node(
                    _jnod, _inum, 2, 1) ;
            _jnod[0] = 
            _jtri->node(_jnod[0]);
            _jnod[1] = 
            _jtri->node(_jnod[1]);
            _jnod[2] = 
            _jtri->node(_jnod[2]);
        
            if (_jnod[2] != _enod[0])
            if (_jnod[2] != _enod[1])
                break ;
        }
        
    /*--------------------------------- push node to mesh */
        typename mesh_type
               ::node_type   _ndat ;
        typename mesh_type
               ::tri3_type   _tdat ;
            
        _ndat.fdim() =   +2;             
        _ndat.feat() = 
            mesh::null_feat;
               
        for(_idim = _dims+0; _idim-- != 0; )
        {
            _ndat.pval(_idim) 
                = (real_type) +0.0 ;
            _ndat.pval(_idim) += 
                 _ibal[_idim] ;
            _ndat.pval(_idim) += 
                 _jbal[_idim] ;
                
            _ndat.pval(_idim)
               /= (real_type) +2.0 ;
        } 
        
        _ndat.pval(_dims)  = (real_type)0. ;
        _ndat.pval(_dims) += 
           _iptr->pval(_dims) ;
        _ndat.pval(_dims) += 
           _jptr->pval(_dims) ;
            
        _ndat.pval(_dims) /= (real_type)2. ;
        
        _ndat.hidx() = 
            size_type::null_hint() ;
         
        iptr_type _nnew = 
            _mesh.push_node(_ndat) ;
            
         auto _nptr  = 
        _mesh._set1.head() + _nnew ;
            
        _hval.set_count(std::max (
            _nnew + 1, 
        (iptr_type)_hval.count())) ;
        
        _hval[_nnew] = (real_type)-1. ;
        
    /*--------------------------------- redo mesh indexes */            
        _tdat.node(0) = _nnew;
        _tdat.node(1) = _inod[  2] ;
        _tdat.node(2) = _inod[  0] ;

        _tnew.push_tail(
        _mesh.push_tri3(_tdat)) ;
        
        _tdat.node(0) = _nnew;
        _tdat.node(1) = _inod[  1] ;
        _tdat.node(2) = _inod[  2] ;
        
        _tnew.push_tail(
        _mesh.push_tri3(_tdat)) ;
        
        _tdat.node(0) = _nnew;
        _tdat.node(1) = _jnod[  1] ;
        _tdat.node(2) = _jnod[  2] ;
        
        _tnew.push_tail(
        _mesh.push_tri3(_tdat)) ;
        
        _tdat.node(0) = _nnew;
        _tdat.node(1) = _jnod[  2] ;
        _tdat.node(2) = _jnod[  0] ;
        
        _tnew.push_tail(
        _mesh.push_tri3(_tdat)) ;
        
    /*--------------------------------- optim. node coord */
        iptr_type static 
            constexpr _INUM = (iptr_type) +4 ;
               
        for (auto _iloc = +0; _iloc != _INUM ; 
                ++_iloc )
        {
            bool_type _move = false ;
        
            _ttmp.set_count(0) ;
            _dsrc.set_count(0) ;
        
            real_type  _minC = 
            loop_tscr( _mesh, _pred , 
                       _tnew, 
                       _ttmp ) ;
                       
            real_type  _minD = 
            loop_dscr( _mesh, _pred , 
                       _tnew, 
                       _dsrc ) ;
            
            move_node( _geom, _mesh ,
                _hfun, _pred, _hval , 
                _opts, _nptr, +1    , 
                _move, _tnew, 
                _ttmp, _tdst,
                _dsrc, _ddst, 
                _minC, _good,
                _minD, _good ) ;
               
            if (_move) continue;
            
            move_node( _geom, _mesh ,
                _hfun, _pred, _hval , 
                _opts, _nptr, +2    , 
                _move, _tnew, 
                _ttmp, _tdst,
                _dsrc, _ddst, 
                _minC, _good,
                _minD, _good ) ;
               
            if (_move) continue;
            
            break ;
        }
        
    /*--------------------------------- compare cost scr. */
        _tdst.set_count(0) ;
    
        loop_tscr(_mesh, _pred, _tnew, 
                  _tdst) ;
    
        move_okay(_tdst, _tsrc, _okay, 
                  _good,
                  _qinc) ;
              
        if (_okay)
        {  
    /*--------------------------------- delete old cavity */   
        for (auto _tria  = _tset.head() ;
                  _tria != _tset.tend() ;
                ++_tria  )
        {              
            _mesh.
                _pop_tri3(* _tria) ;
        }
        }
        else
        {
    /*--------------------------------- delete new cavity */
        for (auto _tria  = _tnew.head() ;
                  _tria != _tnew.tend() ;
                ++_tria  )
        {              
            _mesh.
                _pop_tri3(* _tria) ;
        }
        
            _mesh.
                _pop_node(& _nnew) ;
        }
        
    }
    
    
    
