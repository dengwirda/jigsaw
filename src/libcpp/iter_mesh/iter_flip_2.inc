
    /*
    --------------------------------------------------------
     * ITER-FLIP-2: optim. schemes to "flip" topology.
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
     * Last updated: 12 December, 2017
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
    __static_call
    __normal_call bool_type need_flip_costfunc (
     */
    
    
    __static_call
    __normal_call bool_type need_flip_delaunay (
        geom_type &_geom ,
        mesh_type &_mesh ,
        pred_type &_pred ,
    __const_ptr  (iptr_type) _inod ,
    __const_ptr  (iptr_type) _jnod
        )
    {
        real_type static const _LTOL = 
        std::pow(std::numeric_limits
       <real_type>::epsilon(), +.80) ;

        __unreferenced(_pred) ; // for MSVC...

        iptr_type  _iloc[3] ;
        _iloc[0] = _inod[0] ;
        _iloc[1] = _inod[1] ;
        _iloc[2] = _inod[2] ;
        
        algorithms::isort(
            &_iloc[0], &_iloc[3] , 
                std::less<iptr_type>());
        
        iptr_type  _jloc[3] ;
        _jloc[0] = _jnod[0] ;
        _jloc[1] = _jnod[1] ;
        _jloc[2] = _jnod[2] ;
        
        algorithms::isort(
            &_jloc[0], &_jloc[3] , 
                std::less<iptr_type>());
        
        real_type _ibal[_dims+1] ;
        _pred.circ_ball(_ibal ,
       &_mesh._set1[_iloc[0]].pval(0),
       &_mesh._set1[_iloc[1]].pval(0),
       &_mesh._set1[_iloc[2]].pval(0)
            ) ;
            
        real_type _jbal[_dims+1] ;
        _pred.circ_ball(_jbal ,
       &_mesh._set1[_jloc[0]].pval(0),
       &_mesh._set1[_jloc[1]].pval(0),
       &_mesh._set1[_jloc[2]].pval(0)
            ) ;

        real_type _null[_dims] = {
            (real_type) +0.  } ;

        _pred. proj_node (
            _geom, _null, _ibal) ;
            
        _pred. proj_node (
            _geom, _null, _jbal) ;

        real_type  _ilen = 
            _pred.length_sq(_jbal, 
      &_mesh._set1[_inod[2]].pval(0));
      
        real_type  _jlen = 
            _pred.length_sq(_ibal, 
      &_mesh._set1[_jnod[2]].pval(0));
      
        real_type _btol = _LTOL * 
            std::max(_ilen, _jlen) ;
        
        real_type _idel = 
            _ilen - _jbal[_dims] ;
        real_type _jdel = 
            _jlen - _ibal[_dims] ;
        
        if (_idel >= - _btol ||
                _jdel >= - _btol )
        {
            return false ;
        }
        else
        {
            return  true ;
        }          
    }
    
    __static_call
    __normal_call bool_type need_flip_weighted (
        geom_type &_geom ,
        mesh_type &_mesh ,
        pred_type &_pred ,
    __const_ptr  (iptr_type) _inod ,
    __const_ptr  (iptr_type) _jnod
        )
    {
        real_type static const _LTOL = 
        std::pow(std::numeric_limits
       <real_type>::epsilon(), +.80) ;

        __unreferenced(_pred) ; // for MSVC...

        iptr_type  _iloc[3] ;
        _iloc[0] = _inod[0] ;
        _iloc[1] = _inod[1] ;
        _iloc[2] = _inod[2] ;
        
        algorithms::isort(
            &_iloc[0], &_iloc[3] , 
                std::less<iptr_type>());
        
        iptr_type  _jloc[3] ;
        _jloc[0] = _jnod[0] ;
        _jloc[1] = _jnod[1] ;
        _jloc[2] = _jnod[2] ;
        
        algorithms::isort(
            &_jloc[0], &_jloc[3] , 
                std::less<iptr_type>());
        
        real_type _ibal[_dims+1] ;
        _pred.perp_ball(_ibal ,
       &_mesh._set1[_iloc[0]].pval(0),
       &_mesh._set1[_iloc[1]].pval(0),
       &_mesh._set1[_iloc[2]].pval(0)
            ) ;
            
        real_type _jbal[_dims+1] ;
        _pred.perp_ball(_jbal ,
       &_mesh._set1[_jloc[0]].pval(0),
       &_mesh._set1[_jloc[1]].pval(0),
       &_mesh._set1[_jloc[2]].pval(0)
            ) ;

        real_type _null[_dims] = {
            (real_type) +0.  } ;

        _pred. proj_node (
            _geom, _null, _ibal) ;
            
        _pred. proj_node (
            _geom, _null, _jbal) ;

        real_type  _ilen = 
            _pred.length_sq(_jbal, 
      &_mesh._set1[_inod[2]].pval(   +0));
      
        real_type  _ipwr = 
       _mesh._set1[_inod[2]].pval(_dims) ;
      
        real_type  _jlen = 
            _pred.length_sq(_ibal, 
      &_mesh._set1[_jnod[2]].pval(   +0));
    
        real_type  _jpwr = 
       _mesh._set1[_jnod[2]].pval(_dims) ;
    
        _ilen -=   _ipwr ;
        _jlen -=   _jpwr ;
        
        real_type _btol = _LTOL * 
            std::max(_ilen, _jlen) ;
        
        real_type _idel = 
            _ilen - _jbal[_dims] ;
        real_type _jdel = 
            _jlen - _ibal[_dims] ;
        
        if (_idel >= - _btol ||
                _jdel >= - _btol )
        {
            return false ;
        }
        else
        {
            return  true ;
        }       
    } 
    
     
    /*
    --------------------------------------------------------
     * FLIP-T2T2: 2-simplex topological flip. 
    --------------------------------------------------------
     */
    
    __static_call
    __normal_call void_type flip_t2t2 (
        geom_type &_geom ,
        mesh_type &_mesh ,
        size_type &_hfun ,
        pred_type &_pred ,
        iptr_type  _tria ,
        iptr_type  _edge ,
        iptr_list &_told ,
        iptr_list &_tnew ,
        bool_type &_flip , 
        real_list &_qold ,
        real_list &_qnew
        )
    {
        __unreferenced(  _hfun) ;
        __unreferenced(  _qold) ;
        __unreferenced(  _qnew) ;
    
        _flip = false ;
    
        _told.set_count(0) ;      
        _tnew.set_count(0) ;
    
        iptr_type _enod[3] ;
        mesh_type::tri3_type::
        face_node(_enod, _edge, 2, 1) ;
        _enod[ 0] = _mesh.
        _set3[_tria].node(_enod[0]) ;
        _enod[ 1] = _mesh.
        _set3[_tria].node(_enod[1]) ;
      
        iptr_type _epos = -1 ;
        _mesh.find_edge(_enod, _epos) ;
        
        if (_epos==-1) return;
        
         auto _eptr = 
        _mesh._set2.head() + _epos ;
        
        if (_eptr->self()>=+1) return ;
      
        _mesh.edge_tri3(_enod, _told) ;
    
        if (_told.count()!=+2) return ;
    
        iptr_type _itri = _told[0] ;
        iptr_type _jtri = _told[1] ;
        
         auto _iptr = 
        _mesh._set3.head() + _itri ;
         auto _jptr = 
        _mesh._set3.head() + _jtri ;
    
        iptr_type _inod [3] ;
        iptr_type _jnod [3] ;
        for(auto _inum = 3; _inum-- != 0; )
        {
            mesh_type::tri3_type::
            face_node(_inod, _inum, 2, 1) ;
            _inod[0] = 
                _iptr->node(_inod[0]);
            _inod[1] = 
                _iptr->node(_inod[1]);
            _inod[2] = 
                _iptr->node(_inod[2]);
        
            if (_inod[2] != _enod[0])
            if (_inod[2] != _enod[1])
                break ;
        }
        for(auto _inum = 3; _inum-- != 0; )
        {
            mesh_type::tri3_type::
            face_node(_jnod, _inum, 2, 1) ;
            _jnod[0] = 
                _jptr->node(_jnod[0]);
            _jnod[1] = 
                _jptr->node(_jnod[1]);
            _jnod[2] = 
                _jptr->node(_jnod[2]);
        
            if (_jnod[2] != _enod[0])
            if (_jnod[2] != _enod[1])
                break ;
        }

        __assert(_inod[0]==_jnod[1] &&
                 _inod[1]==_jnod[0] &&
        "FLIP-T2T2: bad orientation!");

        if(!need_flip_weighted (
                _geom, _mesh , 
                _pred, 
                _inod, _jnod)) return ;
        
        _mesh._pop_tri3(_itri) ;
        _mesh._pop_tri3(_jtri) ;
        
        _flip = true  ;
       
        typename mesh_type
               ::tri3_type _tdat ;
        _tdat.node(0) = _inod[0] ;
        _tdat.node(1) = _jnod[2] ;
        _tdat.node(2) = _inod[2] ;
        
        _tnew.push_tail(
            _mesh.push_tri3(_tdat)) ;
        
        _tdat.node(0) = _jnod[0] ;
        _tdat.node(1) = _inod[2] ;
        _tdat.node(2) = _jnod[2] ;
        
        _tnew.push_tail(
            _mesh.push_tri3(_tdat)) ;         
    }
    
    
    
