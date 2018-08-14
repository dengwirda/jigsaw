
    /*
    --------------------------------------------------------
     * RDEL-REFINE-3: refine restricted subfaces in R^3. 
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
     * Last updated: 13 April, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
    // from rdel_mesh_3.hpp
    
    
    /*
    --------------------------------------------------------
     * _BAD-BALL: refine a "bad" ball.
    --------------------------------------------------------
     */
    
    __static_call 
    __normal_call 
    typename rdel_opts::node_kind _bad_ball (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        typename 
    mesh_type::edge_list &_pedg ,
        typename 
    mesh_type::face_list &_pfac ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        ball_heap &_nbpq ,
        edat_list &_eold ,
        edat_list &_ecav ,
        escr_list &_escr ,
        fdat_list &_fold ,
        fdat_list &_fcav ,
        fscr_list &_fscr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_bcav ,
        ball_list &_bscr ,
        char_type &_tdim ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
        class node_pred
            {
        /*----------------- find adj. set of tria-to-node */
            public  :
                typedef typename 
                    mesh_type::
                tria_type tria_type ;
  
            public  :          
                iptr_type _npos;

            public  :
            __inline_call node_pred (
                iptr_type _nsrc
                ) : _npos(_nsrc) {}
        /*----------------- find adj. set of tria-to-node */
            __inline_call 
                bool_type operator()(
                tria_type&_tria,
                iptr_type _tpos,
                iptr_type 
                ) 
            {   return ( 
            _tria.tria(_tpos)->node(+0) 
                == this->_npos ||
            _tria.tria(_tpos)->node(+1) 
                == this->_npos ||
            _tria.tria(_tpos)->node(+2) 
                == this->_npos );
            }
            } ;
    
        typename rdel_opts::node_kind 
        _kind =  rdel_opts::null_kind ;
    
        __unreferenced( _hfun) ;
        __unreferenced( _mode) ;
        __unreferenced( _pfac) ;
        __unreferenced( _eold) ;
        __unreferenced( _fold) ; 
        __unreferenced( _tdim) ;
        __unreferenced( _pass) ;
        __unreferenced( _args) ;

        _nnew.set_count(  +0 ) ;
        _nold.set_count(  +0 ) ;
        _tnew.set_count(  +0 ) ;
        _told.set_count(  +0 ) ;
        
        _escr.set_count(  +0 ) ;
        _ecav.set_count(  +0 ) ;
        _fscr.set_count(  +0 ) ;
        _fcav.set_count(  +0 ) ;
        _tscr.set_count(  +0 ) ;
        _tcav.set_count(  +0 ) ;
        
        _bscr.set_count(  +0 ) ;
        _bcav.set_count(  +0 ) ;

    /*--------------------- pop() leading element from PQ */
        bool_type _find = false;
        ball_data _ball ;
        for ( ; !_nbpq.empty() ; )
        {
    /*--------------------- cycles until valid ball found */
            _nbpq._pop_root(_ball) ;

            typename mesh_type::
                     ball_list::
                 item_type *_bptr = nullptr ;    
            if (_mesh.find_ball(_ball,_bptr))
            {
            if (_bptr->_data._pass ==
                       _ball._pass )
            {
                _find = true; break;
            }
            }
        }
        if (!_find) return (_kind) ;
 
    /*--------------------- check dist. to existing nodes */ 
        real_type static const _HTOL = 
            (real_type) +.33 * .33 ;
       
        real_type static const _STOL = 
            (real_type) +.50 * .50 ;
        
        iptr_list _tset;
        _tset.set_alloc( +32) ;

        iptr_type _nadj = 
            _ball._node [ 0]  ;

        _mesh._tria.walk_node (
            _nadj, 
         node_pred(_nadj), _tset) ;

        real_type _rmax = 
            _ball._ball [ 3]  ;

        for (auto _tpos  = _tset.head() ; 
                  _tpos != _tset.tend() ;
                ++_tpos  )
        {
        iptr_type _npos;
        for (_npos = +3; _npos-- != +0; )
        {
            iptr_type _node = 
                _mesh._tria.
            tria(*_tpos)->node(_npos) ;

            if (_node == _nadj) continue;

            iptr_type _feat = 
                _mesh.
            _tria.node(_node)->feat() ;

            real_type _dsqr = 
                geometry::lensqr_3d (
           &_mesh._tria.
             node(_nadj)->pval(0) ,
           &_mesh._tria.
             node(_node)->pval(0) ) ;

            if (_feat == hard_feat)
            {
                _rmax = std::min (
                _rmax, _HTOL * _dsqr) ;
            }
            else
            {
                _rmax = std::min (
                _rmax, _STOL * _dsqr) ;
            }           
        }    
        }
 
    /*--------------------- otherwise push "collar" nodes */
        _ball._ball[3] = std:: min (
        _ball._ball[3],_rmax) ;

        _ball._ball[3] = std::sqrt (
            _ball._ball [ 3]) ;
          
        iptr_type _hpos = 
            _mesh._tria.
             node(_nadj)->idxh() ;
        iptr_type _ndeg = 
            _mesh._tria.
             node(_nadj)->topo() ;
             
    /*--------------------- find intersections with geom. */       
        mesh::keep_all_3d <
            real_type ,
            iptr_type > _pred ;
        mesh::keep_all_3d <
            real_type ,
            iptr_type > _halo ;
        
        for (auto _iter = 32; _iter-- != 0 ; )
        {         
            _pred.clear() ;
            _halo.clear() ;
                
            real_type _prad = 
           (real_type)+1.0 * _ball._ball[3] ;
           
            real_type _hrad = 
           (real_type)+1.5 * _ball._ball[3] ;
                
            bool_type _okay = true;         
            _okay = _okay & _geom.intersect (
               &_ball._ball[0], 
                    _prad, _pred) ; // actual!
                    
            _okay = _okay & _geom.intersect (
               &_ball._ball[0], 
                    _hrad, _halo) ; // bracket
                                    
            if(!_okay)
            {
                _ball._ball[3] *=   // reduce rad. for geom.
                    (real_type) +.67 ;
            }
            else
            if (_pred._inum !=  _ndeg)
            {
                _ball._ball[3] *=   // reduce rad. for topo.
                    (real_type) +.67 ;
            }
            else { break ; }
        }  
        
    /*--------------------- push finalised "collar" nodes */
        for (auto _iter  = _pred._list.head() ;
                  _iter != _pred._list.tend() ;
                ++_iter  )
        {
            iptr_type _hint = 
                _mesh._tria.
                 node(_nadj)->next() ; 
                   
            iptr_type _node = -1;
            if (_mesh._tria.push_node(
               &_iter->pval(0), _node,
                _hint )  )
            {
                _mesh._tria.node(
                _node)->idxh() = _hpos ;
                _mesh._tria.node(
                    _node)->fdim() = 1 ;
                _mesh._tria.node(
                    _node)->feat() = 0 ;
                _mesh._tria.node(
                    _node)->topo() = 2 ;
            
                edge_data  _edat;
                _edat._node[0] = _nadj ;
                _edat._node[1] = _node ;

                _pedg.push(_edat) ;
            }        
        }
     
        return ( _kind ) ;   
    }
    
    
    
