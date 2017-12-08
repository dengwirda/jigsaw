
    /*
    --------------------------------------------------------
     * RDEL-REFINE-2: refine restricted subfaces in R^2. 
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
     * Last updated: 07 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
    // from rdel_mesh_2.hpp
     
   
    /*
    --------------------------------------------------------
     * ETOP-DISK: return TRUE if 1-manifold disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type etop_disk (
        mesh_type &_mesh,
        iptr_type  _npos,
        iptr_list &_tset,
        edat_list &_eset
        )
    {
        typename 
            mesh_type::edge_list _hedg (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(), 
            .8, _mesh._eset.get_alloc()) ;
        
        _hedg._lptr.set_count (
            _tset.count() * +3 , 
        containers::loose_alloc, nullptr);
        
        iptr_type _feat = 
        _mesh._tria.node(_npos)->feat() ;
        uint_type _topo = 
        _mesh._tria.node(_npos)->topo() ;
 
        _topo = 
       (_feat != null_feat) ? _topo : 2 ;
 
        for (auto _tpos  = _tset.head() ; 
                  _tpos != _tset.tend() ; 
                ++_tpos  )
        { 
    /*--------------------- assemble 1-dim. topo disk */       
        for (auto _epos = +3; _epos-- != +0; )
        {
            iptr_type _enod [ +3];
            mesh_type::tria_type::tria_type::
            face_node(_enod, _epos, 2, 1);
            _enod[0] = _mesh._tria.
             tria(*_tpos)->node(_enod[0]);
            _enod[1] = _mesh._tria.
             tria(*_tpos)->node(_enod[1]);
            
            iptr_type _same = +0 ;
            if (_enod[0] == _npos)
                _same += +1;
            if (_enod[1] == _npos)
                _same += +1;
            if (_same != +1) continue ;
            
            algorithms::isort (
                &_enod[0], &_enod[2], 
                    std::less<iptr_type>()) ;
            
            edge_data _edat;
            _edat._node[0] = _enod[0] ;
            _edat._node[1] = _enod[1] ;
            
            typename mesh_type::
                     edge_list::
                 item_type*_mptr=nullptr;
            typename mesh_type::
                     edge_list::
                 item_type*_eptr=nullptr;
            if(!_mesh.
                 find_edge(_edat, _mptr))
                continue ;
            if( _hedg.find(_edat, _eptr))
                continue ;
            else
            {   _hedg.push(_edat) ;
            }
            
           _eset.push_tail(_mptr->_data ) ;           
        }
        }   
    
        return ( _eset.count() == _topo ) ;
    }

    /*
    --------------------------------------------------------
     * TOPO-NODE: refinement point for bad disk.
    --------------------------------------------------------
     */

    __static_call
    __normal_call 
    typename rdel_opts::node_kind topo_node (
        geom_type &_geom,
        size_type & ,
        mesh_type &_mesh,
        edat_list &_eset,
        real_type *_pmax,
        char_type &_tdim,
        iptr_type &_tadj
        )
    {
        typename rdel_opts::node_kind 
        _kind =  rdel_opts::null_kind ;
    
        _pmax[2] = (real_type)+0.;
    
    /*--------------------- find max. SDB in edge-set */
        for (auto _edge  = _eset.head() ;
                  _edge != _eset.tend() ;
                ++_edge  )
        {   
            char_type _feat, _topo ;
            iptr_type _part;
            real_type _fbal[ +3] ;
            real_type _sbal[ +3] ;
            if (!mesh_pred::edge_ball  (
                _geom, _mesh , 
                _edge->_tadj , 
                _edge->_eadj , 
                _fbal, _sbal , 
                _feat, _topo , _part ) )
            __assert( false && 
                "TOPO-NODE: edge-ball" ) ;
                      
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _kind    = 
                 rdel_opts::circ_kind  ;

                _tdim    = +1;
    
                _tadj =  _edge->_tadj  ;
    
                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
            }
        }
   /*--------------------- return steiner point kind */
        return ( _kind ) ; 
    }

    /*
    --------------------------------------------------------
     * _BAD-ETOP: refine a "bad" edge-disk.
    --------------------------------------------------------
     */

    __static_call 
    __normal_call 
    typename rdel_opts::node_kind _bad_etop (
        geom_type &_geom ,
        size_type &_size ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        typename 
    mesh_type::edge_list &_pedg ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        node_heap &_etpq ,
        typename 
    mesh_type::node_list &_etin ,
        edat_list &_etmp ,
        edat_list &_ecav ,
        escr_list &_escr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        nbal_list &_nbal ,
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

    /*--------------------- pop() leading element from PQ */
        _nnew.set_count(   +0) ;
        _nold.set_count(   +0) ;
        _tnew.set_count(   +0) ;
        _told.set_count(   +0) ;

        _etmp.set_count(   +0) ;
        _escr.set_count(   +0) ;
        _ecav.set_count(   +0) ;
        _tscr.set_count(   +0) ;
        _tcav.set_count(   +0) ;

        _nbal.set_count(   +0) ;
        
        iptr_list _tset;
        _tset.set_alloc(  +16) ;
        edat_list _eset;
        _eset.set_alloc(  +16) ;

    /*--------------------- cycles until a bad disk found */
        bool_type _find = false;
        iptr_type _npos = -1 ;
        for ( ; !_etpq.empty() ;  )
        {
            node_data _ndat, _same ;
            _etpq._pop_root( _ndat);
            _etin._pop(_ndat,_same);

            _npos = _ndat._node[ +0] ;

            _tset.set_count( +0) ;
            _eset.set_count( +0) ;

            _mesh._tria.walk_node( 
                    _npos , 
                node_pred(_npos),_tset) ;

            if(!etop_disk(_mesh, _npos,
                          _tset, _eset) )
            {
                _find = true ; break ;
            }
        }
        if (!_find) return _kind ;
    
    /*------------------------- calc. best node in cavity */
        real_type _pmax [ +3] ;
        char_type _tmax = -1;
        iptr_type _hint = -1;
      
        _kind = topo_node( _geom,
            _size , _mesh ,
            _eset , _pmax , 
            _tmax , _hint ) ;
      
    /*------------------------- push node via constraints */
        _kind = push_node( _geom, 
            _size , _mesh, _mode,
            _tmax , _pmax, _kind,
            _pedg , _nnew, _nold,
            _tnew , _told, _etmp, 
            _ecav , _escr, 
            _tcav , _tscr,
            _nbal , _hint,
            _tdim , _pass, _args) ;
        
        return ( _kind ) ;
    }
 
 

