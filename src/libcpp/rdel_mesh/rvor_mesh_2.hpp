
    /*
    --------------------------------------------------------
     * RVOR-MESH-2: restricted voronoi mesh-gen. in R^2. 
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
     * Last updated: 9 August, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RVOR_MESH_2__
#   define __RVOR_MESH_2__

    namespace mesh {
    
    template <
    typename M ,
    typename P ,
    typename G ,
    typename H ,
    typename A = allocators::basic_alloc
             >
    class rvor_mesh_2d
    { 
    public  : 
    
    /*--------- restricted voronoi mesh-generation in R^2 */  
    
    typedef M                               mesh_type ;
    typedef P                               mesh_pred ;
    typedef G                               geom_type ;
    typedef H                               hfun_type ;
    typedef A                               allocator ;

    typedef typename 
            mesh_type::real_type            real_type ;
    typedef typename 
            mesh_type::iptr_type            iptr_type ;
            
    char_type static constexpr null_ball = +0 ;
    char_type static constexpr feat_ball = +1 ;
    char_type static constexpr voro_ball = +2 ;

    typedef char_type mode_type ;
    
    char_type static constexpr null_mode = +0 ;
    char_type static constexpr feat_mode = +1 ;
    char_type static constexpr edge_mode = +2 ;
    char_type static constexpr topo_mode = +3 ;
    char_type static constexpr ball_mode = +4 ;
    char_type static constexpr tria_mode = +5 ;

    class node_pred ; class ball_pred ; 
    class edge_pred ; class tria_pred ;
    
    class edge_cost : public mesh_pred::edge_data
        {
        public  :
        iptr_type           _node[ +2] ;
        iptr_type           _pass;
        } ;
    class tria_cost : public mesh_pred::tria_data
        {
        public  :
        iptr_type           _node[ +3] ;
        iptr_type           _pass;
        } ;

    typedef typename 
            mesh_type::node_data            node_data ;
    typedef typename 
            mesh_type::ball_data            ball_data ;
    typedef typename 
            mesh_type::edge_data            edge_data ;
    typedef typename 
            mesh_type::tria_data            tria_data ;

    typedef mesh::rdel_params       <
                real_type, 
                iptr_type           >       rdel_opts ;

    typedef containers::array       <
                iptr_type           >       iptr_list ;

/*------------------------------------------ cavity lists */
    typedef containers::array       <
                edge_data           >       edat_list ;
    
    typedef containers::array       <
                edge_cost           >       escr_list ;
    
    typedef containers::array       <
                tria_data           >       tdat_list ;
    
    typedef containers::array       <
                tria_cost           >       tscr_list ;

/*------------------------------------------ refine lists */
    typedef containers::priorityset <
                node_data, 
                node_pred           >       node_heap ;

    typedef containers::priorityset <
                edge_cost, 
                edge_pred           >       edge_heap ;

    typedef containers::priorityset <
                tria_cost, 
                tria_pred           >       tria_heap ;
                
/*------------------------------------------ collar lists */
    typedef containers::array       <
                ball_data           >       ball_list ;
                
    typedef containers::priorityset <
                ball_data, 
                ball_pred           >       ball_heap ;
                                         
    class node_pred
        {
    /*---------------------- "less-than" for node objects */
        public  :
        __inline_call bool_type operator() (
            node_data const& _idat,
            node_data const& _jdat
            ) const
        {   return _idat._pass < _jdat._pass ; 
        }
        } ;
    class ball_pred
        {
    /*---------------------- "less-than" for ball objects */
        public  :
        __inline_call bool_type operator() (
            ball_data const& _idat,
            ball_data const& _jdat
            ) const
        {   return _idat._ball[2] > 
                   _jdat._ball[2] ; 
        }
        } ;
    class edge_pred
        {
    /*---------------------- "less-than" for edge objects */
        public  :
        __inline_call bool_type operator() (
            edge_cost const& _idat,
            edge_cost const& _jdat
            ) const
        {   return mesh_pred
                 ::edge_pred(_idat, _jdat) ;
        }
        } ;
    class tria_pred
        {
    /*---------------------- "less-than" for tria objects */
        public  :
        __inline_call bool_type operator() (
            tria_cost const& _idat,
            tria_cost const& _jdat
            ) const
        {   return mesh_pred
                 ::tria_pred(_idat, _jdat) ;
        }
        } ;
  
        
    /*
    --------------------------------------------------------
     * RVOR-UPDATE: update the restricted-tria. 
    --------------------------------------------------------
     */
     
  //#include "rvor_update_face_2.inc"
    
    
    
    
    /*
    --------------------------------------------------------
     * INIT-RDEL: init. face-set in rDT. 
    --------------------------------------------------------
     */

    __static_call 
    __normal_call void_type init_rdel (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_bnds ,
        iptr_list &_nnew ,
        iptr_list &_tnew ,
        edat_list &_ecav ,
        escr_list &_escr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_nbal ,
        iptr_type  _pass ,
        mode_type  _fdim ,
        rdel_opts &_args
        )
    {
    /*-------------------- mark all existing elem. as new */
        iptr_type _npos = 0, _tpos = 0 ;
        for (auto _iter  = 
            _bnds._tria._nset.head() ; 
                  _iter != 
            _bnds._tria._nset.tend() ; 
                ++_iter  , ++_npos)
        {
            if (_iter->mark() >= +0)
            {
            _nnew. push_tail( _npos) ;
            }
        }
        for (auto _iter  = 
            _bnds._tria._tset.head() ; 
                  _iter != 
            _bnds._tria._tset.tend() ; 
                ++_iter  , ++_tpos)
        {
            if (_iter->mark() >= +0)
            {
            _tnew. push_tail( _tpos) ;
            }
        }
    /*-------------------- init. restricted triangulation */
        /*
        push_rdel( _geom, _hfun, 
            _bnds, _nnew, _tnew, 
            _escr, _ecav,
            _tscr, _tcav, 
            _nbal, _vbal,    
               -1, _pass,
            _fdim, _fdim, _args) ;
         */
    }
    
    /*
    --------------------------------------------------------
     * INIT-RVOR: init. face-set in rVT. 
    --------------------------------------------------------
     */

    __static_call 
    __normal_call void_type init_rvor (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        iptr_list &_NNEW ,
        iptr_list &_TNEW ,
        tdat_list &_TCAV ,
        tscr_list &_TSCR ,
        iptr_type  _pass ,
        mode_type  _fdim ,
        rdel_opts &_args
        )
    {
    /*-------------------- mark all existing elem. as new */
        iptr_type _npos = 0, _tpos = 0 ;
        for (auto _iter  = 
            _mesh._tria._nset.head() ; 
                  _iter != 
            _mesh._tria._nset.tend() ; 
                ++_iter  , ++_npos)
        {
            if (_iter->mark() >= +0)
            {
            _NNEW. push_tail( _npos) ;
            }
        }
        for (auto _iter  = 
            _mesh._tria._tset.head() ; 
                  _iter != 
            _mesh._tria._tset.tend() ; 
                ++_iter  , ++_tpos)
        {
            if (_iter->mark() >= +0)
            {
            _TNEW. push_tail( _tpos) ;
            }
        }
    /*-------------------- init. restricted triangulation */
    
        /*
        push_rvor( _geom, _hfun, 
            _mesh, _NNEW, _TNEW, 
            _TSCR, _TCAV, 
               -1, _pass,
            _fdim, _fdim, _args) ;
         */
    }
    
    /*
    --------------------------------------------------------
     * INIT-BNDS: init. the bounding DT. 
    --------------------------------------------------------
     */

    template <
    typename      init_type
             >
    __static_call
    __normal_call void_type init_bnds (
        geom_type &_geom,
        init_type &_init,
        hfun_type &_hfun,
        mesh_type &_bnds,
        rdel_opts &_opts
        )
    {
        __unreferenced(_init) ;
        __unreferenced(_hfun) ;

    /*------------------------------ initialise mesh bbox */
        real_type _plen[ +2];
        _plen[ 0] = _geom._bmax[ 0] - 
                    _geom._bmin[ 0] ;
        _plen[ 1] = _geom._bmax[ 1] - 
                    _geom._bmin[ 1] ;
        
        _plen[ 0]*= (real_type)+2.0 ;
        _plen[ 1]*= (real_type)+2.0 ;
        
        real_type _pmin[ +2];
        real_type _pmax[ +2];
        _pmin[ 0] = _geom._bmin[ 0] - 
                          _plen[ 0] ;
        _pmin[ 1] = _geom._bmin[ 1] - 
                          _plen[ 1] ;
        
        _pmax[ 0] = _geom._bmax[ 0] + 
                          _plen[ 0] ;
        _pmax[ 1] = _geom._bmax[ 1] + 
                          _plen[ 1] ;
        
        _bnds.
        _tria.push_root(_pmin, _pmax) ;

    /*------------------------------ initialise mesh root */
        _bnds.
        _tria.node(+0)->fdim() = +3 ;
        _bnds.
        _tria.node(+1)->fdim() = +3 ;
        _bnds.
        _tria.node(+2)->fdim() = +3 ;
        
        _bnds.
        _tria.node(+0)->feat() = +0 ;
        _bnds.
        _tria.node(+1)->feat() = +0 ;
        _bnds.
        _tria.node(+2)->feat() = +0 ;
        
        _bnds.
        _tria.node(+0)->topo() = +0 ;
        _bnds.
        _tria.node(+1)->topo() = +0 ;
        _bnds.
        _tria.node(+2)->topo() = +0 ;
    
    /*------------------------------ seed feat from geom. */
        _geom.seed_feat(_bnds, _opts) ;
                        
    /*------------------------------ seed bnds from geom. */
        _geom.seed_mesh(_bnds, _opts) ;            
    }
    
    /*
    --------------------------------------------------------
     * INIT-MESH: init. the bounding DT. 
    --------------------------------------------------------
     */

    template <
    typename      init_type
             >
    __static_call
    __normal_call void_type init_mesh (
        geom_type &_geom,
        init_type &_init,
        hfun_type &_hfun,
        mesh_type &_mesh,
        rdel_opts &_opts
        )
    {
        __unreferenced(_hfun) ;

    /*------------------------------ initialise mesh bbox */
        real_type _plen[ +2];
        _plen[ 0] = _geom._bmax[ 0] - 
                    _geom._bmin[ 0] ;
        _plen[ 1] = _geom._bmax[ 1] - 
                    _geom._bmin[ 1] ;
        
        _plen[ 0]*= (real_type)+2.0 ;
        _plen[ 1]*= (real_type)+2.0 ;
        
        real_type _pmin[ +2];
        real_type _pmax[ +2];
        _pmin[ 0] = _geom._bmin[ 0] - 
                          _plen[ 0] ;
        _pmin[ 1] = _geom._bmin[ 1] - 
                          _plen[ 1] ;
        
        _pmax[ 0] = _geom._bmax[ 0] + 
                          _plen[ 0] ;
        _pmax[ 1] = _geom._bmax[ 1] + 
                          _plen[ 1] ;
        
        _mesh.
        _tria.push_root(_pmin, _pmax) ;

    /*------------------------------ initialise mesh root */
        _mesh.
        _tria.node(+0)->fdim() = +3 ;
        _mesh.
        _tria.node(+1)->fdim() = +3 ;
        _mesh.
        _tria.node(+2)->fdim() = +3 ;
        
        _mesh.
        _tria.node(+0)->feat() = +0 ;
        _mesh.
        _tria.node(+1)->feat() = +0 ;
        _mesh.
        _tria.node(+2)->feat() = +0 ;
        
        _mesh.
        _tria.node(+0)->topo() = +0 ;
        _mesh.
        _tria.node(+1)->topo() = +0 ;
        _mesh.
        _tria.node(+2)->topo() = +0 ;
                   
    /*------------------------------ seed mesh from init. */
        for (auto _node  = 
            _init._mesh._set1.head(); 
                  _node != 
            _init._mesh._set1.tend();
                ++_node  )
        {
            if (_node->mark() >= +0)
            {
            iptr_type _npos = -1 ;
            if (_mesh._tria.push_node (
               &_node->pval(0), _npos))
            {
                _mesh._tria.node
               (_npos)->fdim() = +0 ;
                        
                _mesh._tria.node
               (_npos)->feat() = +0 ;
                        
                _mesh._tria.node
               (_npos)->topo() = +0 ;  
                    
            }     
            }
        }        
    }
    
    /*
    --------------------------------------------------------
     * RVOR-MESH: build an rVT mesh in R^2. 
    --------------------------------------------------------
     */

    template <
    typename      init_type ,
    typename      jlog_file
             >
    __static_call
    __normal_call void_type rvor_mesh (
        geom_type &_geom ,
        init_type &_init ,
        hfun_type &_hfun ,
        mesh_type &_bnds ,
        mesh_type &_mesh ,
        rdel_opts &_args ,
        jlog_file &_dump
        )
    {   
        mode_type _mode = null_mode ;
        
    /*------------------------------ push log-file header */
        _dump.push (
    "#------------------------------------------------------------\n"
    "#    |ITER.|      |DEL-1|      |DEL-2| \n"
    "#------------------------------------------------------------\n"
            ) ;
          
    /*------------------------------ ensure deterministic */  
        std::srand( +1 ) ;
           
    /*------------------------------ init. list workspace */
        iptr_list _nnew, _nold ;
        iptr_list _NNEW, _NOLD ;
        iptr_list _tnew, _told ;
        iptr_list _TNEW, _TOLD ;
        
        escr_list _escr  ;
        tscr_list _TSCR  ;
        
        edat_list _edat, _etmp ;
        tdat_list _TDAT  ;
        ball_list _nbal, _vbal ;

    /*------------------------------ refinement p.-queues */
        edge_heap _eepq  ;
        tria_heap _TTPQ  ;
        node_heap _etpq  ;
        ball_heap _nbpq  ;
        ball_heap _vbpq  ;
        
    /*------------------------------ alloc. for hash obj. */
        _mesh._eset._lptr. set_count (
        _mesh._tria._tset.count()*+3 , 
        containers::loose_alloc, nullptr);
        
        _mesh._tset._lptr. set_count (
        _mesh._tria._tset.count()*+1 , 
        containers::loose_alloc, nullptr);

    /*------------------------------ init. topo hash obj. */
        typename 
            mesh_type::node_list _etin (
        typename mesh_type::node_hash(),
        typename mesh_type::node_pred(), 
            +.8,_mesh._nset.get_alloc()) ;
        
        typename 
            mesh_type::edge_list _pedg (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(), 
            +.8,_mesh._eset.get_alloc()) ;
 
    /*------------------------------ init. point counters */
        containers:: fixed_array<
            iptr_type ,
            rdel_opts::last_kind> _enod;
        
        _enod.fill( +0 ) ;
        
        containers:: fixed_array<
            iptr_type , 
            rdel_opts::last_kind> _tnod;

        _tnod.fill( +0 ) ;

    /*------------------------------ initialise mesh obj. */
        init_bnds( _geom , _init, _hfun, 
            _bnds, _args ) ;
        
        init_mesh( _geom , _init, _hfun, 
            _mesh, _args ) ;

    /*-------------------- init. data for seed ball/nodes */ 
    
        for (auto _node  = 
            _bnds._tria._nset.head() ; 
                  _node != 
            _bnds._tria._nset.tend() ;
                ++_node  )
        {
            if (_node->mark() >= +0)
            {
                _node->idxh()  = 
                    hfun_type::null_hint() ;
            }
        }
        
        for (auto _node  = 
            _mesh._tria._nset.head() ; 
                  _node != 
            _mesh._tria._nset.tend() ;
                ++_node  )
        {
            if (_node->mark() >= +0)
            {
                _node->idxh()  = 
                    hfun_type::null_hint() ;
            }
        }

    /*-------------------- main: refine edges/faces/trias */
    
        iptr_type _pass  =   +0  ;
    
        for(bool_type _done=false; !_done ; )
        {
            iptr_type _trim_freq = +10000 ;
    #       ifdef _DEBUG
            iptr_type _jlog_freq = +250 ;
    #       else
            iptr_type _jlog_freq = +50000 ;
    #       endif

            if(++_pass>_args.iter()) break;
            
        /*------------------------- init. array workspace */
            _nnew.set_count( +0) ;
            _nold.set_count( +0) ;      // for rDT
            _tnew.set_count( +0) ;
            _told.set_count( +0) ;
            
            _etmp.set_count( +0) ;
            _edat.set_count( +0) ;
            _escr.set_count( +0) ;
            
            _nbal.set_count( +0) ;
            _vbal.set_count( +0) ;
        
            _NNEW.set_count( +0) ;      // for rVT
            _NOLD.set_count( +0) ;
            _TNEW.set_count( +0) ;
            _TOLD.set_count( +0) ;
        
            _TDAT.set_count( +0) ;
            _TSCR.set_count( +0) ;
                  
        /*--------- calc. "restricted-ness" incrementally */

            if (_mode == null_mode )
            {
        /*------------------------- init. protecting ball */
                _mode  = feat_mode ;
             
                /*
                init_rdel( _geom, _hfun, 
                    _bnds, _nnew, _tnew, 
                    _edat, _escr, 
                    _tdat, _tscr, 
                    _nbal, _pass, 
                    _mode, _args)  ;
                 */    
            }
       
            if (_mode == feat_mode &&
                    _nbpq. empty() &&
                    _nbal. empty() )
            {
        /*------------------------- init. restricted edge */
                _mode  = edge_mode ;
               
                /*
                init_rdel( _geom, _hfun, 
                    _bnds, _nnew, _tnew, 
                    _edat, _escr, 
                    _tdat, _tscr, 
                    _nbal, _pass, 
                    _mode, _args)  ;
                 */
            }          
            if (_mode == edge_mode && 
                    _eepq. empty() &&
                    _edat. empty() )
            {
        /*------------------------- init. restricted topo */
                _mode  = topo_mode ;
            }
                  
            if (_mode == topo_mode && 
                    _etpq. empty() &&
                    _edat. empty() )
            {
        /*------------------------- init. restricted dual */
                _mode  = tria_mode ;
        
                /*
                init_rvor( _geom, _hfun,    //!! different t/n sets?
                    _mesh, _nnew, _tnew, 
                    _tdat, _tscr, _pass,
                    _mode, _args)  ;
                 */
            }
            
            
            
            
            
        /*--------------------------- update faces in rDT */
                    
            for (auto _npos  = _nold.head() ;
                      _npos != _nold.tend() ; 
                    ++_npos  )
            {
                ball_data  _ball, _same;
                _ball._node[0] = *_npos;
                _ball._kind = feat_ball;                
                _bnds.
                _pop_ball( _ball, _same) ;
                
                _ball._node[0] = *_npos;
                _ball._kind = voro_ball;                
                _bnds.
                _pop_ball( _ball, _same) ;
            }
                    
            for (auto _tpos  = _told.head() ;
                      _tpos != _told.tend() ; 
                    ++_tpos  )
            {
                _pop_edge(_bnds, *_tpos) ;
            }
            for (auto _tpos  = _told.head() ;
                      _tpos != _told.tend() ; 
                    ++_tpos  )
            {
                _bnds.
                _tria._put_tria( *_tpos) ;
            }
            
            for (auto _iter  = _nbal.head() ;
                      _iter != _nbal.tend() ; 
                    ++_iter  )
            {
                _nbpq .push( *_iter ) ;
            }
            for (auto _iter  = _vbal.head() ;
                      _iter != _vbal.tend() ; 
                    ++_iter  )
            {
                _vbpq .push( *_iter ) ;
            }
            for (auto _iter  = _escr.head() ;
                      _iter != _escr.tend() ; 
                    ++_iter  )
            {
                _eepq .push( *_iter ) ;
            }
            
            for (auto _iter  = _nbal.head() ;
                      _iter != _nbal.tend() ; 
                    ++_iter  )
            {
                _bnds.push_ball( *_iter) ;
            }
            for (auto _iter  = _vbal.head() ;
                      _iter != _vbal.tend() ; 
                    ++_iter  )
            {
                _bnds.push_ball( *_iter) ;
            }
            for (auto _iter  = _edat.head() ;
                      _iter != _edat.tend() ; 
                    ++_iter  )
            {
                _bnds.push_edge( *_iter) ;
            }                   
            
        /*--------------------------- update faces in rVT */
            
            for (auto _tpos  = _TOLD.head() ;
                      _tpos != _TOLD.tend() ; 
                    ++_tpos  )
            {
                _pop_tria(_mesh, *_tpos) ;
            }
            for (auto _tpos  = _TOLD.head() ;
                      _tpos != _TOLD.tend() ; 
                    ++_tpos  )
            {
                _mesh.
                _tria._put_tria( *_tpos) ;
            }
            
            for (auto _iter  = _TSCR.head() ;
                      _iter != _TSCR.tend() ; 
                    ++_iter  )
            {
                _TTPQ .push( *_iter ) ;
            }
            for (auto _iter  = _TDAT.head() ;
                      _iter != _TDAT.tend() ; 
                    ++_iter  )
            {
                _mesh.push_tria( *_iter) ;
            }
             
        }

    /*
        if (_args.verb() >= +2 )
        {
       
    //-------------------- push refinement memory metrics *
        
        _dump.push("\n")  ;
        _dump.push("  MEMORY statistics... \n") ;
        _dump.push("\n")  ;
        
        }
         */
        
        if (_args.verb() >= +2 )
        {
    
    /*-------------------- push refinement scheme metrics */
        
        _dump.push("\n")  ;
        _dump.push("  REFINE statistics... \n") ;
        _dump.push("\n")  ;

        _dump.push("  |TYPE-1| (edge) = ");
        _dump.push(std::to_string(
             _enod[rdel_opts::circ_kind]));
        _dump.push("\n")  ;
        _dump.push("  |TYPE-2| (edge) = ");
        _dump.push(std::to_string(
             _enod[rdel_opts::offh_kind]));
        _dump.push("\n")  ;
        _dump.push("  |TYPE-D| (edge) = ");
        _dump.push(std::to_string(
             _enod[rdel_opts::disk_kind]));
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        _dump.push("  |TYPE-1| (tria) = ");
        _dump.push(std::to_string(
             _tnod[rdel_opts::circ_kind]));
        _dump.push("\n")  ;
        _dump.push("  |TYPE-2| (tria) = ");
        _dump.push(std::to_string(
             _tnod[rdel_opts::offh_kind]));
        _dump.push("\n")  ;
        _dump.push("  |TYPE-3| (tria) = ");
        _dump.push(std::to_string(
             _tnod[rdel_opts::offc_kind]));
        _dump.push("\n")  ;
        _dump.push("  |TYPE-4| (tria) = ");
        _dump.push(std::to_string(
             _tnod[rdel_opts::sink_kind]));
        _dump.push("\n")  ;
        
        }
        
        _dump.push("\n")  ;
    }  
    
    } ;

 
    }
    
#   endif   //__RVOR_MESH_2__    
    
    
    
