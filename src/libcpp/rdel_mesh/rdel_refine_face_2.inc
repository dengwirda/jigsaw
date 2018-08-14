
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
     * Last updated: 08 April, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
    // from rdel_mesh_2.hpp
     
   
    /*
    --------------------------------------------------------
     * _BAD-EDGE: refine a "bad" edge.
    --------------------------------------------------------
     */

    __static_call 
    __normal_call 
    typename rdel_opts::node_kind _bad_edge (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        typename 
    mesh_type::edge_list &_pedg ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        edge_heap &_eepq ,
        edat_list &_eold ,
        edat_list &_ecav ,
        escr_list &_escr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_bcav ,
        ball_list &_bscr ,
        char_type &_tdim ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind 
        _kind =  rdel_opts::null_kind ;

        _nnew.set_count( +0) ;
        _nold.set_count( +0) ;
        _tnew.set_count( +0) ;
        _told.set_count( +0) ;
        
        _eold.set_count( +0) ;
        _escr.set_count( +0) ;
        _ecav.set_count( +0) ;
        _tscr.set_count( +0) ;
        _tcav.set_count( +0) ;
        _bscr.set_count( +0) ;
        _bcav.set_count( +0) ;
        
    /*--------------------- pop() leading element from PQ */
        iptr_type _hint= -1;
        real_type _ppos[ +3] ;
        edge_cost _qdat ;
        edge_data _edat ;
        for ( ; !_eepq.empty() ; )
        {
    /*--------------------- cycles until valid face found */
            _eepq._pop_root(_qdat) ;

            _edat._node[0] = _qdat._node[0] ;
            _edat._node[1] = _qdat._node[1] ;
            
            typename mesh_type::
                     edge_list::
                 item_type *_eptr = nullptr ;    
            if (_mesh.find_edge(_edat,_eptr))
            {
                if(_eptr->_data._pass ==
                   _qdat. _pass )
                {
            /*----------------- new steiner vertex coords */
                _kind = 
                 mesh_pred::edge_node (
                    _geom, _hfun , 
                    _mesh, _qdat , 
                    _eptr->_data._tadj, 
                    _eptr->_data._eadj, 
                    _ppos, _args ) ;

                _hint = 
                _eptr->_data._tadj ; 

                if (_kind == rdel_opts::
                             fail_kind)
                {
                    return ( rdel_opts::
                             null_kind);
                }
                else
                if (_kind == rdel_opts::
                             null_kind)
                {
            /*----------------- reject "non-frontal" edge */
                    _eepq.push ( _qdat);
                }
                else { break ; }
                }
            }
        }
        
        if (_kind == rdel_opts::null_kind)
        {
            return _kind ;
        }
 
    /*------------------------- push node via constraints */
        _kind = push_node( _geom , 
            _hfun , _mesh, _mode , 
            +1    , _ppos, _kind ,
            _pedg , _nnew, _nold ,
            _tnew , _told, _eold , 
            _ecav , _escr, 
            _tcav , _tscr,
            _bcav , _bscr, _hint,
            _tdim , _pass, _args ) ;
            
        if (_kind != rdel_opts::null_kind)
        {
            if (_tdim != +1)
            {
    /*------------------------- re-queue face if un-split */
            _escr.push_tail(_qdat) ;
            }
        }
       
        return  _kind ;
    }
    
    /*
    --------------------------------------------------------
     * _BAD-TRIA: refine a "bad" tria.
    --------------------------------------------------------
     */
    
    __static_call 
    __normal_call 
    typename rdel_opts::node_kind _bad_tria (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        typename 
    mesh_type::edge_list &_pedg ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        tria_heap &_ttpq ,
        edat_list &_eold ,
        edat_list &_ecav ,
        escr_list &_escr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_bcav ,
        ball_list &_bscr ,
        char_type &_tdim ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind 
        _kind =  rdel_opts::null_kind ;

        _nnew.set_count( +0) ;
        _nold.set_count( +0) ;
        _tnew.set_count( +0) ;
        _told.set_count( +0) ;
        
        _eold.set_count( +0) ;
        _escr.set_count( +0) ;
        _ecav.set_count( +0) ;
        _tscr.set_count( +0) ;
        _tcav.set_count( +0) ;
        _bscr.set_count( +0) ;
        _bcav.set_count( +0) ;
        
    /*--------------------- pop() leading element from PQ */
        iptr_type _hint= -1;
        real_type _ppos[ +3] ;
        tria_cost _qdat ;
        tria_data _tdat ;
        for ( ; !_ttpq.empty() ; )
        {
    /*--------------------- cycles until valid face found */
            _ttpq._pop_root(_qdat) ;

            _tdat._node[0] = _qdat._node[0] ;
            _tdat._node[1] = _qdat._node[1] ;
            _tdat._node[2] = _qdat._node[2] ;

            typename mesh_type::
                     tria_list::
                 item_type *_tptr = nullptr ;
            if (_mesh.find_tria(_tdat,_tptr))
            {
                if(_tptr->_data._pass ==
                   _qdat. _pass )
                {
            /*----------------- new steiner vertex coords */
                _kind = 
                 mesh_pred::tria_node (
                    _geom, _hfun , 
                    _mesh, _qdat ,
                    _tptr->_data._tadj , 
                    _ppos, _args ) ;

                _hint = 
                _tptr->_data._tadj ;

                if (_kind == rdel_opts::
                             fail_kind)
                {
                    return ( rdel_opts::
                             null_kind);
                }
                if (_kind == rdel_opts::
                             null_kind)
                {
            /*----------------- reject "non-frontal" face */
                    _ttpq.push ( _qdat);
                }
                else { break ; }
                }
            }
        }
        
        if (_kind == rdel_opts::null_kind)
        {
            return _kind ;
        }
 
    /*------------------------- push node via constraints */
        _kind = push_node( _geom , 
            _hfun , _mesh, _mode , 
            +2    , _ppos, _kind ,
            _pedg , _nnew, _nold ,
            _tnew , _told, _eold ,
            _ecav , _escr, 
            _tcav , _tscr,
            _bcav , _bscr, _hint ,
            _tdim , _pass, _args ) ;
            
        if (_kind != rdel_opts::null_kind)
        {
            if (_tdim != +2)
            {
    /*------------------------- re-queue face if un-split */
            _tscr.push_tail(_qdat) ;
            }
        }
       
        return  _kind ;
    }

 

