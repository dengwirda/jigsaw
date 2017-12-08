
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
     * TRIA-SIGN: return sign/part assoc. with tria. 
    --------------------------------------------------------
     */
     
    __static_call
    __normal_call iptr_type tria_sign (
        mesh_type &_mesh,
        iptr_type  _tpos
        )
    {
        iptr_type _sign = -1 ;

        tria_data _tdat;
        _tdat._node[0] = 
        _mesh._tria.tria(_tpos)->node(+0) ;
        _tdat._node[1] = 
        _mesh._tria.tria(_tpos)->node(+1) ;
        _tdat._node[2] = 
        _mesh._tria.tria(_tpos)->node(+2) ;
        
        typename mesh_type::
                 tria_list::
             item_type *_tptr = nullptr ;
        if (_mesh.find_tria(_tdat,_tptr))
        {
            _sign =  _tptr->_data._part ;
        }

        return ( _sign ) ;
    }

    /*
    --------------------------------------------------------
     * PUSH-NODE: try to push a new node onto rDT.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind push_node (
        geom_type &_geom ,
        size_type &_size ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        char_type  _topo ,
        real_type *_ppos ,
        typename 
    rdel_opts::node_kind  _kind ,
        typename 
    mesh_type::edge_list &_pedg ,
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        edat_list &_eold ,
        edat_list &_ecav ,
        escr_list &_escr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        nbal_list &_nbal ,
        iptr_type  _hint ,
        char_type &_tdim ,
        iptr_type  _pass ,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind 
             _kout = _kind;
    
    /*---------------- try lower-dim. faces until success */
        for (_tdim = _topo; _tdim > 0 ; )
        {
            _eold.set_count(  +0) ;
        
            iptr_type _node = -1;
            if(!_mesh._tria.push_node (
                _ppos,  _node, 
                _hint, &_tnew, &_told))
            {
        /*------------------------- bail if DT push fails */
                _nnew.set_count(+0) ;
                _nold.set_count(+0) ;
                _tnew.set_count(+0) ;
                _told.set_count(+0) ;
                
                _eold.set_count(+0) ;
                _escr.set_count(+0) ;
                _ecav.set_count(+0) ;
                _tscr.set_count(+0) ;
                _tcav.set_count(+0) ;
                
                _nbal.set_count(+0) ;
                
                _kout = 
                 rdel_opts ::null_kind;
            
                break ;
            }
           
        /*------------------------- config. new node data */
            _mesh._tria.
                node(_node)->idxh() = 
                    size_type::null_hint() ;

            _mesh._tria.
             node(_node)->fdim() = 0  ;
            _mesh._tria.
             node(_node)->feat() = 0  ;
            _mesh._tria.
             node(_node)->topo() = 0  ;

            iptr_type _sign   = -1 ;
            if (_told.count() > +0 )
            {
                _sign = 
            tria_sign(_mesh, _told[0]);
            }

        /*------------------------- keep "old" rDT cavity */
            find_rdel( _geom, _mesh , 
                       _told, _eold ) ;

        /*------------------------- form "new" rDT cavity */
            mode_type  _dim0=null_mode;
            push_rdel( _geom, _size , 
                _mesh, _nnew, _tnew , 
                _escr, _ecav, 
                _tscr, _tcav,
                _nbal, _sign, _pass , 
                _dim0, _mode, _args ) ;
                
            _mesh._tria.
             node(_node)->fdim()=_tdim;

            real_type  _pnow[3] ;
            _pnow[0] = _ppos[0] ;
            _pnow[1] = _ppos[1] ;

            if (_cav_bnds( _geom, _mesh , 
                    _told, _pedg, _ecav , 
                    _eold) )
            {
        /*-------- new node modifies restricted sub-faces */
                _mesh._tria.roll_back (
                    _tnew, _told) ;

                _tdim = +0 ;

                _kout = 
                 rdel_opts::null_kind ;
                
                _nnew.set_count(  +0) ;
                _nold.set_count(  +0) ;
                _tnew.set_count(  +0) ;
                _told.set_count(  +0) ;
                
                _eold.set_count(  +0) ;
                _escr.set_count(  +0) ;
                _ecav.set_count(  +0) ;
                _tscr.set_count(  +0) ;
                _tcav.set_count(  +0) ;
                
                _nbal.set_count(  +0) ;
            }
            else
            if (_cav_bnds( _geom, _mesh , 
                    _told, _ecav, _pnow , 
                    _ppos, _tdim) )
            {
        /*-------- new node modifies restricted sub-faces */
                _mesh._tria.roll_back (
                    _tnew, _told) ;

                _cav_bnds( _geom, _mesh , 
                    _told, _ecav, 
                    _eold, _pnow, 
                    _ppos, _tdim) ;

                _kout = 
                 rdel_opts::circ_kind ;
                
                _nnew.set_count(  +0) ;
                _nold.set_count(  +0) ;
                _tnew.set_count(  +0) ;
                _told.set_count(  +0) ;
                
                _eold.set_count(  +0) ;
                _escr.set_count(  +0) ;
                _ecav.set_count(  +0) ;
                _tscr.set_count(  +0) ;
                _tcav.set_count(  +0) ;
                
                _nbal.set_count(  +0) ;
            }
        /*------------------- new node acceptable "as-is" */
            else break ;
        }
        
        return   _kout ;
    }

 

