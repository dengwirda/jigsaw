
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
     * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant
     * this code in any way whatsoever.  This code is
     * provided "as-is" to be used at your own risk.
     *
     * THE CONTRIBUTORS include:
     * (a) The University of Sydney
     * (b) The Massachusetts Institute of Technology
     * (c) Columbia University
     * (d) The National Aeronautics & Space Administration
     * (e) Los Alamos National Laboratory
     *
    --------------------------------------------------------
     *
     * Last updated: 01 Sept., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
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
        _tdat._node[0] = _mesh.
            _tria.tria(_tpos)->node(+0) ;
        _tdat._node[1] = _mesh.
            _tria.tria(_tpos)->node(+1) ;
        _tdat._node[2] = _mesh.
            _tria.tria(_tpos)->node(+2) ;

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
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        mode_type  _mode ,
        char_type  _topo ,
        real_type *_ppos ,
        typename
    rdel_opts::node_kind  _kind ,
        typename
    mesh_type::edge_list &_epro ,   // "protected" 1-edge
        iptr_list &_nnew ,
        iptr_list &_nold ,
        iptr_list &_tnew ,
        iptr_list &_told ,
        edat_list &_eold ,
        edat_list &_ecav ,
        escr_list &_escr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_bcav ,
        ball_list &_bscr ,
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

            real_type  _PPOS[3] ;
            _PPOS[0] = _ppos[0] ;
            _PPOS[1] = _ppos[1] ;
            _PPOS[2] = (real_type)+0. ; // laguerre w = 0

            iptr_type _node = -1;
            if(!_mesh._tria.push_node (
                _PPOS,  _node,
                _hint, &_tnew, &_told))
            {
        /*------------------------- bail if DT push fails */
                _kout =
                 rdel_opts ::null_kind;

                _nnew.set_count(+0) ;
                _nold.set_count(+0) ;
                _tnew.set_count(+0) ;
                _told.set_count(+0) ;

                _eold.set_count(+0) ;
                _escr.set_count(+0) ;
                _ecav.set_count(+0) ;
                _tscr.set_count(+0) ;
                _tcav.set_count(+0) ;
                _bscr.set_count(+0) ;
                _bcav.set_count(+0) ;

                break ;
            }

        /*------------------------- config. new node data */
            _mesh._tria.
                node(_node)->idxh() =
                    hfun_type::null_hint() ;

            _mesh._tria.
             node(_node)->fdim () = 0 ;
            _mesh._tria.
             node(_node)->feat () = 0 ;
            _mesh._tria.
             node(_node)->topo(0) = 2 ;
            _mesh._tria.
             node(_node)->topo(1) = 2 ;
            _mesh._tria.
             node(_node)->part () = 0 ;

            iptr_type _sign   = -1 ;
            if (_told.count() > +0 )
            {
                _sign =
            tria_sign(_mesh, _told[0]);
            }

        /*------------------------- keep "old" rDT cavity */
            find_rdel(
                _geom, _mesh, _told, _eold
                ) ;

        /*------------------------- form "new" rDT cavity */
            push_rdel( _geom, _hfun,
                _mesh,  true, _nnew, _tnew ,
                _escr, _ecav,
                _tscr, _tcav,
                _bscr, _bcav, _sign, _pass ,
                   null_mode,
                _mode, _args ) ;

            _mesh._tria.node(
                _node)->fdim() =  _tdim ;

        /*------------------------- test "new" rDT cavity */
            real_type  _pnow[2] ;
            _pnow[0] = _ppos[0] ;
            _pnow[1] = _ppos[1] ;

            if (_cav_bnds( _geom, _mesh ,
                    _told, _epro,
                    _ecav, _eold) )
            {
        /*-------- new node modifies constraint sub-faces */
                _mesh.
                _tria.roll_back(_tnew, _told);

                _tdim = +0 ;

                _kout = rdel_opts::null_kind ;

                _nnew.set_count(+0) ;
                _nold.set_count(+0) ;
                _tnew.set_count(+0) ;
                _told.set_count(+0) ;

                _eold.set_count(+0) ;
                _escr.set_count(+0) ;
                _ecav.set_count(+0) ;
                _tscr.set_count(+0) ;
                _tcav.set_count(+0) ;

                _bscr.set_count(+0) ;
                _bcav.set_count(+0) ;
            }
            else
            if (_cav_bnds( _geom, _mesh ,
                    _told, _ecav, _pnow ,
                    _ppos, _tdim) )
            {
        /*-------- new node modifies restricted sub-faces */
                _mesh.
                _tria.roll_back(_tnew, _told);

                _cav_bnds( _geom, _mesh ,
                    _told, _ecav,
                    _eold, _pnow,
                    _ppos, _tdim) ;

                _kout = rdel_opts::circ_kind ;

                _nnew.set_count(+0) ;
                _nold.set_count(+0) ;
                _tnew.set_count(+0) ;
                _told.set_count(+0) ;

                _eold.set_count(+0) ;
                _escr.set_count(+0) ;
                _ecav.set_count(+0) ;
                _tscr.set_count(+0) ;
                _tcav.set_count(+0) ;

                _bscr.set_count(+0) ;
                _bcav.set_count(+0) ;
            }
        /*------------------- new node acceptable "as-is" */
            else break ;
        }

        return   _kout ;
    }



