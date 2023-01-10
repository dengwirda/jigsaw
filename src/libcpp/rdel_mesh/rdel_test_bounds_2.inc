
    /*
    --------------------------------------------------------
     * RDEL-BOUNDS-2: test restricted boundaries in R^2.
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
     * Last updated: 03 Aug., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_2.hpp


    #   define __GEOMBNDS

    /*
    --------------------------------------------------------
     * FIND-RDEL: assemble restricted faces in cavity.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type find_rdel (
        geom_type &_geom ,
        mesh_type &_mesh ,
        iptr_list &_tset ,
        edat_list &_ecav
        )
    {
    /*------------------------- extract restricted k-face */
        __unreferenced(_geom) ;

        for ( auto _tpos  = _tset.head();
                   _tpos != _tset.tend();
                 ++_tpos  )
        {
    /*------------------------- extract restricted 1-face */
        for ( auto _epos = 3; _epos-- != 0; )
        {
            iptr_type _tnod[ +3] ;
            mesh_type::tria_type::tria_type::
            face_node(_tnod, _epos, 2, 1) ;
            _tnod[0] = _mesh._tria.
             tria(*_tpos)->node(_tnod[0]) ;
            _tnod[1] = _mesh._tria.
             tria(*_tpos)->node(_tnod[1]) ;

            if (_mesh._tria.
                 node(_tnod[0])->fdim()>1 ||
                _mesh._tria.
                 node(_tnod[1])->fdim()>1 )
                continue ;

            algorithms::isort (
                &_tnod[0], &_tnod[2],
                    std::less<iptr_type>()) ;

            edge_data _fdat;
            _fdat._node[0] = _tnod[0] ;
            _fdat._node[1] = _tnod[1] ;

            typename mesh_type::
                     edge_list::
                 item_type *_mptr=nullptr ;
            if (_mesh.
                 find_edge(_fdat, _mptr))
                _ecav.
                 push_tail(_mptr->_data);
        }
        }
    }

    /*
    --------------------------------------------------------
     * _CAV-BNDS: return TRUE if encroachment found.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type _cav_bnds (
        geom_type &_geom ,
        mesh_type &_mesh ,
        iptr_list &_told ,
        edat_list &_enew ,
        real_type *_pcur ,
        real_type *_pmax ,
        char_type &_mode
        )
    {
        bool_type _bnds = false;

        __unreferenced( _told );

        real_type static const _pert =
            std::pow(
        std::numeric_limits<real_type>
        ::epsilon(),(real_type)+.8);

        _pmax[2] = (real_type) +.0 ;

        if  (_mode > +1)
        {
        for ( auto _epos  = _enew.head();
                   _epos != _enew.tend();
                 ++_epos  )
        {
        /*----------------- compute restricted surf. ball */
            char_type _feat;
            char_type _topo[ +2] ;
            iptr_type _part;
            real_type _fbal[ +3] ;
            real_type _sbal[ +3] ;
            if (!mesh_pred::edge_ball  (
                _geom, _mesh ,
                _epos->_tadj ,
                _epos->_eadj ,
                _fbal, _sbal ,
                _feat, _topo , _part ) )
            __assert( false &&
                "_NEW-BNDS: edge-ball" ) ;

            real_type _blen =
            geometry::lensqr_2d(_pcur,
                                _sbal) ;

            _blen +=  _pert * _blen;

            typename mesh_type::
                     edge_list::
                 item_type *_mptr=nullptr;
            if(!_mesh.
                 find_edge(*_epos, _mptr))
            {
        /*----------------- keep worst TOPO. encroachment */
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _mode    =  + 1 ;

                _bnds    = true ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
            }
            }

        #   ifdef  __GEOMBNDS
            if (_blen    < _sbal[2])
            {
        /*----------------- keep worst GEOM. encroachment */
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _mode    =  + 1 ;

                _bnds    = true ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
            }
            }
        #   endif//__GEOMBNDS
        }
        }

        return  _bnds ;
    }

    /*
    --------------------------------------------------------
     * _CAV-BNDS: return TRUE if encroachment found.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type _cav_bnds (
        geom_type &_geom ,
        mesh_type &_mesh ,
        iptr_list &_told ,
        edat_list &_enew ,
        edat_list &_eold ,
        real_type *_ppos ,
        real_type *_pmax ,
        char_type &_mode
        )
    {
        real_type static const _pert =
            std::pow(
        std::numeric_limits<real_type>
        ::epsilon(),(real_type)+.80) ;

        bool_type _bnds = false;

        __unreferenced( _told );

        if  (_mode > +0)
        {
    /*------------------------- init. for local hash obj. */
        _mesh._etwo.clear();
        _mesh._etwo.
            set_slots(_enew.count() * +2) ;

         auto& _eset = _mesh. _etwo;

    /*---------------------------------------- test edges */
        for ( auto _epos  = _enew.head();
                   _epos != _enew.tend();
                 ++_epos  )
        {
            _eset.push(*_epos) ;
        }
        for ( auto _epos  = _eold.head();
                   _epos != _eold.tend();
                 ++_epos  )
        {
        /*----------------- compute restricted surf. ball */
            char_type _feat;
            char_type _topo[ +2] ;
            iptr_type _part;
            real_type _fbal[ +3] ;
            real_type _sbal[ +3] ;
            if (!mesh_pred::edge_ball  (
                _geom, _mesh ,
                _epos->_tadj ,
                _epos->_eadj ,
                _fbal, _sbal ,
                _feat, _topo , _part ) )
            __assert( false &&
                "_OLD-BNDS: edge-ball" ) ;

            real_type _blen =
            geometry::lensqr_2d(_ppos,
                                _sbal) ;

            _blen +=  _pert * _blen;

            typename mesh_type::
                     edge_list::
                 item_type *_mptr=nullptr;
            if(!_eset.find(*_epos, _mptr))
            {
        /*----------------- keep worst TOPO. encroachment */
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _mode    =  + 1 ;

                _bnds    = true ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
            }
            }

        #   ifdef  __GEOMBNDS
            if (_blen    < _sbal[2])
            {
        /*----------------- keep worst GEOM. encroachment */
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _mode    =  + 1 ;

                _bnds    = true ;

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
            }
            }
        #   endif//__GEOMBNDS
        }
        }

        return  _bnds  ;
    }

    /*
    --------------------------------------------------------
     * _CAV-BNDS: return TRUE if encroachment found.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type _cav_bnds (
        geom_type &_geom ,
        mesh_type &_mesh ,
        iptr_list &_told ,
        typename
    mesh_type::edge_list &_epro ,   // "protected" edges
        edat_list &_enew ,
        edat_list &_eold
        )
    {
        __unreferenced( _geom ) ;
        __unreferenced( _told ) ;

        if (_epro.count() > +0)
        {
    /*------------------------- init. for local hash obj. */
        _mesh._etwo.clear();
        _mesh._etwo.
            set_slots(_enew.count() * +2) ;

         auto& _eset = _mesh. _etwo;

    /*------------------------- push alloc. for hash obj. */
        for ( auto _epos  = _enew.head();
                   _epos != _enew.tend();
                 ++_epos  )
        {
            _eset.push(*_epos) ;
        }

        for ( auto _epos  = _eold.head();
                   _epos != _eold.tend();
                 ++_epos  )
        {
            typename mesh_type::
                     edge_list::
                 item_type *_mptr=nullptr;
            if(!_eset.find(*_epos, _mptr))
            {
            if( _epro.find(*_epos, _mptr))
            {
        /*----------------- bail on constraint violation! */
                return true ;
            }
            }
        }
        }

        return false ;
    }


    #undef __GEOMBNDS



