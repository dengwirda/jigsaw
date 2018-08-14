
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
     * GEOM-BNDS: full geometrical encroachment tests. 
    --------------------------------------------------------
     */
 
    // Generally, it seems that topology-only type encroa-
    // chment is better ==> meshes are sparser, etc...
 
    // define __geomBNDS
    
      
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
    /*------------------------- init. for local hash obj. */
        typename 
            mesh_type::edge_list _eset (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(), 
           +.8, _mesh._eset.get_alloc()) ;
        
        __unreferenced(_geom) ;

    /*------------------------- push alloc. for hash obj. */
        _eset._lptr.set_count (
            _tset.count() * +3 , 
        containers::loose_alloc, nullptr);
        
        for ( auto _tpos  = _tset.head();
                   _tpos != _tset.tend();
                 ++_tpos  )
        {
    /*------------------------- extract restricted 1-face */   
        iptr_type  _epos ;
        for (_epos = +3; _epos-- != +0; )
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
                 item_type *_mptr=nullptr;
            if (_eset.find(_fdat,_mptr)) 
                continue ;
            if (_mesh.
                 find_edge(_fdat,_mptr)) 
            _ecav.push_tail(_mptr->_data);
            
            _eset.push( _fdat) ;
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
            char_type _feat, _topo ;
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

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
                
                _mode    = +1 ;
                
                _bnds    = true ;
            }
            }
            
        #   ifdef __geomBNDS
            if (_blen    < _sbal[2])
            {
        /*----------------- keep worst GEOM. encroachment */
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
                
                _mode    = +1 ;
                
                _bnds    = true ;
            }
            }
        #   endif
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
        bool_type _bnds = false;
    
    /*------------------------- init. for local hash obj. */
        typename 
            mesh_type::edge_list _eset (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(), 
           +.8, _mesh._eset.get_alloc()) ;
        
    /*------------------------- push alloc. for hash obj. */
        _eset._lptr.set_count (
            _enew.count() * +2 , 
        containers::loose_alloc, nullptr);
            
        __unreferenced( _told );

        real_type static const _pert = 
            std::pow(
        std::numeric_limits<real_type>
        ::epsilon(),(real_type)+.80) ;
    
        if  (_mode > +0)
        {
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
            char_type _feat, _topo ;
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

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
                
                _mode    = +1 ;
                
                _bnds    = true ;
            }           
            }
            
        #   ifdef __geomBNDS
            if (_blen    < _sbal[2])
            {
        /*----------------- keep worst GEOM. encroachment */
            if (_pmax[2] < _sbal[2])
            {
                _pmax[2] = _sbal[2];

                _pmax[0] = _sbal[0];
                _pmax[1] = _sbal[1];
                
                _mode    = +1 ;
                
                _bnds    = true ;
            }
            }
        #   endif
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
    mesh_type::edge_list &_pedg ,
        edat_list &_enew ,
        edat_list &_eold
        )
    {
    /*------------------------- init. for local hash obj. */
        typename 
            mesh_type::edge_list _eset (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(), 
           +.8, _mesh._eset.get_alloc()) ;
        
    /*------------------------- push alloc. for hash obj. */
        _eset._lptr.set_count (
            _enew.count() * +2 , 
        containers::loose_alloc, nullptr);
            
        __unreferenced( _geom );
        __unreferenced( _told );

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
            if( _pedg.find(*_epos, _mptr))
            {
        /*----------------- bail on constraint violation! */
                return true ;
            }
            }
        }
        
        return false ;
    }


    #undef  __geomBNDS
    
    
    
