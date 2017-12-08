
    /*
    --------------------------------------------------------
     * RDEL-PRED-DELFRONT-3: "frontal-DEL" kernel in R^3. 
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
     * Last updated: 04 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
    // from rdel_pred_delfront_3.hpp
    
    
    /*
    --------------------------------------------------------
     * EDGE-COST: calc. edge refinement "cost".
    --------------------------------------------------------
     */    
    
    __static_call
    __normal_call void_type edge_cost (
        geom_type &_geom,
        hfun_type &_hfun,
        mesh_type &_mesh,
        iptr_type  _tadj,
        iptr_type  _eadj,
        rdel_opts &_args,
        edge_data &_edat,
        iptr_type &_part,
        char_type &_hits,
        char_type &_feat,
        char_type &_topo,
        char_type &_kind,
        real_type *_ebal,
        real_type *_pmax
        )
    {
        _kind = mesh::null_item ;
        
        if (_args.dims() < +1) return;

    /*--------------------------- assemble local indexing */
        iptr_type _enod[ +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_enod, _eadj, 3, 1);
        _enod[0] =_mesh._tria.
         tria(_tadj)->node(_enod[ 0]);
        _enod[1] =_mesh._tria.
         tria(_tadj)->node(_enod[ 1]);
  
    /*--------------------------------- calc. circumballs */
        if (!base_type::edge_ball ( 
            _geom, _mesh, _tadj, 
            _eadj, _ebal, 
            _pmax, _hits, 
            _feat, _topo, _part)  )
    /*--------------------------------- is not restricted */
            return  ;
                 
    /*------------------------- calc. refinement priority */
        _edat._mark =        +0  ;
        
        _edat._cost = _pmax[ +3] ;

    /*------------------------- eval. size func. at _tbal */
        real_type _esiz = (real_type)+0.;
        _esiz   += _hfun.eval(
           &_mesh._tria.
             node(_enod[0])->pval(0) , 
            _mesh._tria.
             node(_enod[0])->idxh()) ;
        _esiz   += _hfun.eval(
           &_mesh._tria.
             node(_enod[1])->pval(0) , 
            _mesh._tria.
             node(_enod[1])->idxh()) ;
        
        _esiz /= (real_type)+2.0 ;

        real_type _eave =_pmax [ +3] ;
        _eave /= _esiz * _esiz ;
        _eave *= (real_type)+4.0 ;

    /*------------------------- eval. surface-disc.-error */
        real_type _srat = 
            geometry::lensqr_3d (
                _ebal, _pmax)/(_esiz * _esiz) ;

    /*------------------------- refinement classification */
        if (_eave >= _args.siz1() * 
                     _args.siz1() ||
            _srat >= _args.eps1() *
                     _args.eps1() )
        {
            _kind = mesh::ring_item; return ;
        }
        else
        {
            _kind = mesh::good_item; return ;
        }
    }
   
    /*
    --------------------------------------------------------
     * FACE-COST: calc. face refinement "cost".
    --------------------------------------------------------
     */
   
    __static_call
    __normal_call void_type face_cost (
        geom_type &_geom,
        hfun_type &_hfun,
        mesh_type &_mesh,
        iptr_type  _tadj,
        iptr_type  _fadj,
        rdel_opts &_args,
        face_data &_fdat,
        iptr_type &_part,
        char_type &_feat,
        char_type &_topo,
        char_type &_kind,
        real_type *_tbal,
        real_type *_pmax
        )
    {
        _kind = mesh::null_item ;

        if (_args.dims() < +2) return;

    /*--------------------------- assemble local indexing */
        iptr_type _fnod[ +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_fnod, _fadj, 3, 2);
        _fnod[0] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 0]);
        _fnod[1] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 1]);
        _fnod[2] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 2]);

    /*--------------------------------- calc. circumballs */
        if (!base_type::face_ball ( 
            _geom, _mesh, _tadj, 
            _fadj, _tbal, _pmax, 
            _feat, _topo, _part)  )
    /*--------------------------------- is not restricted */
            return  ;
                             
    /*--------------------------------- find edge lengths */
        real_type _llen[  +3];
        _llen[0] = geometry::lensqr_3d (
           &_mesh._tria.
             node(_fnod[0])->pval(0),
           &_mesh._tria.
             node(_fnod[1])->pval(0)) ;
        _llen[1] = geometry::lensqr_3d (
           &_mesh._tria.
             node(_fnod[1])->pval(0),
           &_mesh._tria.
             node(_fnod[2])->pval(0)) ;
        _llen[2] = geometry::lensqr_3d (
           &_mesh._tria.
             node(_fnod[2])->pval(0),
           &_mesh._tria.
             node(_fnod[0])->pval(0)) ;

    /*--------------------------------- find min/max edge */
        iptr_type _enum ;
	    iptr_type _emin = (iptr_type)+0;
	    iptr_type _emax = (iptr_type)+0;
	    for(_enum = +3; _enum-- != +1; )
	    {
	    if (_llen[_emax] < _llen[_enum]) 
            _emax = _enum ;
	    if (_llen[_emin] > _llen[_enum]) 
            _emin = _enum ;
	    }
	
    /*------------------------- calc. refinement priority */
        _fdat._mark =   +0 ;
    
    #   ifdef __bias_BNDS      
        for(_enum = +3; _enum-- != +0; )
	    {
	        iptr_type  _enod[ +3];
            mesh_type::tria_type::
                tria_type::
            face_node(_enod, _enum, 2, 1);
            _enod[ +0] = _fnod[_enod[+0]];
            _enod[ +1] = _fnod[_enod[+1]];

	        algorithms::isort(
                &_enod[0] , &_enod[2] ,
                    std::less<iptr_type>()) ;
                    
            typename mesh_type::
                     edge_data _edat;
            _edat._node[0] = _enod[0] ;
            _edat._node[1] = _enod[1] ;
                
            typename mesh_type::
                     edge_list::
                 item_type *_eptr = nullptr ;
            if (_mesh.find_edge(_edat,_eptr))
            {
                _fdat._mark -= +1 ;
            }	    
	    }
    #   endif

    /*------------------------- eval. radius--edge ratios */
        real_type _prat = _pmax[    3] / 
                          _llen[_emin] ;
        real_type _erat = _tbal[    3] / 
                          _llen[_emin] ;

        __unreferenced(_prat) ;
        __unreferenced(_erat) ;

    /*------------------------- calc. refinement priority */
        _fdat._cost = _prat * _pmax[3] ;
    
        uint32_t  _hash = hash_ball (_pmax) ;

        real_type static const _seps = 
            std::sqrt(std::
        numeric_limits<real_type>::epsilon()) ;

        real_type _pert =
      ((real_type)(_hash%65536))/65536 ;
      
        _fdat._cost -=  
           _seps * _pert * _fdat._cost ;
         
    /*------------------------- refinement classification */
        if (_prat >= _args.rad2() *
                     _args.rad2() )
        {
            _kind = mesh::ring_item; return ;
        }
    
    /*------------------------- eval. size func. at _fbal */
        real_type _tsiz = (real_type)+0.;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_fnod[0])->pval(0) , 
            _mesh._tria.
             node(_fnod[0])->idxh()) ;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_fnod[1])->pval(0) , 
            _mesh._tria.
             node(_fnod[1])->idxh()) ;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_fnod[2])->pval(0) , 
            _mesh._tria.
             node(_fnod[2])->idxh()) ;

        _tsiz /= (real_type)+3.0 ;

        real_type _eave =_pmax [ +3] ;
        _eave /= _tsiz * _tsiz;
        _eave *= (real_type)+3.0 ;

    /*------------------------- eval. surface-disc.-error */
        real_type _srat = 
            geometry::lensqr_3d(
                _tbal, _pmax)/(_tsiz * _tsiz) ;

    /*------------------------- refinement classification */
        if (_eave >= _args.siz2() * 
                     _args.siz2() ||
            _srat >= _args.eps2() *
                     _args.eps2() )
        {
            _kind = mesh::ring_item; return ;
        }
        else
        {
            _kind = mesh::good_item; return ;
        }
    }
    
    /*
    --------------------------------------------------------
     * TRIA-COST: calc. tria refinement "cost".
    --------------------------------------------------------
     */
 
    __static_call
    __normal_call void_type tria_cost (
        geom_type &_geom,
        hfun_type &_hfun,
        mesh_type &_mesh,
        iptr_type  _tpos,
        rdel_opts &_args,
        tria_data &_tdat,
        iptr_type &_part,
        char_type &_kind
        )
    {
        _kind = mesh::null_item ;

        __unreferenced(_geom);
        __unreferenced(_part);

        if (_args.dims() < +3) return ;

    /*--------------------------------- get nodes in tria */
        iptr_type _tnod[4] = {
            _mesh.
		_tria.tria(_tpos)->node(+0) ,
            _mesh.
		_tria.tria(_tpos)->node(+1) ,
            _mesh.
		_tria.tria(_tpos)->node(+2) ,
            _mesh.
		_tria.tria(_tpos)->node(+3)
            } ;

    /*--------------------------------- tria circumcentre */
        real_type  _tbal[4] ;
        _tbal[0] = _mesh.
            _tria.tria(_tpos)->circ(+0);
        _tbal[1] = _mesh.
            _tria.tria(_tpos)->circ(+1);
        _tbal[2] = _mesh.
            _tria.tria(_tpos)->circ(+2);
  
        _tbal[ 3] = (real_type)+.0 ; 
        _tbal[ 3]+= 
            geometry::lensqr_3d (_tbal, 
               &_mesh._tria.node(
                    _tnod[0])->pval(0)) ;
        _tbal[ 3]+= 
            geometry::lensqr_3d (_tbal, 
               &_mesh._tria.node(
                    _tnod[1])->pval(0)) ;
        _tbal[ 3]+= 
            geometry::lensqr_3d (_tbal, 
               &_mesh._tria.node(
                    _tnod[2])->pval(0)) ;
        _tbal[ 3]+= 
            geometry::lensqr_3d (_tbal, 
               &_mesh._tria.node(
                    _tnod[3])->pval(0)) ;

        _tbal[ 3]/= (real_type)+4. ;

    /*------------------------- evaluate "in--out" status */
        if (_part<=-1)
        if(!_geom.is_inside(&_tbal[+0]))
        {
    /*------------------------- is not a restricted facet */
            return   ;
        }
        
        _part =   +0 ;  //!! this should come from GEOM!!
        
    /*--------------------------------- find edge lengths */
        real_type  _llen[6] ;
        iptr_type  _enum ;
	    for (_enum = +6; _enum-- != +0; )
	    {
            iptr_type _enod[ +4] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_enod, _enum, 3, 1);
            _enod[0] = _mesh._tria.
             tria(_tpos)->node(_enod[ 0]);
            _enod[1] = _mesh._tria.
             tria(_tpos)->node(_enod[ 1]);

            _llen[_enum] = 
                geometry::lensqr_3d (
           &_mesh._tria.
             node(_enod[0])->pval(0),
           &_mesh._tria.
             node(_enod[1])->pval(0)) ;
        }
        
    /*--------------------------------- find 2-face radii */
        real_type  _frad[4] ;
        iptr_type  _fnum ;
        for(_fnum = +4; _fnum-- != +0; )
        {
            iptr_type _fnod [ +4];
            mesh_type::tria_type::
                tria_type::
            face_node(_fnod, _fnum, 3, 2);
            _fnod[ +0] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 0]);
            _fnod[ +1] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 1]);
            _fnod[ +2] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 2]);

            real_type _fbal [ +4];
            geometry::tria_ball_3d(_fbal ,
           &_mesh._tria.
             node(_fnod[0])->pval(0),
           &_mesh._tria.
             node(_fnod[1])->pval(0),
           &_mesh._tria.
             node(_fnod[2])->pval(0)) ;

            _frad[_fnum] = _fbal[ 3]  ;
        }
        
    /*--------------------------------- find min/max edge */
	    iptr_type _emin = (iptr_type)+0;
	    iptr_type _emax = (iptr_type)+0;
	    for(_enum = +6; _enum-- != +1; )
	    {
	    if (_llen[_emax] < _llen[_enum]) 
            _emax = _enum ;
	    if (_llen[_emin] > _llen[_enum]) 
            _emin = _enum ;
	    }

    /*--------------------------------- find min/max face */
	    iptr_type _fmin = (iptr_type)+0;
	    iptr_type _fmax = (iptr_type)+0;
	    for(_fnum = +4; _fnum-- != +1; )
	    {
	    if (_frad[_fmax] < _frad[_fnum]) 
            _fmax = _fnum ;
	    if (_frad[_fmin] > _frad[_fnum]) 
            _fmin = _fnum ;
	    }

    /*------------------------- calc. refinement priority */
        _tdat._mark =   +0 ;
        
    #   ifdef __bias_BNDS
        for(_fnum = +4; _fnum-- != +0; )
	    {
	        iptr_type  _fnod[ +4];
            mesh_type::tria_type::
                tria_type::
            face_node(_fnod, _fnum, 3, 2);
            _fnod[0] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 0]);
            _fnod[1] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 1]);
            _fnod[2] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 2]);

	        algorithms::isort(
                &_fnod[0] , &_fnod[3] ,
                    std::less<iptr_type>()) ;
                    
            typename mesh_type::
                     face_data _fdat;
            _fdat._node[0] = _fnod[0] ;
            _fdat._node[1] = _fnod[1] ;
            _fdat._node[2] = _fnod[2] ;
                
            typename mesh_type::
                     face_list::
                 item_type *_fptr = nullptr ;
            if (_mesh.find_face(_fdat,_fptr))
            {
                _tdat._mark -= +1 ;
            }	    
	    }
    #   endif

    /*------------------------- eval. radius--edge ratios */
        real_type _erat = _tbal[   +3] / 
                          _llen[_emin] ;
                          
        real_type _frat = _frad[_fmax] / 
                          _llen[_emin] ;
                          
        __unreferenced(   _erat) ;
        __unreferenced(   _frat) ;
        
    /*------------------------- calc. refinement priority */
        _tdat._cost = _erat * _tbal[3] ;

        uint32_t  _hash = hash_ball (_tbal) ;

        real_type static const _seps = 
            std::sqrt(std::
        numeric_limits<real_type>::epsilon()) ;

        real_type _pert =
      ((real_type)(_hash%65536))/65536 ;
      
        _tdat._cost -=  
           _seps * _pert * _tdat._cost ;

    /*------------------------- refinement classification */
        if (_erat >= _args.rad3() *
                     _args.rad3() )
        {
            _kind = mesh::ring_item; return ;
        }

    /*------------------------- eval. size func. at _tbal */
        real_type _tsiz = (real_type)+0.;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_tnod[0])->pval(0) , 
            _mesh._tria.
             node(_tnod[0])->idxh()) ;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_tnod[1])->pval(0) , 
            _mesh._tria.
             node(_tnod[1])->idxh()) ;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_tnod[2])->pval(0) , 
            _mesh._tria.
             node(_tnod[2])->idxh()) ;
        _tsiz   += _hfun.eval(
           &_mesh._tria.
             node(_tnod[3])->pval(0) , 
            _mesh._tria.
             node(_tnod[3])->idxh()) ;
        
        _tsiz /= (real_type)+4.0 ;

        real_type _eave = _tbal [+3] ;
        _eave /= (_tsiz*_tsiz) ;
        _eave *= (real_type)+8.0 /
                 (real_type)+3.0 ;
                 
    /*------------------------- refinement classification */
        if (_eave >= _args.siz3() *
                     _args.siz3() )
        {
            _kind = mesh::ring_item; return ;
        }

    /*------------------------- eval. vol.--length metric */
        real_type _vlen = 
        geometry::tria_quality_3d (
           &_mesh._tria.
             node(_tnod[0])->pval(0),
           &_mesh._tria.
             node(_tnod[1])->pval(0),
           &_mesh._tria.
             node(_tnod[2])->pval(0),
           &_mesh._tria.
             node(_tnod[3])->pval(0)) ;

        if (_vlen <= _args.vol3() )
        {
            _kind = mesh::ring_item; return ;
        }
        else
        {
            _kind = mesh::good_item; return ;
        }
    }
    
    
    
