
    /*
    --------------------------------------------------------
     * RDEL-PRED-DELFRONT-2: "frontal-DEL" kernel in R^2.
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
     * Last updated: 21 Apr., 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_pred_delfront_2.hpp


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
        iptr_type _enod[ +3] ;
        mesh_type::tria_type::tria_type::
        face_node(_enod, _eadj, 2, 1);
        _enod[0] =_mesh._tria.
         tria(_tadj)->node(_enod[ 0]);
        _enod[1] =_mesh._tria.
         tria(_tadj)->node(_enod[ 1]);

    /*--------------------------------- calc. circumballs */
        if (!base_type::edge_ball (
            _geom, _mesh, _tadj,
            _eadj, _ebal, _pmax,
            _feat, _topo, _part)  )
    /*--------------------------------- is not restricted */
            return  ;

    /*------------------------- calc. refinement priority */
        _edat._mark =
        _mesh._tria.node(_enod[0])->fdim() +
        _mesh._tria.node(_enod[1])->fdim() ;

        _edat._cost =
                (float)  _pmax[2] ;

    //  scaling via node numbering seems to make slightly
    //  better (more topologically regular) meshes
    //  this implicitly boosts the refinement priority of
    //  "old" cells, which typically produces meshes with
    //  fewer "seams"...

        _edat._cost/=    _enod[0]+1  ;
        _edat._cost/=    _enod[1]+1  ;

    /*------------------------- eval. size func. at _tbal */
        real_type _esiz = (real_type)+0.;
        _esiz += (real_type)+1.*_hfun.eval(
           &_mesh._tria.
             node(_enod[0])->pval(0) ,
            _mesh._tria.
             node(_enod[0])->idxh()) ;
        _esiz += (real_type)+1.*_hfun.eval(
           &_mesh._tria.
             node(_enod[1])->pval(0) ,
            _mesh._tria.
             node(_enod[1])->idxh()) ;

         auto _hint =
        _mesh._tria.node(_enod[ 1])->idxh() ;

        _esiz += (real_type)+1.0 *
                 _hfun.eval(  _pmax, _hint) ;

        _esiz /= (real_type)+3.0 ;

        real_type _eave =_pmax [ +2] ;
        _eave /= _esiz * _esiz ;
        _eave *= (real_type)+4.0 ;

    /*------------------------- eval. surface-disc.-error */
        real_type _srat =
            geometry::lensqr_2d(
                _ebal, _pmax)/(_esiz * _esiz) ;

    /*------------------------- refinement classification */

    //  flag cells "trapped" between boundaries -- having
    //  all nodes on lower dimensional facets and with at
    //  least one boundary face.

        bool_type _trap =
            _mesh._tria.
        node(_enod[0])->fdim() < +1 &&
            _mesh._tria.
        node(_enod[1])->fdim() < +1  ;

        real_type _scal = 1. + _trap * .100 ;

        if (_eave >= _args.siz1() * _scal *
                     _args.siz1() * _scal ||
            _srat >= _args.eps1() * _scal *
                     _args.eps1() * _scal )
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

        if (_args.dims() < +2) return ;

    /*--------------------------------- get nodes in tria */
        iptr_type  _tnod[3] = {
            _mesh.
        _tria.tria(_tpos)->node(+0) ,
            _mesh.
        _tria.tria(_tpos)->node(+1) ,
            _mesh.
        _tria.tria(_tpos)->node(+2)
            } ;

    /*--------------------------------- calc. circumballs */
        real_type  _tbal[3] ;
        if (!base_type::tria_ball (
            _geom, _mesh,
            _tpos, _tbal, _part)  )
    /*--------------------------------- is not restricted */
            return  ;

    /*--------------------------------- find edge lengths */
        real_type _llen[  +3];
        _llen[0] = geometry::lensqr_2d (
           &_mesh._tria.
             node(_tnod[0])->pval(0),
           &_mesh._tria.
             node(_tnod[1])->pval(0)) ;
        _llen[1] = geometry::lensqr_2d (
           &_mesh._tria.
             node(_tnod[1])->pval(0),
           &_mesh._tria.
             node(_tnod[2])->pval(0)) ;
        _llen[2] = geometry::lensqr_2d (
           &_mesh._tria.
             node(_tnod[2])->pval(0),
           &_mesh._tria.
             node(_tnod[0])->pval(0)) ;

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
        iptr_type _ebnd = _emin ;
        real_type _best =
            std::numeric_limits
                <real_type>::infinity () ;

        bool_type _trap = false ;

        _tdat._mark =  (iptr_type) +0;

    #   ifdef __bias_bnds
        for(_enum = +3; _enum-- != +0; )
        {
            iptr_type  _enod[ +3];
            mesh_type::tria_type::
                tria_type::
            face_node(_enod, _enum, 2, 1);
            _enod[ +0] = _tnod[_enod[+0]];
            _enod[ +1] = _tnod[_enod[+1]];

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
                if (_best > _llen[_enum])
                {
                    _trap =  true ;
                    _ebnd = _enum ;
                    _best = _llen[_enum];
                }

                _tdat._mark -= +1 ;
            }
        }
    #   endif

    /*------------------------- eval. radius--edge ratios */
        real_type _erat = _tbal[    2] /
                          _llen[_emin] ;

    /*------------------------- calc. refinement priority */
        _tdat._cost =
            (float) (_tbal[2] * _erat) ;

    //  scaling via node numbering seems to make slightly
    //  better (more topologically regular) meshes
    //  this implicitly boosts the refinement priority of
    //  "old" cells, which typically produces meshes with
    //  fewer "seams"...

        iptr_type  _enod[ +3];
        mesh_type::tria_type::tria_type::
        face_node(_enod, _ebnd, 2, 1);
        _enod[ +0] = _tnod[_enod[+0]];
        _enod[ +1] = _tnod[_enod[+1]];

        _tdat._cost/=     _enod[0] + 1 ;
        _tdat._cost/=     _enod[1] + 1 ;

    /*------------------------- refinement classification */

    //  flag cells "trapped" between boundaries -- having
    //  all nodes on lower dimensional facets and with at
    //  least one boundary face.

        _trap = _trap &&
            _mesh._tria.
        node(_tnod[0])->fdim() < +2 &&
            _mesh._tria.
        node(_tnod[1])->fdim() < +2 &&
            _mesh._tria.
        node(_tnod[2])->fdim() < +2  ;

        real_type _scal = 1. ;

        if (_erat >= _args.rad2() * _scal *
                     _args.rad2() * _scal )
        {
            _kind = mesh::ring_item; return ;
        }

    /*------------------------- eval. size func. at _tbal */
        real_type _tsiz = (real_type)+0.;
        _tsiz += (real_type)+1.*_hfun.eval(
           &_mesh._tria.
             node(_tnod[0])->pval(0) ,
            _mesh._tria.
             node(_tnod[0])->idxh()) ;
        _tsiz += (real_type)+1.*_hfun.eval(
           &_mesh._tria.
             node(_tnod[1])->pval(0) ,
            _mesh._tria.
             node(_tnod[1])->idxh()) ;
        _tsiz += (real_type)+1.*_hfun.eval(
           &_mesh._tria.
             node(_tnod[2])->pval(0) ,
            _mesh._tria.
             node(_tnod[2])->idxh()) ;

         auto _hint =
        _mesh._tria.node(_tnod[ 2])->idxh() ;

        _tsiz += (real_type)+1.0 *
                 _hfun.eval(  _tbal, _hint) ;

        _tsiz /= (real_type)+4.0 ;

        real_type _eave = _tbal[ +2] ;
        _eave /= _tsiz * _tsiz;
        _eave *= (real_type)+3.0 ;

    /*------------------------- refinement classification */
        _scal += (real_type)(_trap * 0.100) ;

        if (_eave >= _args.siz2() * _scal *
                     _args.siz2() * _scal )
        {
            _kind = mesh::ring_item; return ;
        }
        else
        {
            _kind = mesh::good_item; return ;
        }
    }



