
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
     * Last updated: 22 Apr., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     *
     * This class defines the behaviour of the "frontal"
     * restricted delaunay refinement scheme for meshes
     * in R^3. Routines are provided to assess the
     * "restricted-ness" of faces, to compute refinement
     * "costs" (priorities), and to insert new vertices.
     *
     * My implementation is described in:
     *
     * D. Engwirda and D. Ivers, (2016): Off-centre Steiner
     * points for Delaunay-refinement on curved surfaces,
     * Computer-Aided Design, 72, pp. 157-171,
     * http://dx.doi.org/10.1016/j.cad.2015.10.007
     *
     * D. Engwirda, (2016): "Conforming restricted Delaunay
     * mesh generation for piecewise smooth complexes",
     * Procedia Engineering, 163, pp. 84-96,
     * http://dx.doi.org/10.1016/j.proeng.2016.11.024
     *
     * D. Engwirda, (2014): "Locally-optimal Delaunay-
     * refinement and optimisation-based mesh generation",
     * Ph.D. Thesis, School of Mathematics and Statistics,
     * Univ. of Sydney.
     * http://hdl.handle.net/2123/13148
     *
     * which is based on various previous works, including
     * (primarily):
     *
     * J.D. Boissonnat, S. Oudot, (2005): "Provably Good
     * Sampling and Meshing of Surfaces", Graphical Models,
     * 67, pp. 405-451,
     * https://doi.org/10.1016/j.gmod.2005.01.004
     *
     * C. Jamin, P. Alliez, M. Yvinec, and J.D. Boissonnat,
     * (2015): "CGALmesh: a generic framework for Delaunay
     * mesh generation", ACM Transactions on Mathematical
     * Software (TOMS), 41, pp. 23
     * https://doi.org/10.1145/2699463
     *
     * S.W. Cheng, T.K. Dey, E.A. Ramos, (2010): "Delaunay
     * Refinement for Piecewise Smooth Complexes",
     * Discrete & Computational Geometry, 43, pp. 121-166,
     * https://doi.org/10.1007/s00454-008-9109-3
     *
     * The "frontal" insertion scheme is based on the use
     * of "off-centre" refinement rules - a generalisation
     * of ideas introduced in:
     *
     * S. Rebay, (1993): "Efficient Unstructured Mesh
     * Generation by Means of Delaunay Triangulation and
     * the Bowyer-Watson Algorithm", J. Comp. Phys., 106,
     * pp. 125-138
     * https://doi.org/10.1006/jcph.1993.1097
     *
     * H. Erten, A. �ng�r, (2009): "Quality Triangulations
     * with Locally Optimal Steiner Points", SIAM J. Sci.
     * Comp., 31, pp. 2103-2130,
     * https://doi.org/10.1137/080716748
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_PRED_DELFRONT_3__
#   define __RDEL_PRED_DELFRONT_3__

    namespace mesh {

#   define __bias_bnds          // preference "bnd" faces

    template <
    typename G ,
    typename H ,
    typename M
             >
    class rdel_delfront_3d :
        public rdel_pred_base_3<G, M>
    {
    public  :

    /*-------------- "frontal" delaunay refinement in R^3 */

    typedef G                           geom_type ;
    typedef H                           hfun_type ;
    typedef M                           mesh_type ;

    typedef typename
            mesh_type::real_type        real_type ;
    typedef typename
            mesh_type::iptr_type        iptr_type ;

    typedef mesh::mesh_params  <
                real_type,
                iptr_type      >        rdel_opts ;

    typedef mesh::rdel_pred_base_3      <
                geom_type ,
                mesh_type      >        base_type ;

    typedef mesh::keep_minmax_length_3d <
                real_type,
                iptr_type      >    minmax_intersect ;
    typedef mesh::keep_minmax_cosine_3d <
                real_type,
                iptr_type      >    cosine_intersect ;

    /*------------------------ max. frontal "redo" steps */

    iptr_type static constexpr _REDO = 32 ;
    iptr_type static constexpr _STEP =  8 ;
    real_type static constexpr
        _OKAY = (real_type) +5./6. ;

    /*------------------------ refinement priority types */

    class edge_data
        {
        public  :
        iptr_type       _mark = +0 ;
        float           _cost ;
        } ;

    class face_data
        {
        public  :
        iptr_type       _mark = +0 ;
        float           _cost ;
        } ;

    class tria_data
        {
        public  :
        iptr_type       _mark = +0 ;
        float           _cost ;
        } ;

    /*------------------------ refinement priority pred. */

    __static_call
    __inline_call bool_type edge_pred (
        edge_data const&_idat ,
        edge_data const&_jdat
        )
    {
        if(_idat._mark == _jdat._mark )
        return _idat._cost >
               _jdat._cost ;
        else
        return _idat._mark <
               _jdat._mark ;
    }

    __static_call
    __inline_call bool_type face_pred (
        face_data const&_idat ,
        face_data const&_jdat
        )
    {
        if(_idat._mark == _jdat._mark )
        return _idat._cost >
               _jdat._cost ;
        else
        return _idat._mark <
               _jdat._mark ;
    }

    __static_call
    __inline_call bool_type tria_pred (
        tria_data const&_idat ,
        tria_data const&_jdat
        )
    {
        if(_idat._mark == _jdat._mark )
        return _idat._cost >
               _jdat._cost ;
        else
        return _idat._mark <
               _jdat._mark ;
    }

    /*------------------------ deterministic coord. hash */

    __static_call
    __inline_call uint32_t  hash_ball (
    __const_ptr(real_type) _ball
        )
    {
        uint32_t  _rsiz =
            sizeof(real_type) * +3 ;
        uint32_t  _usiz =
            sizeof(uint32_t ) * +1 ;

        uint32_t  _hash ;
        _hash = hash::hashword (
       (uint32_t*)_ball, _rsiz
                       / _usiz, +13)  ;

        return (  _hash ) ;
    }

    /*
    --------------------------------------------------------
     * FACE-COST: calc. face refinement "cost".
    --------------------------------------------------------
     */

    #include "rdel_cost_delfront_3.inc"


    /*
    --------------------------------------------------------
     * FACE-OFFH: form "size"-based off-centre.
    --------------------------------------------------------
     */

    #include "rdel_offh_delfront_3.inc"


    /*
    --------------------------------------------------------
     * FACE-SINK: form "sink"-based off-centre.
    --------------------------------------------------------
     */

    #include "rdel_sink_delfront_3.inc"


    /*
    --------------------------------------------------------
     * BASE-EDGE: TRUE if edge is "frontal".
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type base_edge (
        mesh_type &_mesh,
        iptr_type *_FNOD,
        iptr_type  _tadj,
        iptr_type  _fadj
        )
    {
        containers::array<iptr_type> _tset;

        _tset.set_alloc(8) ;

        iptr_type _epos;
        for (_epos = +3; _epos-- != +0; )
        {
    /*--------------------- is any face of FPOS "frontal" */
        iptr_type  _enod[ +3];
        mesh_type::tria_type::tria_type::
        face_node(_enod, _epos, 2, 1);
        _enod[0] = _FNOD[_enod[ 0]];
        _enod[1] = _FNOD[_enod[ 1]];

        _tset.set_count(0) ;

        base_type::edge_loop (_mesh ,
            _enod, _tadj ,
            _fadj, _tset ) ;

        algorithms::isort (
            &_enod[0], &_enod[2],
                std::less<iptr_type>());

        typename mesh_type::
                 edge_data _edat;
        _edat._node[0] = _enod[0] ;
        _edat._node[1] = _enod[1] ;

        typename mesh_type::
                 edge_list::
            item_type *_eptr = nullptr ;
        if(_mesh.find_edge(_edat,_eptr))
        {
            return (  true ) ;
        }

        iptr_type _fpos ;
        for (auto _tpos  = _tset.head();
                  _tpos != _tset.tend();
                ++_tpos )
        {
        for (_fpos = +4 ; _fpos-- != +0; )
        {
            iptr_type _fnod [ +4];
            mesh_type::tria_type::
                tria_type::
            face_node(_fnod, _fpos, 3, 2);
            _fnod[0] = _mesh._tria.
             tria(*_tpos)->node(_fnod[0]);
            _fnod[1] = _mesh._tria.
             tria(*_tpos)->node(_fnod[1]);
            _fnod[2] = _mesh._tria.
             tria(*_tpos)->node(_fnod[2]);

            iptr_type _same = +0 ;
            if (_fnod[0] == _enod[0] ||
                _fnod[0] == _enod[1] )
                _same += +1;
            if (_fnod[1] == _enod[0] ||
                _fnod[1] == _enod[1] )
                _same += +1;
            if (_fnod[2] == _enod[0] ||
                _fnod[2] == _enod[1] )
                _same += +1;
            if (_same != +2) continue ;

            if (_mesh._tria.
            node(_fnod[0])->fdim() > +2 ||
                _mesh._tria.
            node(_fnod[1])->fdim() > +2 ||
                _mesh._tria.
            node(_fnod[2])->fdim() > +2 )
                _same  = +0;
            if (_same != +2) continue ;

            algorithms::isort (
                &_fnod[0], &_fnod[3],
                    std::less<iptr_type> ());

            typename mesh_type::
                     face_data _fdat ;
            _fdat._node[0] = _fnod[0];
            _fdat._node[1] = _fnod[1];
            _fdat._node[2] = _fnod[2];

            typename mesh_type::
                     face_list::
                 item_type *_fptr = nullptr ;
            if (_mesh.find_face(_fdat,_fptr))
            {
                if (_fptr->_data._kind
                        ==  mesh::good_item )
                return ( true ) ;
            }
        }
        }

        }

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * BASE-FACE: TRUE if face is "frontal".
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type base_face (
        mesh_type &_mesh,
        iptr_type  _tpos,
        iptr_type/*_fpos*/
        )
    {
        iptr_type _fpos;
        for (_fpos = +4; _fpos-- != +0; )
        {
    /*--------------------- is any face of TPOS "frontal" */
        iptr_type _fnod[  +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_fnod, _fpos, +3, +2);
        _fnod[0] = _mesh._tria.
            tria(_tpos)->node(_fnod[0]);
        _fnod[1] = _mesh._tria.
            tria(_tpos)->node(_fnod[1]);
        _fnod[2] = _mesh._tria.
            tria(_tpos)->node(_fnod[2]);

        algorithms::isort (
            &_fnod[0], &_fnod[3],
                std::less<iptr_type>());

        typename mesh_type::
                 face_data _fdat ;
        _fdat._node[0] = _fnod[0];
        _fdat._node[1] = _fnod[1];
        _fdat._node[2] = _fnod[2];

        typename mesh_type::
                 face_list::
            item_type *_fptr = nullptr ;
        if(_mesh.find_face(_fdat,_fptr))
        {
            return (  true ) ;
        }

        iptr_type _tadj = +0 , _fadj = +0 ,
                  _tmrk = +0 ;
        _mesh. _tria.find_pair (
            _tpos, _tadj,
                _fpos, _fadj, _tmrk) ;

        if (_tadj != _mesh._tria.null_flag())
        {
            typename mesh_type::
                     tria_data _tdat ;
            _tdat._node[0] =
            _mesh._tria.tria(_tadj)->node(0);
            _tdat._node[1] =
            _mesh._tria.tria(_tadj)->node(1);
            _tdat._node[2] =
            _mesh._tria.tria(_tadj)->node(2);
            _tdat._node[3] =
            _mesh._tria.tria(_tadj)->node(3);

            typename mesh_type::
                     tria_list::
                item_type *_tptr = nullptr ;
            if(_mesh.find_tria(_tdat,_tptr))
            {
                if (_tptr->_data._kind
                        ==  mesh::good_item)
                return ( true ) ;
            }
        }
        else
        {
        return ( true  ) ;
        }

        }

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * EDGE-NODE: compute edge-based steiner point.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind edge_node (
        geom_type &_geom,
        hfun_type &_size,
        mesh_type &_mesh,
        edge_data &_edat,
        iptr_type  _tadj,
        iptr_type  _eadj,
        real_type *_ppos,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind
        _kind  = rdel_opts::null_kind ;

        __unreferenced(_edat);

    /*--------------------------- assemble local indexing */
        iptr_type _enod[ +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_enod, _eadj, 3, 1) ;
        _enod[0] =_mesh._tria.
         tria(_tadj)->node(_enod[ 0]) ;
        _enod[1] =_mesh._tria.
         tria(_tadj)->node(_enod[ 1]) ;

    /*--------------------------------- calc. circumballs */
        char_type _feat, _topo ;
        iptr_type _part;
        real_type _ebal[ +4] ;
        real_type _pmax[ +4] ;
        if (!base_type::edge_ball (
            _geom, _mesh, _tadj, _eadj,
            _ebal, _pmax,
            _feat, _topo, _part)  )
        {
    /*--------------------------------- is not restricted */
        __assert( false &&
            "EDGE-NODE: interior facet!") ;

        return ( _kind ) ;
        }

        if (_kind == rdel_opts::null_kind )
        {
    /*----------------------- attempt offcentre placement */
            if(__chkbit(_args.rule(),
                    rdel_opts
                  ::offH_kind) )
            _kind = edge_offh( _geom,
                _size , _mesh, _enod,
                _pmax , _ppos, _args) ;
        }
        if (_kind == rdel_opts::null_kind ||
            _kind == rdel_opts::circ_kind )
        {
            if (_edat._mark++ <= 0)
            {
        /*---------------------- reject as "void" element */
            _kind =  rdel_opts::null_kind ;
            }
            else
            {
        /*---------------------- take circumcentre direct */
            _ppos[ 0] = _pmax[ 0] ;
            _ppos[ 1] = _pmax[ 1] ;
            _ppos[ 2] = _pmax[ 2] ;

            _kind =  rdel_opts::circ_kind ;
            }
        }

    /*----------------------- report point-placement kind */
        return ( _kind ) ;
    }

    /*
    --------------------------------------------------------
     * FACE-NODE: compute face-based steiner point.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind face_node (
        geom_type &_geom,
        hfun_type &_size,
        mesh_type &_mesh,
        face_data &_fdat,
        iptr_type  _tadj,
        iptr_type  _fadj,
        real_type *_ppos,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind
        _kind  = rdel_opts::null_kind ;

    /*--------------------------- assemble local indexing */
        iptr_type _fnod[ +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_fnod, _fadj, 3, 2) ;
        _fnod[0] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 0]) ;
        _fnod[1] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 1]) ;
        _fnod[2] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 2]) ;

    /*--------------------------------- calc. circumballs */
        char_type _feat, _topo ;
        iptr_type _part;
        real_type _fbal[ +4] ;
        real_type _pmax[ +4] ;
        if (!base_type::face_ball (
            _geom, _mesh, _tadj, _fadj,
            _fbal, _pmax,
            _feat, _topo, _part)  )
        {
    /*--------------------------------- is not restricted */
        __assert( false &&
            "FACE-NODE: interior facet!") ;

        return ( _kind ) ;
        }

    /*--------------------------------- find edge lengths */
        real_type _llen[3] ;
        iptr_type _enum ;
        for(_enum = +3; _enum-- != +0; )
        {
            iptr_type _enod[ +3] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_enod, _enum, 2, 1) ;
            _enod[ 0] = _fnod[_enod[0]] ;
            _enod[ 1] = _fnod[_enod[1]] ;

            _llen[_enum] =
                geometry::lensqr_3d (
           &_mesh._tria.
             node(_enod[0])->pval(0),
           &_mesh._tria.
             node(_enod[1])->pval(0)) ;
        }

    /*--------------------------------- find min/max edge */
        iptr_type _emin = +0;
        iptr_type _emax = +0;
        for(_enum = +3; _enum-- != +1; )
        {
        if (_llen[_emax] < _llen[_enum])
            _emax = _enum ;
        if (_llen[_emin] > _llen[_enum])
            _emin = _enum ;
        }

    /*-------------------------- ask for "frontal" status */
        if(!base_edge(_mesh, _fnod, _tadj ,
                      _fadj) )
        {
            if (_fdat._mark <= _REDO)   // finite cycles!
            {
        /*---------------------- reject as "void" element */
            uint32_t _push =
            (hash_ball(_pmax) % _STEP) + 1;

            _fdat._mark += _push;

            return rdel_opts::null_kind;
            }
        /*----------------------------- just fall through */
        }

    /*--------------------------------- eval edge metrics */
        iptr_type  _enod[ +3];
        mesh_type::tria_type::tria_type::
        face_node(_enod, _emin, 2, 1) ;
        _enod[0] = _fnod[_enod[ 0]] ;
        _enod[1] = _fnod[_enod[ 1]] ;
        _enod[2] = _fnod[_enod[ 2]] ;

        real_type  _evec[4] = {
            _mesh.
        _tria.node(_enod[1])->pval(0) -
            _mesh.
        _tria.node(_enod[0])->pval(0) ,
            _mesh.
        _tria.node(_enod[1])->pval(1) -
            _mesh.
        _tria.node(_enod[0])->pval(1) ,
            _mesh.
        _tria.node(_enod[1])->pval(2) -
            _mesh.
        _tria.node(_enod[0])->pval(2)
            } ;

        real_type  _pmid[3] = {
            _mesh.
        _tria.node(_enod[0])->pval(0) +
            _mesh.
        _tria.node(_enod[1])->pval(0) ,
            _mesh.
        _tria.node(_enod[0])->pval(1) +
            _mesh.
        _tria.node(_enod[1])->pval(1) ,
            _mesh.
        _tria.node(_enod[0])->pval(2) +
            _mesh.
        _tria.node(_enod[1])->pval(2)
            } ;

        _pmid[ 0]*= (real_type) +.5 ;
        _pmid[ 1]*= (real_type) +.5 ;
        _pmid[ 2]*= (real_type) +.5 ;

        real_type _dvec[ +5] ;
        geometry::
        vector_3d( _pmid, _pmax, _dvec) ;

        _evec[ 3] =
        geometry:: normalise_3d (_evec) ;

        _dvec[ 3] =
        geometry:: normalise_3d (_dvec) ;

    /*-------------------------- off-centre, a la "ungor" */
        _dvec[4] =
            +std::numeric_limits
                <real_type>::infinity() ;

        real_type _rtri =
            _evec[3] * _args.off2() ;
        real_type _rfac =
            _evec[3] * (real_type)+.5 ;

        if (_rtri > _rfac)
        _dvec[ 4] = _rtri + std::sqrt (
            _rtri * _rtri -
            _rfac * _rfac);

        if (_kind == rdel_opts::null_kind )
        {
    /*----------------------- attempt offcentre placement */
            if(__chkbit(_args.rule(),
                    rdel_opts
                  ::offH_kind) )
            _kind = face_offh( _geom,
                _size , _mesh, _enod,
                _evec , _dvec, _pmax,
                _ppos , _args) ;
        }
        if (_kind == rdel_opts::null_kind ||
            _kind == rdel_opts::circ_kind )
        {
    /*----------------------- attempt sink-type placement */
            if(__chkbit(_args.rule(),
                    rdel_opts
                  ::sink_kind) )
            _kind = face_sink( _geom,
                _size , _mesh, _tadj,
                _enod , _pmax,
                _ppos , _args) ;
        }
        if (_kind == rdel_opts::null_kind ||
            _kind == rdel_opts::circ_kind )
        {
            real_type _alen =
                geometry::tria_quality_3d (
                &_mesh.
            _tria.node(_fnod[0])->pval(0) ,
                &_mesh.
            _tria.node(_fnod[1])->pval(0) ,
                &_mesh.
            _tria.node(_fnod[2])->pval(0)
                ) ;

            if (_alen <= _OKAY &&
                    _fdat._mark <= _REDO)
            {
        /*---------------------- reject as "void" element */
            uint32_t _push =
            (hash_ball(_pmax) % _STEP) + 1;

            _fdat._mark += _push ;

            _kind =  rdel_opts::null_kind ;
            }
            else
            {
        /*---------------------- take circumcentre direct */
            _ppos[ 0] = _pmax[ 0] ;
            _ppos[ 1] = _pmax[ 1] ;
            _ppos[ 2] = _pmax[ 2] ;

            _kind =  rdel_opts::circ_kind ;
            }
        }

    /*----------------------- report point-placement kind */
        return ( _kind ) ;
    }

    /*
    --------------------------------------------------------
     * TRIA-NODE: compute tria-based steiner point.
    --------------------------------------------------------
     */

    __static_call
    __normal_call
    typename rdel_opts::node_kind tria_node (
        geom_type &_geom,
        hfun_type &_size,
        mesh_type &_mesh,
        tria_data &_tdat,
        iptr_type  _tpos,
        real_type *_ppos,
        rdel_opts &_args
        )
    {
        typename rdel_opts::node_kind
        _kind  = rdel_opts::null_kind ;

        __unreferenced(_geom) ;
        __unreferenced(_size) ;
        __unreferenced(_args) ;

    /*--------------------------------- pop tria indexing */
        iptr_type _tnod[4] =  {
            _mesh.
        _tria.tria(_tpos)->node( +0),
            _mesh.
        _tria.tria(_tpos)->node( +1),
            _mesh.
        _tria.tria(_tpos)->node( +2),
            _mesh.
        _tria.tria(_tpos)->node( +3)
            } ;

    /*--------------------------------- calc. 3-tria ball */
        real_type  _tbal[4];
        _tbal[0] = _mesh.
            _tria.tria(_tpos)->circ(0);
        _tbal[1] = _mesh.
            _tria.tria(_tpos)->circ(1);
        _tbal[2] = _mesh.
            _tria.tria(_tpos)->circ(2);

        _tbal[3] = (real_type)+0. ;
        _tbal[3]+=
        geometry::lensqr_3d (_tbal,
           &_mesh._tria.node(
                _tnod[0])->pval( 0)) ;
        _tbal[3]+=
        geometry::lensqr_3d (_tbal,
           &_mesh._tria.node(
                _tnod[1])->pval( 0)) ;
        _tbal[3]+=
        geometry::lensqr_3d (_tbal,
           &_mesh._tria.node(
                _tnod[2])->pval( 0)) ;
        _tbal[3]+=
        geometry::lensqr_3d (_tbal,
           &_mesh._tria.node(
                _tnod[3])->pval( 0)) ;

        _tbal[3]/= (real_type)+4. ;

    /*--------------------------------- find edge lengths */
        real_type _llen[6] ;
        iptr_type _enum ;
        for(_enum = +6; _enum-- != +0; )
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

    /*--------------------------------- find min/max edge */
        iptr_type _emin = +0;
        iptr_type _emax = +0;
        for(_enum = +6; _enum-- != +1; )
        {
        if (_llen[_emax] < _llen[_enum])
            _emax = _enum ;
        if (_llen[_emin] > _llen[_enum])
            _emin = _enum ;
        }

    /*--------------------------------- find 2-face radii */
        real_type _frad[4] ;
        iptr_type _fpos ;
        for(_fpos = +4; _fpos-- != +0; )
        {
            iptr_type _fnod [ +4];
            mesh_type::tria_type::
                tria_type::
            face_node(_fnod, _fpos, 3, 2);
            _fnod[ +0] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 0]);
            _fnod[ +1] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 1]);
            _fnod[ +2] = _mesh._tria.
             tria(_tpos)->node(_fnod[ 2]);

            real_type _fbal [ +4];
            geometry::circ_ball_3d(_fbal ,
           &_mesh._tria.
             node(_fnod[0])->pval(0),
           &_mesh._tria.
             node(_fnod[1])->pval(0),
           &_mesh._tria.
             node(_fnod[2])->pval(0)) ;

            _frad[_fpos] =  _fbal[3];
        }

    /*--------------------------------- find min/max face */
        iptr_type _fmin = +0;
        iptr_type _fmax = +0;
        for(_fpos = +4; _fpos-- != +1; )
        {
        if (_frad[_fmax] < _frad[_fpos])
            _fmax = _fpos ;
        if (_frad[_fmin] > _frad[_fpos])
            _fmin = _fpos ;
        }

    /*-------------------------- ask for "frontal" status */
        if(!base_face(_mesh, _tpos, _fmin))
        {
            if (_tdat._mark <= _REDO)   // finite cycles!
            {
        /*---------------------- reject as "void" element */
            uint32_t _push =
            (hash_ball(_tbal) % _STEP) + 1;

            _tdat._mark += _push;

            return rdel_opts::null_kind;
            }
        /*----------------------------- just fall through */
        }

    /*--------------------------------- get face indexing */
        iptr_type _enod [ +4];
        mesh_type::tria_type::tria_type::
        face_node(_enod, _emin, 3, 1) ;
        _enod[0] = _mesh._tria.
         tria(_tpos)->node(_enod[+0]) ;
        _enod[1] = _mesh._tria.
         tria(_tpos)->node(_enod[+1]) ;

        iptr_type _fnod [ +4];
        mesh_type::tria_type::tria_type::
        face_node(_fnod, _fmin, 3, 2) ;
        _fnod[0] = _mesh._tria.
         tria(_tpos)->node(_fnod[+0]) ;
        _fnod[1] = _mesh._tria.
         tria(_tpos)->node(_fnod[+1]) ;
        _fnod[2] = _mesh._tria.
         tria(_tpos)->node(_fnod[+2]) ;

    /*----------------------------------- calc. edge-ball */
        real_type _elen ;
        _elen = geometry::lensqr_3d (
           &_mesh.
        _tria.node(_enod[0])->pval(0) ,
           &_mesh.
        _tria.node(_enod[1])->pval(0)
            ) ;

    /*----------------------------------- calc. face-ball */
        real_type _fbal [ +4];
        geometry::circ_ball_3d(_fbal,
           &_mesh.
        _tria.node(_fnod[0])->pval(0) ,
           &_mesh.
        _tria.node(_fnod[1])->pval(0) ,
           &_mesh.
        _tria.node(_fnod[2])->pval(0)
            ) ;

    /*------------------------- assemble "frontal" vector */
        real_type  _dvec[5] ;
        geometry::
        vector_3d( _fbal, _tbal, _dvec) ;

        _dvec[3] =
        geometry:: normalise_3d (_dvec) ;

    /*------------------------- off-centre - a'la "ungor" */
        _dvec[4] =
            std::numeric_limits
                <real_type>::infinity() ;

        real_type _rfac =
            std::sqrt(_fbal[ +3]) ;
        real_type _rtri =
        std::sqrt(_elen) * _args.off3() ;

        if (_rtri > _rfac)
        _dvec[4] =  _rtri + std::sqrt (
            _rtri * _rtri -
            _rfac * _rfac);

        if (_kind == rdel_opts::null_kind )
        {
    /*----------------------- attempt offcentre placement */
            if(__chkbit(_args.rule(),
                    rdel_opts
                  ::offH_kind) )
            _kind = tria_offh( _geom,
                _size , _mesh, _fnod,
                _fbal , _tbal, _dvec,
                _ppos , _args) ;
        }
        if (_kind == rdel_opts::null_kind ||
            _kind == rdel_opts::circ_kind )
        {
    /*----------------------- attempt sink-type placement */
            if(__chkbit(_args.rule(),
                    rdel_opts
                  ::sink_kind) )
            _kind = tria_sink( _geom,
                _size , _mesh,
                _tpos , _tnod, _tbal,
                _ppos , _args) ;
        }
        if (_kind == rdel_opts::null_kind ||
            _kind == rdel_opts::circ_kind )
        {
            real_type _vlen =
                geometry::tria_quality_3d (
                &_mesh.
            _tria.node(_tnod[0])->pval(0) ,
                &_mesh.
            _tria.node(_tnod[1])->pval(0) ,
                &_mesh.
            _tria.node(_tnod[2])->pval(0) ,
                &_mesh.
            _tria.node(_tnod[3])->pval(0)
                ) ;

            if (_vlen <= _OKAY &&
                    _tdat._mark <= _REDO)
            {
        /*---------------------- reject as "void" element */
            uint32_t _push =
            (hash_ball(_tbal) % _STEP) + 1;

            _tdat._mark += _push ;

            _kind =  rdel_opts::null_kind ;
            }
            else
            {
        /*---------------------- take circumcentre direct */
            _ppos[ 0] = _tbal[ 0] ;
            _ppos[ 1] = _tbal[ 1] ;
            _ppos[ 2] = _tbal[ 2] ;

            _kind =  rdel_opts::circ_kind ;
            }
        }

    /*----------------------- report point-placement kind */
        return ( _kind ) ;
    }

    } ;

#   undef  __bias_bnds


    }

#   endif   //__RDEL_PRED_DELFRONT_3__



