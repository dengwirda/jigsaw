
    /*
    --------------------------------------------------------
     * RDEL-PRED-DELAUNAY-3: "classical" kernel in R^3.
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
     * Last updated: 28 Dec., 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     *
     * This class defines the behaviour of the "standard"
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
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_PRED_DELAUNAY_3__
#   define __RDEL_PRED_DELAUNAY_3__

    namespace mesh {

    template <
    typename G ,
    typename H ,
    typename M
             >
    class rdel_delaunay_3d :
        public rdel_pred_base_3<G, M>
    {
    public  :

    /*------------ "classical" delaunay refinement in R^3 */

    typedef G                           geom_type ;
    typedef H                           hfun_type ;
    typedef M                           mesh_type ;

    typedef typename
            mesh_type::real_type        real_type ;
    typedef typename
            mesh_type::iptr_type        iptr_type ;

    typedef mesh::mesh_params  <
                real_type ,
                iptr_type      >        rdel_opts ;

    typedef mesh::rdel_pred_base_3      <
                geom_type ,
                mesh_type      >        base_type ;

    typedef mesh::keep_minmax_length_3d <
                real_type ,
                iptr_type      >    minmax_intersect ;

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
        return _idat._cost >
               _jdat._cost ;
    }

    __static_call
    __inline_call bool_type face_pred (
        face_data const&_idat ,
        face_data const&_jdat
        )
    {
        return _idat._cost >
               _jdat._cost ;
    }

    __static_call
    __inline_call bool_type tria_pred (
        tria_data const&_idat ,
        tria_data const&_jdat
        )
    {
        return _idat._cost >
               _jdat._cost ;
    }

    /*
    --------------------------------------------------------
     * FACE-COST: calc. face refinement "cost".
    --------------------------------------------------------
     */

    #include "rdel_cost_delaunay_3.inc"


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
        __unreferenced(_size);
        __unreferenced(_edat);
        __unreferenced(_args);

    /*--------------------------------- calc. circumballs */
        char_type _feat;        
        char_type _topo[ +2] ;
        iptr_type _part;
        real_type _ebal[ +4] ;
        real_type _pmax[ +4] ;
        if (!base_type::edge_ball (
            _geom, _mesh, _tadj,
            _eadj, _ebal, _pmax,
            _feat, _topo, _part)  )
        {
    /*--------------------------------- is not restricted */
        __assert( false &&
        "EDGE-NODE: interior facet!") ;

        return rdel_opts:: null_kind;
        }

    /*--------------------------------- push circumcentre */
        _ppos[ 0] = _pmax[ 0] ;
        _ppos[ 1] = _pmax[ 1] ;
        _ppos[ 2] = _pmax[ 2] ;

        return rdel_opts:: circ_kind;
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
        __unreferenced(_size);
        __unreferenced(_fdat);
        __unreferenced(_args);

    /*--------------------------------- calc. circumballs */
        char_type _feat;        
        char_type _topo[ +2] ;
        iptr_type _part;
        real_type _fbal[ +4] ;
        real_type _pmax[ +4] ;
        if (!base_type::face_ball (
            _geom, _mesh, _tadj,
            _fadj, _fbal, _pmax,
            _feat, _topo, _part)  )
        {
    /*--------------------------------- is not restricted */
        __assert( false &&
        "FACE-NODE: interior facet!") ;

        return rdel_opts:: null_kind;
        }

    /*--------------------------------- push circumcentre */
        _ppos[ 0] = _pmax[ 0] ;
        _ppos[ 1] = _pmax[ 1] ;
        _ppos[ 2] = _pmax[ 2] ;

        return rdel_opts:: circ_kind;
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
        __unreferenced(_geom) ;
        __unreferenced(_size) ;
        __unreferenced(_tdat) ;
        __unreferenced(_args) ;

    /*--------------------------------- tria circumcentre */
        real_type  _tbal[4];
        _tbal[0] = _mesh.
        _tria.tria(_tpos)->circ(+0);
        _tbal[1] = _mesh.
        _tria.tria(_tpos)->circ(+1);
        _tbal[2] = _mesh.
        _tria.tria(_tpos)->circ(+2);

    /*--------------------------------- push circumcentre */
        _ppos[ 0] = _tbal[ 0] ;
        _ppos[ 1] = _tbal[ 1] ;
        _ppos[ 2] = _tbal[ 2] ;

        return rdel_opts:: circ_kind ;
    }

    } ;


    }

#   endif   //__RDEL_PRED_DELAUNAY_3__



