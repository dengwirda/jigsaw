
    /*
    --------------------------------------------------------
     * RDEL-PRED-DELAUNAY-2: "classical" kernel in R^2. 
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

#   pragma once

#   ifndef __RDEL_PRED_DELAUNAY_2__
#   define __RDEL_PRED_DELAUNAY_2__

    namespace mesh {

    template <
    typename G ,
    typename H ,
    typename M
             >
    class rdel_delaunay_2d :
        public rdel_pred_base_2<G, M>
    {
    public  :
    
    /*------------ "classical" delaunay refinement in R^2 */
    
    typedef G						    geom_type ;
    typedef H						    hfun_type ;
    typedef M						    mesh_type ;

    typedef typename 
            mesh_type::real_type        real_type ;
    typedef typename 
            mesh_type::iptr_type        iptr_type ;

    typedef mesh::rdel_params  <
                real_type , 
                iptr_type      >        rdel_opts ;

    typedef mesh::rdel_pred_base_2      <
                geom_type , 
                mesh_type      >        base_type ;

    typedef mesh::keep_minmax_length_2d <
                real_type , 
                iptr_type      >    minmax_intersect ;
    
    /*------------------------ refinement priority types */

    class edge_data
        {
        public  :
        real_type               _cost ;
        } ;

    class tria_data
        {
        public  :
        real_type               _cost ;
        } ;

    /*------------------------ refinement priority pred. */
    
    __static_call
    __inline_call bool_type edge_pred (
        edge_data const&_idat ,
        edge_data const&_jdat
        )
    {   return ( _idat. _cost > 
                 _jdat. _cost ) ;
    }
    
    __static_call
    __inline_call bool_type tria_pred (
        tria_data const&_idat ,
        tria_data const&_jdat
        )
    {   return ( _idat. _cost > 
                 _jdat. _cost ) ;
    }
    
    /*
    --------------------------------------------------------
     * FACE-COST: calc. face refinement "cost".
    --------------------------------------------------------
     */
     
    #include "rdel_cost_delaunay_2.inc"
    

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
        char_type _feat, _topo ;
        iptr_type _part;
        real_type _ebal[ +3];
        real_type _pmax[ +3];
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
        real_type  _tbal[3];
        _tbal[0] = _mesh.
        _tria.tria(_tpos)->circ(+0);
        _tbal[1] = _mesh.
        _tria.tria(_tpos)->circ(+1);

    /*--------------------------------- push circumcentre */
        _ppos[ 0] = _tbal[ 0] ;
        _ppos[ 1] = _tbal[ 1] ;

        return rdel_opts:: circ_kind ;
    }

    } ;
    
    
    }
    
#   endif   //__RDEL_PRED_DELAUNAY_2__



