
    /*
    --------------------------------------------------------
     * RDEL-PRED-BASE-2: base-class for refinement in R^2. 
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
     * Last updated: 10 December, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_PRED_BASE_2__
#   define __RDEL_PRED_BASE_2__

    namespace mesh {
    
    template <
    typename G ,
    typename M
             >
    class rdel_pred_base_2
    {
    public  :
    
    /*------------ base-class for refinement kern. in R^2 */
    
    typedef G                           geom_type ;
    typedef M                           mesh_type ;

    typedef typename 
            mesh_type::real_type        real_type ;
    typedef typename 
            mesh_type::iptr_type        iptr_type ;

    /*
    --------------------------------------------------------
     * EDGE-BALL: calc. edge-based circumballs.
    --------------------------------------------------------
     */
     
    __static_call 
    __normal_call bool_type edge_ball (
        geom_type &_geom,
        mesh_type &_mesh,
        iptr_type  _tadj,
        iptr_type  _eadj,
        real_type *_ebal,
        real_type *_sbal,
        char_type &_feat,
        char_type &_topo,
        iptr_type &_part
        )
    {
    /*--------------------------- assemble local indexing */
        iptr_type _enod[ +3] ;
        mesh_type::tria_type::tria_type::
        face_node(_enod, _eadj, 2, 1);
        _enod[0] =_mesh._tria.
         tria(_tadj)->node(_enod[ 0]);
        _enod[1] =_mesh._tria.
         tria(_tadj)->node(_enod[ 1]);
        
    /*--------------------------- init. output balls = 0. */
        _ebal[0] = (real_type) +0. ;
        _ebal[1] = (real_type) +0. ;
        _ebal[2] = (real_type) +0. ;
        
        _sbal[0] = (real_type) +0. ;
        _sbal[1] = (real_type) +0. ;
        _sbal[2] = (real_type) +0. ;
        
    /*--------------------------- get local neighbourhood */
        iptr_type _topp = +0 ;
        iptr_type _eopp = +0 ;
        iptr_type _tmrk = +0 ;
        _mesh. _tria.find_pair (
        _tadj, _topp,
        _eadj, _eopp, _tmrk) ;

    /*--------------------------- skip faces adj. to hull */
        if (_topp == 
            _mesh._tria.null_flag())
            return ( false ) ;

        if (_topp < _tadj)
            std::swap(_tadj, _topp);

#       ifdef __testdual
    /*--------------------------- test dual-edge validity */
        real_type _nvec[ +3] ;
            geometry::vector_2d (
       &_mesh._tria.
         node(_enod[ 0])->pval(0) ,
       &_mesh._tria.
         node(_enod[ 1])->pval(0) ,
        _nvec ) ;
        
        real_type _cvec[2] = {
        _mesh._tria.
            tria( _tadj)->circ(0) - 
        _mesh._tria.
            tria( _topp)->circ(0) ,
        _mesh._tria.
            tria( _tadj)->circ(1) - 
        _mesh._tria.
            tria( _topp)->circ(1)
            } ;

        _nvec[2] =  _nvec[1] ;
        _nvec[1] = +_nvec[0] ;
        _nvec[0] = -_nvec[2] ;

        real_type _vdot = 
        geometry::dot_2d(_nvec, _cvec) ;
        
        if (_vdot < (real_type)0.)
        {
            std::cout << 
        "EDGE-BALL: bad dual-edge! \n" ;
        }
#       endif

    /*------------------------------ bipolar voronoi edge */
        real_type  _ibal[ 2] = {
        _mesh.
        _tria.tria(_tadj)->circ( 0),
        _mesh.
        _tria.tria(_tadj)->circ( 1)
            } ;
        
        real_type  _jbal[ 2] = {
        _mesh.
        _tria.tria(_topp)->circ( 0),
        _mesh.
        _tria.tria(_topp)->circ( 1)
            } ;
 
    /*------------------------------ calc. diametric ball */
        geometry::circ_ball_2d(_ebal , 
       &_mesh._tria.
            node(_enod[0])->pval(0),
       &_mesh._tria.
            node(_enod[1])->pval(0)) ;

    /*------------------------- find surface intersection */
        mesh::keep_minmax_length_2d<
            real_type, 
            iptr_type> _pred(_ebal ) ;
        if (!_geom.intersect (
               &_ibal[+0] ,
               &_jbal[+0] ,  _pred))
    /*------------------------- is not a restricted facet */
        return ( _pred._find ) ;
   
        _sbal[ 0] 
            = _pred._pmax.pval(0);
        _sbal[ 1] 
            = _pred._pmax.pval(1);
   
        _feat=_pred._pmax.feat ();
        _topo=_pred._pmax.topo ();
        _part=_pred._pmax.itag ();
   
    /*------------------------- eval. surface ball radius */
        _sbal[ 2]+= 
        geometry::lensqr_2d (_sbal, 
       &_mesh._tria.
            node(_enod[0])->pval(0)) ;
        _sbal[ 2]+= 
        geometry::lensqr_2d (_sbal, 
       &_mesh._tria.
            node(_enod[1])->pval(0)) ;
        
        _sbal[ 2]/= (real_type)+2.;

    /*----------------------- return TRUE if "restricted" */
        return ( _pred._find ) ;
    }

    } ;

    
    }
    
#   endif   //__RDEL_PRED_BASE_2__



