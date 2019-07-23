
    /*
    --------------------------------------------------------
     * GEOM-BASE-3: base class for geom. in R^3.
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
     * Last updated: 16 January, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __GEOM_BASE_3__
#   define __GEOM_BASE_3__

    namespace mesh {

    template <
    typename R ,
    typename I 
             >
    class geom_base_3d
    {
    public  :
    
    /*------------------------- base geometry type in R^3 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;

    /*------------------------- types for intersect calls */

    struct line_type
        {
        real_type  _ipos[3] ;   // line endpoints
        real_type  _jpos[3] ;
        } ;

    struct flat_type
        {
        real_type  _ppos[3] ;   // point on flat
        real_type  _nvec[3] ;   // norm. to flat
        real_type  _rmin[3] ;   // bounding-aabb
        real_type  _rmax[3] ;
        } ;

    struct disc_type
        {
        real_type  _pmid[3] ;   // point on disc
        real_type  _nvec[3] ;   // norm. to disc
        real_type  _rrad ;      // radius
        } ;

    struct ball_type
        {
        real_type  _pmid[3] ;   // centre of ball
        real_type  _rrad ;      // radius
        } ;

    /*
    --------------------------------------------------------
     * INIT-GEOM: init. geometry data structures.
    --------------------------------------------------------
     */

    template <
        typename  geom_opts
             >
    __normal_call void_type init_geom (
        geom_opts &_opts
        )
    {
        __unreferenced(_opts) ;
    }

    /*
    --------------------------------------------------------
     * HAVE-FEAT: TRUE if has k-dim. FEAT.
    --------------------------------------------------------
     */

    __inline_call bool_type have_feat (
        iptr_type _fdim
        )
    {
        __unreferenced(_fdim) ;

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * SEED-MESH: init. "seed" vertex set on geom.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  geom_opts
             >
    __normal_call void_type seed_mesh (
        mesh_type &_mesh,
        geom_opts &_opts
        )
    {
        __unreferenced(_mesh) ;
        __unreferenced(_opts) ;
    }

    /*
    --------------------------------------------------------
     * SEED-FEAT: init. "seed" vertex set on geom.
    --------------------------------------------------------
     */

    template <
        typename  mesh_type ,
        typename  geom_opts
             >
    __normal_call void_type seed_feat (
        mesh_type &_mesh,
        geom_opts &_opts
        )
    {
        __unreferenced(_mesh) ;
        __unreferenced(_opts) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: find LINE/2-GEOM. intersections.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        line_type &_line ,
        hits_func &_hfun
        )
    {   
        __unreferenced(_line) ;
        __unreferenced(_hfun) ;

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: find FLAT/1-GEOM. intersections.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        flat_type &_flat ,
        hits_func &_hfun
        )
    {   
        __unreferenced(_flat) ;
        __unreferenced(_hfun) ;

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: find DISC/2-GEOM. intersections.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        disc_type &_disc ,
        real_type *_sbal ,
        hits_func &_hfun
        )
    {   
        __unreferenced(_disc) ;
        __unreferenced(_sbal) ;
        __unreferenced(_hfun) ;

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * INTERSECT: find BALL/1-GEOM. intersections.
    --------------------------------------------------------
     */

    template <
        typename  hits_func
             >
    __normal_call bool_type intersect (
        ball_type &_ball ,
        hits_func &_hfun
        )
    {   
        __unreferenced(_ball) ;
        __unreferenced(_hfun) ;

        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * IS-INSIDE: find NODE/3-GEOM. intersections.
    --------------------------------------------------------
     */

    __normal_call iptr_type is_inside (
        real_type *_ppos
        )
    {
        __unreferenced(_ppos) ;

        return (iptr_type) -1 ;
    }


    } ;


    }

#   endif//__GEOM_BASE_3__


