
    /*
    --------------------------------------------------------
     * MSH-CULL: cull MESH data via BFS from SEEDS.
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
     * Last updated: 08 June, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_CULL__
#   define __MSH_CULL__

    /*
    --------------------------------------------------------
     * CULL-MESH: cull mesh-complex via BFS from SEEDS.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type cull_mesh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom ,
        mesh_data &_mesh
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;
        __unreferenced(_jlog) ;

        if (_mesh._ndim == +2 &&
            _geom._kind ==
            jmsh_kind::euclidean_mesh)
        {
    /*---------------------------------- cull 2-dim. mesh */

          //find_near(_mesh, _aabb, +2, _geom._seed, _near);


          //build aabb for mesh
          //find nearest cell for each seed (find_near)
          //find connected cells for nearest to each seed (find_conn)



        }
        else
        if (_mesh._ndim == +3 &&
            _geom._kind ==
            jmsh_kind::euclidean_mesh)
        {
    /*---------------------------------- cull 3-dim. mesh */


        }
        else
        if (_mesh._ndim == +3 &&
            _geom._kind ==
            jmsh_kind::ellipsoid_mesh)
        {
    /*---------------------------------- cull 2-sph. mesh */


        }

        return ( _errv ) ;
    }


#   endif   //__MSH_CULL__



