
/* 
------------------------------------------------------------
 * data-structures/algorithms for mesh generation. 
------------------------------------------------------------
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
------------------------------------------------------------
 *
 * Last updated: 10 September, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __GEOMMESH__
#   define __GEOMMESH__
 
#   include "containers.hpp"
#   include "algorithms.hpp"

    namespace mesh
    {
/*-------------------------- classification of mesh items */
    char_type null_item = 0;
    char_type void_item = 1;
    char_type ring_item = 2;
    char_type good_item = 3;
        
/*-------------------------- classification of mesh feat. */   
    char_type null_feat = 0;
    char_type user_feat = 1;
    char_type soft_feat = 2;
    char_type hard_feat = 3; 
    }

#   include  "hashfunc.hpp"
#   include  "geometry.hpp"

/*------------------------------------------ AABB indexes */

#   include "geom_mesh/aabb_type_k.hpp"
#   include "geom_mesh/aabb_tree_k.hpp"

#   include "geom_mesh/aabb_pred_k.hpp"

#   include "geom_mesh/aabb_mesh_k.hpp"

/*------------------------------------------ tria-complex */

#   include "geom_mesh/tria_complex_type_k.hpp"

#   include "geom_mesh/tria_complex_1.hpp"
#   include "geom_mesh/tria_complex_2.hpp"
#   include "geom_mesh/tria_complex_3.hpp"

/*------------------------------------------ delaunay-tri */

#   include "geom_mesh/delaunay_tri_type_2.hpp"
#   include "geom_mesh/delaunay_tri_type_3.hpp"

#   include "geom_mesh/delaunay_tri_k.hpp"

#   include "geom_mesh/del_tri_euclidean_2.hpp"
#   include "geom_mesh/del_tri_euclidean_3.hpp"
    
/*-------------------------------------- GEOM kern. types */

#   include "geom_mesh/geom_mesh_euclidean_2.hpp"
#   include "geom_mesh/geom_mesh_euclidean_3.hpp"

#   include "geom_mesh/geom_mesh_ellipsoid_3.hpp"

/*-------------------------------------- HFUN kern. types */

#   include "geom_mesh/hfun_base_k.hpp"

#   include "geom_mesh/hfun_constant_value_k.hpp"

#   include "geom_mesh/hfun_mesh_euclidean_2.hpp"
#   include "geom_mesh/hfun_mesh_euclidean_3.hpp"

#   include "geom_mesh/hfun_grid_euclidean_2.hpp"
#   include "geom_mesh/hfun_grid_euclidean_3.hpp"

#   include "geom_mesh/hfun_grid_ellipsoid_3.hpp"

/*-------------------------------------- RDEL tria. types */

#   include "geom_mesh/rdel_complex_2.hpp"
#   include "geom_mesh/rdel_complex_3.hpp"

/*-------------------------- restricted delaunay mesh-gen */

#   include "geom_mesh/rdel_params.hpp"

#   include "geom_mesh/rdel_filt_k.hpp"

#   include "geom_mesh/rdel_base_2.hpp"
#   include "geom_mesh/rdel_base_3.hpp"

#   include "geom_mesh/rdel_pred_delaunay_2.hpp"
#   include "geom_mesh/rdel_pred_delaunay_3.hpp"

#   include "geom_mesh/rdel_pred_delfront_2.hpp"
#   include "geom_mesh/rdel_pred_delfront_3.hpp"

#   include "geom_mesh/rdel_mesh_2.hpp"
#   include "geom_mesh/rdel_mesh_3.hpp"

/*-------------------------- hill-climb mesh optimisation */

#   include "geom_mesh/iter_params.hpp"

#   include "geom_mesh/iter_mesh_euclidean_2.hpp"
#   include "geom_mesh/iter_pred_euclidean_2.hpp"

#   include "geom_mesh/iter_mesh_euclidean_3.hpp"
//  include "geom_mesh/iter_pred_euclidean_3.hpp"

//  include "geom_mesh/iter_mesh_ellipsoid_3.hpp"
#   include "geom_mesh/iter_pred_ellipsoid_3.hpp"

#   include "geom_mesh/iter_mesh_2.hpp"
//  include "geom_mesh/iter_mesh_3.hpp"


#   endif   //  __GEOMMESH__


