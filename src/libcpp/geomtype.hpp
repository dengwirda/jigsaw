
/* 
------------------------------------------------------------
 * data-structures for geometry rep. in R^d. 
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
 * Last updated: 13 August, 2018
 *
 * Copyright 2013-2018
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __GEOMTYPE__
#   define __GEOMTYPE__
 
    namespace mesh
    {        
/*-------------------------- classification of mesh feat. */   
    char_type null_feat = 0;
    char_type user_feat = 1;
    char_type soft_feat = 2;
    char_type hard_feat = 3; 
    }
 
#   include "containers.hpp"

#   include "geometry.hpp"

#   include "aabbtree.hpp"
#   include "meshtype.hpp"

#   include "geom_type/geom_mesh_euclidean_2.hpp"
#   include "geom_type/geom_mesh_euclidean_3.hpp"

#   include "geom_type/geom_mesh_ellipsoid_3.hpp"
    
#   endif//__GEOMTYPE__



