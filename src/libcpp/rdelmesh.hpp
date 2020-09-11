
/*
------------------------------------------------------------
 * restricted delaunay mesh generation in R^d.
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
 * Last updated: 20 February, 2019
 *
 * Copyright 2013-2019
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __RDELMESH__
#   define __RDELMESH__

    namespace mesh
    {
/*-------------------------- classification of mesh cells */
    char_type constexpr null_item = 0 ;
    char_type constexpr void_item = 1 ;
    char_type constexpr ring_item = 2 ;
    char_type constexpr good_item = 3 ;

    char_type constexpr null_ball = 0 ;
    char_type constexpr feat_ball = 1 ;
    }

#   include "containers.hpp"
#   include "algorithms.hpp"

#   include "geometry.hpp"
#   include "hashfunc.hpp"

#   include "tessellate.hpp"

#   include "geomtype.hpp"
#   include "meshfunc.hpp"
#   include "meshtype.hpp"

#   include "rdel_mesh/rdel_complex_2.hpp"
#   include "rdel_mesh/rdel_complex_3.hpp"

#   include "rdel_mesh/rdel_timers.hpp"

#   include "rdel_mesh/rdel_filt_k.hpp"

#   include "rdel_mesh/rdel_base_2.hpp"
#   include "rdel_mesh/rdel_base_3.hpp"

#   include "rdel_mesh/rdel_pred_delaunay_2.hpp"
#   include "rdel_mesh/rdel_pred_delaunay_3.hpp"

#   include "rdel_mesh/rdel_pred_delfront_2.hpp"
#   include "rdel_mesh/rdel_pred_delfront_3.hpp"

#   include "rdel_mesh/rdel_mesh_2.hpp"
#   include "rdel_mesh/rdel_mesh_3.hpp"

#   include "rdel_mesh/rdel_make_2.hpp"
#   include "rdel_mesh/rdel_make_3.hpp"

//  include "rdel_mesh/rvor_mesh_2.hpp"
//  include "rdel_mesh/rvor_mesh_3.hpp"

#   endif//__RDELMESH__



