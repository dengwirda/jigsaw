
/*
------------------------------------------------------------
 * data-structures for discrete complexes in R^d.
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
 * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant 
 * this code in any way whatsoever.  This code is 
 * provided "as-is" to be used at your own risk.
 *
 * THE CONTRIBUTORS include:
 * (a) The University of Sydney
 * (b) The Massachusetts Institute of Technology
 * (c) Columbia University
 * (d) The National Aeronautics & Space Administration
 * (e) Los Alamos National Laboratory
 *
------------------------------------------------------------
 *
 * Last updated: 01 Sept., 2021
 *
 * Copyright 2013-2021
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __MESHREPS__
#   define __MESHREPS__

    namespace mesh
    {
/*-------------------------- classification of mesh cells */
    char_type constexpr EMPTY_tag = +0 ;
    char_type constexpr
        POINT_tag = JIGSAW_POINT_TAG ;
    char_type constexpr
        EDGE2_tag = JIGSAW_EDGE2_TAG ;
    char_type constexpr
        TRIA3_tag = JIGSAW_TRIA3_TAG ;
    char_type constexpr
        QUAD4_tag = JIGSAW_QUAD4_TAG ;
    char_type constexpr
        TRIA4_tag = JIGSAW_TRIA4_TAG ;
    char_type constexpr
        HEXA8_tag = JIGSAW_HEXA8_TAG ;
    char_type constexpr
        WEDG6_tag = JIGSAW_WEDG6_TAG ;
    char_type constexpr
        PYRA5_tag = JIGSAW_PYRA5_TAG ;
    }

#   include "containers.hpp"

#   include "hashfunc.hpp"

#   include "mesh_reps/mesh_complex_type_k.hpp"

#   include "mesh_reps/hash_complex_k.hpp"
#   include "mesh_reps/base_complex_k.hpp"

#   include "mesh_reps/tria_complex_1.hpp"
#   include "mesh_reps/tria_complex_2.hpp"
#   include "mesh_reps/tria_complex_3.hpp"

#   include "mesh_reps/mesh_complex_1.hpp"
#   include "mesh_reps/mesh_complex_2.hpp"
#   include "mesh_reps/mesh_complex_3.hpp"

#   endif//__MESHREPS__



