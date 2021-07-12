
/*
------------------------------------------------------------
 * algorithms for delaunay tessellation in R^d.
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
 * Last updated: 21 Apr., 2021
 *
 * Copyright 2013-2021
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __TESSELLATE__
#   define __TESSELLATE__

#   include "containers.hpp"
#   include "algorithms.hpp"

#   include "hashfunc.hpp"
#   include "geometry.hpp"

#   include "tessellate/delaunay_tri_type_2.hpp"
#   include "tessellate/delaunay_tri_type_3.hpp"
//  include "tessellate/delaunay_tri_type_4.hpp"

#   include "tessellate/delaunay_tri_k.hpp"

#   include "tessellate/delaunay_tri_pred_2.hpp"
#   include "tessellate/delaunay_tri_pred_3.hpp"
//  include "tessellate/delaunay_tri_pred_4.hpp"

#   include "tessellate/laguerre_tri_pred_2.hpp"
#   include "tessellate/laguerre_tri_pred_3.hpp"
//  include "tessellate/laguerre_tri_pred_4.hpp"

#   endif//__TESSELLATE__



