
/* 
------------------------------------------------------------
 * std::library-like containers, but not quite...
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
 * Last updated: 04 May, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __CONTAINERS__
#   define __CONTAINERS__

#   include "allocators.hpp"

    namespace containers 
    {
/*----- buffer reallocation options */
    enum alloc_types
    {    loose_alloc, // zero manip.
         tight_alloc  // full manip.
    } ;
    } //     containers

#   include "containers/fixed_array.hpp"
#   include "containers/array.hpp"
#   include "containers/block_array.hpp"

#   include "containers/single_list.hpp"
#   include "containers/double_list.hpp"

#   include "containers/arraylist.hpp"
#   include "containers/hashtable.hpp"

#   include "containers/priorityset.hpp"
#   include "containers/prioritymap.hpp"


#   endif//__CONTAINERS__



