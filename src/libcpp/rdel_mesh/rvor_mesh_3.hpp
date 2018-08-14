
    /*
    --------------------------------------------------------
     * RVOR-MESH-3: restricted voronoi mesh-gen. in R^3. 
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
     * Last updated: 06 April, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RVOR_MESH_3__
#   define __RVOR_MESH_3__

    namespace mesh {
    
    template <
    typename M ,
    typename P ,
    typename G ,
    typename H ,
    typename A = allocators::basic_alloc
             >
    class rvor_mesh_3d
    { 
    public  : 
    
    /*--------- restricted voronoi mesh-generation in R^3 */  
    
    typedef M                               mesh_type ;
    typedef P                               mesh_pred ;
    typedef G                               geom_type ;
    typedef H                               hfun_type ;
    typedef A                               allocator ;

    typedef typename 
            mesh_type::real_type            real_type ;
    typedef typename 
            mesh_type::iptr_type            iptr_type ;
            
            
    } ;

 
    }
    
#   endif   //__RVOR_MESH_3__               
            
            
            
