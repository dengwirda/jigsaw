
    /*
    --------------------------------------------------------
     * ITER-MESH-3: mesh-optimisation for 3-complexes.
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
     * Last updated: 9 August, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_MESH_3__
#   define __ITER_MESH_3__

    namespace  mesh {
    
    /*
    --------------------------------------------------------
     * ITER-MESH-3: hill-climbing surf. iter.
    --------------------------------------------------------
     */

    template <
    typename G ,
    typename M ,
    typename H ,
    typename P
             >
    class iter_mesh_3
    {
    public  :  
    typedef  M                          mesh_type ;
    typedef  G                          geom_type ;
    typedef  H                          size_type ;
    typedef  P                          pred_type ;

    typedef typename 
            mesh_type::real_type        real_type ;
    typedef typename 
            mesh_type::iptr_type        iptr_type ;
    
    iptr_type static 
        constexpr _dims = pred_type::_dims ;
    
    typedef mesh::iter_params  <
            real_type ,
            iptr_type          >        iter_opts ;
            
    typedef mesh::iter_timers  <
            real_type ,
            iptr_type          >        iter_stat ;
    
    typedef containers
            ::array< iptr_type >        iptr_list ;        
    typedef containers
            ::array< real_type >        real_list ;
    
    public  :
    
    
    
    } ;
    
    }
    
#   endif   //__ITER_MESH_3__
    
    
    
