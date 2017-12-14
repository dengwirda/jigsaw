
    /*
    --------------------------------------------------------
     * ITER-PRED-EUCLIDEAN-2: predicates for MESH-ITER-2.
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
     * Last updated: 12 December, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_PRED_EUCLIDEAN_2__
#   define __ITER_PRED_EUCLIDEAN_2__

    namespace mesh {

    //_pred.vnrm();
    //_pred.vcos();
    
    
    template <
    typename R ,
    typename I
             >
    class iter_pred_euclidean_2d
        {
        public  :
        typedef R                   real_type ;
        typedef I                   iptr_type ;
        
        iptr_type static constexpr _dims = +2 ; 
         
        public  :
        
        __static_call
        __inline_call real_type mass_tria (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {   return geometry
                ::tria_area_2d (
                   _ipos, _jpos, _kpos) ;
        }
        
        __static_call
        __inline_call void_type circ_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
                      bool_type  _bind = false
            )
        {   return geometry
                ::circ_ball_2d (
            _ball, _ipos, _jpos, _kpos, _bind) ;
        }
        
        __static_call
        __inline_call void_type perp_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
                      bool_type  _bind = false
            )
        {   return geometry
                ::perp_ball_2d (
            _ball, _ipos, _jpos, _kpos, _bind) ;
        }
        
        __static_call
        __inline_call real_type cost_tria (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {   return geometry
                ::tria_quality_2d (
                   _ipos, _jpos, _kpos) ;
        }
        
        __static_call
        __inline_call real_type cost_dual (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {   return geometry
                ::dual_quality_2d (
                   _ipos, _jpos, _kpos) ;
        }
        
        __static_call
        __inline_call real_type length_sq (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos
            )
        {   return geometry::
                lensqr_2d(_ipos, _jpos) ;
        }
        
        __static_call
        __inline_call real_type length_sq (
          __const_ptr(real_type) _vvec
            )
        {   return 
            geometry::lensqr_2d (_vvec) ;
        }
      
        template <
        typename      geom_type
                 >
        __static_call
        __inline_call void_type proj_node (
            geom_type &_geom ,
          __const_ptr(real_type) _save ,
          __write_ptr(real_type) _proj
            )
        {   // R^2: do nothing!
            __unreferenced(_geom) ;
            __unreferenced(_save) ;
            __unreferenced(_proj) ;
        }
        
        } ;

    }


#   endif   //__ITER_PRED_EUCLIDEAN_2__



