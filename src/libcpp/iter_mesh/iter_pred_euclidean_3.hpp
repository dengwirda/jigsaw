
    /*
    --------------------------------------------------------
     * ITER-PRED-EUCLIDEAN-3: predicates for MESH-ITER-3.
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
     * Last updated: 12 August, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_PRED_EUCLIDEAN_3__
#   define __ITER_PRED_EUCLIDEAN_3__

    namespace mesh {

    //_pred.vnrm();
    //_pred.vcos();
    
    
    template <
    typename R ,
    typename I
             >
    class iter_pred_euclidean_3d
        {
        public  :
        typedef R                   real_type ;
        typedef I                   iptr_type ;
        
        iptr_type static constexpr _dims = +3 ; 
         
        public  :
        
        __static_call
        __inline_call real_type mass_tria (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {   return geometry
                ::tria_area_3d (
                   _ipos, _jpos, _kpos) ;
        }
        
        __static_call
        __inline_call real_type mass_tria (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos
            )
        {   return geometry
                ::tetra_vol_3d (
                   _ipos, _jpos, _kpos, _lpos) ;
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
                ::circ_ball_3d (
            _ball, _ipos, _jpos, _kpos, _bind) ;
        }
        
        __static_call
        __inline_call void_type circ_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos ,
                      bool_type  _bind = false
            )
        {   return geometry
                ::circ_ball_3d ( _ball, 
            _ipos, _jpos, _kpos, _lpos, _bind) ;
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
                ::perp_ball_3d (
            _ball, _ipos, _jpos, _kpos, _bind) ;
        }
        
        __static_call
        __inline_call void_type perp_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos ,
                      bool_type  _bind = false
            )
        {   return geometry
                ::perp_ball_3d ( _ball, 
            _ipos, _jpos, _kpos, _lpos, _bind) ;
        }
        
        __static_call
        __inline_call real_type cost_tria (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {   return geometry
                ::tria_quality_3d (
                   _ipos, _jpos, _kpos) ;
        }
        
        __static_call
        __inline_call real_type cost_tria (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos
            )
        {   return geometry
                ::tria_quality_3d (
                   _ipos, _jpos, _kpos, _lpos) ;
        }
        
        __static_call
        __inline_call real_type cost_dual (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {   return geometry
                ::dual_quality_3d (
                   _ipos, _jpos, _kpos) ;
        }
        
        __static_call
        __inline_call real_type cost_dual (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos
            )
        {   return geometry
                ::dual_quality_3d (
                   _ipos, _jpos, _kpos, _lpos) ;
        }
        
        __static_call
        __inline_call real_type length_sq (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos
            )
        {   return geometry::
                lensqr_3d(_ipos, _jpos) ;
        }
        
        __static_call
        __inline_call real_type length_sq (
          __const_ptr(real_type) _vvec
            )
        {   return 
            geometry::lensqr_3d (_vvec) ;
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
        {   // R^3: to-do!
            __unreferenced(_geom) ;
            __unreferenced(_save) ;
            __unreferenced(_proj) ;
        }
        
        } ;

    }


#   endif   //__ITER_PRED_EUCLIDEAN_3__



