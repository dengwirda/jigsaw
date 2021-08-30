
    /*
    --------------------------------------------------------
     * ITER-PRED-ELLIPSOID-3: predicates for ITER-MESH-k.
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
     * Last updated: 28 Aug., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_PRED_ELLIPSOID_3__
#   define __ITER_PRED_ELLIPSOID_3__

    namespace mesh {

    template <
    typename G  ,
    size_t   C  ,
    typename M
             >
    class iter_pred_ellipsoid_3d
        {
/*------------------------ optimisation predicates in S^2 */
        public  :

        typedef G                   geom_type ;
        typedef M                   mesh_type ;

        typedef typename
            mesh_type::real_type    real_type ;
        typedef typename
            mesh_type::iptr_type    iptr_type ;

        struct cell_kind {};
        struct dual_kind {};

        size_t    static constexpr _cost =  C ;

        iptr_type static
            constexpr  topo_dims        =  +2 ;
        iptr_type static
            constexpr  geom_dims        =  +3 ;
        iptr_type static
            constexpr  real_dims        =  +4 ;

        public  :

        __static_call
        __inline_call real_type innerprod (
          __const_ptr(real_type) _avec ,
          __const_ptr(real_type) _bvec
            )
        {   return geometry
                ::dot_3d (_avec, _bvec) ;
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

        __static_call
        __inline_call void_type proj_node (
            geom_type &_geom ,
            real_type *_save ,
            real_type *_proj
            )
        {
            __unreferenced(_save);

            real_type _zero[3] ;
            _zero[0] = (real_type) +.0 ;
            _zero[1] = (real_type) +.0 ;
            _zero[2] = (real_type) +.0 ;

            real_type _ttaa, _ttbb ;
            if (_geom.line_surf(
                   _zero, _proj, _ttaa, _ttbb) )
            {
            real_type _pmid[3] = {
            _proj[0] * (real_type) +.5 +
            _zero[0] * (real_type) +.5 ,
            _proj[1] * (real_type) +.5 +
            _zero[1] * (real_type) +.5 ,
            _proj[2] * (real_type) +.5 +
            _zero[2] * (real_type) +.5
                } ;
            real_type _pdel[3] = {
            _proj[0] * (real_type) +.5 -
            _zero[0] * (real_type) +.5 ,
            _proj[1] * (real_type) +.5 -
            _zero[1] * (real_type) +.5 ,
            _proj[2] * (real_type) +.5 -
            _zero[2] * (real_type) +.5
                } ;

            if (_ttaa > (real_type)-1.)
            {
            _proj[0] =
            _pmid[0] + _ttaa*_pdel[ 0] ;
            _proj[1] =
            _pmid[1] + _ttaa*_pdel[ 1] ;
            _proj[2] =
            _pmid[2] + _ttaa*_pdel[ 2] ;
            }
            else
            if (_ttbb > (real_type)-1.)
            {
            _proj[0] =
            _pmid[0] + _ttbb*_pdel[ 0] ;
            _proj[1] =
            _pmid[1] + _ttbb*_pdel[ 1] ;
            _proj[2] =
            _pmid[2] + _ttbb*_pdel[ 2] ;
            }
            }
        }

        __static_call
        __inline_call real_type tri3_mass (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos
            )
        {
            return geometry::tria_area_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ) ;
        }

        __static_call
        __inline_call real_type quad_mass (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos
            )
        {
            return geometry::quad_area_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ,
               &_lpos[0] ) ;
        }

        __static_call
        __inline_call void_type edge_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
                      bool_type  _bind = false
            )
        {
            return geometry::perp_ball_3d (
               &_ball[0],
               &_ipos[0],
               &_jpos[0], _bind) ;
        }

        __static_call
        __inline_call void_type tri3_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
                      bool_type  _bind = false
            )
        {
            return geometry::perp_ball_3d (
               &_ball[0],
               &_ipos[0],
               &_jpos[0],
               &_kpos[0], _bind) ;
        }

        __static_call
        __inline_call void_type quad_ball (
          __write_ptr(real_type) _ball ,
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos ,
                      bool_type  _bind = false
            )
        {
            return geometry::quad_ball_3d (
               &_ball[0],
               &_ipos[0],
               &_jpos[0],
               &_kpos[0],
               &_lpos[0], _bind) ;
        }

        __static_call
        __inline_call real_type tri3_cost (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
            cell_kind const&
            )
        {
            if constexpr (
                _cost == JIGSAW_KERN_AREA_LEN)
            {
            return geometry::tria_quality_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ) ;
            }

            if constexpr (
                _cost == JIGSAW_KERN_SKEW_COS)
            {
            return geometry::tria_skewcos_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ) ;
            }
        }

        __static_call
        __inline_call real_type tri3_cost (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
            dual_kind const&
            )
        {
            return geometry::tria_duality_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ) ;
        }

        __static_call
        __inline_call real_type quad_cost (
          __const_ptr(real_type) _ipos ,
          __const_ptr(real_type) _jpos ,
          __const_ptr(real_type) _kpos ,
          __const_ptr(real_type) _lpos ,
            cell_kind const&
            )
        {
            return geometry::quad_quality_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ,
               &_lpos[0] ) ;
        }

        __static_call
        __inline_call real_type quad_cost (
          __const_ptr(real_type) /*_ipos*/ ,
          __const_ptr(real_type) /*_jpos*/ ,
          __const_ptr(real_type) /*_kpos*/ ,
          __const_ptr(real_type) /*_lpos*/ ,
            dual_kind const&
            )
        {   /*
            return geometry::quad_duality_3d (
               &_ipos[0] ,
               &_jpos[0] ,
               &_kpos[0] ,
               &_lpos[0] ) ;
            */
            return (real_type) +1. ;
        }

        } ;


    }

#   endif   //__ITER_PRED_ELLIPSOID_3__



