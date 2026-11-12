
    /*
    --------------------------------------------------------
     * MESH-PARAMS: parameters for meshing kernels.
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
    --------------------------------------------------------
     *
     * Last updated: 11 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MESH_PARAMS__
#   define __MESH_PARAMS__

#   include "../../liblib/init_jig_t.hpp"

    namespace mesh {

    /*
    --------------------------------------------------------
     * MESH-PARAMS: user-param's for meshing kernels
    --------------------------------------------------------
     */

    class mesh_params : public jigsaw_jig_t
        {
        public  :

        typedef  real_t         real_type ;
        typedef  indx_t         iptr_type ;

        enum node_kind {        // "off-centre" descriptors
        null_kind = +0 ,
        fail_kind = +1 ,
        circ_kind = +2 ,        // "circ"-type refinement
        sink_kind = JIGSAW_MESH_RULE_SINK ,
                                // "sink"-type off-centre
        offH_kind = JIGSAW_MESH_RULE_OFFH ,
                                // "size"-type off-centre
        offC_kind = JIGSAW_MESH_RULE_OFFC ,
                                // "circ"-type off-centre
        offE_kind = JIGSAW_MESH_RULE_OFFE ,
                                // "err."-type off-centre
        offT_kind = JIGSAW_MESH_RULE_OFFT ,
                                // "topo"-type off-centre
        last_kind }    ;

        public  :

    /*-------------------------- construct default param. */
        __inline_call mesh_params  (
            )
        {   // load default values
            jigsaw_null_jig_t(this) ;
        }

    /*------------------------------------ "write" access */
        __inline_call iptr_type      & verb (
            )
        {   return  this->_verbosity ;
        }
        __inline_call iptr_type      & nprt (
            )
        {   return  this->_numthread ;
        }

        __inline_call iptr_type      & rule (
            )
        {   return  this->_mesh_rule ;
        }

        __inline_call iptr_type      & iter (
            )
        {   return  this->_mesh_iter ;
        }

        __inline_call iptr_type      & seed (
            )
        {   return  this->_geom_seed ;
        }

        __inline_call real_type      & phi1 (
            )
        {   return  this->_mesh_phi1 ;
        }
        __inline_call real_type      & phi2 (
            )
        {   return  this->_mesh_phi2 ;
        }

        __inline_call real_type      & eta1 (
            )
        {   return  this->_geom_eta1 ;
        }
        __inline_call real_type      & eta2 (
            )
        {   return  this->_geom_eta2 ;
        }

        __inline_call iptr_type      & feat (
            )
        {   return  this->_geom_feat ;
        }

        __inline_call real_type      & near (
            )
        {   return  this->_init_near ;
        }

        __inline_call iptr_type      & dims (
            )
        {   return  this->_mesh_dims ;
        }

        __inline_call real_type      & siz1 (
            )
        {   return  this->_mesh_siz1 ;
        }
        __inline_call real_type      & siz2 (
            )
        {   return  this->_mesh_siz2 ;
        }
        __inline_call real_type      & siz3 (
            )
        {   return  this->_mesh_siz3 ;
        }

        __inline_call real_type      & eps1 (
            )
        {   return  this->_mesh_eps1 ;
        }
        __inline_call real_type      & eps2 (
            )
        {   return  this->_mesh_eps2 ;
        }

        __inline_call real_type      & rad2 (
            )
        {   return  this->_mesh_rad2 ;
        }
        __inline_call real_type      & rad3 (
            )
        {   return  this->_mesh_rad3 ;
        }

        __inline_call real_type      & off2 (
            )
        {   return  this->_mesh_off2 ;
        }
        __inline_call real_type      & off3 (
            )
        {   return  this->_mesh_off3 ;
        }

        __inline_call real_type      & snk2 (
            )
        {   return  this->_mesh_snk2 ;
        }
        __inline_call real_type      & snk3 (
            )
        {   return  this->_mesh_snk3 ;
        }

        __inline_call real_type      & vol3 (
            )
        {   return  this->_mesh_vol3 ;
        }

        __inline_call iptr_type      & orph (
            )
        {   return  this->_mesh_orph ;
        }

        __inline_call iptr_type      & lock (
            )
        {   return  this->_mesh_lock ;
        }

        __inline_call iptr_type      & top1 (
            )
        {   return  this->_mesh_top1 ;
        }
        __inline_call iptr_type      & top2 (
            )
        {   return  this->_mesh_top2 ;
        }

    /*------------------------------------ "const" access */
        __inline_call iptr_type const& verb (
            ) const
        {   return  this->_verbosity ;
        }
        __inline_call iptr_type const& nprt (
            ) const
        {   return  this->_numthread ;
        }

        __inline_call iptr_type const& rule (
            ) const
        {   return  this->_mesh_rule ;
        }

        __inline_call iptr_type const& iter (
            ) const
        {   return  this->_mesh_iter ;
        }

        __inline_call iptr_type const& seed (
            ) const
        {   return  this->_geom_seed ;
        }

        __inline_call real_type const& phi1 (
            ) const
        {   return  this->_mesh_phi1 ;
        }
        __inline_call real_type const& phi2 (
            ) const
        {   return  this->_mesh_phi2 ;
        }

        __inline_call real_type const& eta1 (
            ) const
        {   return  this->_geom_eta1 ;
        }
        __inline_call real_type const& eta2 (
            ) const
        {   return  this->_geom_eta2 ;
        }

        __inline_call iptr_type const& feat (
            ) const
        {   return  this->_geom_feat ;
        }

        __inline_call real_type const& near (
            ) const
        {   return  this->_init_near ;
        }

        __inline_call iptr_type const& dims (
            ) const
        {   return  this->_mesh_dims ;
        }

        __inline_call real_type const& siz1 (
            ) const
        {   return  this->_mesh_siz1 ;
        }
        __inline_call real_type const& siz2 (
            ) const
        {   return  this->_mesh_siz2 ;
        }
        __inline_call real_type const& siz3 (
            ) const
        {   return  this->_mesh_siz3 ;
        }

        __inline_call real_type const& eps1 (
            ) const
        {   return  this->_mesh_eps1 ;
        }
        __inline_call real_type const& eps2 (
            ) const
        {   return  this->_mesh_eps2 ;
        }

        __inline_call real_type const& rad2 (
            ) const
        {   return  this->_mesh_rad2 ;
        }
        __inline_call real_type const& rad3 (
            ) const
        {   return  this->_mesh_rad3 ;
        }

        __inline_call real_type const& off2 (
            ) const
        {   return  this->_mesh_off2 ;
        }
        __inline_call real_type const& off3 (
            ) const
        {   return  this->_mesh_off3 ;
        }

        __inline_call real_type const& snk2 (
            ) const
        {   return  this->_mesh_snk2 ;
        }
        __inline_call real_type const& snk3 (
            ) const
        {   return  this->_mesh_snk3 ;
        }

        __inline_call real_type const& vol3 (
            ) const
        {   return  this->_mesh_vol3 ;
        }

        __inline_call iptr_type const& orph (
            ) const
        {   return  this->_mesh_orph ;
        }

        __inline_call iptr_type const& lock (
            ) const
        {   return  this->_mesh_lock ;
        }

        __inline_call iptr_type const& top1 (
            ) const
        {   return  this->_mesh_top1 ;
        }
        __inline_call iptr_type const& top2 (
            ) const
        {   return  this->_mesh_top2 ;
        }

        } ;

    }

#   endif   // __MESH_PARAMS__

