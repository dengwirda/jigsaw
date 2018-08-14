
    /*
    --------------------------------------------------------
     * INIT-JIG_T: init. jig_t for lib_jigsaw.
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
     * Last updated: 31 July, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __INIT_JIG_T__
#   define __INIT_JIG_T__

#   include <limits>

    /*
    --------------------------------------------------------
     * INIT-JIG_T: init. jig_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    void jigsaw_init_jig_t (                // init. jig_t
        jigsaw_jig_t *_jcfg
        )
    {
    /*------------------------------------- MISC keywords */
        _jcfg->_verbosity = (indx_t) +0 ;
    
    /*------------------------------------- BNDS keywords */
        _jcfg->_bnds_kern = JIGSAW_BNDS_TRIACELL ;
    
    /*------------------------------------- GEOM keywords */
        _jcfg->_geom_seed = (indx_t) +8 ;
        
        _jcfg->_geom_feat = (indx_t) +0 ;
        _jcfg->_geom_eta1 = (real_t) +45.0 ;
        _jcfg->_geom_eta2 = (real_t) +45.0 ;
        
    /*------------------------------------- HFUN keywords */
        _jcfg->_hfun_scal = JIGSAW_HFUN_RELATIVE ;
        _jcfg->_hfun_hmax = (real_t) +0.02 ;
        _jcfg->_hfun_hmin = (real_t) +0.00 ;
        
    /*------------------------------------- MESH keywords */
        _jcfg->_mesh_kern = JIGSAW_KERN_DELFRONT ;
        _jcfg->_mesh_dims = (indx_t) +3 ;
        _jcfg->_mesh_iter = (indx_t) 
            std::numeric_limits<iptr_type>::max();
        
        _jcfg->_mesh_top1 = (indx_t) +0 ;
        _jcfg->_mesh_top2 = (indx_t) +0 ;
        
        _jcfg->_mesh_rad2 = (real_t) +1.05 ;
        _jcfg->_mesh_rad3 = (real_t) +2.05 ;
        _jcfg->_mesh_off2 = (real_t) +0.90 ;
        _jcfg->_mesh_off3 = (real_t) +1.10 ;
     
        real_t _SIZ1 = +1. *   4./3. ;
        real_t _SIZ2 = +.5 *  (4./3. +
            2. / (1.+std::sqrt(1./3.))) ;
        real_t _SIZ3 = +.5 *  (4./3. +
            2. / (1.+std::sqrt(3./8.))) ;
        
        _jcfg->_mesh_siz1 = (real_t) _SIZ1 ;
        _jcfg->_mesh_siz2 = (real_t) _SIZ2 ;
        _jcfg->_mesh_siz3 = (real_t) _SIZ3 ;
        
        _jcfg->_mesh_snk2 = (real_t) +0.20 ;
        _jcfg->_mesh_snk3 = (real_t) +0.33 ;
        
        _jcfg->_mesh_eps1 = (real_t) +0.33 ;
        _jcfg->_mesh_eps2 = (real_t) +0.33 ;
        
        _jcfg->_mesh_vol3 = (real_t) +0.00 ;
        
    /*------------------------------------- OPTM keywords */
        _jcfg->_optm_iter = (indx_t) +16;
        
        _jcfg->_optm_qtol = (real_t) +1.E-04 ;
        _jcfg->_optm_qlim = (real_t) +0.9375 ; 
        
        _jcfg->_optm_tria = (indx_t) +1 ;
        _jcfg->_optm_dual = (indx_t) +0 ;    
        _jcfg->_optm_div_ = (indx_t) +1 ;
        _jcfg->_optm_zip_ = (indx_t) +1 ; 
    }
    
#   endif//__lib_jigsaw

#   undef  SIZ1
#   undef  SIZ2
#   undef  SIZ3

#   endif
    
    
    
