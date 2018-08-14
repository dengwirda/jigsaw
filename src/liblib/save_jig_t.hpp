
    /*
    --------------------------------------------------------
     * SAVE-JIG_T: save jig_t for lib_jigsaw.
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
     * Last updated: 30 July, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __SAVE_JIG_T__
#   define __SAVE_JIG_T__

    /*
    --------------------------------------------------------
     * SAVE-JIG_T: save jig_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    indx_t jigsaw_save_jig_t (              // save jig_t
        char *_fchr,
        jigsaw_jig_t *_jcfg
        )
    {
        iptr_type _errv = __no_error ;
    
        try
        {
            std::ofstream  _file ;

            std::string _fstr(_fchr) ;
            std::string _path ;
            std::string _name ;
            std::string _fext ; 
            file_part(
                _fstr, _path, _name, _fext) ;

            _file.open(
                _fstr, std::ofstream::out ) ;
                     
            if (_file.is_open())
            {
                _file << "# " << _name << ".jig"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                
        /*--------------------------------- MISC keywords */
                _file << "VERBOSITY = " << 
                    _jcfg->_verbosity << "\n" ;
    
        /*--------------------------------- GEOM keywords */
                _file << "GEOM_SEED = " << 
                    _jcfg->_geom_seed << "\n" ;
                
                _file << "GEOM_FEAT = " << 
                    _jcfg->_geom_feat << "\n" ;
                    
                _file << "GEOM_ETA1 = " << 
                    _jcfg->_geom_eta1 << "\n" ;
                _file << "GEOM_ETA2 = " << 
                    _jcfg->_geom_eta2 << "\n" ;
   
        /*--------------------------------- HFUN keywords */
                if (_jcfg->_hfun_scal == 
                        JIGSAW_HFUN_RELATIVE)
                _file << "HFUN_SCAL = " 
                      << "RELATIVE\n" ;
                else
                if (_jcfg->_hfun_scal == 
                        JIGSAW_HFUN_ABSOLUTE)
                _file << "HFUN_SCAL = " 
                      << "ABSOLUTE\n" ;
                    
                _file << "HFUN_HMIN = " << 
                    _jcfg->_hfun_hmax << "\n" ;
                _file << "HFUN_HMAX = " << 
                    _jcfg->_hfun_hmin << "\n" ;
        
        /*--------------------------------- BNDS keywords */
                if (_jcfg->_bnds_kern == 
                        JIGSAW_BNDS_TRIACELL)
                _file << "BNDS_KERN = " 
                      << "BND-TRIA\n" ;
                else
                if (_jcfg->_bnds_kern == 
                        JIGSAW_BNDS_DUALCELL)
                _file << "BNDS_KERN = " 
                      << "BND-DUAL\n" ;
        
        /*--------------------------------- MESH keywords */
                if (_jcfg->_mesh_kern == 
                        JIGSAW_KERN_DELFRONT)
                _file << "MESH_KERN = " 
                      << "DELFRONT\n" ;
                else
                if (_jcfg->_mesh_kern == 
                        JIGSAW_KERN_DELAUNAY)
                _file << "MESH_KERN = " 
                      << "DELAUNAY\n" ;
        
        /*--------------------------------- MESH keywords */
                _file << "MESH_DIMS = " << 
                    _jcfg->_mesh_dims << "\n" ;
                
                _file << "MESH_ITER = " << 
                    _jcfg->_mesh_iter << "\n" ;
            
                _file << "MESH_SIZ1 = " << 
                    _jcfg->_mesh_siz1 << "\n" ;
                _file << "MESH_SIZ2 = " << 
                    _jcfg->_mesh_siz2 << "\n" ;
                _file << "MESH_SIZ3 = " << 
                    _jcfg->_mesh_siz3 << "\n" ;
            
                _file << "MESH_TOP1 = " << 
                    _jcfg->_mesh_top1 << "\n" ;
                _file << "MESH_TOP2 = " << 
                    _jcfg->_mesh_top2 << "\n" ;
            
                _file << "MESH_RAD2 = " << 
                    _jcfg->_mesh_rad2 << "\n" ;
                _file << "MESH_RAD3 = " << 
                    _jcfg->_mesh_rad3 << "\n" ;
            
                _file << "MESH_OFF2 = " << 
                    _jcfg->_mesh_off2 << "\n" ;
                _file << "MESH_OFF3 = " <<
                    _jcfg->_mesh_off3 << "\n" ;
            
                _file << "MESH_SNK2 = " << 
                    _jcfg->_mesh_snk2 << "\n" ;
                _file << "MESH_SNK3 = " << 
                    _jcfg->_mesh_snk3 << "\n" ;
            
                _file << "MESH_EPS1 = " << 
                    _jcfg->_mesh_eps1 << "\n" ;
                _file << "MESH_EPS2 = " << 
                    _jcfg->_mesh_eps2 << "\n" ;
            
                _file << "MESH_VOL3 = " <<
                    _jcfg->_mesh_vol3 << "\n" ;
        
        /*--------------------------------- OPTM keywords */
                _file << "OPTM_ITER = " << 
                    _jcfg->_optm_iter << "\n" ;
            
                _file << "OPTM_QTOL = " << 
                    _jcfg->_optm_qtol << "\n" ;
                _file << "OPTM_QLIM = " << 
                    _jcfg->_optm_qlim << "\n" ; 
            
                _file << "OPTM_TRIA = " << 
                    _jcfg->_optm_tria << "\n" ;
                _file << "OPTM_DUAL = " << 
                    _jcfg->_optm_dual << "\n" ;
                _file << "OPTM_DIV_ = " << 
                    _jcfg->_optm_div_ << "\n" ;
                _file << "OPTM_ZIP_ = " << 
                    _jcfg->_optm_zip_ << "\n" ;
            }
            else
            {
                _errv = __file_not_created ;
            }
            
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }
    
        return  _errv ;
    }

#   endif//__lib_jigsaw

#   endif



