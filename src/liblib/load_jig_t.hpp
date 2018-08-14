
    /*
    --------------------------------------------------------
     * LOAD-JIG_T: load jig_t for lib_jigsaw.
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

#   ifndef __LOAD_JIG_T__
#   define __LOAD_JIG_T__

    /*
    --------------------------------------------------------
     * LOAD-JIG_T: load jig_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

#   include  "../jig_read.hpp"

    indx_t jigsaw_load_jig_t (              // load jig_t
        char *_fchr,
        jigsaw_jig_t *_jcfg
        )
    {
        class jcfg_loader: 
            public jcfg_reader_base
        {
        public  :
            jigsaw_jig_t      *_jjig;
        
        public  :
    /*------------------------------------- create loader */
        __normal_call jcfg_loader (
            jigsaw_jig_t *_jsrc = nullptr
            ) : _jjig( _jsrc ) {}
        
    /*------------------------------------- MISC keywords */
        __normal_call void_type push_verbosity (
            std::int32_t  _verb
            ) 
        { 
            this->
           _jjig->_verbosity = _verb;       
        }          
        
    /*------------------------------------- GEOM keywords */       
        __normal_call void_type push_geom_seed (
            std::int32_t  _seed
            ) 
        {
            this->
           _jjig->_geom_seed = _seed; 
        }
        __normal_call void_type push_geom_feat (
            bool          _feat
            ) 
        {
            this->
           _jjig->_geom_feat = _feat; 
        }
        __normal_call void_type push_geom_eta1 (
            double        _eta1
            ) 
        {
            this->
           _jjig->_geom_eta1 = _eta1; 
        }
        __normal_call void_type push_geom_eta2 (
            double        _eta2
            ) 
        {
            this->
           _jjig->_geom_eta2 = _eta2; 
        }
        
    /*------------------------------------- HFUN keywords */
        __normal_call void_type push_hfun_scal (
            std::int32_t  _scal
            ) 
        {
            this->
           _jjig->_hfun_scal = _scal; 
        }
        __normal_call void_type push_hfun_hmax (
            double        _hmax
            ) 
        {
            this->
           _jjig->_hfun_hmax = _hmax; 
        }
        __normal_call void_type push_hfun_hmin (
            double        _hmin
            ) 
        {
            this->
           _jjig->_hfun_hmin = _hmin; 
        }
    
    /*------------------------------------- KERN keywords */
        __normal_call void_type push_bnds_kern (
            std::int32_t  _kern
            ) 
        { 
            this->
           _jjig->_bnds_kern = _kern;
        }
        
        __normal_call void_type push_mesh_kern (
            std::int32_t  _kern
            ) 
        {
            this->
           _jjig->_mesh_kern = _kern;
        }
        
    /*------------------------------------- MESH keywords */
        __normal_call void_type push_mesh_dims (
            std::int32_t  _dims
            ) 
        { 
            this->
           _jjig->_mesh_dims = _dims;
        }
        __normal_call void_type push_mesh_iter (
            std::int32_t  _iter
            ) 
        { 
            this->
           _jjig->_mesh_iter = _iter;
        }
        __normal_call void_type push_mesh_siz1 (
            double        _siz1
            ) 
        { 
            this->
           _jjig->_mesh_siz1 = _siz1;
        }
        __normal_call void_type push_mesh_siz2 (
            double        _siz2
            ) 
        { 
            this->
           _jjig->_mesh_siz2 = _siz2;
        }
        __normal_call void_type push_mesh_siz3 (
            double        _siz3
            ) 
        {
            this->
           _jjig->_mesh_siz3 = _siz3; 
        }      
        __normal_call void_type push_mesh_top1 (
            bool          _top1
            ) 
        { 
            this->
           _jjig->_mesh_top1 = _top1;
        }
        __normal_call void_type push_mesh_top2 (
            bool          _top2
            ) 
        {
            this->
           _jjig->_mesh_top2 = _top2; 
        }
        __normal_call void_type push_mesh_rad2 (
            double        _rad2
            ) 
        {
            this->
           _jjig->_mesh_rad2 = _rad2; 
        }
        __normal_call void_type push_mesh_rad3 (
            double        _rad3
            ) 
        {
            this->
           _jjig->_mesh_rad3 = _rad3; 
        }
        __normal_call void_type push_mesh_off2 (
            double        _off2
            ) 
        { 
            this->
           _jjig->_mesh_off2 = _off2;
        }
        __normal_call void_type push_mesh_off3 (
            double        _off3
            ) 
        {
            this->
           _jjig->_mesh_off3 = _off3; 
        }
        __normal_call void_type push_mesh_snk2 (
            double        _snk2
            ) 
        {
            this->
           _jjig->_mesh_snk2 = _snk2; 
        }
        __normal_call void_type push_mesh_snk3 (
            double        _snk3
            ) 
        {
            this->
           _jjig->_mesh_snk3 = _snk3; 
        }
        __normal_call void_type push_mesh_eps1 (
            double        _eps1
            ) 
        {
            this->
           _jjig->_mesh_eps1 = _eps1; 
        }
        __normal_call void_type push_mesh_eps2 (
            double        _eps2
            ) 
        {
            this->
           _jjig->_mesh_eps2 = _eps2; 
        }
        __normal_call void_type push_mesh_vol3 (
            double        _vol3
            ) 
        {
            this->
           _jjig->_mesh_vol3 = _vol3; 
        }
        
    /*------------------------------------- OPTM keywords */
        __normal_call void_type push_optm_iter (
            std::int32_t  _iter
            ) 
        {
            this->
           _jjig->_optm_iter = _iter;  
        }
        __normal_call void_type push_optm_qtol (
            double        _qtol
            ) 
        {
            this->
           _jjig->_optm_qtol = _qtol; 
        }
        __normal_call void_type push_optm_qlim (
            double        _qlim
            ) 
        {
            this->
           _jjig->_optm_qlim = _qlim; 
        }
        __normal_call void_type push_optm_tria (
            bool          _flag
            ) 
        {
            this->
           _jjig->_optm_tria = _flag; 
        }
        __normal_call void_type push_optm_dual (
            bool          _flag
            ) 
        {
            this->
           _jjig->_optm_dual = _flag; 
        }
        __normal_call void_type push_optm_div_ (
            bool          _flag
            ) 
        {
            this->
           _jjig->_optm_div_ = _flag; 
        }
        __normal_call void_type push_optm_zip_ (
            bool          _flag
            ) 
        {
            this->
           _jjig->_optm_zip_ = _flag; 
        }
        
        } ;
    
    /*---------------------------------- parse JIG_T file */
        iptr_type _errv  = __no_error ;
    
        try
        {
            jcfg_reader   _read ;
            std::ifstream _file ;
            
            std::string _fstr(_fchr) ;
            std::string _path ;
            std::string _name ;
            std::string _fext ;
            file_part(_fstr, 
                _path, _name, _fext) ;
             
            _file.open(
                _fstr, std::ifstream:: in) ;

            if (_file.is_open() )
            {
                _read.read_file (
                _file, jcfg_loader(_jcfg)) ;
            }
            else
            {           
                _errv = __file_not_located ;
            }
            
            _file.close ();

            if(!_read._errs.empty())
            {
                _errv = __invalid_argument ;
            }
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

    /*---------------------------------- return read flag */
        return  _errv ;
    
    }

#   endif//__lib_jigsaw

#   endif



