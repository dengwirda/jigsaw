
    /*
    --------------------------------------------------------
     * JIG-READ: parse *.JIG file into config. data.
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

#   ifndef __JIG_READ__
#   define __JIG_READ__

    class jcfg_reader_base
    {
/*------------------------ inherit from base to read file */
    public  :
    __normal_call void_type push_verbosity (
        std::int32_t /*_verb*/
        ) { }
        
    __normal_call void_type push_geom_file (
        std::string  /*_file*/
        ) { }
    __normal_call void_type push_hfun_file (
        std::string  /*_file*/
        ) { }
    __normal_call void_type push_init_file (
        std::string  /*_file*/
        ) { }
    __normal_call void_type push_mesh_file (
        std::string  /*_file*/
        ) { }
    __normal_call void_type push_tria_file (
        std::string  /*_file*/
        ) { }
    __normal_call void_type push_bnds_file (
        std::string  /*_file*/
        ) { }
        
    __normal_call void_type push_geom_seed (
        std::int32_t /*_seed*/
        ) { }
    __normal_call void_type push_geom_feat (
        bool         /*_feat*/
        ) { }
    __normal_call void_type push_geom_phi1 (
        double       /*_phi1*/
        ) { }
    __normal_call void_type push_geom_phi2 (
        double       /*_phi2*/
        ) { }
    __normal_call void_type push_geom_eta1 (
        double       /*_eta1*/
        ) { }
    __normal_call void_type push_geom_eta2 (
        double       /*_eta2*/
        ) { }
        
    __normal_call void_type push_hfun_scal (
        std::int32_t /*_scal*/
        ) { }
    __normal_call void_type push_hfun_hmax (
        double       /*_hmax*/
        ) { }
    __normal_call void_type push_hfun_hmin (
        double       /*_hmin*/
        ) { }
        
    __normal_call void_type push_mesh_kern (
        std::int32_t /*_kern*/
        ) { }
    __normal_call void_type push_bnds_kern (
        std::int32_t /*_kern*/
        ) { }
    __normal_call void_type push_mesh_dims (
        std::int32_t /*_dims*/
        ) { }
    __normal_call void_type push_mesh_iter (
        std::int32_t /*_iter*/
        ) { }
    __normal_call void_type push_mesh_siz1 (
        double       /*_siz1*/
        ) { }
    __normal_call void_type push_mesh_siz2 (
        double       /*_siz2*/
        ) { }
    __normal_call void_type push_mesh_siz3 (
        double       /*_siz3*/
        ) { }      
    __normal_call void_type push_mesh_top1 (
        bool         /*_top1*/
        ) { }
    __normal_call void_type push_mesh_top2 (
        bool         /*_top2*/
        ) { }
    __normal_call void_type push_mesh_rad2 (
        double       /*_rad2*/
        ) { }
    __normal_call void_type push_mesh_rad3 (
        double       /*_rad3*/
        ) { }
    __normal_call void_type push_mesh_off2 (
        double       /*_off2*/
        ) { }
    __normal_call void_type push_mesh_off3 (
        double       /*_off3*/
        ) { }
    __normal_call void_type push_mesh_snk2 (
        double       /*_snk2*/
        ) { }
    __normal_call void_type push_mesh_snk3 (
        double       /*_snk3*/
        ) { }
    __normal_call void_type push_mesh_eps1 (
        double       /*_eps1*/
        ) { }
    __normal_call void_type push_mesh_eps2 (
        double       /*_eps2*/
        ) { }
    __normal_call void_type push_mesh_vol3 (
        double       /*_vol3*/
        ) { }
        
    __normal_call void_type push_optm_iter (
        std::int32_t /*_iter*/
        ) { }
    __normal_call void_type push_optm_qtol (
        double       /*_qtol*/
        ) { }
    __normal_call void_type push_optm_qlim (
        double       /*_qlim*/
        ) { }
    __normal_call void_type push_optm_tria (
        bool         /*_flag*/
        ) { }
    __normal_call void_type push_optm_dual (
        bool         /*_flag*/
        ) { }
    __normal_call void_type push_optm_div_ (
        bool         /*_flag*/
        ) { }
    __normal_call void_type push_optm_zip_ (
        bool         /*_flag*/
        ) { }
        
    } ;
    
    /*
    --------------------------------------------------------
     * JCFG-READER: read *.JIG mesh files
    --------------------------------------------------------
     */

    class jcfg_reader
    {
    public  :
    
    typedef containers::array <
            std::string > string_tokens;
    
    string_tokens           _errs ;
    
    public  :
    
    /*
    --------------------------------------------------------
     * READ-FILE: read *.MSH file into MESH
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_file (
        std::ifstream&_ffid ,
        dest_type   &&_dest
        )
    {
    
    /*---------------------------------- flip to up.-case */
        #define __toUPPER(__str)        \
            std::transform(__str.begin(),   \
                           __str.  end(),   \
                           __str.begin(),   \
                ::toupper ) ;

    /*---------------------------------- read "file" data */
        #define __putFILE(__fun, __tok)     \
            if (__tok.count() == +2)    \
           _dest.__fun(                 \
                 trim(__tok[ 1])) ;     \
            else                        \
           _errs.push_tail(_line) ;
           
    /*---------------------------------- read "MESH" pred */
        #define __putMESH(__fun, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find("DELAUNAY")!= \
                    std::string::npos ) \
                _dest.__fun (           \
                    jcfg_data           \
               ::mesh_pred::delaunay) ; \
            else                        \
            if (__str[1].find("DELFRONT")!= \
                    std::string::npos ) \
                _dest.__fun (           \
                    jcfg_data           \
               ::mesh_pred::delfront) ; \
            else                        \
           _errs.push_tail(_line) ;     \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
           
    /*---------------------------------- read "BNDS" pred */
        #define __putBNDS(__fun, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find("BND-TRIA")!= \
                    std::string::npos ) \
                _dest.__fun (           \
                    jcfg_data           \
               ::bnds_pred::bnd_tria) ; \
            else                        \
            if (__str[1].find("BND-DUAL")!= \
                    std::string::npos ) \
                _dest.__fun (           \
                    jcfg_data           \
               ::bnds_pred::bnd_dual) ; \
            else                        \
           _errs.push_tail(_line) ;     \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
    
    /*---------------------------------- read "SCAL" data */
        #define __putSCAL(__fun, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find("ABSOLUTE")!= \
                    std::string::npos ) \
                _dest.__fun (           \
                    jcfg_data           \
               ::hfun_scal::absolute) ; \
            else                        \
            if (__str[1].find("RELATIVE")!= \
                    std::string::npos ) \
                _dest.__fun (           \
                    jcfg_data           \
               ::hfun_scal::relative) ; \
            else                        \
           _errs.push_tail(_line) ;     \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
    
    /*---------------------------------- read "real" data */
        #define __putREAL(__fun, __tok)     \
            if (__tok.count() == +2)    \
            {                           \
           _dest.__fun(                 \
                std::stod(__tok[1])) ;  \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
           
    /*---------------------------------- read "ints" data */
        #define __putINTS(__fun, __tok)     \
            if (__tok.count() == +2)    \
            {                           \
           _dest.__fun(                 \
                std::stol(__tok[1])) ;  \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
           
    /*---------------------------------- read "bool" data */
        #define __putBOOL(__fun, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find( "TRUE") !=   \
                    std::string::npos)  \
           _dest.__fun( true);          \
            else                        \
            if (__str[1].find("FALSE") !=   \
                    std::string::npos)  \
           _dest.__fun(false);          \
            else                        \
           _errs.push_tail(_line) ;     \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
  
    
        std::string _line;
        while (std::getline(_ffid, _line))
        {
            _line = trim(_line) ;
            
            if (_line.size() <= 0) continue ;
            if (_line[ +0] == '#') continue ;
        
            try
            {
            containers::
                array<std::string> _stok ;
            
            find_toks(_line, "=;", _stok);
                     
            for (auto _iter  = _stok.head() ;
                      _iter != _stok.tend() ;
                    ++_iter  )
        /*---------------------------- trim on each token */
            *_iter = trim( *_iter ) ;
        
            std::transform(_stok[0].begin() , 
                           _stok[0].  end() , 
                           _stok[0].begin() , 
                         ::toupper) ;
        
        /*---------------------------- read MISC keywords */    
            if (_stok[0] == "VERBOSITY")
                {
            __putINTS(push_verbosity, _stok) ;
                }
            else
        /*---------------------------- read GEOM keywords */
            if (_stok[0] == "GEOM_FILE")
                {
            __putFILE(push_geom_file, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_SEED")
                {
            __putINTS(push_geom_seed, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_PHI1")
                {
            __putREAL(push_geom_phi1, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_PHI2")
                {
            __putREAL(push_geom_phi2, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_ETA1")
                {
            __putREAL(push_geom_eta1, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_ETA2")
                {
            __putREAL(push_geom_eta2, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_FEAT")
                {
            __putBOOL(push_geom_feat, _stok) ;
                }
            else
        /*---------------------------- read INIT keywords */
            if (_stok[0] == "INIT_FILE")
                {
            __putFILE(push_init_file, _stok) ;
                }
            else
        /*---------------------------- read HFUN keywords */
            if (_stok[0] == "HFUN_FILE")
                {
            __putFILE(push_hfun_file, _stok) ;
                }
            else
            if (_stok[0] == "HFUN_SCAL")
                {
            __putSCAL(push_hfun_scal, _stok) ;
                }
            else
            if (_stok[0] == "HFUN_HMAX")
                {
            __putREAL(push_hfun_hmax, _stok) ;
                }
            else
            if (_stok[0] == "HFUN_HMIN")
                {
            __putREAL(push_hfun_hmin, _stok) ;
                }          
            else
        /*---------------------------- read MESH keywords */
            if (_stok[0] == "TRIA_FILE")
                {
            __putFILE(push_tria_file, _stok) ;
                }
            else
            if (_stok[0] == "BNDS_FILE")
                {
            __putFILE(push_bnds_file, _stok) ;
                }
            else
        /*---------------------------- read MESH keywords */
            if (_stok[0] == "MESH_FILE")
                {
            __putFILE(push_mesh_file, _stok) ;
                }
            else
            if (_stok[0] == "MESH_KERN")
                {
            __putMESH(push_mesh_kern, _stok) ;
                }
            else
            if (_stok[0] == "BNDS_KERN")
                {
            __putBNDS(push_bnds_kern, _stok) ;
                }
            else
            if (_stok[0] == "MESH_DIMS")
                {
            __putINTS(push_mesh_dims, _stok) ;
                }
            else
            if (_stok[0] == "MESH_SIZ1")
                {
            __putREAL(push_mesh_siz1, _stok) ;
                }
            else
            if (_stok[0] == "MESH_SIZ2")
                {
            __putREAL(push_mesh_siz2, _stok) ;
                }
            else
            if (_stok[0] == "MESH_SIZ3")
                {
            __putREAL(push_mesh_siz3, _stok) ;
                }
            else
            if (_stok[0] == "MESH_EPS1")
                {
            __putREAL(push_mesh_eps1, _stok) ;
                }
            else
            if (_stok[0] == "MESH_EPS2")
                {
            __putREAL(push_mesh_eps2, _stok) ;
                }
            else
            if (_stok[0] == "MESH_RAD2")
                {
            __putREAL(push_mesh_rad2, _stok) ;
                }
            else
            if (_stok[0] == "MESH_RAD3")
                {
            __putREAL(push_mesh_rad3, _stok) ;
                }
            else
            if (_stok[0] == "MESH_OFF2")
                {
            __putREAL(push_mesh_off2, _stok) ;
                }          
            else
            if (_stok[0] == "MESH_OFF3")
                {
            __putREAL(push_mesh_off3, _stok) ;
                }
            else
            if (_stok[0] == "MESH_SNK2")
                {
            __putREAL(push_mesh_snk2, _stok) ;
                }
            else
            if (_stok[0] == "MESH_SNK3")
                {
            __putREAL(push_mesh_snk3, _stok) ;
                }
            else
            if (_stok[0] == "MESH_VOL3")
                {
            __putREAL(push_mesh_vol3, _stok) ;
                }
            else
            if (_stok[0] == "MESH_TOP1")
                {
            __putBOOL(push_mesh_top1, _stok) ;
                }
            else
            if (_stok[0] == "MESH_TOP2")
                {
            __putBOOL(push_mesh_top2, _stok) ;
                }
            else
            if (_stok[0] == "MESH_ITER")
                {
            __putINTS(push_mesh_iter, _stok) ;
                }
            else
        /*---------------------------- read OPTM keywords */
            if (_stok[0] == "OPTM_ITER")
                {
            __putINTS(push_optm_iter, _stok) ;
                }
            else
            if (_stok[0] == "OPTM_QTOL")
                {
            __putREAL(push_optm_qtol, _stok) ;
                }
            if (_stok[0] == "OPTM_QLIM")
                {
            __putREAL(push_optm_qlim, _stok) ;
                }
            else
            if (_stok[0] == "OPTM_ZIP_")
                {
            __putBOOL(push_optm_zip_, _stok) ;
                }
            else
            if (_stok[0] == "OPTM_DIV_")
                {
            __putBOOL(push_optm_div_, _stok) ;
                }
            else
            if (_stok[0] == "OPTM_TRIA")
                {
            __putBOOL(push_optm_tria, _stok) ;
                }
            else
            if (_stok[0] == "OPTM_DUAL")
                {
            __putBOOL(push_optm_dual, _stok) ;
                }
  
            }
            catch (...)
            {
                this->
               _errs.push_tail (_line) ;
            }        
        }
    
        #undef  __toUPPER
    
        #undef  __putFILE
        #undef  __putSCAL
        #undef  __putMESH
        #undef  __putBNDS
        #undef  __putREAL
        #undef  __putINTS 
        #undef  __putBOOL
                
    }    
    
    } ;
      

#   endif   //__JIG_READ__
        
        
        
