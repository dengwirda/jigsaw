
    /*
    --------------------------------------------------------
     * JIG-LOAD: parse *.JIG file into JCFG data.
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
     * Last updated: 14 December, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __JIG_LOAD__
#   define __JIG_LOAD__

    /*
    --------------------------------------------------------
     * JCFG-READER: read *.JIG data files
    --------------------------------------------------------
     */

    class jcfg_reader
    {
    public  :
    
    typedef containers::array <
            std::string > string_tokens ;
    
    string_tokens           _errs ;
    
    public  :

    /*
    --------------------------------------------------------
     * READ-FILE: read *.MSH file into MESH
    --------------------------------------------------------
     */

    __normal_call void_type read_file (
        std::ifstream&_ffid ,
        jcfg_data    &_jcfg
        )
    {
    /*---------------------------------- flip to up.-case */
        #define __toUPPER(__str)        \
            std::transform(__str.begin(),   \
                           __str.  end(),   \
                           __str.begin(),   \
                ::toupper ) ;

    /*---------------------------------- read "file" data */
        #define __putFILE(__var, __tok)     \
            if (__tok.count() == +2)    \
            _jcfg.__var                 \
                = trim(__tok[1]) ;      \
            else                        \
            _errs.push_tail(_line) ;

    /*---------------------------------- read "MESH" pred */
        #define __putMESH(__var, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find("DELAUNAY")!= \
                    std::string::npos ) \
                _jcfg.__var =           \
            jcfg_data::mesh_pred::delaunay; \
            else                        \
            if (__str[1].find("DELFRONT")!= \
                    std::string::npos ) \
                _jcfg.__var =           \
            jcfg_data::mesh_pred::delfront; \
            else                        \
           _errs.push_tail(_line) ;     \
            }                           \
            else                        \
           _errs.push_tail(_line) ;

    /*---------------------------------- read "SCAL" data */
        #define __putSCAL(__var, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find("ABSOLUTE")!= \
                    std::string::npos ) \
                _jcfg.__var =           \
            jcfg_data::hfun_scal::absolute; \
            else                        \
            if (__str[1].find("RELATIVE")!= \
                    std::string::npos ) \
                _jcfg.__var =           \
            jcfg_data::hfun_scal::relative; \
            else                        \
           _errs.push_tail(_line) ;     \
            }                           \
            else                        \
           _errs.push_tail(_line) ;

    /*---------------------------------- read "real" data */
        #define __putREAL(__var, __tok)     \
            if (__tok.count() == +2)    \
            {                           \
            _jcfg.__var =               \
                std::stod(__tok [1]) ;  \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
           
    /*---------------------------------- read "ints" data */
        #define __putINTS(__var, __tok)     \
            if (__tok.count() == +2)    \
            {                           \
            _jcfg.__var =               \
                std::stol(__tok [1]) ;  \
            }                           \
            else                        \
           _errs.push_tail(_line) ;
           
    /*---------------------------------- read "bool" data */
        #define __putBOOL(__var, __str)     \
            if (__str.count() == 2 )    \
            {                           \
                __toUPPER(__str [1])    \
            if (__str[1].find( "TRUE") !=   \
                    std::string::npos)  \
           _jcfg.__var = true ;         \
            else                        \
            if (__str[1].find("FALSE") !=   \
                    std::string::npos)  \
           _jcfg.__var = false;         \
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
            __putINTS(_verbosity, _stok) ;
            __putINTS(_rdel_opts.verb (), _stok) ;
            __putINTS(_iter_opts.verb (), _stok) ;
                }
            else
        /*---------------------------- read GEOM keywords */
            if (_stok[0] == "GEOM_FILE")
                {
            __putFILE(_geom_file, _stok) ;
                }
            else
            if (_stok[0] == "GEOM_SEED")
                {
            __putINTS(_rdel_opts.seed (), _stok) ;
                }
            else
            if (_stok[0] == "GEOM_PHI1")
                {
            __putREAL(_rdel_opts.phi1 (), _stok) ;
                }
            else
            if (_stok[0] == "GEOM_PHI2")
                {
            __putREAL(_rdel_opts.phi2 (), _stok) ;
                }
            else
            if (_stok[0] == "GEOM_ETA1")
                {
            __putREAL(_rdel_opts.eta1 (), _stok) ;
                }
            else
            if (_stok[0] == "GEOM_ETA2")
                {
            __putREAL(_rdel_opts.eta2 (), _stok) ;
                }
            else
            if (_stok[0] == "GEOM_FEAT")
                {
            __putBOOL(_rdel_opts.feat (), _stok) ;
                }
            else
        /*---------------------------- read HFUN keywords */
            if (_stok[0] == "HFUN_FILE")
                {
            __putFILE(_hfun_file, _stok) ;
                }
            else
            if (_stok[0] == "HFUN_SCAL")
                {
            __putSCAL(_hfun_scal, _stok) ;
                }
            else
            if (_stok[0] == "HFUN_HMAX")
                {
            __putREAL(_hfun_hmax, _stok) ;
                }
            else
            if (_stok[0] == "HFUN_HMIN")
                {
            __putREAL(_hfun_hmin, _stok) ;
                }          
            else
        /*---------------------------- read MESH keywords */
            if (_stok[0] == "MESH_FILE")
                {
            __putFILE(_mesh_file, _stok) ;
                }
            else
            if (_stok[0] == "MESH_KERN")
                {
            __putMESH(_mesh_pred, _stok) ;
                }
            else
            if (_stok[0] == "MESH_DIMS")
                {
            __putINTS(_rdel_opts.dims (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_SIZ1")
                {
            __putREAL(_rdel_opts.siz1 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_SIZ2")
                {
            __putREAL(_rdel_opts.siz2 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_SIZ3")
                {
            __putREAL(_rdel_opts.siz3 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_EPS1")
                {
            __putREAL(_rdel_opts.eps1 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_EPS2")
                {
            __putREAL(_rdel_opts.eps2 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_RAD2")
                {
            __putREAL(_rdel_opts.rad2 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_RAD3")
                {
            __putREAL(_rdel_opts.rad3 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_OFF2")
                {
            __putREAL(_rdel_opts.off2 (), _stok) ;
                }          
            else
            if (_stok[0] == "MESH_OFF3")
                {
            __putREAL(_rdel_opts.off3 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_SNK2")
                {
            __putREAL(_rdel_opts.snk2 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_SNK3")
                {
            __putREAL(_rdel_opts.snk3 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_VOL3")
                {
            __putREAL(_rdel_opts.vol3 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_TOP1")
                {
            __putBOOL(_rdel_opts.top1 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_TOP2")
                {
            __putBOOL(_rdel_opts.top2 (), _stok) ;
                }
            else
            if (_stok[0] == "MESH_ITER")
                {
            __putINTS(_rdel_opts.iter (), _stok) ;
                }
            else
        /*---------------------------- read OPTM keywords */
            if (_stok[0] == "OPTM_ITER")
                {
            __putINTS(_iter_opts.iter (), _stok) ;
                }
            else
            if (_stok[0] == "OPTM_QTOL")
                {
            __putREAL(_iter_opts.qtol (), _stok) ;
                }
            if (_stok[0] == "OPTM_QLIM")
                {
            __putREAL(_iter_opts.qlim (), _stok) ;
                }
            else
            if (_stok[0] == "OPTM_ZIP_")
                {
            __putBOOL(_iter_opts.zip_ (), _stok) ;
                }
            else
            if (_stok[0] == "OPTM_DIV_")
                {
            __putBOOL(_iter_opts.div_ (), _stok) ;
                }
            else
            if (_stok[0] == "OPTM_TRIA")
                {
            __putBOOL(_iter_opts.tria (), _stok) ;
                }
            else
            if (_stok[0] == "OPTM_DUAL")
                {
            __putBOOL(_iter_opts.dual (), _stok) ;
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
        #undef  __putMESH
        #undef  __putSCAL
        #undef  __putREAL
        #undef  __putINTS 
        #undef  __putBOOL
                    
    }

    } ;
    
    /*
    --------------------------------------------------------
     * READ-JCFG: read *.JCFG input file.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type read_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog
        )
    {
        iptr_type _errv  = __no_error ;
        try
        {
            jcfg_reader   _read;
            std::ifstream _file; 
            _file. open(
            _jcfg._jcfg_file, std::ifstream::in);

            if (_file.is_open())
            {
                _read.
                 read_file(_file, _jcfg);
            }
            else
            {           
                _errv = __file_not_located ;
            }
            _file.close ();

            for (auto _iter  = 
                      _read._errs.head();
                      _iter != 
                      _read._errs.tend(); 
                    ++_iter  )
            {
                _jlog.push(
            "  **parse error: " + *_iter + "\n");
            }        
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }
    
    /*
    --------------------------------------------------------
     * COPY-JCFG: read *.JCFG input data.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type copy_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        jigsaw_jig_t const&_jjig
        )
    {
        iptr_type _errv  = __no_error ;
        try
        {
    /*------------------------------------- MISC keywords */
            _jcfg._rdel_opts.
                verb() = _jjig._verbosity ;
            _jcfg._iter_opts.
                verb() = _jjig._verbosity ;
    
    /*------------------------------------- GEOM keywords */
            _jcfg._rdel_opts.
                seed() = _jjig._geom_seed ;
            _jcfg._rdel_opts.
                feat() = _jjig._geom_feat ;
            _jcfg._rdel_opts.
                eta1() = _jjig._geom_eta1 ;
            _jcfg._rdel_opts.
                eta2() = _jjig._geom_eta2 ;
            
    /*------------------------------------- HFUN keywords */
            if (_jjig._hfun_scal == 
                    JIGSAW_HFUN_RELATIVE)
            _jcfg._hfun_scal = 
                jcfg_data::hfun_scal::relative ;
            else
            if (_jjig._hfun_scal ==
                    JIGSAW_HFUN_ABSOLUTE)
            _jcfg._hfun_scal = 
                jcfg_data::hfun_scal::absolute ;
            
            _jcfg.
            _hfun_hmax = _jjig._hfun_hmax ;
            _jcfg.
            _hfun_hmin = _jjig._hfun_hmin ;
            
    /*------------------------------------- RDEL keywords */
            if (_jjig._hfun_scal == 
                    JIGSAW_KERN_DELFRONT)
            _jcfg._mesh_pred = 
                jcfg_data::mesh_pred::delfront ;
            else
            if (_jjig._hfun_scal ==
                    JIGSAW_KERN_DELAUNAY)
            _jcfg._mesh_pred = 
                jcfg_data::mesh_pred::delaunay ;
            
            _jcfg._rdel_opts.
                dims() = _jjig._mesh_dims ;
            _jcfg._rdel_opts.
                iter() = _jjig._mesh_iter ;
                
            _jcfg._rdel_opts.
                top1() = _jjig._mesh_top1 ;
            _jcfg._rdel_opts.
                top2() = _jjig._mesh_top2 ;
            
            _jcfg._rdel_opts.
                rad2() = _jjig._mesh_rad2 ;
            _jcfg._rdel_opts.
                rad3() = _jjig._mesh_rad3 ;
            _jcfg._rdel_opts.
                off2() = _jjig._mesh_off2 ;
            _jcfg._rdel_opts.
                off3() = _jjig._mesh_off3 ;
            _jcfg._rdel_opts.
                snk2() = _jjig._mesh_snk2 ;
            _jcfg._rdel_opts.
                snk3() = _jjig._mesh_snk3 ;
            
            _jcfg._rdel_opts.
                eps1() = _jjig._mesh_eps1 ;
            _jcfg._rdel_opts.
                eps2() = _jjig._mesh_eps2 ;
            
            _jcfg._rdel_opts.
                vol3() = _jjig._mesh_vol3 ;
                
    /*------------------------------------- ITER keywords */
            _jcfg._iter_opts.
                iter() = _jjig._optm_iter ;
            
            _jcfg._iter_opts.
                qtol() = _jjig._optm_qtol ;
            _jcfg._iter_opts.
                qlim() = _jjig._optm_qlim ;
                
            _jcfg._iter_opts.
                tria() = _jjig._optm_tria ;
            _jcfg._iter_opts.
                dual() = _jjig._optm_dual ;
            _jcfg._iter_opts.
                div_() = _jjig._optm_div_ ;
            _jcfg._iter_opts.
                zip_() = _jjig._optm_zip_ ;
                
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * TEST-JCFG: test the JCFG data for validity.
    --------------------------------------------------------
     */
 
    template <
    typename      jlog_data
             >
    __normal_call iptr_type test_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog
        )
    {
        std::stringstream  _sstr;

    /*---------------------------------- read "ints" data */
        #define __testINTS(__tag, __var, __vlo, __vhi) \
            if ( (__var < __vlo) ||             \
                 (__var > __vhi) )              \
            {                                   \
                _errv = __invalid_argument ;    \
                _sstr.str("");                  \
                _sstr.clear();                  \
                _sstr <<                        \
                "  **input error: " __tag "="   \
                      << std::setw(+9)          \
                      << std::setfill(' ')      \
                      << __var << "\n";         \
                _jlog.push(_sstr.str()) ;       \
            }                                   \
        
    /*---------------------------------- test "real" data */
        #define __testREAL(__tag, __var, __vlo, __vhi) \
            if ( (__var < __vlo) ||             \
                 (__var > __vhi) )              \
            {                                   \
                _errv = __invalid_argument ;    \
                _sstr.str("");                  \
                _sstr.clear();                  \
                _sstr <<                        \
                "  **input error: " __tag "="   \
                      << std::scientific        \
                      << std::setprecision(2)   \
                      << __var << "\n";         \
                _jlog.push(_sstr.str()) ;       \
            }                                   \
        
    /*---------------------------------- warn "real" data */
        #define __warnREAL(__tag, __var, __vlo, __vhi) \
            if ( (__var < __vlo) ||             \
                 (__var > __vhi) )              \
            {                                   \
                _sstr.str("");                  \
                _sstr.clear();                  \
                _sstr <<                        \
              "  **input warning: " __tag "="   \
                      << std::scientific        \
                      << std::setprecision(2)   \
                      << __var << "\n";         \
                _jlog.push(_sstr.str()) ;       \
            }                                   \


        iptr_type _errv = __no_error ;

    /*---------------------------- test GEOM keywords */
        __testINTS("GEOM_SEED", 
            _jcfg._rdel_opts.seed(), 
            (iptr_type) + 0, 
        std::numeric_limits<iptr_type>::     max())

        __testREAL("GEOM_PHI1", 
            _jcfg._rdel_opts.phi1(), 
            (real_type)  0., 
            (real_type)180.)
        __testREAL("GEOM_PHI2", 
            _jcfg._rdel_opts.phi2(), 
            (real_type)  0., 
            (real_type)180.)
            
        __testREAL("GEOM_ETA1", 
            _jcfg._rdel_opts.eta1(), 
            (real_type)  0., 
            (real_type)180.)
        __testREAL("GEOM_ETA2", 
            _jcfg._rdel_opts.eta2(), 
            (real_type)  0., 
            (real_type)180.)

    /*---------------------------- test HFUN keywords */
        __testREAL("HFUN_HMAX", 
            _jcfg ._hfun_hmax , 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("HFUN_HMIN", 
            _jcfg ._hfun_hmin , 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        
    /*---------------------------- test MESH keywords */
        __testINTS("MESH_ITER", 
            _jcfg._rdel_opts.iter(), 
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())

        __testINTS("MESH_DIMS", 
            _jcfg._rdel_opts.dims(), 
            (iptr_type) + 1,
            (iptr_type) + 3)

        __testREAL("MESH_SIZ1", 
            _jcfg._rdel_opts.siz1(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH_SIZ2", 
            _jcfg._rdel_opts.siz2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH_SIZ3", 
            _jcfg._rdel_opts.siz3(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH_EPS1", 
            _jcfg._rdel_opts.eps1(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH_EPS2", 
            _jcfg._rdel_opts.eps2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH_RAD2", 
            _jcfg._rdel_opts.rad2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH_RAD3", 
            _jcfg._rdel_opts.rad3(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        
        __warnREAL("MESH_RAD2", 
            _jcfg._rdel_opts.rad2(), 
            (real_type)  1., 
        std::numeric_limits<real_type>::infinity())
        __warnREAL("MESH_RAD3", 
            _jcfg._rdel_opts.rad3(), 
            (real_type)  2., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH_OFF2", 
            _jcfg._rdel_opts.off2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH_OFF3", 
            _jcfg._rdel_opts.off3(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH_SNK2", 
            _jcfg._rdel_opts.snk2(),
            (real_type)  0., 
            (real_type)  1.)
        __testREAL("MESH_SNK3", 
            _jcfg._rdel_opts.snk3(),
            (real_type)  0., 
            (real_type)  1.)

        __testREAL("MESH_VOL3", 
            _jcfg._rdel_opts.vol3(),
            (real_type)  0., 
            (real_type)  1.)            
        __warnREAL("MESH_VOL3", 
            _jcfg._rdel_opts.vol3(),
            (real_type)  0., 
            (real_type)  .3)

    /*---------------------------- test OPTM keywords */
        __testINTS("OPTM_ITER", 
            _jcfg._iter_opts.iter(), 
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())
        
        __testREAL("OPTM_QTOL", 
            _jcfg._iter_opts.qtol(),
            (real_type)  0., 
            (real_type)  1.)
            
        __testREAL("OPTM_QLIM", 
            _jcfg._iter_opts.qlim(),
            (real_type)  0., 
            (real_type)  1.)
            

        #undef  __testINTS
        #undef  __testREAL
        #undef  __warnREAL
      
        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * ECHO-JCFG: print a summary of the JCFG data.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type echo_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog
        )
    {
        std::stringstream  _sstr ;

    /*---------------------------------- push "file" data */
        #define __dumpFILE( __tag, __var)   \
            _jlog.push("  " __tag " = " +   \
                _jcfg.__var + " \n" ) ;

    /*---------------------------------- push "ints" data */
        #define __pushIVAL(__var)           \
           (__var !=                        \
        std::numeric_limits<iptr_type>::max() ? \
            std::to_string(__var) : "MAXINT")
        
        #define __dumpINTS(__tag,__var)     \
            _sstr.str("");                  \
            _sstr.clear();                  \
            _sstr << "  " __tag " = "       \
                  << __pushIVAL(__var)      \
                  << " \n" ;                \
            _jlog.push(_sstr.str()) ;       \

    /*---------------------------------- push "bool" data */
        #define __pushBVAL(__var)           \
           ( __var ? "TRUE" : "FALSE" )
        
        #define __dumpBOOL(__tag,__var)     \
            _sstr.str("");                  \
            _sstr.clear();                  \
            _sstr << "  " __tag " = "       \
                  << __pushBVAL(__var)      \
                  << " \n" ;                \
            _jlog.push(_sstr.str()) ;       \
        
    /*---------------------------------- push "real" data */
        #define __dumpREAL(__tag,__var)     \
            _sstr.str("");                  \
            _sstr.clear();                  \
            _sstr << "  " __tag " = "       \
                  << std::scientific        \
                  << std::setprecision(2)   \
                  << __var                  \
                  << " \n" ;                \
            _jlog.push(_sstr.str()) ;       \
        
        
        iptr_type _errv  = __no_error;
        
        __dumpFILE(
            "GEOM_FILE", _geom_file)
        __dumpFILE(
            "MESH_FILE", _mesh_file)
        __dumpFILE(
            "HFUN_FILE", _hfun_file)

        _jlog.push("\n") ;

        if (_jcfg._verbosity > +0)
        {
    /*---------------------------- push GEOM keywords */
        __dumpINTS("GEOM_SEED", 
            _jcfg._rdel_opts.seed())
        
        __dumpREAL("GEOM_PHI1", 
            _jcfg._rdel_opts.phi1())
        __dumpREAL("GEOM_PHI2", 
            _jcfg._rdel_opts.phi2())
        
        __dumpREAL("GEOM_ETA1", 
            _jcfg._rdel_opts.eta1())
        __dumpREAL("GEOM_ETA2", 
            _jcfg._rdel_opts.eta2())
        
        __dumpBOOL("GEOM_FEAT", 
            _jcfg._rdel_opts.feat())

        _jlog.push("\n") ;

    /*---------------------------- push HFUN keywords */
        if(_jcfg._hfun_scal ==
         jcfg_data::hfun_scal::absolute)
        _jlog.push (
            "  HFUN_SCAL = ABSOLUTE \n") ;
        else
        if(_jcfg._hfun_scal ==
         jcfg_data::hfun_scal::relative)
        _jlog.push (
            "  HFUN_SCAL = RELATIVE \n") ;

        __dumpREAL(
            "HFUN_HMAX", _jcfg._hfun_hmax)
        __dumpREAL(
            "HFUN_HMIN", _jcfg._hfun_hmin)

        _jlog.push("\n") ;
        
    /*---------------------------- push MESH keywords */
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::delaunay)
        _jlog.push (
            "  MESH_KERN = DELAUNAY \n") ;
        else
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::delfront)
        _jlog.push (
            "  MESH_KERN = DELFRONT \n") ;

        __dumpBOOL("MESH_TOP1", 
            _jcfg._rdel_opts.top1())
        __dumpBOOL("MESH_TOP2", 
            _jcfg._rdel_opts.top2())

        __dumpINTS("MESH_ITER", 
            _jcfg._rdel_opts.iter())

        __dumpINTS("MESH_DIMS", 
            _jcfg._rdel_opts.dims())

        __dumpREAL("MESH_SIZ1", 
            _jcfg._rdel_opts.siz1())
        __dumpREAL("MESH_SIZ2", 
            _jcfg._rdel_opts.siz2())
        __dumpREAL("MESH_SIZ3", 
            _jcfg._rdel_opts.siz3())

        __dumpREAL("MESH_EPS1", 
            _jcfg._rdel_opts.eps1())
        __dumpREAL("MESH_EPS2", 
            _jcfg._rdel_opts.eps2())

        __dumpREAL("MESH_RAD2", 
            _jcfg._rdel_opts.rad2())
        __dumpREAL("MESH_RAD3", 
            _jcfg._rdel_opts.rad3())

        __dumpREAL("MESH_OFF2", 
            _jcfg._rdel_opts.off2())
        __dumpREAL("MESH_OFF3", 
            _jcfg._rdel_opts.off3())

        __dumpREAL("MESH_SNK2", 
            _jcfg._rdel_opts.snk2())
        __dumpREAL("MESH_SNK3", 
            _jcfg._rdel_opts.snk3())

        __dumpREAL("MESH_VOL3", 
            _jcfg._rdel_opts.vol3())

        _jlog.push("\n") ;
        
    /*---------------------------- push OPTM keywords */
        __dumpINTS("OPTM_ITER", 
            _jcfg._iter_opts.iter())

        __dumpREAL("OPTM_QTOL", 
            _jcfg._iter_opts.qtol())
        __dumpREAL("OPTM_QLIM", 
            _jcfg._iter_opts.qlim())
            
        __dumpBOOL("OPTM_ZIP_", 
            _jcfg._iter_opts.zip_())
        __dumpBOOL("OPTM_DIV_", 
            _jcfg._iter_opts.div_())
        __dumpBOOL("OPTM_TRIA", 
            _jcfg._iter_opts.tria())
        __dumpBOOL("OPTM_DUAL", 
            _jcfg._iter_opts.dual())

        _jlog.push("\n") ;
        
        }

        #undef  __dumpFILE
        #undef  __pushIVAL
        #undef  __pushBVAL
        #undef  __dumpINTS
        #undef  __dumpBOOL
        #undef  __dumpREAL

        return (  _errv) ;
    }
     
     
#   endif   //__JIG_LOAD__



