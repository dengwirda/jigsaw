
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

#   include "jig_read.hpp"

#   ifndef __JIG_LOAD__
#   define __JIG_LOAD__

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
        class jcfg_loader: 
            public jcfg_reader_base
        {
        public  :
        jcfg_data             *_jjig;
        
        public  :
    /*------------------------------------- create loader */
        __normal_call jcfg_loader (
            jcfg_data *_jsrc = nullptr
            ) : _jjig( _jsrc ) {}
        
    /*------------------------------------- MISC keywords */
        __normal_call void_type push_verbosity (
            std::int32_t  _verb
            ) 
        { 
            this->_jjig->
           _rdel_opts.verb() = _verb;
            this->_jjig->
           _iter_opts.verb() = _verb;       
        }          
        
    /*------------------------------------- GEOM keywords */       
        __normal_call void_type push_geom_file (
            std::string   _file
            ) 
        { 
            this->
           _jjig->_geom_file = _file;
        }
        __normal_call void_type push_geom_seed (
            std::int32_t  _seed
            ) 
        {
            this->_jjig->
           _rdel_opts.seed() = _seed; 
        }
        __normal_call void_type push_geom_feat (
            bool          _feat
            ) 
        {
            this->_jjig->
           _rdel_opts.feat() = _feat; 
        }
        __normal_call void_type push_geom_phi1 (
            double        _phi1
            ) 
        {
            this->_jjig->
           _rdel_opts.phi1() = _phi1; 
        }
        __normal_call void_type push_geom_phi2 (
            double        _phi2
            ) 
        {
            this->_jjig->
           _rdel_opts.phi2() = _phi2; 
        }
        __normal_call void_type push_geom_eta1 (
            double        _eta1
            ) 
        {
            this->_jjig->
           _rdel_opts.eta1() = _eta1; 
        }
        __normal_call void_type push_geom_eta2 (
            double        _eta2
            ) 
        {
            this->_jjig->
           _rdel_opts.eta2() = _eta2; 
        }
        
    /*------------------------------------- HFUN keywords */
        __normal_call void_type push_hfun_file (
            std::string   _file
            ) 
        { 
            this->
           _jjig->_hfun_file = _file;
        }
        __normal_call void_type push_hfun_scal (
            std::int32_t  _scal
            ) 
        {
            this->_jjig->_hfun_scal = 
        (jcfg_data::hfun_scal::enum_data)_scal ; 
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
    
    /*------------------------------------- INIT keywords */    
        __normal_call void_type push_init_file (
            std::string   _file
            ) 
        { 
            this->
           _jjig->_init_file = _file;
        }
        
    /*------------------------------------- KERN keywords */
        __normal_call void_type push_bnds_kern (
            std::int32_t  _kern
            ) 
        { 
            this->_jjig->_bnds_pred = 
        (jcfg_data::bnds_pred::enum_data)_kern ;
        }
        
        __normal_call void_type push_mesh_kern (
            std::int32_t  _kern
            ) 
        {
            this->_jjig->_mesh_pred = 
        (jcfg_data::mesh_pred::enum_data)_kern ;
        }
        
    /*------------------------------------- MESH keywords */
        __normal_call void_type push_mesh_file (
            std::string   _file
            ) 
        { 
            this->
           _jjig->_mesh_file = _file;
        }
        __normal_call void_type push_tria_file (
            std::string   _file
            ) 
        { 
            this->
           _jjig->_tria_file = _file;
        }
        __normal_call void_type push_bnds_file (
            std::string   _file
            ) 
        { 
            this->
           _jjig->_bnds_file = _file;
        }
        
        __normal_call void_type push_mesh_dims (
            std::int32_t  _dims
            ) 
        { 
            this->_jjig->
           _rdel_opts.dims() = _dims;
        }
        __normal_call void_type push_mesh_iter (
            std::int32_t  _iter
            ) 
        { 
            this->_jjig->
           _rdel_opts.iter() = _iter;
        }
        __normal_call void_type push_mesh_siz1 (
            double        _siz1
            ) 
        { 
            this->_jjig->
           _rdel_opts.siz1() = _siz1;
        }
        __normal_call void_type push_mesh_siz2 (
            double        _siz2
            ) 
        { 
            this->_jjig->
           _rdel_opts.siz2() = _siz2;
        }
        __normal_call void_type push_mesh_siz3 (
            double        _siz3
            ) 
        {
            this->_jjig->
           _rdel_opts.siz3() = _siz3; 
        }      
        __normal_call void_type push_mesh_top1 (
            bool          _top1
            ) 
        { 
            this->_jjig->
           _rdel_opts.top1() = _top1;
        }
        __normal_call void_type push_mesh_top2 (
            bool          _top2
            ) 
        {
            this->_jjig->
           _rdel_opts.top2() = _top2; 
        }
        __normal_call void_type push_mesh_rad2 (
            double        _rad2
            ) 
        {
            this->_jjig->
           _rdel_opts.rad2() = _rad2; 
        }
        __normal_call void_type push_mesh_rad3 (
            double        _rad3
            ) 
        {
            this->_jjig->
           _rdel_opts.rad3() = _rad3; 
        }
        __normal_call void_type push_mesh_off2 (
            double        _off2
            ) 
        { 
            this->_jjig->
           _rdel_opts.off2() = _off2;
        }
        __normal_call void_type push_mesh_off3 (
            double        _off3
            ) 
        {
            this->_jjig->
           _rdel_opts.off3() = _off3; 
        }
        __normal_call void_type push_mesh_snk2 (
            double        _snk2
            ) 
        {
            this->_jjig->
           _rdel_opts.snk2() = _snk2; 
        }
        __normal_call void_type push_mesh_snk3 (
            double        _snk3
            ) 
        {
            this->_jjig->
           _rdel_opts.snk3() = _snk3; 
        }
        __normal_call void_type push_mesh_eps1 (
            double        _eps1
            ) 
        {
            this->_jjig->
           _rdel_opts.eps1() = _eps1; 
        }
        __normal_call void_type push_mesh_eps2 (
            double        _eps2
            ) 
        {
            this->_jjig->
           _rdel_opts.eps2() = _eps2; 
        }
        __normal_call void_type push_mesh_vol3 (
            double        _vol3
            ) 
        {
            this->_jjig->
           _rdel_opts.vol3() = _vol3; 
        }
        
    /*------------------------------------- OPTM keywords */
        __normal_call void_type push_optm_iter (
            std::int32_t  _iter
            ) 
        {
            this->_jjig->
           _iter_opts.iter() = _iter;  
        }
        __normal_call void_type push_optm_qtol (
            double        _qtol
            ) 
        {
            this->_jjig->
           _iter_opts.qtol() = _qtol; 
        }
        __normal_call void_type push_optm_qlim (
            double        _qlim
            ) 
        {
            this->_jjig->
           _iter_opts.qlim() = _qlim; 
        }
        __normal_call void_type push_optm_tria (
            bool          _flag
            ) 
        {
            this->_jjig->
           _iter_opts.tria() = _flag; 
        }
        __normal_call void_type push_optm_dual (
            bool          _flag
            ) 
        {
            this->_jjig->
           _iter_opts.dual() = _flag; 
        }
        __normal_call void_type push_optm_div_ (
            bool          _flag
            ) 
        {
            this->_jjig->
           _iter_opts.div_() = _flag; 
        }
        __normal_call void_type push_optm_zip_ (
            bool          _flag
            ) 
        {
            this->_jjig->
           _iter_opts.zip_() = _flag; 
        }
        
        } ;
    
    /*---------------------------------- parse JCFG. file */ 
        iptr_type _errv  = __no_error ;
        
        try
        {
            jcfg_reader   _read;
            std::ifstream _file; 
            _file. open(
            _jcfg._jcfg_file, std::ifstream::in);

            if (_file.is_open())
            {
                _read.read_file (
                _file, jcfg_loader(&_jcfg));
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
            "**parse error: " + * _iter + "\n") ;
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
    
    /*------------------------------------- BNDS keywords */
            
            if (_jjig._bnds_kern == 
                    JIGSAW_BNDS_TRIACELL)
            _jcfg._bnds_pred = 
                jcfg_data::bnds_pred::bnd_tria ;
            else
            if (_jjig._bnds_kern ==
                    JIGSAW_BNDS_DUALCELL)
            _jcfg._bnds_pred = 
                jcfg_data::bnds_pred::bnd_dual ;
    
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
            if (_jjig._mesh_kern == 
                    JIGSAW_KERN_DELFRONT)
            _jcfg._mesh_pred = 
                jcfg_data::mesh_pred::delfront ;
            else
            if (_jjig._mesh_kern ==
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
        __testINTS("GEOM-SEED", 
            _jcfg._rdel_opts.seed(), 
            (iptr_type) + 0, 
        std::numeric_limits<iptr_type>::     max())

        __testREAL("GEOM-PHI1", 
            _jcfg._rdel_opts.phi1(), 
            (real_type)  0., 
            (real_type)180.)
        __testREAL("GEOM-PHI2", 
            _jcfg._rdel_opts.phi2(), 
            (real_type)  0., 
            (real_type)180.)
            
        __testREAL("GEOM-ETA1", 
            _jcfg._rdel_opts.eta1(), 
            (real_type)  0., 
            (real_type)180.)
        __testREAL("GEOM-ETA2", 
            _jcfg._rdel_opts.eta2(), 
            (real_type)  0., 
            (real_type)180.)

    /*---------------------------- test HFUN keywords */
        __testREAL("HFUN-HMAX", 
            _jcfg ._hfun_hmax , 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("HFUN-HMIN", 
            _jcfg ._hfun_hmin , 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        
    /*---------------------------- test MESH keywords */
        __testINTS("MESH-ITER", 
            _jcfg._rdel_opts.iter(), 
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())

        __testINTS("MESH-DIMS", 
            _jcfg._rdel_opts.dims(), 
            (iptr_type) + 1,
            (iptr_type) + 3)

        __testREAL("MESH-SIZ1", 
            _jcfg._rdel_opts.siz1(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-SIZ2", 
            _jcfg._rdel_opts.siz2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-SIZ3", 
            _jcfg._rdel_opts.siz3(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-EPS1", 
            _jcfg._rdel_opts.eps1(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-EPS2", 
            _jcfg._rdel_opts.eps2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-RAD2", 
            _jcfg._rdel_opts.rad2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-RAD3", 
            _jcfg._rdel_opts.rad3(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        
        __warnREAL("MESH-RAD2", 
            _jcfg._rdel_opts.rad2(), 
            (real_type)  1., 
        std::numeric_limits<real_type>::infinity())
        __warnREAL("MESH-RAD3", 
            _jcfg._rdel_opts.rad3(), 
            (real_type)  2., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-OFF2", 
            _jcfg._rdel_opts.off2(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-OFF3", 
            _jcfg._rdel_opts.off3(), 
            (real_type)  0., 
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-SNK2", 
            _jcfg._rdel_opts.snk2(),
            (real_type)  0., 
            (real_type)  1.)
        __testREAL("MESH-SNK3", 
            _jcfg._rdel_opts.snk3(),
            (real_type)  0., 
            (real_type)  1.)

        __testREAL("MESH-VOL3", 
            _jcfg._rdel_opts.vol3(),
            (real_type)  0., 
            (real_type)  1.)            
        __warnREAL("MESH-VOL3", 
            _jcfg._rdel_opts.vol3(),
            (real_type)  0., 
            (real_type)  .3)

    /*---------------------------- test OPTM keywords */
        __testINTS("OPTM-ITER", 
            _jcfg._iter_opts.iter(), 
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())
        
        __testREAL("OPTM-QTOL", 
            _jcfg._iter_opts.qtol(),
            (real_type)  0., 
            (real_type)  1.)
            
        __testREAL("OPTM-QLIM", 
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
            "GEOM-FILE", _geom_file)
        __dumpFILE(
            "MESH-FILE", _mesh_file)
        __dumpFILE(
            "HFUN-FILE", _hfun_file)
        __dumpFILE(
            "INIT-FILE", _init_file)
        __dumpFILE(
            "TRIA-FILE", _tria_file)
        __dumpFILE(
            "BNDS-FILE", _bnds_file)

        _jlog.push("\n") ;

        if (_jcfg._verbosity > +0)
        {
    /*---------------------------- push GEOM keywords */
        __dumpINTS("GEOM-SEED", 
            _jcfg._rdel_opts.seed())
        
        __dumpREAL("GEOM-PHI1", 
            _jcfg._rdel_opts.phi1())
        __dumpREAL("GEOM-PHI2", 
            _jcfg._rdel_opts.phi2())
        
        __dumpREAL("GEOM-ETA1", 
            _jcfg._rdel_opts.eta1())
        __dumpREAL("GEOM-ETA2", 
            _jcfg._rdel_opts.eta2())
        
        __dumpBOOL("GEOM-FEAT", 
            _jcfg._rdel_opts.feat())

        _jlog.push("\n") ;

    /*---------------------------- push HFUN keywords */
        if(_jcfg._hfun_scal ==
         jcfg_data::hfun_scal::absolute)
        _jlog.push (
            "  HFUN-SCAL = ABSOLUTE \n") ;
        else
        if(_jcfg._hfun_scal ==
         jcfg_data::hfun_scal::relative)
        _jlog.push (
            "  HFUN-SCAL = RELATIVE \n") ;

        __dumpREAL(
            "HFUN-HMAX", _jcfg._hfun_hmax)
        __dumpREAL(
            "HFUN-HMIN", _jcfg._hfun_hmin)

        _jlog.push("\n") ;
        
    /*---------------------------- push MESH keywords */
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::delaunay)
        _jlog.push (
            "  MESH-KERN = DELAUNAY \n") ;
        else
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::delfront)
        _jlog.push (
            "  MESH-KERN = DELFRONT \n") ;
            
        if(_jcfg._bnds_pred ==
         jcfg_data::bnds_pred::bnd_tria)
        _jlog.push (
            "  BNDS-KERN = BND-TRIA \n") ;
        else
        if(_jcfg._bnds_pred ==
         jcfg_data::bnds_pred::bnd_dual)
        _jlog.push (
            "  BNDS-KERN = BND-DUAL \n") ;

        __dumpBOOL("MESH-TOP1", 
            _jcfg._rdel_opts.top1())
        __dumpBOOL("MESH-TOP2", 
            _jcfg._rdel_opts.top2())

        __dumpINTS("MESH-ITER", 
            _jcfg._rdel_opts.iter())

        __dumpINTS("MESH-DIMS", 
            _jcfg._rdel_opts.dims())

        __dumpREAL("MESH-SIZ1", 
            _jcfg._rdel_opts.siz1())
        __dumpREAL("MESH-SIZ2", 
            _jcfg._rdel_opts.siz2())
        __dumpREAL("MESH-SIZ3", 
            _jcfg._rdel_opts.siz3())

        __dumpREAL("MESH-EPS1", 
            _jcfg._rdel_opts.eps1())
        __dumpREAL("MESH-EPS2", 
            _jcfg._rdel_opts.eps2())

        __dumpREAL("MESH-RAD2", 
            _jcfg._rdel_opts.rad2())
        __dumpREAL("MESH-RAD3", 
            _jcfg._rdel_opts.rad3())

        __dumpREAL("MESH-OFF2", 
            _jcfg._rdel_opts.off2())
        __dumpREAL("MESH-OFF3", 
            _jcfg._rdel_opts.off3())

        __dumpREAL("MESH-SNK2", 
            _jcfg._rdel_opts.snk2())
        __dumpREAL("MESH-SNK3", 
            _jcfg._rdel_opts.snk3())

        __dumpREAL("MESH-VOL3", 
            _jcfg._rdel_opts.vol3())

        _jlog.push("\n") ;
        
    /*---------------------------- push OPTM keywords */
        __dumpINTS("OPTM-ITER", 
            _jcfg._iter_opts.iter())

        __dumpREAL("OPTM-QTOL", 
            _jcfg._iter_opts.qtol())
        __dumpREAL("OPTM-QLIM", 
            _jcfg._iter_opts.qlim())
            
        __dumpBOOL("OPTM-ZIP_", 
            _jcfg._iter_opts.zip_())
        __dumpBOOL("OPTM-DIV_", 
            _jcfg._iter_opts.div_())
        __dumpBOOL("OPTM-TRIA", 
            _jcfg._iter_opts.tria())
        __dumpBOOL("OPTM-DUAL", 
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



