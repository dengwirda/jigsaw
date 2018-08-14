
    //
    // for cmd-marche:
    //
    // g++ -std=c++11 -pedantic -Wall -s -O3 -flto -D NDEBUG 
    // -static-libstdc++ marche.cpp -o marche64r
    //
    //
    // for lib-marche:
    //
    // g++ -std=c++11 -pedantic -Wall -O3 -flto -fPIC 
    // -D NDEBUG -static-libstdc++ marche.cpp -shared 
    // -o libmarche64r.so
    //

    /*
    --------------------------------------------------------
     *                                   `8
     *                                   |8                
     * `~-88e~-88e  /~~~8e  `8-~\ e88~~\ |8-~88e  e88~~8e  
     * |  888  888      88b |8   d888    |8  888 d888  88b 
     * |  888  888 e88~-888 |8   8888    |8  888 8888__888 
     * |  888  888 888  888 |8   Y888    |8  888 Y888    , 
     * ;  888  888 "88_-888 ;8    "88__/ ;8  888  "88___/ 
     *
     *
    --------------------------------------------------------
     * MARCHE: a 'fast-marching' eikonal equations solver.
    --------------------------------------------------------
     *
     * MARCHE release 0.9.0.x
     * Last updated: 04 August, 2018
     *
     * Copyright 2018
     * Darren Engwirda
     * darren.engwirda@columbia.edu
     * https://github.com/dengwirda
     *
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
     * University of Sydney, nor the National Aeronautics
     * and Space Administration warrant this code in any 
     * way whatsoever.  This code is provided "as-is" to be 
     * used at your own risk.
     *
    --------------------------------------------------------
     */
     
#   define __MRCHVSTR "MARCHE VERSION 0.9.0"


    /*---------------------------------- for i/o on files */
   
#   include <iomanip>
#   include <fstream>
#   include <sstream>
#   include <iostream>

    /*---------------------------------- for many things! */

#   include <cmath>

    /*---------------------------------- for ascii string */
    
#   include <string>

    /*---------------------------------- to do cpu timing */

#   define __use_timers

#   ifdef  __use_timers
#   include <chrono>
#   endif//__use_timers

    /*---------------------------------- MARCHE's backend */

#   include "libcpp/libbasic.hpp"
#   include "libcpp/libparse.hpp"

#   include "libcpp/meshfunc.hpp"

    /*---------------------------------- MARCHE's userlib */

//  define __lib_marche               // define makes lib

    extern  "C" 
    {
#   include "../inc/lib_marche.h"
    }

    typedef real_t real_type ;        // double-precision
    typedef indx_t iptr_type ;        // 32bit signed int

    /*---------------------------------- JIGSAW mesh kind */        
        
    struct jmsh_kind {   
        enum enum_data {
            euclidean_mesh                = +0 ,
            euclidean_grid ,
            euclidean_dual ,
            ellipsoid_mesh ,
            ellipsoid_grid ,
            ellipsoid_dual
        } ;
        } ;

    /*---------------------------------- run-time errors! */
    
    iptr_type static constexpr 
        __unknown_error                   = -1 ;
    
    iptr_type static constexpr 
        __no_error                        = +0 ;
    
    iptr_type static constexpr 
        __file_not_located                = +2 ;
    iptr_type static constexpr 
        __file_not_created                = +3 ;
    
    iptr_type static constexpr 
        __invalid_argument                = +4 ;
        
    /*
    --------------------------------------------------------
     * JCFG-DATA: MARCHE's config data.
    --------------------------------------------------------
     */

    class jcfg_data
        {
        public  :       
    /*------------------------------- "top-level" config. */
        std::string             _file_path ;
        std::string             _file_name ;

        std::string             _jcfg_file ;
        std::string             _geom_file ;
        std::string             _init_file ;
        std::string             _hfun_file ;
        std::string             _tria_file ;
        std::string             _mesh_file ;
        std::string             _bnds_file ;
        
        iptr_type               _verbosity = 0 ;
    
    /*--------------------------------- geom-bnd. kernels */
        struct bnds_pred {
            enum enum_data {
            nullkern ,    
            bnd_tria = JIGSAW_BNDS_TRIACELL,
            bnd_dual = JIGSAW_BNDS_DUALCELL
            } ;
            } ;

        bnds_pred::enum_data    
            _bnds_pred = bnds_pred::bnd_tria ;
        
    /*--------------------------------- mesh-gen. kernels */
        struct mesh_pred {
            enum enum_data {
            nullkern ,    
            delfront = JIGSAW_KERN_DELFRONT,
            delaunay = JIGSAW_KERN_DELAUNAY
            } ;
            } ;

        mesh_pred::enum_data    
            _mesh_pred = mesh_pred::delfront ;
        
    /*--------------------------------- H(x) fun. scaling */
        struct hfun_scal { 
            enum enum_data {
            nullscal ,
            absolute = JIGSAW_HFUN_ABSOLUTE,
            relative = JIGSAW_HFUN_RELATIVE
            } ;
            } ;
        
        hfun_scal::enum_data    
            _hfun_scal = hfun_scal::relative ;

        real_type _hfun_hmax = 
            (real_type) +2.00E-02 ;
        real_type _hfun_hmin = 
            (real_type) +0.00E+00 ;

    /*------------------------------- "low-level" config. */        
        typedef mesh::rdel_params <
                real_type ,
                iptr_type >      rdel_opts ;
        
        rdel_opts               _rdel_opts ;
        
        typedef mesh::iter_params <
                real_type ,
                iptr_type >      iter_opts ;
        
        iter_opts               _iter_opts ;

        } ;

    /*
    --------------------------------------------------------
     * aggregated MESH data containers.
    --------------------------------------------------------
     */
        
    class mesh_data         // holds the MESH obj.
        {
        public  :
        
        /*
        typedef mesh ::grad_constant_value_kd  <
                    iptr_type,
                    real_type>   constant_value_kd ;
        
        typedef mesh ::grad_mesh_euclidean_2d  <
                    real_type,
                    iptr_type>   euclidean_mesh_2d ;
                    
        typedef mesh ::grad_mesh_euclidean_3d  <
                    real_type,
                    iptr_type>   euclidean_mesh_3d ;
                    
        typedef mesh ::grad_grid_euclidean_2d  <
                    real_type,
                    iptr_type>   euclidean_grid_2d ;
                    
        typedef mesh ::grad_grid_euclidean_3d  <
                    real_type,
                    iptr_type>   euclidean_grid_3d ;
          
        typedef mesh ::grad_grid_ellipsoid_3d  <
                    iptr_type,
                    real_type>   ellipsoid_grid_3d ;
         
        iptr_type               _ndim = +0;
        
        jmsh_kind ::
        enum_data               _kind ;
      
        constant_value_kd       _constant_value_kd ;
      
        euclidean_mesh_2d       _euclidean_mesh_2d ;
        euclidean_mesh_3d       _euclidean_mesh_3d ;
        
        euclidean_grid_2d       _euclidean_grid_2d ;
        euclidean_grid_3d       _euclidean_grid_3d ;
        
        ellipsoid_grid_3d       _ellipsoid_grid_3d ; 
         */
         
        public  :
        
    /*------------------------- helper: init. everything! */
    
        __normal_call void_type init_mesh (
            jcfg_data &_jcfg
            )
        {
            __unreferenced(_jcfg) ;
    
           /*    
            this->
           _constant_value_kd.init() ;
        
            this->
           _euclidean_mesh_2d.init() ;
            this->
           _euclidean_mesh_3d.init() ;
            
            this->
           _euclidean_grid_2d.init() ;
            this->
           _euclidean_grid_3d.init() ;
            this->
           _ellipsoid_grid_3d.init() ;
            */
        }
        
        } ;
      
    /*
    --------------------------------------------------------
     * MLOG-DATA: log-file writer for MARCHE.
    --------------------------------------------------------
     */

    class mlog_text
        {
    /*-------------------------- a "real" log-file writer */
            std::ofstream    _file ;   
        
        public  :
        
        __inline_call  mlog_text (
            jcfg_data const& _jcfg
            )
        {
            if (_jcfg._file_path.length() == +0)
            {
                this->_file.open (
                _jcfg._file_name + ".log",
        std::ofstream::out|std::ofstream::trunc) ;
            }
            else
            {
                this->_file.open (
                _jcfg._file_path + "/" +
                _jcfg._file_name + ".log",
        std::ofstream::out|std::ofstream::trunc) ;
            }
        }
        
        __inline_call ~mlog_text (
            )
        {   this->_file.close () ;
        }

    /*-------------------------- echo onto logfile/stdout */

        template <
        typename      data_type
                 >
        __inline_call void_type push (
            data_type const&_data
            )
        {
            std :: cout <<  _data ;
            this->_file <<  _data ;
        }
        
        } ;

    class mlog_null
        {
        public  :
    /*-------------------------- a "null" log-file writer */
        
            iptr_type  _verbosity ;
            
        public  :
    
        __inline_call  mlog_null (
            jcfg_data const& _jcfg
            )
        {
    /*-------------------------- def. no: for lib_marche! */
            this->_verbosity = 
                _jcfg._verbosity ;
        }
    
        template <
        typename      data_type
                 >
        __inline_call void_type push (
            data_type const&_data
            )
        {
    /*-------------------------- def. no: for lib_marche! */
            if (this->_verbosity > +0)
            {
            std::cout << _data ;
            }
        }
    
        } ;
        
    /*
    --------------------------------------------------------
     * READ-JCFG: load JCFG data from file.
    --------------------------------------------------------
     */
    
    #   include "jig_load.hpp"
    
    
    
    
    
    
    
        
        
        
    /*
    --------------------------------------------------------
     * MLOG-HEAD: write header for *.LOG file.
    --------------------------------------------------------
     */
 
    template <
    typename      mlog_data
             >
    __normal_call void_type mlog_head (
        mlog_data &_mlog
        )
    {
    /*-- NB: silliness re. escape sequences */
        _mlog.push (
    " \n"
    "#------------------------------------------------------------\n"
    "#                                   `,                  \n"
    "#                                   |8                  \n"
    "# `~-88e~-88e  /~~~8e  `8-~\\ e88~~\\ |8-~88e  e88~~8e  \n"
    "# |  888  888      88b |8   d888    |8  888 d888  88b   \n"
    "# |  888  888 e88~-888 |8   8888    |8  888 8888__888   \n"
    "# |  888  888 888  888 |8   Y888    |8  888 Y888    ,   \n"
    "# ;  888  888 \"88_-888 Y8    \"88__/ Y8  888  \"88___/ \n"
    "#\n"
    "#\n"
    "#------------------------------------------------------------\n"
    "# MARCHE: a 'fast-marching' eikonal equations solver.   \n"
    "#------------------------------------------------------------\n"
    " \n"
    "  " __MRCHVSTR "\n\n"
        ) ;      
    }
    
        #define __mloglndv              \
    "#------------------------------------------------------------\n"
 
    /*
    --------------------------------------------------------
     * TIME-SPAN: elapsed sec. between markers.
    --------------------------------------------------------
     */

#   ifdef  __use_timers

    __inline_call double      time_span (
        typename std::
        chrono::high_resolution_clock
            ::time_point const& _ttic ,
        typename std::
        chrono::high_resolution_clock
            ::time_point const& _ttoc
        )
    {
        return (double)(
        std::chrono::duration_cast<
        std::chrono::microseconds >
        (_ttoc-_ttic).count()) / 1.E+06 ;
    }
    
#   endif//__use_timers
    
    /*
    --------------------------------------------------------
     * DUMP-TIME: "dump" time-span to a string.
    --------------------------------------------------------
     */
    
#   ifdef  __use_timers

    __inline_call std::string dump_time (
        typename std::
        chrono::high_resolution_clock
            ::time_point const& _ttic ,
        typename std::
        chrono::high_resolution_clock
            ::time_point const& _ttoc
        )
    {
        std::stringstream _sstr;
        _sstr << "  Done. ("
              << std::scientific
              << std::setprecision(2)
              << time_span(_ttic, _ttoc )
              << "sec)\n\n" ;

        return _sstr.str () ;     
    }

#   endif//__use_timers

    /*
    --------------------------------------------------------
     * MAIN: jumping-off point for MARCHE!
    --------------------------------------------------------
     */
 
#   ifdef __lib_marche

#   include "liblib/init_jig_t.hpp"
#   include "liblib/init_msh_t.hpp"

#   include "liblib/load_jig_t.hpp"
#   include "liblib/load_msh_t.hpp"

#   include "liblib/save_jig_t.hpp"
#   include "liblib/save_msh_t.hpp"
 
    __normal_call
        iptr_type marche_iter_mesh (    // lib-marche
        jigsaw_jig_t *_jjig ,
        jigsaw_msh_t *_mmsh
        )
    {
        iptr_type _retv = +0;
    
        mesh_data _dest ;               // MESH data
        jcfg_data _mcfg ;
    
#       ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttic ;
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time;
        
        __unreferenced(_time) ;
#       endif//__use_timers
    
    /*--------------------------------- setup *.JLOG data */
        if (_jjig != nullptr )
        {
            _jcfg._verbosity = _jjig->_verbosity ;
        }
        mlog_null _mlog(_jcfg) ;
        mlog_head(_mlog) ;
        
    /*--------------------------------- parse *.JCFG data */
        if (_jjig != nullptr )
        {
            _mlog.push (
                "  Reading CFG. data...\n\n" ) ;

#           ifdef  __use_timers        
            _ttic   = _time.now();
#           endif//__use_timers
        
            if ((_retv = copy_jcfg (
                 _jcfg, 
                 _mlog,*_jjig)) != __no_error)
            {
                return  _retv ;
            }    
            
            if ((_retv = test_jcfg (
                 _jcfg, _mlog)) != __no_error)
            {
                return  _retv ;
            }

            _mlog.push (
                "  CFG. data summary...\n\n" ) ;

            if ((_retv = echo_jcfg (
                 _jcfg, _mlog)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers            
            _ttoc   = _time.now();
            _mlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
  
  
        //!! TODO
 
 
  
  
    /*-------------------------- success, if we got here! */

        return ( _retv ) ;
    }
    
#   else
 
    __normal_call iptr_type   main (    // cmd-marche
        int           _argc , 
        char        **_argv
        )
    {
        mesh_data _mesh ;               // MESH data
        
#       ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttic ;
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time;
        
        __unreferenced(_time) ;
#       endif//__use_timers
        
    /*-------------------------- find *.JFCG file in args */    
        iptr_type _retv = -1  ;
        jcfg_data _jcfg ;
        for (; _argc-- != +0; )
        {
            std::string _ssrc(_argv[_argc]) ;
            
            std::string _path ;
            std::string _name ;
            std::string _fext ;
            file_part ( _ssrc ,
                _path , _name , _fext)  ;

            if (_ssrc.find("-whoami") == 0)
            {
                _retv = -2 ;
                
                std::cout << __MRCHVSTR ;
                std::cout <<  std::endl ;
                
                break ;
            }

            if (_fext.find("jig") == 0)
            {
                _retv = +0 ;

                _jcfg._jcfg_file =_ssrc ;

                _jcfg._file_path =_path ;
                _jcfg._file_name =_name ;

                break ;
            }
        }
        if (_retv != +0) return ( _retv ) ;

    /*--------------------------------- setup *.JLOG file */
        mlog_text _mlog(_jcfg) ;
        mlog_head(_mlog) ;
      
        if(!_jcfg._jcfg_file.empty())
        {
    /*--------------------------------- parse *.JCFG file */
            _mlog.push (
                "  Reading CFG. file...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers
        
            if ((_retv = read_jcfg (
                 _jcfg, _mlog)) != __no_error)
            {
                return  _retv ;
            }           
            
            if ((_retv = test_jcfg (
                 _jcfg, _mlog)) != __no_error)
            {
                return  _retv ;
            }

            _mlog.push (
                "  CFG. data summary...\n\n" ) ;

            if ((_retv = echo_jcfg (
                 _jcfg, _mlog)) != __no_error)
            {
                return  _retv ;
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _mlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        
        //!! TODO
        
        
        
        
    /*-------------------------- success, if we got here! */

        return ( _retv ) ;
    }

        
#   endif   //__lib_marche 
    
    
    
    
