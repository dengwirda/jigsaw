
    //
    // for cmd-jigsaw:
    //
    // g++ -std=c++11 -pedantic -Wall -s -O3 -flto -D NDEBUG 
    // -static-libstdc++ jigsaw.cpp -o jigsaw64r
    //
    //
    // for lib-jigsaw:
    //
    // g++ -std=c++11 -pedantic -Wall -O3 -flto -fPIC 
    // -D NDEBUG -static-libstdc++ jigsaw.cpp -shared 
    // -o libjigsaw64r.so
    //

    /*
    --------------------------------------------------------
     *
     *   ,o, ,o,       /                                
     *    `   `  e88~88e  d88~\   /~~~8e Y88b    e    /
     *   888 888 88   88 C888         88b Y88b  d8b  / 
     *   888 888 "8b_d8"  Y88b   e88~-888  Y888/Y88b/  
     *   888 888  /        888D C88   888   Y8/  Y8/   
     *   88P 888 Cb      \_88P   "8b_-888    Y    Y     
     * \_8"       Y8""8D                                
     *
    --------------------------------------------------------
     * JIGSAW: an unstructured mesh generation package.
    --------------------------------------------------------
     *
     * JIGSAW release 0.9.7.x
     * Last updated: 13 August, 2018
     *
     * Copyright 2013 -- 2018
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
     *
     * JIGSAW is a collection of unstructured triangle- and
     * tetrahedron-based meshing algorithms, designed to 
     * produce very high quality Delaunay-based grids for 
     * computational simulation. JIGSAW includes both 
     * Delaunay 'refinement' based algorithms for the 
     * construction of new meshes, as well as optimisation 
     * driven methods for the 'improvement' of existing 
     * grids. JIGSAW supports both two- and 
     * three-dimensional operations, catering to a variety 
     * of planar, surface and volumetric configurations.
     *
     * JIGSAW has grown out of my Ph.D. research, in which
     * I explored early versions of both the refinement and
     * optimisation-based algorithms described here:
     *
     * D. Engwirda, (2014): "Locally-optimal Delaunay-
     * refinement and optimisation-based mesh generation", 
     * Ph.D. Thesis, School of Mathematics and Statistics, 
     * Univ. of Sydney. 
     * http://hdl.handle.net/2123/13148
     *
    --------------------------------------------------------
     *
     * JIGSAW builds on a variety of "standard" algorithms,
     * building upon technqiues described in, for instance:
     *
     * S. Cheng, T. Dey & J. Shewchuk, (2012): "Delaunay 
     * mesh generation", CRC Press.
     *
     * Additional references are provided inline throughout 
     * the JIGSAW src.
     *
    --------------------------------------------------------
     */


#   define __JGSWVSTR "JIGSAW VERSION 0.9.7"


    /*---------------------------------- for i/o on files */
   
#   include <iomanip>
#   include <fstream>
#   include <sstream>
#   include <iostream>

    /*---------------------------------- for ascii string */
    
#   include <string>

    /*---------------------------------- for many things! */

#   include <cmath>

    /*---------------------------------- to do cpu timing */

#   define __use_timers

#   ifdef  __use_timers
#   include <chrono>
#   endif//__use_timers

    /*---------------------------------- JIGSAW's backend */

#   include "libcpp/libbasic.hpp"
#   include "libcpp/libparse.hpp"

#   include "libcpp/rdelmesh.hpp"
#   include "libcpp/itermesh.hpp"

    /*---------------------------------- JIGSAW's userlib */

//  define __lib_jigsaw               // define makes lib

    extern  "C" 
    {
#   include "../inc/lib_jigsaw.h"
    }

    typedef real_t real_type ;        // double-precision
    typedef indx_t iptr_type ;        // 32bit signed int

    /*---------------------------------- JIGSAW mesh kind */        
        
    struct jmsh_kind {   
        enum enum_data {
            euclidean_mesh      = +0 ,
            euclidean_grid ,
            euclidean_dual ,
            ellipsoid_mesh ,
            ellipsoid_grid ,
            ellipsoid_dual
        } ;
        } ;

    /*---------------------------------- run-time errors! */
    
    iptr_type static constexpr 
        __unknown_error         = -1 ;
    
    iptr_type static constexpr 
        __no_error              = +0 ;
    
    iptr_type static constexpr 
        __file_not_located      = +2 ;
    iptr_type static constexpr 
        __file_not_created      = +3 ;
    
    iptr_type static constexpr 
        __invalid_argument      = +4 ;

 
    /*
    --------------------------------------------------------
     * JCFG-DATA: JIGSAW's config data.
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
     * aggregated GEOM data containers.
    --------------------------------------------------------
     */

    class geom_data         // holds the GEOM obj.
        {
        public  :
        
        typedef mesh ::geom_mesh_euclidean_2d  <
                    real_type,
                    iptr_type>   euclidean_mesh_2d ;
                    
        typedef mesh ::geom_mesh_euclidean_3d  <
                    real_type,
                    iptr_type>   euclidean_mesh_3d ;
        
        typedef mesh ::geom_mesh_ellipsoid_3d  <
                    real_type,
                    iptr_type>   ellipsoid_mesh_3d ;
        
        iptr_type               _ndim = +0;

        jmsh_kind ::
        enum_data               _kind ;

        euclidean_mesh_2d       _euclidean_mesh_2d ;
        euclidean_mesh_3d       _euclidean_mesh_3d ;
        
        ellipsoid_mesh_3d       _ellipsoid_mesh_3d ;
        
        public  :
        
    /*------------------------- helper: init. everything! */
        
        __normal_call void_type init_geom (
            jcfg_data &_jcfg
            )
        {
            this->_euclidean_mesh_2d.
                _tria.make_ptrs() ;
            this->_euclidean_mesh_2d.
                init_geom(_jcfg._rdel_opts) ;
                
            this->_euclidean_mesh_3d.
                _tria.make_ptrs() ;
            this->_euclidean_mesh_3d.
                init_geom(_jcfg._rdel_opts) ;
                
            this->_ellipsoid_mesh_3d.
                init_geom(_jcfg._rdel_opts) ;
        }
        
        } ;

    /*
    --------------------------------------------------------
     * aggregated HFUN data containers.
    --------------------------------------------------------
     */

    class hfun_data         // holds the HFUN obj.
        {
        public  :
        
        typedef mesh ::hfun_constant_value_kd  <
                    iptr_type,
                    real_type>   constant_value_kd ;
        
        typedef mesh ::hfun_mesh_euclidean_2d  <
                    real_type,
                    iptr_type>   euclidean_mesh_2d ;
                    
        typedef mesh ::hfun_mesh_euclidean_3d  <
                    real_type,
                    iptr_type>   euclidean_mesh_3d ;
                    
        typedef mesh ::hfun_grid_euclidean_2d  <
                    real_type,
                    iptr_type>   euclidean_grid_2d ;
                    
        typedef mesh ::hfun_grid_euclidean_3d  <
                    real_type,
                    iptr_type>   euclidean_grid_3d ;
          
        typedef mesh ::hfun_grid_ellipsoid_3d  <
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
        
        public  :
        
    /*------------------------- helper: init. everything! */
    
        __normal_call void_type init_hfun (
            jcfg_data &_jcfg
            )
        {
            __unreferenced(_jcfg) ;
        
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
        }
        
        } ;
        
    /*
    --------------------------------------------------------
     * aggregated RDEL data containers.
    --------------------------------------------------------
     */
    
    class rdel_data         // holds the restrict-del obj.
        {
        public  :
        
        typedef mesh::rdel_complex_2d <
                real_type ,
                iptr_type >      euclidean_rdel_2d ;
                    
        typedef mesh::rdel_complex_3d <
                real_type ,
                iptr_type >      euclidean_rdel_3d ;
           
        iptr_type               _ndim = +0;

        jmsh_kind ::
        enum_data               _kind ;
                
        euclidean_rdel_2d       _euclidean_rdel_2d ;
        euclidean_rdel_3d       _euclidean_rdel_3d ;
        
        euclidean_rdel_2d       _euclidean_rvor_2d ;
        euclidean_rdel_3d       _euclidean_rvor_3d ;
        
        } ;
    
    /*
    --------------------------------------------------------
     * aggregated MESH data containers.
    --------------------------------------------------------
     */
        
    class mesh_data         // holds the tria-complex obj.
        {
        public  :
        
        typedef mesh::iter_mesh_euclidean_2d <
                real_type ,
                iptr_type >      euclidean_mesh_2d ;
                    
        typedef mesh::iter_mesh_euclidean_3d <
                real_type ,
                iptr_type >      euclidean_mesh_3d ;
           
        iptr_type               _ndim = +0;

        jmsh_kind ::
        enum_data               _kind ;
                
        euclidean_mesh_2d       _euclidean_mesh_2d ;
        euclidean_mesh_3d       _euclidean_mesh_3d ;
        
        } ;

    /*
    --------------------------------------------------------
     * JLOG-DATA: log-file writer for JIGSAW.
    --------------------------------------------------------
     */

    class jlog_text
        {
        public  :
    /*-------------------------- a "real" log-file writer */
            std::ofstream    _file ;   
        
        public  :
        
        __inline_call  jlog_text (
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
        
        __inline_call ~jlog_text (
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

    class jlog_null
        {
        public  :
    /*-------------------------- a "null" log-file writer */
        
            iptr_type  _verbosity ;
            
        public  :
    
        __inline_call  jlog_null (
            jcfg_data const& _jcfg
            )
        {
    /*-------------------------- def. no: for lib_jigsaw! */
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
    /*-------------------------- def. no: for lib_jigsaw! */
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
     * READ-DATA: load MESH data from file.
    --------------------------------------------------------
     */

    #   include "geo_load.hpp"
    #   include "ini_load.hpp"
    #   include "hfn_load.hpp"
    
    
    /*
    --------------------------------------------------------
     * INIT-DATA: initialise mesh pointers.
    --------------------------------------------------------
     */

    #   include "msh_init.hpp"
    #   include "hfn_init.hpp"
    
    
    /*
    --------------------------------------------------------
     * SAVE-MESH: push MESH data into file.
    --------------------------------------------------------
     */

    #   include "msh_save.hpp"
    
    
    /*
    --------------------------------------------------------
     * COPY-MESH: copy r-DT to tri-complex.
    --------------------------------------------------------
     */

    #   include "msh_copy.hpp"


    /*
    --------------------------------------------------------
     * CALC-MESH: call mesh generator.
    --------------------------------------------------------
     */
    
    #   include "run_mesh.hpp"
    #   include "run_iter.hpp"


    /*
    --------------------------------------------------------
     * JLOG-HEAD: write header for *.LOG file.
    --------------------------------------------------------
     */
 
    template <
    typename      jlog_data
             >
    __normal_call void_type jlog_head (
        jlog_data &_jlog
        )
    {
    /*-- NB: silliness re. escape sequences */
        _jlog.push (
    " \n"
    "#------------------------------------------------------------\n"
    "#\n"
    "#   ,o, ,o,       /                                 \n"
    "#    `   `  e88~88e  d88~\\   /~~~8e Y88b    e    / \n"
    "#   888 888 88   88 C888         88b Y88b  d8b  /   \n"
    "#   888 888 \"8b_d8\"  Y88b   e88~-888  Y888/Y88b/  \n"
    "#   888 888  /        888D C88   888   Y8/  Y8/     \n"
    "#   88P 888 Cb      \\_88P   \"8b_-888    Y    Y    \n"
    "# \\_8\"       Y8\"\"8D                             \n"
    "#\n"
    "#------------------------------------------------------------\n"
    "# JIGSAW: an unstructured mesh generation package.  \n"
    "#------------------------------------------------------------\n"
    " \n"
    "  " __JGSWVSTR "\n\n"
        ) ;
    }
    
        #define __jloglndv              \
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
     * MAIN: jumping-off point for JIGSAW!
    --------------------------------------------------------
     */
 
#   ifdef __lib_jigsaw

#   include "liblib/init_jig_t.hpp"
#   include "liblib/init_msh_t.hpp"

#   include "liblib/load_jig_t.hpp"
#   include "liblib/load_msh_t.hpp"

#   include "liblib/save_jig_t.hpp"
#   include "liblib/save_msh_t.hpp"
 
    __normal_call
        iptr_type jigsaw_make_mesh (    // lib-jigsaw
        jigsaw_jig_t *_jjig ,
        jigsaw_msh_t *_gmsh ,
        jigsaw_msh_t *_imsh ,
        jigsaw_msh_t *_hmsh ,
        jigsaw_msh_t *_mmsh
        )
    {
        iptr_type _retv = +0;
    
        hfun_data _hfun ;               // HFUN data
        geom_data _geom ;               // GEOM data
        rdel_data _rdel ;               // TRIA data
        mesh_data _mesh ;               // MESH data
        jcfg_data _jcfg ;
    
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
    
    /*--------------------------------- init. output data */    
        jigsaw_init_msh_t(_mmsh) ;
    
    /*--------------------------------- setup *.JLOG data */
        if (_jjig != nullptr )
        {
            _jcfg._verbosity = _jjig->_verbosity ;
        }
        jlog_null _jlog(_jcfg) ;
        jlog_head(_jlog) ;
        
    /*--------------------------------- parse *.JCFG data */
        if (_jjig != nullptr )
        {
            _jlog.push (
                "  Reading CFG. data...\n\n" ) ;

#           ifdef  __use_timers        
            _ttic   = _time.now();
#           endif//__use_timers
        
            if ((_retv = copy_jcfg (
                 _jcfg, 
                 _jlog,*_jjig)) != __no_error)
            {
                return  _retv ;
            }    
            
            if ((_retv = test_jcfg (
                 _jcfg, _jlog)) != __no_error)
            {
                return  _retv ;
            }

            _jlog.push (
                "  CFG. data summary...\n\n" ) ;

            if ((_retv = echo_jcfg (
                 _jcfg, _jlog)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers            
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if (_gmsh != nullptr )
        {
    /*--------------------------------- parse *.GEOM data */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Reading GEOM data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = copy_geom (
                 _jcfg, _jlog , 
                 _geom,*_gmsh)) != __no_error)
            {
                return  _retv ;
            }
            
            if ((_retv = test_geom (
                 _jcfg, 
                 _jlog, _geom)) != __no_error)
            {
                return  _retv ;
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if (_gmsh != nullptr )
        {
    /*--------------------------------- parse *.GEOM data */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming GEOM data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            _geom.init_geom(_jcfg) ;
            
            if (_jcfg._verbosity > 0 )
            {

            _jlog.push (
                "  GEOM data summary...\n\n" ) ;

            if ((_retv = echo_geom (
                 _jcfg, 
                 _jlog, _geom)) != __no_error)
            {
                return  _retv ;
            }
            
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }

        if (_imsh != nullptr )
        {
    /*--------------------------------- parse *.INIT data */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Reading INIT data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = copy_init (
                 _jcfg, _jlog, 
                 _mesh,*_mmsh)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv = test_init (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if (_imsh != nullptr )
        {
    /*--------------------------------- assemble init-con */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming INIT data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            _mesh._euclidean_mesh_2d.
                _mesh.make_ptrs();
            _mesh._euclidean_mesh_3d.
                _mesh.make_ptrs();

            if (_jcfg._verbosity > 0 )
            {

            _jlog.push (
                "  INIT data summary...\n\n" ) ;

            if ((_retv = echo_init (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }
            
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if (_hmsh != nullptr )
        {
    /*--------------------------------- parse *.HFUN data */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Reading HFUN data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = copy_hfun (
                 _jcfg, _jlog , 
                 _hfun,*_hmsh)) != __no_error)
            {
                return  _retv ;
            }
            
            if ((_retv = test_hfun (
                 _jcfg, 
                 _jlog, _hfun)) != __no_error)
            {
                return  _retv ;
            }
            
#           ifdef  __use_timers             
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if (_gmsh != nullptr )
        {
    /*--------------------------------- assemble size-fun */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming HFUN data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = init_hfun (
                 _jcfg, _jlog ,
                 _geom, _hfun)) != __no_error)
            {
                return  _retv ;
            }

            _hfun.init_hfun (_jcfg);
            
            if (_jcfg._verbosity > 0 )
            {

            _jlog.push (
                "  HFUN data summary...\n\n" ) ;
            
            if ((_retv = echo_hfun (
                 _jcfg, 
                 _jlog, _hfun)) != __no_error)
            {
                return  _retv ;
            } 
                       
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if (_gmsh != nullptr )
        {
            if(_jcfg._rdel_opts.iter() != +0 )
            {
    /*--------------------------------- call mesh routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Generate rDT MESH...\n\n" ) ;
                
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = make_mesh (
                 _jcfg, _jlog ,
                 _geom, _mesh , 
                 _hfun, _rdel)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers         
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
            }
        }
        
        if (_gmsh != nullptr )
        {
    /*--------------------------------- call copy routine */
            if(_jcfg._rdel_opts.iter() != +0 &&
               _jcfg._iter_opts.iter() != +0 )
            {
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Pushing MESH data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = copy_mesh (
                 _jcfg, _jlog ,
                 _rdel, _mesh)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
            }
        }
        
        if (_gmsh != nullptr )
        {
            if(_jcfg._iter_opts.iter() != +0 )
            {
    /*--------------------------------- call iter routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  MESH optimisation...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = init_mesh (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv = iter_mesh (
                 _jcfg, _jlog ,
                 _geom, 
                 _hfun, _mesh)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
            }
        }
        
        if (_gmsh != nullptr )
        {
    /*--------------------------------- dump mesh to data */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Writing MESH file...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if (_jcfg._rdel_opts.iter() != +0 &&
                _jcfg._iter_opts.iter() == +0 )
            {

            if ((_retv = save_msht (
                 _jcfg, _jlog , 
                 _rdel,*_mmsh)) != __no_error)
            {
                return  _retv ;
            }
        
            }    
            else
            {

            if ((_retv = save_msht (
                 _jcfg, _jlog , 
                 _mesh,*_mmsh)) != __no_error)
            {
                return  _retv ;
            }
        
            }

#           ifdef  __use_timers         
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
    /*-------------------------- success, if we got here! */

        return ( _retv ) ;
    }
        
#   else
 
    __normal_call iptr_type   main (    // cmd-jigsaw
        int           _argc , 
        char        **_argv
        )
    {
        hfun_data _hfun ;               // HFUN data
        geom_data _geom ;               // GEOM data
        rdel_data _rdel ;               // TRIA data
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
                
                std::cout << __JGSWVSTR ;
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
        jlog_text _jlog(_jcfg) ;
        jlog_head(_jlog) ;
      
        if(!_jcfg._jcfg_file.empty())
        {
    /*--------------------------------- parse *.JCFG file */
            _jlog.push (
                "  Reading CFG. file...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers
        
            if ((_retv = read_jcfg (
                 _jcfg, _jlog)) != __no_error)
            {
                return  _retv ;
            }           
            
            if ((_retv = test_jcfg (
                 _jcfg, _jlog)) != __no_error)
            {
                return  _retv ;
            }

            _jlog.push (
                "  CFG. data summary...\n\n" ) ;

            if ((_retv = echo_jcfg (
                 _jcfg, _jlog)) != __no_error)
            {
                return  _retv ;
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }

        if(!_jcfg._geom_file.empty())
        {
    /*--------------------------------- parse *.GEOM file */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Reading GEOM file...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = read_geom (
                 _jcfg, 
                 _jlog, _geom)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv = test_geom (
                 _jcfg, 
                 _jlog, _geom)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }

        if(!_jcfg._geom_file.empty())
        {
    /*--------------------------------- assemble geometry */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming GEOM data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            _geom.init_geom(_jcfg) ;
            
            if (_jcfg._verbosity > 0 )
            {

            _jlog.push (
                "  GEOM data summary...\n\n" ) ;

            if ((_retv = echo_geom (
                 _jcfg, 
                 _jlog, _geom)) != __no_error)
            {
                return  _retv ;
            }
            
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if(!_jcfg._init_file.empty())
        {
    /*--------------------------------- parse *.INIT file */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Reading INIT file...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = read_init (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv = test_init (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if(!_jcfg._init_file.empty())
        {
    /*--------------------------------- assemble init-con */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming INIT data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            _mesh._euclidean_mesh_2d.
                _mesh.make_ptrs();
            _mesh._euclidean_mesh_3d.
                _mesh.make_ptrs();

            if (_jcfg._verbosity > 0 )
            {

            _jlog.push (
                "  INIT data summary...\n\n" ) ;

            if ((_retv = echo_init (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }
            
            }
            
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }

        if(!_jcfg._hfun_file.empty())
        {
    /*--------------------------------- parse *.HFUN file */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Reading HFUN file...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = read_hfun (
                 _jcfg, 
                 _jlog, _hfun)) != __no_error)
            {
                return  _retv ;
            }
            
            if ((_retv = test_hfun (
                 _jcfg, 
                 _jlog, _hfun)) != __no_error)
            {
                return  _retv ;
            }
             
#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if(!_jcfg._geom_file.empty())
        {
    /*--------------------------------- assemble size-fun */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming HFUN data...\n\n" ) ;
        
#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = init_hfun (
                 _jcfg, _jlog ,
                 _geom, _hfun)) != __no_error)
            {
                return  _retv ;
            }

            _hfun.init_hfun (_jcfg);

            if (_jcfg._verbosity > 0 )
            {

            _jlog.push (
                "  HFUN data summary...\n\n" ) ;
            
            if ((_retv = echo_hfun (
                 _jcfg, 
                 _jlog, _hfun)) != __no_error)
            {
                return  _retv ;
            } 
                       
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }
        
        if(!_jcfg._geom_file.empty())
        {
            if(_jcfg._rdel_opts.iter() != +0 )
            {
    /*--------------------------------- call mesh routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Generate rDT MESH...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = make_mesh (
                 _jcfg, _jlog ,
                 _geom, _mesh ,
                 _hfun, _rdel)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers         
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
            }
        }
        
        if(!_jcfg._geom_file.empty() &&
           !_jcfg._tria_file.empty() )
        {
    /*--------------------------------- dump tria to file */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Writing TRIA file...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = save_tria (
                 _jcfg, 
                 _jlog, _rdel)) != __no_error)
            {
                return  _retv ;
            }
        
#           ifdef  __use_timers         
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }

        if(!_jcfg._geom_file.empty())
        {
            if(_jcfg._rdel_opts.iter() != +0 &&
               _jcfg._iter_opts.iter() != +0 )
            {
    /*--------------------------------- call copy routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Pushing MESH data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = copy_mesh (
                 _jcfg, _jlog ,
                 _rdel, _mesh)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
            }
        }
        
        if(!_jcfg._geom_file.empty())
        {
            if(_jcfg._iter_opts.iter() != +0 )
            {
    /*--------------------------------- call iter routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  MESH optimisation...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = init_mesh (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv = iter_mesh (
                 _jcfg, _jlog ,
                 _geom, 
                 _hfun, _mesh)) != __no_error)
            {
                return  _retv ;
            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
            }
        }
        
        if(!_jcfg._geom_file.empty() &&
           !_jcfg._mesh_file.empty() )
        {
    /*--------------------------------- dump mesh to file */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Writing MESH file...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if (_jcfg._rdel_opts.iter() != +0 &&
                _jcfg._iter_opts.iter() == +0 )
            {

            if ((_retv = save_jmsh (
                 _jcfg, 
                 _jlog, _rdel)) != __no_error)
            {
                return  _retv ;
            }
        
            }    
            else
            {

            if ((_retv = save_jmsh (
                 _jcfg, 
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }
        
            }

#           ifdef  __use_timers         
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
        }

    /*-------------------------- success, if we got here! */

        return ( _retv ) ;
    }

        
#   endif   //__lib_jigsaw

    
    
