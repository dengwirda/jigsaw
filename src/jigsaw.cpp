
    //
    // for cmd-jigsaw:
    //
    // g++ -std=c++11 -pedantic -Wall -s -O3 -flto -D NDEBUG 
    // -D __cmd_jigsaw jigsaw.cpp -o jigsaw
    //
    // g++ -std=c++11 -pedantic -Wall -s -O3 -flto -D NDEBUG 
    // -D __cmd_tripod jigsaw.cpp -o tripod
    //
    //
    // for lib-jigsaw:
    //
    // g++ -std=c++11 -pedantic -Wall -O3 -flto -fPIC 
    // -D NDEBUG -D __lib_jigsaw jigsaw.cpp -shared -o libjigsaw.so
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
     * JIGSAW: an unstructured mesh generation library.
    --------------------------------------------------------
     *
     * JIGSAW release 0.9.10.x
     * Last updated: 10 April, 2019
     *
     * Copyright 2013 -- 2019
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
     * JIGSAW's 'frontal' Delaunay-refinement algorithms
     * are described here:
     *
     * D. Engwirda, D. Ivers, (2016): "Off-centre Steiner 
     * points for Delaunay-refinement on curved surfaces", 
     * Computer-Aided Design, 72, pp. 157-171, 
     * http://dx.doi.org/10.1016/j.cad.2015.10.007
     *
     * D. Engwirda, (2016): "Conforming restricted Delaunay 
     * mesh generation for piecewise smooth complexes", 
     * Procedia Engineering, 163, pp. 84-96, 
     * http://dx.doi.org/10.1016/j.proeng.2016.11.024
     *
     * D. Engwirda, (2015): "Voronoi-based point-placement 
     * for three-dimensional Delaunay-refinement", 
     * Procedia Engineering, 124, pp. 330-342, 
     * http://dx.doi.org/10.1016/j.proeng.2015.10.143
     *
     * JIGSAW's hybrid, optimisation-based mesh improvement
     * schemes are discussed here: 
     *
     * D. Engwirda, (2018): "Generalised primal-dual grids 
     * for unstructured co-volume schemes, J. Comp. Phys., 
     * 375, pp. 155-176, 
     * https://doi.org/10.1016/j.jcp.2018.07.025 
     *
     * JIGSAW originally grew out of my Ph.D. research, in 
     * which I explored initial versions of the refinement 
     * and optimisation-based algorithms:
     *
     * D. Engwirda, (2014): "Locally-optimal Delaunay-
     * refinement and optimisation-based mesh generation", 
     * Ph.D. Thesis, School of Mathematics and Statistics, 
     * Univ. of Sydney. 
     * http://hdl.handle.net/2123/13148
     *
    --------------------------------------------------------
     *
     * JIGSAW's "restricted" Delaunay refinement strategies
     * are generalisations of the methods developed in:
     *
     * J.D. Boissonnat, S. Oudot, (2005): "Provably Good 
     * Sampling and Meshing of Surfaces", Graphical Models, 
     * 67, pp. 405-451,
     * https://doi.org/10.1016/j.gmod.2005.01.004
     *
     * L. Rineau, M. Yvinec, (2008): "Meshing 3D Domains 
     * Bounded by Piecewise Smooth Surfaces", Proc. of the 
     * 16th International Meshing Roundtable, pp. 443-460,
     * https://doi.org/10.1007/978-3-540-75103-8_25
     *
     * C. Jamin, P. Alliez, M. Yvinec, and J.D. Boissonnat, 
     * (2015): "CGALmesh: a generic framework for Delaunay 
     * mesh generation", ACM Transactions on Mathematical 
     * Software (TOMS), 41, pp. 23
     * https://doi.org/10.1145/2699463
     *
     * S.W. Cheng, T.K. Dey, E.A. Ramos, (2010): "Delaunay 
     * Refinement for Piecewise Smooth Complexes", 
     * Discrete & Computational Geometry, 43, pp. 121-166,
     * https://doi.org/10.1007/s00454-008-9109-3
     *
     * JIGSAW employs a "hybrid" mesh-optimisation approach
     * based on a combination of ODT techniques and direct
     * gradient-based optimisation:
     *
     * L. Chen, J.C. Xu, (2004): "Optimal Delaunay 
     * triangulations, J. Comp. Math., 22, pp. 299-308,
     * https://www.jstor.org/stable/43693155
     *
     * L.A. Freitag, C. Ollivier-Gooch, (1997): "Tetrahedral 
     * mesh improvement using swapping and smoothing", 
     * International Journal for Numerical Methods in 
     * Engineering 40 (21), pp. 3979-4002,
     * https://doi.org/10.1002/(SICI)1097-0207
     * (19971115)40:21<3979::AID-NME251>3.0.CO;2-9 
     *
     * B.M. Klingner, J.R. Shewchuk, (2008)" "Aggressive 
     * Tetrahedral Mesh Improvement", Proc. of the 16th 
     * International Meshing Roundtable, pp. 3-23,
     * https://doi.org/10.1007/978-3-540-75103-8_1
     *
     * P. Mullen, P. Memari, F. de Goes, M. Desbrun, (2011): 
     * "HOT: Hodge-optimized triangulations", ACM 
     * Transactions on Graphics (TOG), 30 (4) pp. 103,
     * https://doi.org/10.1145/2010324.1964998
     *
     * Core theory and techniques for Delaunay tessellation
     * and refinement can be found (for example) here:
     *
     * S. Cheng, T. Dey and J. Shewchuk, (2012): "Delaunay 
     * mesh generation", CRC Press.
     *
     * Additional references are provided inline throughout 
     * the JIGSAW src.
     *
    --------------------------------------------------------
     */


#   define __jloglndv   \
"#------------------------------------------------------------\n"

//  define __cmd_jigsaw               // the cmd-ln exe's
//  define __cmd_tripod
//  define __cmd_marche

//  define __lib_jigsaw               // a shared library

#   if !defined(__cmd_jigsaw) && \
       !defined(__cmd_tripod) && \
       !defined(__cmd_marche) && \
       !defined(__lib_jigsaw)

    /*---------------------------------- build by default */
#       define  __cmd_jigsaw

#   endif

#   define __JGSWVSTR "JIGSAW VERSION 0.9.10"

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

    extern  "C" 
    {
#   include "../inc/lib_jigsaw.h"
    }

    typedef real_t real_type ;        // double-precision
    typedef indx_t iptr_type ;        // 32bit signed int

    /*---------------------------------- JIGSAW mesh kind */        
        
    struct jmsh_kind {   
        enum enum_data {
            null_mesh_kind      = +0 ,
            euclidean_mesh ,
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
        
        std::size_t             _ndim = +0;

        jmsh_kind ::
        enum_data               _kind = 
                            jmsh_kind::null_mesh_kind ;

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
         
		std::size_t             _ndim = +0;
        
        jmsh_kind ::
        enum_data               _kind = 
                            jmsh_kind::null_mesh_kind ;
      
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
           
		std::size_t             _ndim = +0;

        jmsh_kind ::
        enum_data               _kind = 
                            jmsh_kind::null_mesh_kind ;
                
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
           
		std::size_t             _ndim = +0;

        jmsh_kind ::
        enum_data               _kind = 
                            jmsh_kind::null_mesh_kind ;
                
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
    #   include "run_tria.hpp"


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
     * Jumping-off points for CMD + LIB JIGSAW!
    --------------------------------------------------------
     */
 
    #   include "jigsaw.hpp"
    #   include "tripod.hpp"
    
    #   include "marche.hpp"

//  #   include "stitch.hpp"
    

    
