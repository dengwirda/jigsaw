
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
     * Last updated: 28 December, 2018
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
     */
     
    template <
        typename  jlog_data
             >
    __normal_call void_type jigsaw_banner (
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
    "# JIGSAW: an unstructured mesh generation library.  \n"
    "#------------------------------------------------------------\n"
    " \n"
    "  " __JGSWVSTR "\n\n"
        ) ;
    }
     
#   ifdef __lib_jigsaw

#   include "liblib/init_jig_t.hpp"
#   include "liblib/init_msh_t.hpp"

#   include "liblib/load_jig_t.hpp"
#   include "liblib/load_msh_t.hpp"

#   include "liblib/save_jig_t.hpp"
#   include "liblib/save_msh_t.hpp"
 
    __normal_call iptr_type jigsaw (    // lib-jigsaw
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
        jigsaw_banner  (_jlog) ;
        
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

#   if defined(__cmd_jigsaw)

    __normal_call iptr_type main (      // cmd-jigsaw
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
        jigsaw_banner  (_jlog) ;
      
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
    
#   endif   //__cmd_jigsaw
    
#   endif   //__lib_jigsaw
     
     
    
