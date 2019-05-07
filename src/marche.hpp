
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
     * MARCHE: "fast-marching" eikonal eqn. solver.
    --------------------------------------------------------
     *
     * Last updated: 22 January, 2019
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
     */
     
    template <
        typename  jlog_data
             >
    __normal_call void_type marche_banner (
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
    "# MARCHE: \"fast-marching\" eikonal eqn. solver.  \n"
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
 
    __normal_call iptr_type marche (    // lib-jigsaw
        jigsaw_jig_t *_jjig ,
        jigsaw_msh_t *_fmsh
        )
    {
        iptr_type _retv = +0;
        
        __unreferenced(_fmsh) ;

        hfun_data _ffun ;               // FUNC data
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
    
    /*--------------------------------- setup *.JLOG data */
        if (_jjig != nullptr )
        {
            _jcfg._verbosity = _jjig->_verbosity ;
        }
        jlog_null _jlog(_jcfg) ;
        marche_banner  (_jlog) ;
        
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
     
    /*-------------------------- success, if we got here! */

        return ( _retv ) ;
    }

#   else

#   if defined(__cmd_marche)

    __normal_call iptr_type main (      // cmd-marche
        int           _argc , 
        char        **_argv
        )
    {
        hfun_data _ffun ;               // FUNC data
        
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
        marche_banner  (_jlog) ;
      
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
        
/*-------------------------- success, if we got here! */

        return ( _retv ) ;        
    }
    
#   endif   //__cmd_marche

#   endif   //__lib_jigsaw 
    
    
    
