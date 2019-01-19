
    
    /*
    --------------------------------------------------------
     * RUN-TRIA: do the restricted tessellation step. 
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
     * Last updated: 28 December, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RUN_TRIA__
#   define __RUN_TRIA__

    /*
    --------------------------------------------------------
     * Call the 2-dimensional tessellator.
    --------------------------------------------------------
     */
 
    template <
    typename      geom_type ,
    typename      init_type ,
    typename      rdel_type ,
    typename      jlog_data
             >
    __normal_call void_type rdel_euclidean_2d (
        geom_type &_geom,
        init_type &_init,
        rdel_type &_rdel,
        jcfg_data &_args,
        jlog_data &_jlog
        )
    {
    /*------------------------------ call rDEL kernel */
        typedef mesh::rdel_make_2d  <
            rdel_type , 
            geom_type >         rdel_func ;

        typedef 
        jcfg_data::rdel_opts    rdel_opts ;

        rdel_opts  *_opts =  
           &_args._rdel_opts ;

        rdel_func::rdel_make ( 
            _geom,  _init ,
            _rdel, *_opts , _jlog ) ;
    }
    
    /*
    --------------------------------------------------------
     * Call the 3-dimensional tessellator.
    --------------------------------------------------------
     */
 
    template <
    typename      geom_type ,
    typename      init_type ,
    typename      rdel_type ,
    typename      jlog_data
             >
    __normal_call void_type rdel_euclidean_3d (
        geom_type &_geom,
        init_type &_init,
        rdel_type &_rdel,
        jcfg_data &_args,
        jlog_data &_jlog
        )
    {
    /*------------------------------ call rDEL kernel */
        typedef mesh::rdel_make_3d  <
            rdel_type , 
            geom_type >         rdel_func ;

        typedef 
        jcfg_data::rdel_opts    rdel_opts ;

        rdel_opts  *_opts =  
           &_args._rdel_opts ;

        rdel_func::rdel_make ( 
            _geom,  _init ,
            _rdel, *_opts , _jlog ) ;
    }
    
    /*
    --------------------------------------------------------
     * MAKE-RDEL: call the tessellator.
    --------------------------------------------------------
     */
 
    template <
    typename      jlog_data
             >
    __normal_call iptr_type make_rdel (
        jcfg_data &_args,
        jlog_data &_jlog,
        mesh_data &_init,
        geom_data &_geom,
        rdel_data &_rdel
        )
    {
        iptr_type _errv = __no_error ;

        try 
        {
            if (_geom._ndim == +2 &&
                _geom._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- have euclidean-mesh GEOM kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +2 ;
                
                rdel_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _rdel._euclidean_rdel_2d,
                _args, _jlog) ;
            }
            else           
            if (_geom._ndim == +3 &&
                _geom._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- have euclidean-mesh GEOM kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                rdel_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_geom._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
        /*----------- have ellipsoid-mesh GEOM kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                rdel_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }


#   endif   // __RUN_TRIA__



