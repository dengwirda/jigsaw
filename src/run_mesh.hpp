
    /*
    --------------------------------------------------------
     * RUN-MESH: do the mesh-generation step. 
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
     * Last updated: 14 March, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RUN_MESH__
#   define __RUN_MESH__

    /*
    --------------------------------------------------------
     * Call the 2-dimensional mesh generator.
    --------------------------------------------------------
     */
 
    template <
    typename      geom_type ,
    typename      init_type ,
    typename      hfun_type ,
    typename      mesh_type ,
    typename      jlog_data
             >
    __normal_call void_type mesh_euclidean_2d (
        geom_type &_geom,
        init_type &_init,
        hfun_type &_hfun,
        mesh_type &_mesh,
        jcfg_data &_args,
        jlog_data &_jlog
        )
    {
        if (_args._mesh_pred == 
                jcfg_data::mesh_pred::delaunay)
        {
    /*-------------------------- call DELAUNAY kernel */
            typedef 
                mesh::rdel_delaunay_2d  <
                geom_type ,
                hfun_type ,
                mesh_type >         rdel_pred ;
                
            typedef mesh::rdel_mesh_2d  <
                mesh_type , 
                rdel_pred ,
                geom_type ,
                hfun_type >         rdel_func ;

            typedef 
            jcfg_data::rdel_opts    rdel_opts ;
 
            rdel_opts *_opts =  
               &_args._rdel_opts;

            rdel_func::rdel_mesh( 
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )  ;
        }
        else
        if (_args._mesh_pred ==
                jcfg_data::mesh_pred::delfront)
        {
    /*-------------------------- call DELFRONT kernel */
            typedef 
                mesh::rdel_delfront_2d  <
                geom_type ,
                hfun_type ,
                mesh_type >         rdel_pred ;
                
            typedef mesh::rdel_mesh_2d  <
                mesh_type , 
                rdel_pred ,
                geom_type ,
                hfun_type >         rdel_func ;

            typedef 
            jcfg_data::rdel_opts    rdel_opts ;
 
            rdel_opts *_opts =  
               &_args._rdel_opts;

            rdel_func::rdel_mesh( 
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )  ;
        }
    }

    /*
    --------------------------------------------------------
     * Call the 3-dimensional mesh generator.
    --------------------------------------------------------
     */
     
    template <
    typename      geom_type ,
    typename      init_type ,
    typename      hfun_type ,
    typename      mesh_type ,
    typename      jlog_data
             >
    __normal_call void_type mesh_euclidean_3d (
        geom_type &_geom,
        init_type &_init,
        hfun_type &_hfun,
        mesh_type &_mesh,
        jcfg_data &_args,
        jlog_data &_jlog
        )
    {    
        if (_args._mesh_pred == 
                jcfg_data::mesh_pred::delaunay)
        {
    /*-------------------------- call DELAUNAY kernel */
            typedef 
                mesh::rdel_delaunay_3d  <
                geom_type ,
                hfun_type ,
                mesh_type >         rdel_pred ;
                
            typedef mesh::rdel_mesh_3d  <
                mesh_type , 
                rdel_pred ,
                geom_type ,
                hfun_type >         rdel_func ;

            typedef 
            jcfg_data::rdel_opts    rdel_opts ;
 
            rdel_opts *_opts =  
               &_args._rdel_opts;

            rdel_func::rdel_mesh( 
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )  ;
        }
        else
        if (_args._mesh_pred ==
                jcfg_data::mesh_pred::delfront)
        {
    /*-------------------------- call DELFRONT kernel */
            typedef 
                mesh::rdel_delfront_3d  <
                geom_type ,
                hfun_type ,
                mesh_type >         rdel_pred ;
                
            typedef mesh::rdel_mesh_3d  <
                mesh_type , 
                rdel_pred ,
                geom_type ,
                hfun_type >         rdel_func ;

            typedef 
            jcfg_data::rdel_opts    rdel_opts ;
 
            rdel_opts *_opts =  
               &_args._rdel_opts;

            rdel_func::rdel_mesh( 
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )  ;
        }
    }

    /*
    --------------------------------------------------------
     * MAKE-MESH: call mesh generator.
    --------------------------------------------------------
     */
 
    template <
    typename      jlog_data
             >
    __normal_call iptr_type make_mesh (
        jcfg_data &_args,
        jlog_data &_jlog,
        geom_data &_geom,
        mesh_data &_init,
        hfun_data &_hfun,
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
            
            if (_hfun._ndim == +0 )
            {
        /*----------- with constant-value HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +2 ;
                
                mesh_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _hfun._constant_value_kd,
                _rdel._euclidean_rdel_2d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +2 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- with euclidean-mesh HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +2 ;
                
                mesh_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _hfun._euclidean_mesh_2d,
                _rdel._euclidean_rdel_2d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +2 &&
                _hfun._kind ==
                 jmsh_kind::euclidean_grid)
            {
        /*----------- with euclidean-grid HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +2 ;
                
                mesh_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _hfun._euclidean_grid_2d,
                _rdel._euclidean_rdel_2d,
                _args, _jlog) ;
            }
            
            }
            else           
            if (_geom._ndim == +3 &&
                _geom._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- have euclidean-mesh GEOM kernel */
            
            if (_hfun._ndim == +0 )
            {
        /*----------- with constant-value HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                mesh_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._constant_value_kd,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- with euclidean-mesh HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                mesh_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._euclidean_mesh_3d,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_grid)
            {
        /*----------- with euclidean-grid HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                mesh_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._euclidean_grid_3d,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            
            }
            else
            if (_geom._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
        /*----------- have ellipsoid-mesh GEOM kernel */
            
            if (_hfun._ndim == +0 )
            {
        /*----------- with constant-value HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._constant_value_kd,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
        /*----------- with ellipsoid-grid HFUN kernel */
                _rdel._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _rdel._ndim  = +3 ;
                
                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._ellipsoid_grid_3d,
                _rdel._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            
            }
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }


#   endif   //__RUN_MESH__



