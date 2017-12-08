
    /*
    --------------------------------------------------------
     * RUN-ITER: do the mesh-optimisation step. 
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
     * Last updated: 24 August, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RUN_ITER__
#   define __RUN_ITER__

    /*
    --------------------------------------------------------
     * Call the 2-dimensional mesh optimiser.
    --------------------------------------------------------
     */
 
    template <
    typename      geom_type ,
    typename      hfun_type ,
    typename      mesh_type ,
    typename      jlog_data
             >
    __normal_call void_type iter_euclidean_2d (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        jcfg_data &_args ,
        jlog_data &_jlog
        )
    {
        if (true) /* kernel? */
        {
            typedef mesh::
                iter_pred_euclidean_2d <
                real_type ,
                iptr_type >         pred_type ;
            
            typedef mesh::iter_mesh_2  <
                geom_type ,
                typename
                mesh_type::
                mesh_type ,
                hfun_type ,
                pred_type >         iter_func ;
            
            typedef 
            jcfg_data::iter_opts    iter_opts ;
 
            iter_opts *_opts =  
               &_args._iter_opts;
            
            pred_type  _pred ;
            iter_func::iter_mesh(
                _geom, _hfun , 
                _mesh. _mesh, 
                _pred,
               *_opts, _jlog )  ;
        }
    }
    
    template <
    typename      geom_type ,
    typename      hfun_type ,
    typename      mesh_type ,
    typename      jlog_data
             >
    __normal_call void_type iter_ellipsoid_3d (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        jcfg_data &_args ,
        jlog_data &_jlog
        )
    {
        if (true) /* kernel? */
        {
            typedef mesh::
                iter_pred_ellipsoid_3d <
                real_type ,
                iptr_type >         pred_type ;
            
            typedef mesh::iter_mesh_2  <
                geom_type ,
                typename
                mesh_type::
                mesh_type ,
                hfun_type ,
                pred_type >         iter_func ;
            
            typedef 
            jcfg_data::iter_opts    iter_opts ;
 
            iter_opts *_opts =  
               &_args._iter_opts;
            
            pred_type  _pred ;
            iter_func::iter_mesh(
                _geom, _hfun , 
                _mesh. _mesh, 
                _pred,
               *_opts, _jlog )  ;
        }
    }
    
    /*
    --------------------------------------------------------
     * Call the 3-dimensional mesh optimiser.
    --------------------------------------------------------
     */
 
    template <
    typename      geom_type ,
    typename      hfun_type ,
    typename      mesh_type ,
    typename      jlog_data
             >
    __normal_call void_type iter_euclidean_3d (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        jcfg_data &_args ,
        jlog_data &_jlog
        )
    {
    
        //!! yup...
    
    }    

    /*
    --------------------------------------------------------
     * ITER-MESH: call mesh optimiser.
    --------------------------------------------------------
     */
 
    template <
    typename      jlog_data
             >
    __normal_call iptr_type iter_mesh (
        jcfg_data &_args,
        jlog_data &_jlog,
        geom_data &_geom,
        hfun_data &_hfun,
        mesh_data &_mesh
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
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +2 ;
        
                iter_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _hfun._constant_value_kd,
                _mesh._euclidean_mesh_2d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +2 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- with euclidean-mesh HFUN kernel */
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +2 ;
        
                iter_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _hfun._euclidean_mesh_2d,
                _mesh._euclidean_mesh_2d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +2 &&
                _hfun._kind ==
                jmsh_kind::euclidean_grid)
            {
        /*----------- with euclidean-grid HFUN kernel */
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +2 ;
                
                iter_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _hfun._euclidean_grid_2d,
                _mesh._euclidean_mesh_2d,
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
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +3 ;
        
                /*
                iter_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _hfun._constant_value_kd,
                _mesh._euclidean_mesh_3d,
                _args, _jlog) ;
                 */
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*----------- with euclidean-mesh HFUN kernel */
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +3 ;
        
                /*
                iter_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _hfun._euclidean_mesh_3d,
                _mesh._euclidean_mesh_3d,
                _args, _jlog) ;
                 */
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_grid)
            {
        /*----------- with euclidean-grid HFUN kernel */
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +3 ;
                
                /*
                iter_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _hfun._euclidean_grid_3d,
                _mesh._euclidean_mesh_3d,
                _args, _jlog) ;
                 */
            }
            
            }
            if (_geom._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
        /*----------- have ellipsoid-mesh GEOM kernel */
            
            if (_hfun._ndim == +0 )
            {
        /*----------- with constant-value HFUN kernel */
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +3 ;
        
                iter_ellipsoid_3d (
                _geom._ellipsoid_mesh_3d,
                _hfun._constant_value_kd,
                _mesh._euclidean_mesh_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
        /*----------- with ellipsoid-grid HFUN kernel */
                _mesh._kind  = 
                jmsh_kind::euclidean_mesh;
                
                _mesh._ndim  = +3 ;
                
                iter_ellipsoid_3d (
                _geom._ellipsoid_mesh_3d,
                _hfun._ellipsoid_grid_3d,
                _mesh._euclidean_mesh_3d,
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


#   endif   //__RUN_ITER__



