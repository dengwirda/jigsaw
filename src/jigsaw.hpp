
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
     * Last updated: 05 March, 2020
     *
     * Copyright 2013 -- 2020
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

    namespace JIGSAW {

    std::string asciibanner =
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
    "  " __JGSWVSTR "\n\n"  ;

    /*
    --------------------------------------------------------
     * Call the 2-dimensional mesh generator.
    --------------------------------------------------------
     */

    template <
        typename  geom_type ,
        typename  init_type ,
        typename  hfun_type ,
        typename  mesh_type ,
        typename  jlog_data
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
                hfun_type ,
                init_type >         mesh_func ;

            typedef
            jcfg_data::mesh_opts    mesh_opts ;

            mesh_opts *_opts =
               &_args._mesh_opts ;

            mesh_func::rdel_mesh (
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )   ;
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
                hfun_type ,
                init_type >         mesh_func ;

            typedef
            jcfg_data::mesh_opts    mesh_opts ;

            mesh_opts *_opts =
               &_args._mesh_opts ;

            mesh_func::rdel_mesh (
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )   ;
        }
        else
        if (_args._mesh_pred ==
                jcfg_data::mesh_pred::bisector)
        {
    /*-------------------------- call BISECTOR kernel */

            /*
            typedef
                mesh::tree_bisector_2d  <
                geom_type ,
                hfun_type ,
                mesh_type >         tree_pred ;

            typedef mesh::tree_mesh_2d  <
                mesh_type ,
                tree_pred ,
                geom_type ,
                hfun_type ,
                init_type >         mesh_func ;

            typedef
            jcfg_data::mesh_opts    mesh_opts ;

            mesh_opts *_opts =
               &_args._mesh_opts ;

            mesh_func::tree_mesh (
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )   ;
            */

        }
    }

    /*
    --------------------------------------------------------
     * Call the 3-dimensional mesh generator.
    --------------------------------------------------------
     */

    template <
        typename  geom_type ,
        typename  init_type ,
        typename  hfun_type ,
        typename  mesh_type ,
        typename  jlog_data
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
                hfun_type ,
                init_type >         mesh_func ;

            typedef
            jcfg_data::mesh_opts    mesh_opts ;

            mesh_opts *_opts =
               &_args._mesh_opts ;

            mesh_func::rdel_mesh (
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )   ;
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
                hfun_type ,
                init_type >         mesh_func ;

            typedef
            jcfg_data::mesh_opts    mesh_opts ;

            mesh_opts *_opts =
               &_args._mesh_opts ;

            mesh_func::rdel_mesh (
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )   ;
        }
        else
        if (_args._mesh_pred ==
                jcfg_data::mesh_pred::bisector)
        {
    /*-------------------------- call BISECTOR kernel */

            /*
            typedef
                mesh::tree_bisector_3d  <
                geom_type ,
                hfun_type ,
                mesh_type >         tree_pred ;

            typedef mesh::tree_mesh_3d  <
                mesh_type ,
                tree_pred ,
                geom_type ,
                hfun_type ,
                init_type >         mesh_func ;

            typedef
            jcfg_data::mesh_opts    mesh_opts ;

            mesh_opts *_opts =
               &_args._mesh_opts ;

            mesh_func::tree_mesh (
                _geom, _init ,
                _hfun, _mesh ,
               *_opts, _jlog )   ;
            */

        }
    }

    /*
    --------------------------------------------------------
     * Call the k-dimensional mesh generator.
    --------------------------------------------------------
     */

    template <
        typename  jlog_data
             >
    __normal_call iptr_type mesh_core (
        jcfg_data &_args,
        jlog_data &_jlog,
        geom_data &_geom,
        mesh_data &_init,
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
        /*--------------- have euclidean-mesh GEOM kernel */

            if (_hfun._ndim == +0 )
            {
        /*--------------- with constant-value HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +2 ;

                mesh_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _hfun._constant_value_kd,
                _mesh._euclidean_rdel_2d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +2 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*--------------- with euclidean-mesh HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +2 ;

                mesh_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _hfun._euclidean_mesh_2d,
                _mesh._euclidean_rdel_2d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +2 &&
                _hfun._kind ==
                 jmsh_kind::euclidean_grid)
            {
        /*--------------- with euclidean-grid HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +2 ;

                mesh_euclidean_2d (
                _geom._euclidean_mesh_2d,
                _init._euclidean_mesh_2d,
                _hfun._euclidean_grid_2d,
                _mesh._euclidean_rdel_2d,
                _args, _jlog) ;
            }

            }
            else
            if (_geom._ndim == +3 &&
                _geom._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*--------------- have euclidean-mesh GEOM kernel */

            if (_hfun._ndim == +0 )
            {
        /*--------------- with constant-value HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._constant_value_kd,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*--------------- with euclidean-mesh HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._euclidean_mesh_3d,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_grid)
            {
        /*--------------- with euclidean-grid HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._euclidean_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._euclidean_grid_3d,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }

            }
            else
            if (_geom._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
        /*--------------- have ellipsoid-mesh GEOM kernel */

            if (_hfun._ndim == +0 )
            {
        /*--------------- with constant-value HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._constant_value_kd,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
        /*--------------- with ellipsoid-grid HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._ellipsoid_grid_3d,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
        /*--------------- with ellipsoid-mesh HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._ellipsoid_mesh_3d,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_mesh)
            {
        /*--------------- with euclidean-mesh HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._euclidean_mesh_3d,
                _mesh._euclidean_rdel_3d,
                _args, _jlog) ;
            }
            else
            if (_hfun._ndim == +3 &&
                _hfun._kind ==
                jmsh_kind::euclidean_grid)
            {
        /*--------------- with euclidean-grid HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                mesh_euclidean_3d (
                _geom._ellipsoid_mesh_3d,
                _init._euclidean_mesh_3d,
                _hfun._euclidean_grid_3d,
                _mesh._euclidean_rdel_3d,
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

    /*
    --------------------------------------------------------
     * Call the 2-dimensional mesh optimiser.
    --------------------------------------------------------
     */

    template <
        typename  geom_type ,
        typename  hfun_type ,
        typename  mesh_type ,
        typename  jlog_data
             >
    __normal_call void_type iter_euclidean_2d (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        jcfg_data &_args ,
        jlog_data &_jlog
        )
    {
        if (_args._iter_pred ==
                jcfg_data::iter_pred::odt_dqdx)
        {
    /*-------------------------- call ODT-ITER kernel */
            typedef mesh::
                iter_pred_euclidean_2d <
                geom_type ,
                typename
                mesh_type::
                mesh_type >         pred_type ;

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
               &_args._iter_opts ;

            iter_func::iter_mesh (
                _geom, _hfun ,
                _mesh. _mesh ,
            iter_func::_odt_kern ,
               *_opts, _jlog )   ;
        }
        else
        if (_args._iter_pred ==
                jcfg_data::iter_pred::cvt_dqdx)
        {
    /*-------------------------- call CVT-ITER kernel */
            typedef mesh::
                iter_pred_euclidean_2d <
                geom_type ,
                typename
                mesh_type::
                mesh_type >         pred_type ;

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
               &_args._iter_opts ;

            iter_func::iter_mesh (
                _geom, _hfun ,
                _mesh. _mesh ,
            iter_func::_cvt_kern ,
               *_opts, _jlog )   ;
        }
    }

    template <
        typename  geom_type ,
        typename  hfun_type ,
        typename  mesh_type ,
        typename  jlog_data
             >
    __normal_call void_type iter_ellipsoid_3d (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        jcfg_data &_args ,
        jlog_data &_jlog
        )
    {
        if (_args._iter_pred ==
                jcfg_data::iter_pred::odt_dqdx)
        {
    /*-------------------------- call ODT-ITER kernel */
            typedef mesh::
                iter_pred_ellipsoid_3d <
                geom_type ,
                typename
                mesh_type::
                mesh_type >         pred_type ;

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
               &_args._iter_opts ;

            iter_func::iter_mesh (
                _geom, _hfun ,
                _mesh. _mesh ,
            iter_func::_odt_kern ,
               *_opts, _jlog )   ;
        }
        else
        if (_args._iter_pred ==
                jcfg_data::iter_pred::cvt_dqdx)
        {
    /*-------------------------- call CVT-ITER kernel */
            typedef mesh::
                iter_pred_ellipsoid_3d <
                geom_type ,
                typename
                mesh_type::
                mesh_type >         pred_type ;

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
               &_args._iter_opts ;

            iter_func::iter_mesh (
                _geom, _hfun ,
                _mesh. _mesh ,
            iter_func::_cvt_kern ,
               *_opts, _jlog )   ;
        }
    }

    /*
    --------------------------------------------------------
     * Call the 3-dimensional mesh optimiser.
    --------------------------------------------------------
     */

    template <
        typename  geom_type ,
        typename  hfun_type ,
        typename  mesh_type ,
        typename  jlog_data
             >
    __normal_call void_type iter_euclidean_3d (
        geom_type &/*_geom*/,
        hfun_type &/*_hfun*/,
        mesh_type &/*_mesh*/,
        jcfg_data &/*_args*/,
        jlog_data &/*_jlog*/
        )
    {   //!! yes, indeed...

    }

    /*
    --------------------------------------------------------
     * Call the k-dimensional mesh optimiser.
    --------------------------------------------------------
     */

    template <
        typename  jlog_data
             >
    __normal_call iptr_type iter_core (
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
            else
            if (_hfun._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
        /*----------- with ellipsoid-mesh HFUN kernel */
                _mesh._kind  =
                jmsh_kind::euclidean_mesh;

                _mesh._ndim  = +3 ;

                iter_ellipsoid_3d (
                _geom._ellipsoid_mesh_3d,
                _hfun._ellipsoid_mesh_3d,
                _mesh._euclidean_mesh_3d,
                _args, _jlog) ;
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

                iter_ellipsoid_3d (
                _geom._ellipsoid_mesh_3d,
                _hfun._euclidean_mesh_3d,
                _mesh._euclidean_mesh_3d,
                _args, _jlog) ;
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

                iter_ellipsoid_3d (
                _geom._ellipsoid_mesh_3d,
                _hfun._euclidean_grid_3d,
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

    /*--------------------------------- init. geo. kernel */
        mp_float::exactinit() ;

    /*--------------------------------- init. output data */
        jigsaw_init_msh_t(_mmsh) ;

    /*--------------------------------- setup *.JLOG data */
        if (_jjig != nullptr )
        {
            _jcfg._verbosity = _jjig->_verbosity ;
        }
         jlog_null _jlog(_jcfg) ;
        _jlog.push(JIGSAW::
                   asciibanner) ;

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
                 _mesh,*_imsh)) != __no_error)
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
    /*--------------------------------- initialise i.c.'s */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming INIT data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            _mesh._euclidean_mesh_2d.
                _mesh.make_link();
            _mesh._euclidean_mesh_3d.
                _mesh.make_link();

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
            if(_jcfg._mesh_opts.iter() != +0 )
            {
    /*--------------------------------- call mesh routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Generate rDT MESH...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv =
                JIGSAW ::mesh_core (
                 _jcfg, _jlog ,
                 _geom, _mesh ,
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
    /*--------------------------------- call copy routine */
            if(_jcfg._mesh_opts.iter() != +0 &&
               _jcfg._iter_opts.iter() != +0 )
            {
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Pushing MESH data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv = copy_mesh (
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
                 _jcfg, _jlog ,
                 _geom, _mesh)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv =
                JIGSAW ::iter_core (
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
                "  Writing MESH data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if (_jcfg._mesh_opts.iter() != +0 &&
                _jcfg._iter_opts.iter() == +0 )
            {

            if ((_retv = save_rdel (
                 _jcfg, _jlog ,
                 _mesh,*_mmsh)) != __no_error)
            {
                return  _retv ;
            }

            }
            else
            {

            if ((_retv = save_mesh (
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

    /*--------------------------------- init. geo. kernel */
        mp_float::exactinit() ;

    /*-------------------------- find *.JFCG file in args */
        iptr_type _retv = -1  ;
        jcfg_data _jcfg ;
        for (; _argc-- != +0; )
        {
            std::string _ssrc(_argv[_argc]) ;

            if (_ssrc.find("-h") == 0 ||
                _ssrc.find(
                       "--help") == 0 )
            {
                _retv = -2 ;

                std::cout <<
                "run jigsaw jigname.jig";
                std::cout <<  std::endl ;

                break ;
            }

            if (_ssrc.find("-v") == 0 ||
                _ssrc.find(
                    "--version") == 0 ||
                _ssrc.find(
                      "-whoami") == 0 )
            {
                _retv = -2 ;

                std::cout << __JGSWVSTR ;
                std::cout <<  std::endl ;

                break ;
            }

            std::string _path ;
            std::string _name ;
            std::string _fext ;
            file_part ( _ssrc ,
                _path , _name , _fext ) ;

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
        _jlog.push(JIGSAW::
                   asciibanner) ;

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
    /*--------------------------------- initialise i.c.'s */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Forming INIT data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            _mesh._euclidean_mesh_2d.
                _mesh.make_link();
            _mesh._euclidean_mesh_3d.
                _mesh.make_link();

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
            if(_jcfg._mesh_opts.iter() != +0 )
            {
    /*--------------------------------- call mesh routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Generate rDT MESH...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if ((_retv =
                JIGSAW ::mesh_core (
                 _jcfg, _jlog ,
                 _geom, _mesh ,
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
                 _jlog, _mesh)) != __no_error)
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
            if(_jcfg._mesh_opts.iter() != +0 &&
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
        }

        if(!_jcfg._geom_file.empty())
        {
    /*--------------------------------- call cull routine */
            _jlog.push (  __jloglndv    "\n" ) ;
            _jlog.push (
                "  Culling MESH data...\n\n" ) ;

#           ifdef  __use_timers
            _ttic   = _time.now();
#           endif//__use_timers

            if(_jcfg._mesh_opts.iter() != +0 &&
               _jcfg._iter_opts.iter() == +0 )
            {

            /*
            if ((_retv = cull_rdel (
                 _jcfg, _jlog ,
                 _geom, _mesh)) != __no_error)
            {
                return  _retv ;
            }
            */

            }
            else
            {

            /*
            if ((_retv = cull_mesh (
                 _jcfg, _jlog ,
                 _geom, _mesh)) != __no_error)
            {
                return  _retv ;
            }
            */

            }

#           ifdef  __use_timers
            _ttoc   = _time.now();
            _jlog.push(dump_time(_ttic, _ttoc));
#           endif//__use_timers
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
                 _jcfg, _jlog ,
                 _geom, _mesh)) != __no_error)
            {
                return  _retv ;
            }

            if ((_retv =
                JIGSAW ::iter_core (
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

            if (_jcfg._mesh_opts.iter() != +0 &&
                _jcfg._iter_opts.iter() == +0 )
            {

            if ((_retv = save_rdel (
                 _jcfg,
                 _jlog, _mesh)) != __no_error)
            {
                return  _retv ;
            }

            }
            else
            {

            if ((_retv = save_mesh (
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



