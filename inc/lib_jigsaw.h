
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
     * JIGSAW: Interface to the JIGSAW meshing library.
    --------------------------------------------------------
     *
     * Last updated: 19 June, 2019
     *
     * Copyright 2013 -- 2019
     * Darren Engwirda
     * d.engwirda@gmail.com
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
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __LIB_JIGSAW__
#   define __LIB_JIGSAW__

#   ifdef  __cplusplus
    extern "C" {
#   endif

#   ifdef  _WIN32

#   ifdef  __lib_jigsaw
#       define SHARED __declspec(dllexport)
#   else
#       define SHARED __declspec(dllimport)
#   endif

#   else
#       define SHARED
#   endif

#   include "stdint.h"
#   include "stddef.h"

    typedef int32_t         indx_t ;    // yes, 32 bit ints,
                                        // but, 64 bit ptrs!
    typedef double          real_t ;

#   include "jigsaw_const.h"

#   include "jigsaw_jig_t.h"
#   include "jigsaw_msh_t.h"

    /*
    --------------------------------------------------------
     * generate mesh via JIGSAW.
    --------------------------------------------------------
     */

#   define jigsaw_make_mesh jigsaw

    SHARED indx_t jigsaw (

    /* JCFG (REQUIRED): settings obj. definition.
     */
        jigsaw_jig_t   *_jcfg,

    /* GEOM (REQUIRED): geometry obj. definition.
     */
        jigsaw_msh_t   *_geom,

    /* INIT (OPTIONAL): initial tria. definition.
     * => NULL for empty IC's.
     */
        jigsaw_msh_t   *_init,

    /* HFUN (OPTIONAL): mesh-spacing function h(x).
     * => NULL for empty h(x).
     */
        jigsaw_msh_t   *_hfun,

    /* MESH (REQUIRED): output mesh data-structure.
     */
        jigsaw_msh_t   *_mesh
        ) ;

    /*
    --------------------------------------------------------
     * compute rDT's via TRIPOD.
    --------------------------------------------------------
     */

    SHARED indx_t tripod (

    /* JCFG (REQUIRED): settings obj. definition.
     */
        jigsaw_jig_t   *_jcfg,

    /* INIT (REQUIRED): initial point definition.
     */
        jigsaw_msh_t   *_init,

    /* GEOM (OPTIONAL): geometry obj. definition.
     * => NULL for empty GEOM.
     */
        jigsaw_msh_t   *_geom,

    /* MESH (REQUIRED): output mesh data-structure.
     */
        jigsaw_msh_t   *_mesh
        ) ;

    /*
    --------------------------------------------------------
     * limit |df/dx| via MARCHE.
    --------------------------------------------------------
     */

    SHARED indx_t marche (

    /* JCFG (REQUIRED): settings obj. definition.
     */
        jigsaw_jig_t   *_jcfg,

    /* FFUN (REQUIRED): apply limiter to |df/dx|.
     */
        jigsaw_msh_t   *_ffun
        ) ;

    /*
    --------------------------------------------------------
     * setup objects for JIGSAW.
    --------------------------------------------------------
     */

    SHARED void   jigsaw_init_msh_t (
        jigsaw_msh_t   *_mesh
        ) ;

    SHARED void   jigsaw_init_jig_t (
        jigsaw_jig_t   *_jjig
        ) ;

    /*
    --------------------------------------------------------
     * parse-to-file for JIGSAW.
    --------------------------------------------------------
     */

    /*
    SHARED indx_t jigsaw_save_msh_t (
        char *_file,
        jigsaw_msh_t   *_mesh
        ) ;
     */

    SHARED indx_t jigsaw_save_jig_t (
        char *_file,
        jigsaw_jig_t   *_jjig
        ) ;

    SHARED indx_t jigsaw_load_msh_t (
        char *_file,
        jigsaw_msh_t   *_mesh
        ) ;

    SHARED indx_t jigsaw_load_jig_t (
        char *_file,
        jigsaw_jig_t   *_jjig
        ) ;

    /*
    --------------------------------------------------------
     * alloc objects for JIGSAW.
    --------------------------------------------------------
     */

    SHARED void  jigsaw_alloc_vert2 (
        jigsaw_VERT2_array_t *_xsrc ,
        size_t _size
        ) ;

    SHARED void  jigsaw_alloc_vert3 (
        jigsaw_VERT3_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_edge2 (
        jigsaw_EDGE2_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_tria3 (
        jigsaw_TRIA3_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_quad4 (
        jigsaw_QUAD4_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_tria4 (
        jigsaw_TRIA4_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_hexa8 (
        jigsaw_HEXA8_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_wedg6 (
        jigsaw_WEDG6_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_pyra5 (
        jigsaw_PYRA5_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_bound (
        jigsaw_BOUND_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_index (
        jigsaw_INDEX_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void  jigsaw_alloc_reals (
        jigsaw_REALS_array_t *_xsrc ,
		size_t _size
        ) ;

    SHARED void   jigsaw_free_msh_t (
        jigsaw_msh_t   *_mesh
        ) ;

    SHARED void   jigsaw_free_vert2 (
        jigsaw_VERT2_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_vert3 (
        jigsaw_VERT3_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_edge2 (
        jigsaw_EDGE2_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_tria3 (
        jigsaw_TRIA3_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_quad4 (
        jigsaw_QUAD4_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_tria4 (
        jigsaw_TRIA4_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_hexa8 (
        jigsaw_HEXA8_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_wedg6 (
        jigsaw_WEDG6_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_pyra5 (
        jigsaw_PYRA5_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_bound (
        jigsaw_BOUND_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_index (
        jigsaw_INDEX_array_t *_xsrc
        ) ;

    SHARED void   jigsaw_free_reals (
        jigsaw_REALS_array_t *_xsrc
        ) ;

#   ifdef  __cplusplus
    }
#   endif

#   undef   SHARED

#   endif   // __LIB_JIGSAW__



