
    /*
    --------------------------------------------------------
     * Interface to the JIGSAW meshing library.
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
     * Last updated: 31 July, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __LIB_JIGSAW__
#   define __LIB_JIGSAW__

#   include "stdint.h"

    typedef int32_t         indx_t ;    // yes, 32 bit ints,
                                        // but, 64 bit ptrs!
    typedef double          real_t ;

#   include "jigsaw_jig_t.h"
#   include "jigsaw_msh_t.h"

    /* 
    --------------------------------------------------------
     * return codes for JIGSAW.
    --------------------------------------------------------
     */ 

#   define JIGSAW_UNKNOWN_ERROR    -1

#   define JIGSAW_NO_ERROR         +0

#   define JIGSAW_FILE_NOT_LOCATED +2
#   define JIGSAW_FILE_NOT_CREATED +3

#   define JIGSAW_INVALID_ARGUMENT +4

    /* 
    --------------------------------------------------------
     * constants for JIGSAW.
    --------------------------------------------------------
     */ 

#   define JIGSAW_NULL_FLAG      -100

#   define JIGSAW_EUCLIDEAN_MESH +100
#   define JIGSAW_EUCLIDEAN_GRID +101
#   define JIGSAW_EUCLIDEAN_DUAL +102
    
#   define JIGSAW_ELLIPSOID_MESH +200
#   define JIGSAW_ELLIPSOID_GRID +201
#   define JIGSAW_ELLIPSOID_DUAL +202

#   define JIGSAW_HFUN_RELATIVE  +300
#   define JIGSAW_HFUN_ABSOLUTE  +301
    
#   define JIGSAW_KERN_DELFRONT  +400
#   define JIGSAW_KERN_DELAUNAY  +401

#   define JIGSAW_BNDS_TRIACELL  +402
#   define JIGSAW_BNDS_DUALCELL  +403

    /* 
    --------------------------------------------------------
     * generate mesh via JIGSAW.
    --------------------------------------------------------
     */ 
    
    extern indx_t jigsaw_make_mesh (    
    
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
        
    /* HFUN (OPTIONAL): mesh-spacing function H(x). 
     * => NULL for empty H(x). 
     */
        jigsaw_msh_t   *_hfun,
        
    /* MESH (REQUIRED): mesh data-structure.
     * Defines initial complex upon entry, returns 
     * ouput at function exit.
     */
        jigsaw_msh_t   *_mesh
        ) ;

    /* 
    --------------------------------------------------------
     * setup objects for JIGSAW.
    --------------------------------------------------------
     */
        
    extern void   jigsaw_init_msh_t (
        jigsaw_msh_t   *_mesh
        ) ;

    extern void   jigsaw_init_jig_t (
        jigsaw_jig_t   *_jjig
        ) ;

    /* 
    --------------------------------------------------------
     * alloc objects for JIGSAW.
    --------------------------------------------------------
     */
    
    extern void  jigsaw_alloc_vert2 (
        jigsaw_VERT2_array_t *_xsrc ,
        indx_t _size
        ) ;
        
    extern void  jigsaw_alloc_vert3 (
        jigsaw_VERT3_array_t *_xsrc ,
        indx_t _size
        ) ;
        
    extern void  jigsaw_alloc_edge2 (
        jigsaw_EDGE2_array_t *_xsrc ,
        indx_t _size
        ) ;
        
    extern void  jigsaw_alloc_tria3 (
        jigsaw_TRIA3_array_t *_xsrc ,
        indx_t _size
        ) ;
        
    extern void  jigsaw_alloc_tria4 (
        jigsaw_TRIA4_array_t *_xsrc ,
        indx_t _size
        ) ;
        
    extern void  jigsaw_alloc_reals (
        jigsaw_REALS_array_t *_xsrc ,
        indx_t _size
        ) ;
        
    extern void   jigsaw_free_msh_t (
        jigsaw_msh_t   *_mesh
        ) ;
    
    /* 
    --------------------------------------------------------
     * parse-to-file for JIGSAW.
    --------------------------------------------------------
     */
      
    /*
    extern indx_t jigsaw_save_msh_t (
        char *_file,
        jigsaw_msh_t   *_mesh
        ) ;
     */

    extern indx_t jigsaw_save_jig_t (
        char *_file,
        jigsaw_jig_t   *_jjig
        ) ;
        
    extern indx_t jigsaw_load_msh_t (
        char *_file,
        jigsaw_msh_t   *_mesh
        ) ;

    extern indx_t jigsaw_load_jig_t (
        char *_file,
        jigsaw_jig_t   *_jjig
        ) ;
     
     
#   endif   // __LIB_JIGSAW__



