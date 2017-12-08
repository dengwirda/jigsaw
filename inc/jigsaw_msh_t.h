
    /*
    --------------------------------------------------------
     * JIGSAW's "msh" mesh type definition.
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
     * Last updated: 14 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __JIGSAW_MSH_T__
#   define __JIGSAW_MSH_T__

/*------------------------------------------- geom. types */

    typedef struct 
        {
        real_t                  _ppos [2] ;
        indx_t                  _itag ;
        } jigsaw_VERT2_t ;
        
    typedef struct
        {
        real_t                  _ppos [3] ;
        indx_t                  _itag ;
        } jigsaw_VERT3_t ;
    
    typedef struct
        {
        indx_t                  _node [2] ;
        indx_t                  _itag ;
        } jigsaw_EDGE2_t ;
    
    typedef struct
        {
        indx_t                  _node [3] ;
        indx_t                  _itag ;
        } jigsaw_TRIA3_t ;
        
    typedef struct
        {
        indx_t                  _node [4] ;
        indx_t                  _itag ;
        } jigsaw_TRIA4_t ;

/*------------------------------------------- array types */    
    
    typedef struct
        {
        indx_t                  _size ;
        jigsaw_VERT2_t         *_data ;
        } jigsaw_VERT2_array_t ;
        
    typedef struct
        {
        indx_t                  _size ;
        jigsaw_VERT3_t         *_data ;
        } jigsaw_VERT3_array_t ;
        
    typedef struct
        {
        indx_t                  _size ;
        jigsaw_EDGE2_t         *_data ;
        } jigsaw_EDGE2_array_t ;
        
    typedef struct
        {
        indx_t                  _size ;
        jigsaw_TRIA3_t         *_data ;
        } jigsaw_TRIA3_array_t ;
        
    typedef struct
        {
        indx_t                  _size ;
        jigsaw_TRIA4_t         *_data ;
        } jigsaw_TRIA4_array_t ;
        
    typedef struct
        {
        indx_t                  _size ;
        real_t                 *_data ;
        } jigsaw_REALS_array_t ;
    
/*------------------------------------------- "msh" class */

    typedef struct 
        {
    
        indx_t                  _flags;
    
    /* if (_flags == EUCLIDEAN_MESH) */
    
        jigsaw_VERT2_array_t    _vert2;
        jigsaw_VERT3_array_t    _vert3;
    
        jigsaw_EDGE2_array_t    _edge2;
        jigsaw_TRIA3_array_t    _tria3;
        jigsaw_TRIA4_array_t    _tria4;
        
        
    /* if (_flags == ELLIPSOID_MESH) */
        
        jigsaw_REALS_array_t    _radii;
        
       
    /* if (_flags == EUCLIDEAN_GRID) */
       
        jigsaw_REALS_array_t    _xgrid;
        jigsaw_REALS_array_t    _ygrid;
        jigsaw_REALS_array_t    _zgrid;


    /* if (_flags == EUCLIDEAN_MESH) */
    /* OR (_flags == EUCLIDEAN_GRID) */
    
        jigsaw_REALS_array_t    _value;
       
       
    /* if (_flags == EUCLIDEAN_DUAL) */
    
    /*!! todo: dual-mesh data here!! */

    
        } jigsaw_msh_t ;
   
          
#   endif   //__JIGSAW_MSH_T__      
        
        

