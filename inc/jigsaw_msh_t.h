
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
     * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant 
     * this code in any way whatsoever.  This code is 
     * provided "as-is" to be used at your own risk.
     *
     * THE CONTRIBUTORS include:
     * (a) The University of Sydney
     * (b) The Massachusetts Institute of Technology
     * (c) Columbia University
     * (d) The National Aeronautics & Space Administration
     * (e) Los Alamos National Laboratory
     *
    --------------------------------------------------------
     *
     * Last updated: 08 Feb., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __JIGSAW_MSH_T__
#   define __JIGSAW_MSH_T__

/*------------------------------------------- geom. types */

    typedef struct
        {
        real_t                  _ppos [2] ; // coord.'s
        indx_t                  _itag ;     // ID tag
        } jigsaw_VERT2_t ;

    typedef struct
        {
        real_t                  _ppos [3] ; // coord.'s
        indx_t                  _itag ;     // ID tag
        } jigsaw_VERT3_t ;

    typedef struct
        {
        indx_t                  _node [2] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_EDGE2_t ;

    typedef struct
        {
        indx_t                  _node [3] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_TRIA3_t ;

    typedef struct
        {
        indx_t                  _node [4] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_QUAD4_t ;

    typedef struct
        {
        indx_t                  _node [4] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_TRIA4_t ;

    typedef struct
        {
        indx_t                  _node [8] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_HEXA8_t ;

    typedef struct
        {
        indx_t                  _node [6] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_WEDG6_t ;

    typedef struct
        {
        indx_t                  _node [5] ; // indexing
        indx_t                  _itag ;     // ID tag
        } jigsaw_PYRA5_t ;

    typedef struct
        {
        indx_t                  _itag ;     // ID tag
        indx_t                  _indx ;     // MSH num.
        indx_t                  _kind ;     // MSH obj.
        } jigsaw_BOUND_t ;

/*------------------------------------------- array types */

    typedef struct
        {
        size_t                  _size ;
        jigsaw_VERT2_t         *_data ;
        } jigsaw_VERT2_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_VERT3_t         *_data ;
        } jigsaw_VERT3_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_EDGE2_t         *_data ;
        } jigsaw_EDGE2_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_TRIA3_t         *_data ;
        } jigsaw_TRIA3_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_QUAD4_t         *_data ;
        } jigsaw_QUAD4_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_TRIA4_t         *_data ;
        } jigsaw_TRIA4_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_HEXA8_t         *_data ;
        } jigsaw_HEXA8_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_WEDG6_t         *_data ;
        } jigsaw_WEDG6_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_PYRA5_t         *_data ;
        } jigsaw_PYRA5_array_t ;

    typedef struct
        {
        size_t                  _size ;
        jigsaw_BOUND_t         *_data ;
        } jigsaw_BOUND_array_t ;

    typedef struct
        {
        size_t                  _size ;
        indx_t                 *_data ;
        } jigsaw_INDEX_array_t ;

    typedef struct
        {
        size_t                  _size ;
        real_t                 *_data ;
        } jigsaw_REALS_array_t ;

    typedef struct
        {
        size_t                  _size ;
        fp32_t                 *_data ;
        } jigsaw_FLT32_array_t ;

/*------------------------------------------- "msh" class */

    typedef struct
        {

        indx_t                  _flags;

    /* if (_flags == EUCLIDEAN_MESH) */

        jigsaw_VERT2_array_t    _vert2;
        jigsaw_VERT3_array_t    _vert3;

        jigsaw_REALS_array_t    _power;

        jigsaw_EDGE2_array_t    _edge2;

        jigsaw_TRIA3_array_t    _tria3;
        jigsaw_QUAD4_array_t    _quad4;

        jigsaw_TRIA4_array_t    _tria4;
        jigsaw_HEXA8_array_t    _hexa8;
        jigsaw_WEDG6_array_t    _wedg6;
        jigsaw_PYRA5_array_t    _pyra5;

        jigsaw_BOUND_array_t    _bound;

        jigsaw_VERT2_array_t    _seed2;
        jigsaw_VERT3_array_t    _seed3;

    /* if (_flags == ELLIPSOID_MESH) */

        jigsaw_REALS_array_t    _radii;

    /* if (_flags == EUCLIDEAN_GRID) */
    /* OR (_flags == ELLIPSOID_GRID) */

        jigsaw_REALS_array_t    _xgrid;
        jigsaw_REALS_array_t    _ygrid;
        jigsaw_REALS_array_t    _zgrid;

    /* if (_flags == EUCLIDEAN_MESH) */
    /* OR (_flags == EUCLIDEAN_GRID) */
    /* OR (_flags == ELLIPSOID_MESH) */
    /* OR (_flags == ELLIPSOID_GRID) */

        jigsaw_FLT32_array_t    _value;

    /* for |dfdx| limiting in MARCHE */

        jigsaw_FLT32_array_t    _slope;

    /* if (_flags == EUCLIDEAN_DUAL) */

    /*!! todo: dual-mesh data here!! */


        } jigsaw_msh_t ;


#   endif   //__JIGSAW_MSH_T__



