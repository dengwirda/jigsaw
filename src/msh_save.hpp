
    /*
    --------------------------------------------------------
     * MSH-SAVE: parse MESH data into *.MSH / MSH_t.
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
     * Last updated: 28 May, 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_SAVE__
#   define __MSH_SAVE__

    class mshfile_tag {} ;
    class _netcdf_tag {} ;

    /*
    --------------------------------------------------------
     * for C-style sprintf to C++ stream.
    --------------------------------------------------------
     */

#   define PRINTCHUNK   (+32768)

#   define VERT2CHUNK   (+ 512 )
#   define VERT3CHUNK   (+ 320 )

#   define VALUECHUNK   (+1024 )

#   define EDGE2CHUNK   (+ 832 )
#   define TRIA3CHUNK   (+ 640 )
#   define QUAD4CHUNK   (+ 512 )
#   define TRIA4CHUNK   (+ 512 )
#   define HEXA8CHUNK   (+ 288 )
#   define WEDG6CHUNK   (+ 352 )
#   define PYRA5CHUNK   (+ 416 )

#   define BOUNDCHUNK   (+ 832 )

#   define PRINTFINAL                   \
        {                           \
            _fbuf[_next] ='\0';     \
            _file << _fbuf;         \
            _roll = _next = 0 ;     \
        }

#   define PRINTCHARS(_char, _size) \
        {                           \
            _next +=  _char;        \
                                    \
            if ((_roll++) >= _size)     \
            {                       \
                PRINTFINAL ;        \
            }                       \
        }

    /*
    --------------------------------------------------------
     * inc. low-level file write kernels.
    --------------------------------------------------------
     */

#   include  "rdt_save.inc"
#   include  "tri_save.inc"
#   include  "msh_save.inc"
#   include  "hfn_save.inc"


#   endif   //__MSH_SAVE__



