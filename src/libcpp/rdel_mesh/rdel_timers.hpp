
    /*
    --------------------------------------------------------
     * RDEL-TIMERS: CPU timers for RDEL-MESH-K.
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
     * Last updated: 20 Jul., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_TIMERS__
#   define __RDEL_TIMERS__

    namespace mesh {

    /*
    --------------------------------------------------------
     * RDEL-TIMERS: cpu timers for RDEL-MESH-K
    --------------------------------------------------------
     */

    class rdel_timers
        {
        public  :
        double      _mesh_seed = (double   )  +0. ;

        double      _node_init = (double   )  +0. ;
        double      _node_rule = (double   )  +0. ;
        double      _edge_init = (double   )  +0. ;
        double      _edge_rule = (double   )  +0. ;
        double      _face_init = (double   )  +0. ;
        double      _face_rule = (double   )  +0. ;
        double      _tria_init = (double   )  +0. ;
        double      _tria_rule = (double   )  +0. ;

        double      _list_trim = (double   )  +0. ;

        double      _topo_init = (double   )  +0. ;

        public  :

    /*-------------------------------------- elapsed time */

    #   ifdef  __use_timers

        __static_call
        __inline_call double time_span (
            typename std::
                chrono::high_resolution_clock
                    ::time_point const& _ttic,
            typename std::
                chrono::high_resolution_clock
                    ::time_point const& _ttoc
            )
        {
            return (double)(
                std::chrono::duration_cast<
                std::chrono::microseconds >
                (_ttoc-_ttic).count()) / +1.0E+06 ;
        }

        __static_call
        __inline_call double nano_span (
            typename std::
                chrono::high_resolution_clock
                    ::time_point const& _ttic,
            typename std::
                chrono::high_resolution_clock
                    ::time_point const& _ttoc
            )
        {
            return (double)(
                std::chrono::duration_cast<
                std::chrono:: nanoseconds >
                (_ttoc-_ttic).count()) / +1.0E+09 ;
        }

    #   endif//__use_timers

        } ;

    }

#   endif // __RDEL_TIMERS__



