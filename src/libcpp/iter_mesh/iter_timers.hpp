
    /*
    --------------------------------------------------------
     * ITER-TIMERS: timers for ITER-MESH-K.
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
     * Last updated: 20 Jan., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_TIMERS__
#   define __ITER_TIMERS__

    namespace mesh {

    /*
    --------------------------------------------------------
     * ITER-TIMERS: cpu timers for ITER-MESH-K
    --------------------------------------------------------
     */

    template <
    typename R ,
    typename I
             >
    class iter_timers
        {
        public  :

        typedef R                       real_type ;
        typedef I                       iptr_type ;

        typedef iter_timers<R, I>       self_type ;

        double      _iter_full = (double   )  +0. ;

        double      _move_node = (double   )  +0. ;
        double      _init_node = (double   )  +0. ;
        double      _core_node = (double   )  +0. ;
        double      _ldir_node = (double   )  +0. ;
        double      _lopt_node = (double   )  +0. ;

        double      _move_dual = (double   )  +0. ;
        double      _init_dual = (double   )  +0. ;
        double      _core_dual = (double   )  +0. ;
        double      _ldir_dual = (double   )  +0. ;
        double      _lopt_dual = (double   )  +0. ;

        double      _topo_flip = (double   )  +0. ;
        double      _init_flip = (double   )  +0. ;
        double      _core_flip = (double   )  +0. ;

        double      _topo_zips = (double   )  +0. ;
        double      _init_zips = (double   )  +0. ;
        double      _core_zips = (double   )  +0. ;
        double      _core_divs = (double   )  +0. ;

        public  :

    /*-------------------------------------- elapsed time */

    #   ifdef  __use_timers

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


#   endif   //__ITER_TIMERS__



