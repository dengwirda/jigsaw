
    /*
    --------------------------------------------------------
     * ITER-PARAMS: parameters for ITER-MESH-K.
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
     * Last updated: 11 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_PARAMS__
#   define __ITER_PARAMS__

#   include "../../liblib/init_jig_t.hpp"

    namespace mesh {

    /*
    --------------------------------------------------------
     * ITER-PARAMS: user-param's for ITER-MESH-K
    --------------------------------------------------------
     */

    template <
    typename R ,
    typename I
             >
    class iter_params : public jigsaw_jig_t
        {
        public  :

        typedef R                       real_type ;
        typedef I                       iptr_type ;

        typedef iter_params<R, I>       self_type ;

        public  :

    /*-------------------------- construct default param. */
        __inline_call iter_params (
            )
        {   // load default values
            jigsaw_null_jig_t(this) ;
        }

    /*------------------------------------ "write" access */
        __inline_call iptr_type      & verb (
            )
        {   return  this->_verbosity ;
        }
        __inline_call iptr_type      & nprt (
            )
        {   return  this->_numthread ;
        }

        __inline_call iptr_type      & iter (
            )
        {   return  this->_optm_iter ;
        }

        __inline_call real_type      & beta (
            )
        {   return  this->_optm_beta ;
        }
        __inline_call real_type      & zeta (
            )
        {   return  this->_optm_zeta ;
        }

        __inline_call real_type      & qtol (
            )
        {   return  this->_optm_qtol ;
        }
        __inline_call real_type      & qlim (
            )
        {   return  this->_optm_qlim ;
        }

        __inline_call real_type      & wmin (
            )
        {   return  this->_optm_wmin ;
        }
        __inline_call real_type      & wmax (
            )
        {   return  this->_optm_wmax ;
        }

        __inline_call iptr_type      & zip_ (
            )
        {   return  this->_optm_zip_ ;
        }
        __inline_call iptr_type      & div_ (
            )
        {   return  this->_optm_div_ ;
        }
        __inline_call iptr_type      & tria (
            )
        {   return  this->_optm_tria ;
        }
        __inline_call iptr_type      & dual (
            )
        {   return  this->_optm_dual ;
        }

    /*------------------------------------ "const" access */
        __inline_call iptr_type const& verb (
            ) const
        {   return  this->_verbosity ;
        }
        __inline_call iptr_type const& nprt (
            ) const
        {   return  this->_numthread ;
        }

        __inline_call iptr_type const& iter (
            ) const
        {   return  this->_optm_iter ;
        }

        __inline_call real_type const& beta (
            ) const
        {   return  this->_optm_beta ;
        }
        __inline_call real_type const& zeta (
            ) const
        {   return  this->_optm_zeta ;
        }

        __inline_call real_type const& qtol (
            ) const
        {   return  this->_optm_qtol ;
        }
        __inline_call real_type const& qlim (
            ) const
        {   return  this->_optm_qlim ;
        }

        __inline_call real_type const& wmin (
            ) const
        {   return  this->_optm_wmin ;
        }
        __inline_call real_type const& wmax (
            ) const
        {   return  this->_optm_wmax ;
        }

        __inline_call iptr_type const& zip_ (
            ) const
        {   return  this->_optm_zip_ ;
        }
        __inline_call iptr_type const& div_ (
            ) const
        {   return  this->_optm_div_ ;
        }
        __inline_call iptr_type const& tria (
            ) const
        {   return  this->_optm_tria ;
        }
        __inline_call iptr_type const& dual (
            ) const
        {   return  this->_optm_dual ;
        }

        } ;

    }

#   endif   //__ITER_PARAMS__



