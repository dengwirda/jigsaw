
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 10 Jul., 2021
     *
     * Copyright 2013-2021
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
    class iter_params
        {
        public  :

        typedef R                       real_type ;
        typedef I                       iptr_type ;

        typedef iter_params<R, I>       self_type ;

        iptr_type        _verb ;

        iptr_type        _iter ;

        real_type        _qtol ;
        real_type        _qlim ;

        bool_type        _zip_ ;
        bool_type        _div_ ;
        bool_type        _tria ;
        bool_type        _dual ;

        public  :

    /*-------------------------- construct default param. */
        __inline_call iter_params (
            )
        {   // load default values
            jigsaw_jig_t _jjig ;
            jigsaw_null_jig_t( &_jjig ) ;

            this->_verb = 
            iptr_type(_jjig._verbosity) ;
            
            this->_iter = 
            iptr_type(_jjig._optm_iter) ;

            this->_qtol =
            real_type(_jjig._optm_qtol) ;
            this->_qlim =
            real_type(_jjig._optm_qlim) ;

            this->_zip_ =
            bool_type(_jjig._optm_zip_) ;
            this->_div_ =
            bool_type(_jjig._optm_div_) ;
            this->_tria =
            bool_type(_jjig._optm_tria) ;
            this->_dual =
            bool_type(_jjig._optm_dual) ;
        }

    /*------------------------------------ "write" access */
        __inline_call iptr_type      & verb (
            )
        {   return  this->_verb ;
        }

        __inline_call iptr_type      & iter (
            )
        {   return  this->_iter ;
        }

        __inline_call real_type      & qtol (
            )
        {   return  this->_qtol ;
        }
        __inline_call real_type      & qlim (
            )
        {   return  this->_qlim ;
        }

        __inline_call bool_type      & zip_ (
            )
        {   return  this->_zip_ ;
        }
        __inline_call bool_type      & div_ (
            )
        {   return  this->_div_ ;
        }
        __inline_call bool_type      & tria (
            )
        {   return  this->_tria ;
        }
        __inline_call bool_type      & dual (
            )
        {   return  this->_dual ;
        }

    /*------------------------------------ "const" access */
        __inline_call iptr_type const& verb (
            ) const
        {   return  this->_verb ;
        }

        __inline_call iptr_type const& iter (
            ) const
        {   return  this->_iter ;
        }

        __inline_call real_type const& qtol (
            ) const
        {   return  this->_qtol ;
        }
        __inline_call real_type const& qlim (
            ) const
        {   return  this->_qlim ;
        }

        __inline_call bool_type const& zip_ (
            ) const
        {   return  this->_zip_ ;
        }
        __inline_call bool_type const& div_ (
            ) const
        {   return  this->_div_ ;
        }
        __inline_call bool_type const& tria (
            ) const
        {   return  this->_tria ;
        }
        __inline_call bool_type const& dual (
            ) const
        {   return  this->_dual ;
        }

        } ;

    }

#   endif   //__ITER_PARAMS__



