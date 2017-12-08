
    /*
    --------------------------------------------------------
     * HFUN-BASE-K: base class for H(x) functions. 
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
     * Last updated: 03 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
#   pragma once

#   ifndef __HFUN_BASE_K__
#   define __HFUN_BASE_K__

    namespace mesh {

    template <
    typename I ,
    typename R 
             >
    class hfun_base_kd
    {
    public :
     
    /*-------------------- base class for size-fun in R^k */
    
    typedef R               real_type ;
    typedef I               iptr_type ;
    
    typedef iptr_type       hint_type ;
    
    __static_call
    __inline_call hint_type null_hint (
        )
    {   return ( (hint_type) -1 );
    }
    
    public  :

    /*-------------------- forward dec.'s for sub-classes */
    
    __normal_call void_type init (
        )  {}
        
    __normal_call void_type eval (
        real_type *_ppos ,
        hint_type &_hint
        )  {}
        
    __inline_call real_type eval (
        real_type *_ppos ,
        hint_type&&_hint = null_hint ()
        )
    {   return eval(_ppos, _hint);
    }
    
    } ;
    
    
    }
    
#   endif   //__HFUN_BASE_K__



