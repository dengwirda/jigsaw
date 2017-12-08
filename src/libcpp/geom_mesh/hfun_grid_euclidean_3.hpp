
    /*
    --------------------------------------------------------
     * HFUN-GRID-EUCLIDEAN-kD: structured H(X) in R^k.
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
     * Last updated: 23 August, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_GRID_EUCLIDEAN_3__
#   define __HFUN_GRID_EUCLIDEAN_3__

    namespace mesh {

    template <
    typename R ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_grid_euclidean_3d 
        : public hfun_base_kd <I, R>
    {
    public  :
    
    /*---------------------- "grid"-based size-fun in R^3 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef hfun_grid_euclidean_3d  <
            real_type ,
            iptr_type >             hfun_type ;

    typedef typename  hfun_base_kd  <
            iptr_type , 
            real_type >::hint_type  hint_type ;
      
      
    containers::array <
        real_type, allocator>      _xpos;
    containers::array <
        real_type, allocator>      _ypos;
    containers::array <
        real_type, allocator>      _zpos;
    
    containers::array <
        real_type, allocator>      _hmat; 
        
    bool_type                      _xvar;
    bool_type                      _yvar;
    bool_type                      _zvar;    
    
    public  :
    
    __inline_call void_type indx_from_subs (
        iptr_type _ipos,
        iptr_type _jpos,
        iptr_type&_indx
        )
    {
    }
    
    /*
    --------------------------------------------------------
     * EVAL: eval. size-fun. value.
    --------------------------------------------------------
     */
    
    __inline_call real_type eval (
        real_type *_ppos,
        hint_type &_hint
        )
    {
        real_type _hval = 
    +std::numeric_limits<real_type>::infinity();
    
        __unreferenced (_hint) ;
    
        if (this->_xpos.count() == +0)
            return _hval ;
        
        if (_ppos[0] < *this->_xpos.head() )
            return _hval ;
        if (_ppos[0] > *this->_xpos.tail() )
            return _hval ;
        
        if (this->_ypos.count() == +0)
            return _hval ;
            
        if (_ppos[1] < *this->_ypos.head() )
            return _hval ;
        if (_ppos[1] > *this->_ypos.tail() )
            return _hval ;
            
        if (this->_zpos.count() == +0)
            return _hval ;
            
        if (_ppos[2] < *this->_zpos.head() )
            return _hval ;
        if (_ppos[2] > *this->_zpos.tail() )
            return _hval ;
    
    
    
        return _hval ;
    }    
    
    
    } ;
    
    
    }
    
#   endif   //__HFUN_GRID_EUCLIDEAN_3__



