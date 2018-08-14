
/* 
------------------------------------------------------------
 * AABB-PRED-K: search pred. for AABB-trees in R^k. 
------------------------------------------------------------
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
------------------------------------------------------------
 *
 * Last updated: 10 September, 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __AABB_PRED_K__
#   define __AABB_PRED_K__

    namespace geom_tree {
    

    template <
    typename R ,
    typename I ,
    size_t   K
             >
    class aabb_pred_node_k
    {
/*---------------------- node-tree intersection predicate */
    public  :
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    
    iptr_type static constexpr     _dims = K;

    private :
    
    real_type                      _ppos [ K] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_node_k (
    __write_ptr  (real_type) _psrc
        ) 
    {   
    /*----------------------------------- make local copy */
        for (auto _idim = _dims; _idim-- != +0; )
        {
            this->_ppos[_idim] = _psrc[_idim] ;
        }
    }
    
/*----------------------- TRUE for node-aabb intersection */
    __normal_call bool_type operator() (
    __const_ptr  (real_type) _bmin,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*----------------------- test bounds along each axis */
        for (auto _idim = _dims; _idim-- != +0; )
        {
        if (_bmin[_idim] > this->_ppos[_idim] ||
            _bmax[_idim] < this->_ppos[_idim] )
        return ( false ) ;
        }
        
    /*----------------------- intersecting if we got here */
        return (  true ) ;
    }
    
    } ;

    
    template <
    typename R ,
    typename I ,
    size_t   K
             >
    class aabb_pred_rect_k
    {
/*---------------------- rect-tree intersection predicate */
    public  :
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    
    iptr_type static constexpr     _dims = K;

    private :
    
    real_type                      _rmin [ K] ;
    real_type                      _rmax [ K] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_rect_k (
    __write_ptr  (real_type) _lsrc ,
    __write_ptr  (real_type) _rsrc
        ) 
    {   
    /*----------------------------------- make local copy */
        for (auto _idim = _dims; _idim-- != +0; )
        {
            this->_rmin[_idim] = _lsrc[_idim] ;
            this->_rmax[_idim] = _rsrc[_idim] ;
        }
    }
    
/*----------------------- TRUE for rect-aabb intersection */
    __normal_call bool_type operator() (
    __const_ptr  (real_type) _bmin,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*----------------------- test bounds along each axis */
        for (auto _idim = _dims; _idim-- != +0; )
        {
        if (_bmin[_idim] > this->_rmax[_idim] ||
            _bmax[_idim] < this->_rmin[_idim] )
        return ( false ) ;
        }
        
    /*----------------------- intersecting if we got here */
        return (  true ) ;
    }
    
    } ;
    
    
    template <
    typename R ,
    typename I ,
    size_t   K
             >
    class aabb_pred_line_k
    {
/*---------------------- line-tree intersection predicate */
    public  :
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    
    iptr_type static constexpr     _dims = K;

    private :
    
    real_type                      _ipos [ K] ;
    real_type                      _jpos [ K] ;

    real_type                      _xden [ K] ;

    real_type                      _rmin [ K] ;
    real_type                      _rmax [ K] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_line_k (
    __write_ptr  (real_type) _isrc ,
    __write_ptr  (real_type) _jsrc
        ) 
    {   
    /*------------------------------ init bbox at -+ inf. */
        for (auto _idim = _dims; _idim-- != +0; )
        {
            _rmin[_idim] = 
               +std::numeric_limits
                    <real_type>::infinity() ;
                    
            _rmax[_idim] = 
                -std::numeric_limits
                    <real_type>::infinity() ;
        }

    /*------------------------------ copy and calc. bbox. */
        for (auto _idim = _dims; _idim-- != +0; )
        {
            this->_ipos[_idim] = _isrc[_idim] ;
            this->_jpos[_idim] = _jsrc[_idim] ;

            this->_xden[_idim] = _jsrc[_idim] - 
                                 _isrc[_idim] ;

            real_type _xmin = std::min (
                _isrc[_idim],
                _jsrc[_idim]) ;
                
            real_type _xmax = std::max (
                _isrc[_idim],
                _jsrc[_idim]) ;

            if (this->_rmin[_idim] > _xmin)
                this->_rmin[_idim] = _xmin;

            if (this->_rmax[_idim] < _xmax)
                this->_rmax[_idim] = _xmax;
        }
    }

/*----------------------- TRUE for line-aabb intersection */
    __normal_call bool_type operator() (
    __const_ptr  (real_type) _bmin,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*-------------------------------- test bbox. overlap */
        for (auto _idim = _dims; _idim-- != +0; )
        {
        if (_bmin[_idim] > this->_rmax[_idim] ||
            _bmax[_idim] < this->_rmin[_idim] )
    /*----------------------------------- can't intersect */
            return ( false )  ;
        }

    /*-------------------------------- endpoints enclosed */
        bool_type _ends[+2];
        _ends[+0] =  true;
        _ends[+1] =  true;
        
        for (auto _idim = _dims; _idim-- != +0; )
        {
    /*------------------------------------- test i-th end */
        if (_bmin[_idim] > this->_ipos[_idim] ||
            _bmax[_idim] < this->_ipos[_idim] )
            _ends[+0] = false ;
            
    /*------------------------------------- test j-th end */
        if (_bmin[_idim] > this->_jpos[_idim] ||
            _bmax[_idim] < this->_jpos[_idim] )
            _ends[+1] = false ;
        }
        
        if (_ends[+0] || _ends[+1])
        {
    /*------------------------------------- quick return! */
            return (  true )  ;
        }

    /*------------------ test for line-face intersections */
        real_type _tmin = 
            -std::numeric_limits
                <real_type>::infinity() ;
                
        real_type _tmax = 
            +std::numeric_limits
                <real_type>::infinity() ;

        for (auto _idim = _dims; _idim-- != +0; )
        {
            real_type _ival=(_bmin[_idim] - 
                       this->_ipos[_idim])/
                       this->_xden[_idim] ; 
                         
            real_type _jval=(_bmax[_idim] - 
                       this->_ipos[_idim])/
                       this->_xden[_idim] ;

            _tmin = std::max(_tmin, 
                    std::min(_ival, _jval)) ;
                    
            _tmax = std::min(_tmax, 
                    std::max(_ival, _jval)) ; 
        }
        
        return ( _tmax >= _tmin ) ;
    }
    
    } ;

    
    }
    
#   endif   //__AABB_PRED_K__



