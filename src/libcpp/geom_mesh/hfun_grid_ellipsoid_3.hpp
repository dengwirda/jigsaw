
    /*
    --------------------------------------------------------
     * HFUN-GRID-ELLIPSOID-3D: ellipsoidal H(x) in R^3.
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
     * Last updated: 16 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
#   pragma once

#   ifndef __HFUN_GRID_ELLIPSOID_3__
#   define __HFUN_GRID_ELLIPSOID_3__

    namespace mesh {

    template <
    typename I ,
    typename R ,
    typename A = allocators::basic_alloc
             >
    class hfun_grid_ellipsoid_3d
        : public hfun_base_kd <I, R>
    {
    public :
     
    /*-------------- lat.-lon. spheroidal size-fun in R^3 */
    
    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;
    
    typedef hfun_grid_ellipsoid_3d  <
            iptr_type ,
            real_type >             hfun_type ;
                
    typedef typename  hfun_base_kd  <
            iptr_type , 
            real_type >::hint_type  hint_type ;
        
    public  :
    
    containers::array <
        real_type, allocator>      _xpos;
    containers::array <
        real_type, allocator>      _ypos;
    
    containers::array <
        real_type, allocator>      _hmat; 
        
    bool_type                      _xvar;
    bool_type                      _yvar;
               
    public  :

    __inline_call void_type indx_from_subs (
        iptr_type _ipos,
        iptr_type _jpos,
        iptr_type&_indx
        )
    {
        iptr_type _ynum = 
       (iptr_type)this->_ypos.count() ;

        _indx = _jpos * _ynum + _ipos ;
    }
    
    /*
    --------------------------------------------------------
     * INIT: init. size-fun. class.
    --------------------------------------------------------
     */
    
    __inline_call void_type init (
        )
    {
        real_type static const _FTOL = 
            std::pow(std::numeric_limits
       <real_type>::epsilon(), (real_type).8);
    
        this->_xvar = false ;
        this->_yvar = false ;
    
        if (this->_xpos.empty()) return ;
        if (this->_ypos.empty()) return ;
   
        real_type _xbar, _ybar ;
        _xbar = *this->_xpos.tail() - 
                *this->_xpos.head() ;
        
        _xbar /=(this->_xpos.count () - 1) ;
        
        real_type _xmin, _xmax ;
        _xmin = _xbar - _FTOL * _xbar ;
        _xmax = _xbar + _FTOL * _xbar ;
        
        _ybar = *this->_ypos.tail() - 
                *this->_ypos.head() ;

        _ybar /=(this->_ypos.count () - 1) ;
        
        real_type _ymin, _ymax ;
        _ymin = _ybar - _FTOL * _ybar ;
        _ymax = _ybar + _FTOL * _ybar ;
        
        for (auto _iter  = this->_xpos.head();
                  _iter != this->_xpos.tail();
                ++_iter  )
        {
            real_type _xdel =*(_iter+1) - 
                             *(_iter+0) ;
        
            if (_xdel < _xmin && 
                _xdel > _xmax ) 
            {
                _xvar =  true ; break ;
            }
        }
         
        for (auto _iter  = this->_ypos.head();
                  _iter != this->_ypos.tail();
                ++_iter  )
        {
            real_type _ydel =*(_iter+1) - 
                             *(_iter+0) ;
        
            if (_ydel < _ymin && 
                _ydel > _ymax ) 
            {
                _yvar =  true ; break ;
            }
        }       
    
    }

    /*
    --------------------------------------------------------
     * EVAL: eval. size-fun. value.
    --------------------------------------------------------
     */
    
    __normal_call real_type eval (
        real_type *_ppos ,
        hint_type &_hint
        )
    {   
        real_type static const FT = 
            std::pow (
        +std::numeric_limits <real_type>
            ::epsilon(),(real_type)+.8);
    
        real_type _hval = 
            +std::numeric_limits
                <real_type>::infinity();
    
        __unreferenced(_hint) ;
    
        if (this->_xpos.count() == +0)
            return _hval ;
            
        if (this->_ypos.count() == +0)
            return _hval ;
    
    /*---------------------------- compute xyz to lat-lon */
        real_type _radius = std::sqrt (
            _ppos[0]*_ppos[0] + 
            _ppos[1]*_ppos[1] + 
            _ppos[2]*_ppos[2] ) ; 
         
        if (_radius < FT)
        _radius = _radius + FT;
         
        real_type _alat = 
        std::asin (_ppos[2]/ _radius ) ;
        
        real_type _alon = 
        std::atan2(_ppos[1], _ppos[0]) ;

        real_type static const PI =
       (real_type)std::atan(+1.0) * 4. ;
       
        real_type static const PI_h = 
       (real_type)+.5 * PI ; 
        real_type static const PI_1 = 
       (real_type)+1. * PI ; 
        real_type static const PI_2 = 
       (real_type)+2. * PI ; 
 
        if (_alat<-PI_h) _alat  =-PI_h ;
        if (_alat>=PI_h) _alat  = PI_h ;
        
        if (_alon<-PI_1) _alon += PI_2 ;
        if (_alon>=PI_1) _alon -= PI_2 ;   
        
        if (_alon < *this->_xpos.head())
            return _hval ;
        if (_alon > *this->_xpos.tail())
            return _hval ;
        
        if (_alat < *this->_ypos.head())
            return _hval ;
        if (_alat > *this->_ypos.tail())
            return _hval ;
        
    /*---------------------------- find enclosing x-range */
        iptr_type _ipos = (iptr_type) -1 ;
        iptr_type _jpos = (iptr_type) -1 ;
           
        if (this->_xvar == true)
        {
            auto _jits = 
            algorithms::upper_bound (
                this->_xpos.head(), 
                this->_xpos.tend(), 
            _alon,std::less<real_type>() ) ;
           
            _jpos = (iptr_type)
           (_jits - this->_xpos.head()- 1) ;
        }
        else
        {
            real_type _xmin, _xmax, _xdel;
            _xmin = *this->_xpos.head();
            _xmax = *this->_xpos.tail();
            
            _xdel = (_xmax - _xmin) /
                (this->_xpos.count () - 1) ;
            
            _jpos = (iptr_type)
                ((_alon-_xmin)/_xdel) ;
        }
        
    /*---------------------------- find enclosing y-range */
        if (this->_yvar == true)
        {
            auto _iits = 
            algorithms::upper_bound (
                this->_ypos.head(), 
                this->_ypos.tend(), 
            _alat,std::less<real_type>() ) ;
           
            _ipos = (iptr_type)
           (_iits - this->_ypos.head()- 1) ;
        }
        else
        {
            real_type _ymin, _ymax, _ydel;
            _ymin = *this->_ypos.head();
            _ymax = *this->_ypos.tail();
            
            _ydel = (_ymax - _ymin) /
                (this->_ypos.count () - 1) ;
            
            _ipos = (iptr_type)
                ((_alat-_ymin)/_ydel) ;
        }
        
        if (_ipos == 
       (iptr_type)this->_ypos.count() - 1)
            _ipos = _ipos - 1 ;
        if (_jpos == 
       (iptr_type)this->_xpos.count() - 1)
            _jpos = _jpos - 1 ;
        
    /*---------------------------- a linear interpolation */
        real_type _xp11 = 
            this->_xpos[_jpos + 0] ;
        real_type _xp22 = 
            this->_xpos[_jpos + 1] ;
            
        real_type _yp11 = 
            this->_ypos[_ipos + 0] ;
        real_type _yp22 = 
            this->_ypos[_ipos + 1] ;

        real_type _xval = _alon ;
        real_type _yval = _alat ;
        
        real_type _aa22 = 
           (_yval-_yp11) * (_xval-_xp11) ;
        real_type _aa21 = 
           (_yval-_yp11) * (_xp22-_xval) ;
        real_type _aa12 = 
           (_yp22-_yval) * (_xval-_xp11) ;
        real_type _aa11 = 
           (_yp22-_yval) * (_xp22-_xval) ;
    
        iptr_type _kk11 ;
        indx_from_subs(
            _ipos + 0, _jpos + 0, _kk11) ;
        iptr_type _kk12 ;
        indx_from_subs(
            _ipos + 0, _jpos + 1, _kk12) ;
        iptr_type _kk21 ;
        indx_from_subs(
            _ipos + 1, _jpos + 0, _kk21) ;
        iptr_type _kk22 ;
        indx_from_subs(
            _ipos + 1, _jpos + 1, _kk22) ;
        
        real_type _hbar = 
          ( _aa11*this->_hmat[_kk11]
          + _aa12*this->_hmat[_kk12]
          + _aa21*this->_hmat[_kk21]
          + _aa22*this->_hmat[_kk22] )
          /(_aa11+_aa12+_aa21+_aa22) ;

        return (  _hbar ) ; 
    }
     
    } ;
     
     
    }

#   endif   //__HFUN_GRID_ELLIPSOID_3__



