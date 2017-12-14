
/*
------------------------------------------------------------
 * (double-double)^1 precision data-type
------------------------------------------------------------
 */

#pragma once

#ifndef __DD_FLT__
#define __DD_FLT__

#include "dd_imp.hpp"

class dd_flt ;

__inline_call dd_flt operator + (
    dd_flt const&, 
    double       ) ;
__inline_call dd_flt operator + (
    double       ,
    dd_flt const&) ;
__inline_call dd_flt operator + (
    dd_flt const&, 
    dd_flt const&) ;
   
__inline_call dd_flt operator - (
    dd_flt const&, 
    double       ) ;
__inline_call dd_flt operator - (
    double       ,
    dd_flt const&) ;
__inline_call dd_flt operator - (
    dd_flt const&, 
    dd_flt const&) ;
    
__inline_call dd_flt operator * (
    dd_flt const&, 
    double       ) ;
__inline_call dd_flt operator * (
    double       ,
    dd_flt const&) ;
__inline_call dd_flt operator * (
    dd_flt const&, 
    dd_flt const&) ;
    
__inline_call dd_flt operator / (
    dd_flt const&, 
    double       ) ;
__inline_call dd_flt operator / (
    double       ,
    dd_flt const&) ;
__inline_call dd_flt operator / (
    dd_flt const&, 
    dd_flt const&) ;

class dd_flt
    {
    private :
    double           _xx[2] ;
    
    public  :
    __inline_call double&       hi (
        )
    {   return this->_xx[0] ;
    }
    __inline_call double&       lo (
        )
    {   return this->_xx[1] ;
    }    
    __inline_call double const& hi (
        ) const
    {   return this->_xx[0] ;
    }
    __inline_call double const& lo (
        ) const
    {   return this->_xx[1] ;
    } 
    
    public  :
    __inline_call dd_flt (
        double _hi = double(+0.0) ,
        double _lo = double(+0.0)
        )
    {   this->_xx[0] = _hi ;
        this->_xx[1] = _lo ;
    }
    
    __inline_call dd_flt (
        double const* _dx
        )
    {   this->_xx[0] = _dx[0] ;
        this->_xx[1] = _dx[1] ;
    }
    
    __inline_call dd_flt (
        int    _ii
        )
    {   this->_xx[0] =(double)_ii ;
        this->_xx[1] = +0. ;
    }
    
    
    __inline_call dd_flt (
        dd_flt const& _aa
        )
    {
        this->_xx[0] = _aa.hi();
        this->_xx[1] = _aa.lo();
    }
    
    
    __inline_call dd_flt& operator = (
        double        _aa
        )
    {   
        this->_xx[0] = _aa ;
        this->_xx[1] = +0. ;
        
        return ( *this ) ;
    }
    __inline_call dd_flt& operator = (
        dd_flt const& _aa
        )
    {
        this->_xx[0] = _aa.hi();
        this->_xx[1] = _aa.lo();
        
        return ( *this ) ;
    } 
    
    
    __inline_call operator double (
        )
    {   return ( (double)hi() ) ;
    }
    __inline_call operator int    (
        )
    {   return ( (int   )hi() ) ;
    }
    
    
    __inline_call dd_flt operator + (
        ) const
    {   return dd_flt(+hi(), +lo()) ;
    }
    __inline_call dd_flt operator - (
        ) const
    {   return dd_flt(-hi(), -lo()) ;
    }
    
    
    __inline_call dd_flt& operator += (
        double _aa
        )
    {   
        dd_flt _tt = *this + _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    __inline_call dd_flt& operator -= (
        double _aa
        )
    {   
        dd_flt _tt = *this - _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    __inline_call dd_flt& operator *= (
        double _aa
        )
    {   
        dd_flt _tt = *this * _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    __inline_call dd_flt& operator /= (
        double _aa
        )
    {   
        dd_flt _tt = *this / _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    
    __inline_call dd_flt& operator += (
        dd_flt const& _aa
        )
    {   
        dd_flt _tt = *this + _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    __inline_call dd_flt& operator -= (
        dd_flt const& _aa
        )
    {   
        dd_flt _tt = *this - _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    __inline_call dd_flt& operator *= (
        dd_flt const& _aa
        )
    {   
        dd_flt _tt = *this * _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    }
    
    __inline_call dd_flt& operator /= (
        dd_flt const& _aa
        )
    {   
        dd_flt _tt = *this / _aa ;

        hi() = _tt.hi();
        lo() = _tt.lo();
    
        return ( *this ) ;
    } 
    
    } ;
    
// a+b
    __inline_call dd_flt operator + (
        dd_flt const& _aa, 
        double        _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_one_add (
            _aa.hi(), 
            _aa.lo(),
            _bb, _y1, _y0 ) ;
                
        return dd_flt(_y1, _y0) ;
    }
    
    __inline_call dd_flt operator + (
        double        _aa,
        dd_flt const& _bb
        )
    {   return ( _bb + _aa ) ;
    }
    
    __inline_call dd_flt operator + (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_two_add (
            _aa.hi(), _aa.lo(),
            _bb.hi(), _bb.lo(),
            _y1, _y0) ;
                
        return dd_flt(_y1, _y0) ;
    }
    
// a-b
    __inline_call dd_flt operator - (
        dd_flt const& _aa, 
        double        _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_one_sub (
            _aa.hi(), 
            _aa.lo(),
            _bb, _y1, _y0 ) ;
                
        return dd_flt(_y1, _y0) ;
    }
    
    __inline_call dd_flt operator - (
        double        _aa,
        dd_flt const& _bb
        )
    {
        double _y0, _y1;
        dd_imp::one_two_sub (_aa,
            _bb.hi(), 
            _bb.lo(),
            _y1, _y0 ) ;
                
        return dd_flt(_y1, _y0) ;
    }
 
    __inline_call dd_flt operator - (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_two_sub (
            _aa.hi(), _aa.lo(),
            _bb.hi(), _bb.lo(),
            _y1, _y0) ;
                
        return dd_flt(_y1, _y0) ;
    }
  
// a * b
    __inline_call dd_flt operator * (
        dd_flt const& _aa, 
        double        _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_one_mul (
            _aa.hi(), 
            _aa.lo(),
            _bb, _y1, _y0 ) ;
                
        return dd_flt(_y1, _y0) ;
    }
    
    __inline_call dd_flt operator * (
        double        _aa,
        dd_flt const& _bb
        )
    {   return ( _bb * _aa ) ;
    }
    
    __inline_call dd_flt operator * (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_two_mul (
            _aa.hi(), _aa.lo(),
            _bb.hi(), _bb.lo(),
            _y1, _y0) ;
                
        return dd_flt(_y1, _y0) ;
    }
    
// a / b
    __inline_call dd_flt operator / (
        dd_flt const& _aa,
        double        _bb
        )
    {   return ( _aa / dd_flt(_bb) ) ;
    }
    
    __inline_call dd_flt operator / (
        double        _aa,
        dd_flt const& _bb
        )
    {   return ( dd_flt(_aa) / _bb ) ;
    }

    __inline_call dd_flt operator / (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {
        double _y0, _y1;
        dd_imp::two_two_div (
            _aa.hi(), _aa.lo(),
            _bb.hi(), _bb.lo(),
            _y1, _y0) ;
                
        return dd_flt(_y1, _y0) ;
    }

  
    __inline_call  bool operator == (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {   return _aa.hi() == _bb.hi() &&
               _aa.lo() == _bb.lo() ;
    }
    
    __inline_call  bool operator != (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {   return _aa.hi() != _bb.hi() ||
               _aa.lo() != _bb.lo() ;
    }
    
    __inline_call  bool operator <  (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {   return _aa.hi() != _bb.hi() ?
               _aa.hi() <  _bb.hi() :
               _aa.lo() <  _bb.lo() ;
    }
    
    __inline_call  bool operator >  (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {   return _aa.hi() != _bb.hi() ?
               _aa.hi() >  _bb.hi() :
               _aa.lo() >  _bb.lo() ;
    }
    
    __inline_call  bool operator <= (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {   return _aa.hi() == _bb.hi() ?
               _aa.lo() <= _bb.lo() :
               _aa.hi() <= _bb.hi() ;
    }
    
    __inline_call  bool operator >= (
        dd_flt const& _aa, 
        dd_flt const& _bb
        )
    {   return _aa.hi() == _bb.hi() ?
               _aa.lo() >= _bb.lo() :
               _aa.hi() >= _bb.hi() ;
    }
  
#endif//__DD_FLT__


