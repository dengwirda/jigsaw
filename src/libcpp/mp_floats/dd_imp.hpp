
/*
------------------------------------------------------------
 double-double precision arithmetic
------------------------------------------------------------
 */

class dd_imp
    {
    public  :
    typedef double word ;
    
    public  :        
    __static_call 
    __inline_call void_type one_one_add_fast (
        word  _aa ,
        word  _bb ,
        word& _yy ,
        word& _er
        )
    {
        _yy = _aa + _bb ;
        _er = _yy - _aa ;
        _er = _bb - _er ;
    }
    
    __static_call 
    __inline_call void_type one_one_sub_fast (
        word  _aa ,
        word  _bb ,
        word& _yy ,
        word& _er
        )
    {
        _yy = _aa - _bb ;
        _er = _aa - _yy ;
        _er = _er - _bb ;
    }

    __static_call 
    __inline_call void_type one_one_add_full (
        word  _aa ,
        word  _bb ,
        word& _yy ,
        word& _er
        )
    {
        _yy = _aa + _bb ;
        
        word  _t0 , _t1 ;
        _t0 = _yy - _aa ;
        _t1 = _yy - _t0 ;
        
        word  _t2 , _t3 ;
        _t2 = _aa - _t1 ;
        _t3 = _bb - _t0 ;
        
        _er = _t2 + _t3 ;
    }
    
    __static_call 
    __inline_call void_type one_one_sub_full (
        word  _aa ,
        word  _bb ,
        word& _yy ,
        word& _er
        )
    {
        _yy = _aa - _bb ;
        
        word  _t0 , _t1 ;
        _t0 = _yy - _aa ;  
        _t1 = _yy - _t0 ;
        
        word  _t2 , _t3 ;
        _t2 = _aa - _t1 ;
        _t3 = _bb + _t0 ;
        
        _er = _t2 - _t3 ;
    }
    
    __static_call
    __inline_call void_type one_split (
        word  _aa ,
        word& _a1 ,
        word& _a0
        )
    {  
        if (_aa > +6.69692879491417E+299 ||  
            _aa < -6.69692879491417E+299 )
        {
        word  _tt , _ss ;
        _aa*= +3.7252902984619140625E-09 ;
        _tt = _aa * 134217729.0 ;
        _ss = _tt - _aa ;
        _a1 = _tt - _ss ;
        _a0 = _aa - _a1 ;
        _a1*=       268435456.0 ;
        _a0*=       268435456.0 ;
        }
        else
        {
        word  _tt , _ss ;
        _tt = _aa * 134217729.0 ;
        _ss = _tt - _aa ;
        _a1 = _tt - _ss ;
        _a0 = _aa - _a1 ;
        }
    }
    
    __static_call
    __inline_call void_type one_one_mul_full (
        word  _aa ,
        word  _bb ,
        word& _y1 ,
        word& _y0
        )
    {    
        word _a1, _a0, _b1, _b0 ;
        one_split(_aa, _a1, _a0);
        one_split(_bb, _b1, _b0);
        
        _y1    = _aa * _bb ;
        
        word _e1, _e2, _e3 ;
        _e1 = _y1 - (_a1 * _b1) ;
        _e2 = _e1 - (_a0 * _b1) ;
        _e3 = _e2 - (_a1 * _b0) ;
        
        _y0 = (_a0 * _b0) - _e3 ;
        
    }
    
    __static_call
    __inline_call void_type two_one_add (
        word  _a1 ,
        word  _a0 ,
        word  _bb ,
        word& _y1 ,
        word& _y0
        )
    {
        word _t0, _t1 ;
        one_one_add_full (
            _a1, _bb, 
            _t1, _t0) ;
        
        _t0 += _a0 ;
        one_one_add_fast (
            _t1, _t0, 
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type two_one_sub (
        word  _a1 ,
        word  _a0 ,
        word  _bb ,
        word& _y1 ,
        word& _y0
        )
    {
        word _t0, _t1 ;
        one_one_sub_full (
            _a1, _bb, 
            _t1, _t0) ;
            
        _t0 += _a0 ;
        one_one_add_fast (
            _t1, _t0, 
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type one_two_sub (
        word  _aa ,
        word  _b1 ,
        word  _b0 ,
        word& _y1 ,
        word& _y0
        )
    {
        word _t0, _t1 ;
        one_one_sub_full (
            _aa, _b1, 
            _t1, _t0) ;
            
        _t0 -= _b0 ;
        one_one_add_fast (
            _t1, _t0,
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type two_two_add (
        word  _a1 ,
        word  _a0 ,
        word  _b1 ,
        word  _b0 ,
        word& _y1 ,
        word& _y0
        )
    {
        word _s0, _s1 ;
        word _t0, _t1 ;
        word _w0, _w1 ;
        one_one_add_full (
            _a1, _b1, 
            _s1, _s0) ;
        one_one_add_full (
            _a0, _b0, 
            _t1, _t0) ;
    
        _s0 += _t1;
        one_one_add_fast (
            _s1, _s0, 
            _w1, _w0) ;
            
        _w0 += _t0;
        one_one_add_fast (
            _w1, _w0, 
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type two_two_sub (
        word  _a1 ,
        word  _a0 ,
        word  _b1 ,
        word  _b0 ,
        word& _y1 ,
        word& _y0
        )
    {
        word _s0, _s1 ;
        word _t0, _t1 ;
        word _w0, _w1 ;
        one_one_sub_full (
            _a1, _b1, 
            _s1, _s0) ;
        one_one_sub_full (
            _a0, _b0, 
            _t1, _t0) ;
            
        _s0 += _t1;
        one_one_add_fast (
            _s1, _s0, 
            _w1, _w0) ;
            
        _w0 += _t0;
        one_one_add_fast (
            _w1, _w0, 
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type two_one_mul (
        word  _a1 ,
        word  _a0 ,
        word  _bb ,
        word& _y1 ,
        word& _y0
        )
    {
        word _t0, _t1, _ss ;
        one_one_mul_full (
            _a1, _bb, 
            _t1, _t0) ;
      
        _ss = _a0 * _bb ;
        _t0 = _t0 + _ss ;
        
        one_one_add_fast (
            _t1, _t0, 
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type two_two_mul (
        word  _a1 ,
        word  _a0 ,
        word  _b1 ,
        word  _b0 ,
        word& _y1 ,
        word& _y0
        )
    {  
        word _t0, _t1, _ss ;
        one_one_mul_full (
            _a1, _b1, 
            _t1, _t0) ;
            
        _ss = _a1 * _b0 +
              _a0 * _b1 ;
        _t0 = _t0 + _ss ;
        
        one_one_add_fast (
            _t1, _t0, 
            _y1, _y0) ;
    }
    
    __static_call
    __inline_call void_type two_two_div (
        word  _a1 ,
        word  _a0 ,
        word  _b1 ,
        word  _b0 ,
        word& _y1 ,
        word& _y0
        )
    {
        word  _t0 , _t1 , _ee ;
        _t1 = _a1 / _b1 ;

        word  _r0 , _r1 ;
        word  _w0 , _w1 ;
        two_one_mul(_b1, _b0, // rr = bb * t1 
                    _t1, 
                    _r1, _r0) ;
                    
        two_two_sub(_a1, _a0, // ww = aa - rr 
                    _r1, _r0,
                    _w1, _w0) ;
        
        _t0 = _w1 / _b1 ;
        
        word  _u0 , _u1 ;
        two_one_mul(_b1, _b0, // rr = bb * t0
                    _t0, 
                    _r1, _r0) ;
        two_two_sub(_w1, _w0, // uu = ww - rr 
                    _r1, _r0,
                    _u1, _u0) ;
    
        _ee = _u1 / _b1 ;
      
        word  _q0 , _q1 ;
        one_one_add_fast(_t1, _t0, _q1, _q0);
        
        two_one_add(_q1, _q0, 
                    _ee, 
                    _t1, _t0) ;
        
        _y1 = _t1; _y0 = _t0;          
    }
    
    } ;
    
    
