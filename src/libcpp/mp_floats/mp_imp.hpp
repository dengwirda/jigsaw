
/*
------------------------------------------------------------
 multi-double precision arithmetic
------------------------------------------------------------
 */

class mp_imp
    {
    public  :
    typedef double word ;

    public  :

    __static_call
    __inline_call void_type one_one_add_full (
        word  _aa ,
        word  _bb ,
        word& _y0 ,
        word& _y1
        )
    {
        _y0 = _aa + _bb ;

        word  _t0 , _t1 ;
        _t0 = _y0 - _aa ;
        _t1 = _y0 - _t0 ;

        word  _t2 , _t3 ;
        _t2 = _aa - _t1 ;
        _t3 = _bb - _t0 ;

        _y1 = _t2 + _t3 ;
    }

    __static_call
    __inline_call void_type one_one_sub_full (
        word  _aa ,
        word  _bb ,
        word& _y0 ,
        word& _y1
        )
    {
        _y0 = _aa - _bb ;

        word  _t0 , _t1 ;
        _t0 = _y0 - _aa ;
        _t1 = _y0 - _t0 ;

        word  _t2 , _t3 ;
        _t2 = _aa - _t1 ;
        _t3 = _bb + _t0 ;

        _y1 = _t2 - _t3 ;
    }

    __static_call
    __inline_call void_type one_split (
        word  _aa ,
        word& _a0 ,
        word& _a1
        )
    {
        if (_aa > +6.69692879491417E+299 ||
            _aa < -6.69692879491417E+299 )
        {
            word  _tt , _ss ;
            _aa*=
             + 3.7252902984619140625E-09 ;
            _tt = _aa * 134217729.0 ;
            _ss = _tt - _aa ;
            _a0 = _tt - _ss ;
            _a1 = _aa - _a0 ;
            _a0*=       268435456.0 ;
            _a1*=       268435456.0 ;
        }
        else
        {
            word  _tt , _ss ;
            _tt = _aa * 134217729.0 ;
            _ss = _tt - _aa ;
            _a0 = _tt - _ss ;
            _a1 = _aa - _a0 ;
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
        one_split(_aa, _a0, _a1);
        one_split(_bb, _b0, _b1);

        _y0    = _aa * _bb ;

        word _e1, _e2, _e3 ;
        _e1 = _y0 - (_a0 * _b0) ;
        _e2 = _e1 - (_a1 * _b0) ;
        _e3 = _e2 - (_a0 * _b1) ;

        _y1 = (_a1 * _b1) - _e3 ;
    }

    } ;


    template <
    size_t   N
             >
    class mp_flt
    {
    public  :

    double              _dx[N] = {0.0} ;
    size_t              _dn    =   +0  ;

    } ;


    template <
    size_t   X
             >
    __inline_call double    mp_est (
        mp_flt<X> const& _xx
        )
    {
        double _tt = +0.0 ;
        size_t _ii ;

        for (_ii = _xx._dn; _ii-- != 0; )
        {
            _tt += _xx._dx [_ii];
        }
        return _tt ;
    }

    template <
    size_t   X
             >
    __inline_call void_type mp_neg (
        mp_flt<X>      & _xx
        )
    {
        size_t _ii ;

        for (_ii = _xx._dn; _ii-- != 0; )
        {
            _xx._dx[_ii] *= -1. ;
        }
    }

    template <
    size_t   X
             >
    __inline_call void_type mp_add (
        double           _aa,
        double           _bb,
        mp_flt<X>      & _xx
        )
    {
        size_t _nn =  +2 ;
        _nn =    std::min (_nn, X) ;

        double _tt ;

        __assert( _nn <= X &&
        "MP-FLT: N out of range" ) ;

        _xx._dn = _nn ;

        if (_nn == 2)
        mp_imp::one_one_add_full (
            _aa, _bb,
        _xx._dx [ 0],
        _xx._dx [ 1]) ;
        else
        if (_nn == 1)
        mp_imp::one_one_add_full (
            _aa, _bb,
        _xx._dx [ 0],
            _tt     ) ;
    }

    template <
    size_t   A ,
    size_t   X
             >
    __inline_call void_type mp_add (
        mp_flt<A> const& _aa,
        double           _bb,
        mp_flt<X>      & _xx
        )
    {
        mp_add (_aa, _bb, +0, _xx) ;
    }

    template <
    size_t   A ,
    size_t   X
             >
    __inline_call void_type mp_add (
        mp_flt<A> const& _aa,
        double           _bb,
        size_t           _ib,
        mp_flt<X>      & _xx
        )
    {
        size_t _nn = _aa._dn+1 ;
        _nn = std::min(_nn, X) ;

        __assert( _nn <= X &&
        "MP-FLT: N out of range" ) ;

        _xx._dn =     _nn ;

        double _tt = +0.0 ;
        size_t _ii ;

        for (_ii = _ib; _ii < _xx._dn-1; ++_ii)
        {
            mp_imp::one_one_add_full (
                _aa._dx[_ii], _bb,
                _xx._dx[_ii], _tt) ;

            _bb = _tt ;
        }
        _xx._dx[_xx._dn - 1]= _tt  ;
    }

    template <
    size_t   A ,
    size_t   B ,
    size_t   X
             >
    __inline_call void_type mp_add (
        mp_flt<A> const& _aa,
        mp_flt<B> const& _bb,
        mp_flt<X>      & _xx
        )
    {
        mp_add(_aa, _bb._dx[ +0],  +0, _xx) ;

        size_t _ii;

        for (_ii = 1; _ii < _bb._dn; ++_ii)
        {
        mp_add(_xx, _bb._dx[_ii], _ii, _xx) ;
        }
    }

    template <
    size_t   X
             >
    __inline_call void_type mp_sub (
        double           _aa,
        double           _bb,
        mp_flt<X>      & _xx
        )
    {
        size_t _nn =  +2 ;
        _nn =    std::min (_nn, X) ;

        double _tt ;

        __assert( _nn <= X &&
        "MP-FLT: N out of range" ) ;

        _xx._dn = _nn ;

        if (_nn == 2)
        mp_imp::one_one_sub_full (
            _aa, _bb,
        _xx._dx [ 0],
        _xx._dx [ 1]) ;
        else
        if (_nn == 1)
        mp_imp::one_one_sub_full (
            _aa, _bb,
        _xx._dx [ 0],
            _tt     ) ;
    }

    template <
    size_t   A ,
    size_t   X
             >
    __inline_call void_type mp_sub (
        mp_flt<A> const& _aa,
        double           _bb,
        mp_flt<X>      & _xx
        )
    {
        mp_sub (_aa, _bb, +0, _xx) ;
    }

    template <
    size_t   A ,
    size_t   X
             >
    __inline_call void_type mp_sub (
        mp_flt<A> const& _aa,
        double           _bb,
        size_t           _ib,
        mp_flt<X>      & _xx
        )
    {
        size_t _nn = _aa._dn+1 ;
        _nn = std::min(_nn, X) ;

        __assert( _nn <= X &&
        "MP-FLT: N out of range" ) ;

        _xx._dn =     _nn ;

        double _tt = +0.0 ;
        size_t _ii ;

        for (_ii = _ib; _ii < _xx._dn-1; ++_ii)
        {
            mp_imp::one_one_sub_full (
                _aa._dx[_ii], _bb,
                _xx._dx[_ii], _tt) ;

            _bb = _tt ;
        }
        _xx._dx[_xx._dn - 1]= _tt  ;
    }

    template <
    size_t   A ,
    size_t   B ,
    size_t   X
             >
    __inline_call void_type mp_sub (
        mp_flt<A> const& _aa,
        mp_flt<B> const& _bb,
        mp_flt<X>      & _xx
        )
    {
        mp_sub(_aa, _bb._dx[ +0],  +0, _xx) ;

        size_t _ii;

        for (_ii = 1; _ii < _bb._dn; ++_ii)
        {
        mp_sub(_xx, _bb._dx[_ii], _ii, _xx) ;
        }
    }


