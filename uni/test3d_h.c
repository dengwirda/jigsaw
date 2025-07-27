
//  gcc -Wall -Wextra test3d_h.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_h

//  Use MARCHE to set "gradient-limits" on mesh-spacing data
//
//  ensures: |GRAD(h)| <= slope_limit(x),
//
//  via a "fast-marching" solver for the Hamilton-Jacobi eq.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test3d_h (int _verb)
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _hfun ;
        jigsaw_init_msh_t(&_hfun) ;

    /*-------------------------------- setup JIGSAW hfun. */

        real_t         _hfun_xgrid[3] = {   // setup hfun.
            0., .5, 1.
            } ;

        real_t         _hfun_ygrid[3] = {
            0., .5, 1.
            } ;

        real_t         _hfun_zgrid[3] = {
            0., .5, 1.
            } ;

        fp32_t         _hfun_value[27]= {
            2.f, 2.f, 2.f, 2.f, 1.f, 2.f, 2.f, 2.f, 2.f,
            2.f, 2.f, 2.f, 2.f, 1.f, 2.f, 2.f, 2.f, 2.f,
            2.f, 2.f, 2.f, 2.f, 1.f, 2.f, 2.f, 2.f, 2.f
            } ;

        fp32_t         _hfun_slope[27]= {
            .4f, .4f, .4f, .4f, .1f, .4f, .4f, .4f, .4f,
            .4f, .4f, .4f, .4f, .1f, .4f, .4f, .4f, .4f,
            .4f, .4f, .4f, .4f, .1f, .4f, .4f, .4f, .4f
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_GRID;

        _hfun._xgrid._data = &_hfun_xgrid[0] ;
        _hfun._xgrid._size = +3 ;

        _hfun._ygrid._data = &_hfun_ygrid[0] ;
        _hfun._ygrid._size = +3 ;

        _hfun._zgrid._data = &_hfun_zgrid[0] ;
        _hfun._zgrid._size = +3 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +27;

        _hfun._slope._data = &_hfun_slope[0] ;
        _hfun._slope._size = +27;

    /*-------------------------------- build MARCHE hfun. */

        _jjig._verbosity = _verb;

        _retv = marche (
            & _jjig ,               // the config. opts
            & _hfun ) ;             // the spacing h(x)

    /*-------------------------------- print MARCHE hfun. */

        if (_verb > 0 )
        {
        printf("\n VALUE: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _hfun._value._size ;
                   ++_ipos )
        {
            printf("%1.4f\n",
                _hfun._value._data[_ipos]
                ) ;
        }
        }

        printf (
       "[3d_h] MARCHE returned code : %d \n", _retv) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test3d_h(1) ; }
#   endif//__SKIP_MAIN__



