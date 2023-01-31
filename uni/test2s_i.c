
//  gcc -Wall -Wextra test2s_i.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2s_i

//  Use MARCHE to set "gradient-limits" on mesh-spacing data
//
//  ensures: |GRAD(h)| <= slope_limit(x),
//
//  via a "fast-marching" solver for the Hamilton-Jacobi eq.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2s_i (int _verb)
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _hfun ;
        jigsaw_init_msh_t(&_hfun) ;

    /*-------------------------------- setup JIGSAW hfun. */

        real_t PI = 3.14159265358979323846 ;

        real_t _radii[3] = {+1., +1., +1.
            } ;

        jigsaw_VERT2_t _hfun_vert2[6] = {
            { {-1.0 * PI, +0.0 * PI}, +0 } ,
            { {-.33 * PI, +0.0 * PI}, +0 } ,
            { {+.33 * PI, +0.0 * PI}, +0 } ,
            { {+1.0 * PI, +0.0 * PI}, +0 } ,
            { {+0.0 * PI, +0.5 * PI}, +0 } ,
            { {+0.0 * PI, -0.5 * PI}, +0 } ,
            } ;

        jigsaw_TRIA3_t _hfun_tria3[6] = {
            { {+0, +1, +4}, +0 } ,
            { {+1, +2, +4}, +0 } ,
            { {+2, +3, +4}, +0 } ,
            { {+1, +0, +5}, +0 } ,
            { {+2, +1, +5}, +0 } ,
            { {+3, +2, +5}, +0 } ,
            } ;

        fp32_t _hfun_value[6] = {
            2.f, 2.f, 2.f, 2.f, 1.f, 2.f
            } ;

        fp32_t _hfun_slope[1] = {
            .1f
            } ;

        _hfun._flags
            = JIGSAW_ELLIPSOID_MESH;

        _hfun._radii._data = &_radii[0] ;
        _hfun._radii._size = +3 ;

        _hfun._vert2._data = &_hfun_vert2[0] ;
        _hfun._vert2._size = +6 ;

        _hfun._tria3._data = &_hfun_tria3[0] ;
        _hfun._tria3._size = +6 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +6 ;

        _hfun._slope._data = &_hfun_slope[0] ;
        _hfun._slope._size = +1 ;


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
       "[2s_i] MARCHE returned code : %d \n", _retv) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2s_i(1) ; }
#   endif//__SKIP_MAIN__



