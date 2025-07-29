
//  gcc -Wall -Wextra test3d_i.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_i

//  Use MARCHE to set "gradient-limits" on mesh-spacing data
//
//  ensures: |GRAD(h)| <= slope_limit(x),
//
//  via a "fast-marching" solver for the Hamilton-Jacobi eq.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test3d_i (int _verb)
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _hfun ;
        jigsaw_init_msh_t(&_hfun) ;

    /*-------------------------------- setup JIGSAW hfun. */

        jigsaw_VERT3_t _hfun_vert3[9] = {  // setup hfun.
            { {.5, .5, 0.}, +0 } ,
            { {0., 0., 0.}, +0 } ,
            { {1., 0., 0.}, +0 } ,
            { {1., 1., 0.}, +0 } ,
            { {0., 1., 0.}, +0 } ,
            { {0., 0., 1.}, +0 } ,
            { {1., 0., 1.}, +0 } ,
            { {1., 1., 1.}, +0 } ,
            { {0., 1., 1.}, +0 } ,
            } ;

        jigsaw_TRIA4_t _hfun_tria4[10] = {
            { {0, 5, 2, 1}, +0 } ,
            { {0, 6, 5, 2}, +0 } ,
            { {7, 0, 6, 2}, +0 } ,
            { {7, 0, 3, 2}, +0 } ,
            { {8, 0, 5, 1}, +0 } ,
            { {8, 0, 4, 1}, +0 } ,
            { {8, 0, 6, 5}, +0 } ,
            { {8, 7, 0, 6}, +0 } ,
            { {8, 0, 4, 3}, +0 } ,
            { {8, 7, 0, 3}, +0 } ,
            } ;

        fp32_t         _hfun_value[9]= {
            .5f, .5f, .5f, .5f, .5f,
            1.f, 1.f, 1.f, 1.f
            } ;

        fp32_t _hfun_slope[1] = {
            .1f
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _hfun._vert3._data = &_hfun_vert3[0] ;
        _hfun._vert3._size = +9 ;

        _hfun._tria4._data = &_hfun_tria4[0] ;
        _hfun._tria4._size = +10;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +9 ;

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
       "[3d_i] MARCHE returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test3d_i(1) ; }
#   endif//__SKIP_MAIN__



