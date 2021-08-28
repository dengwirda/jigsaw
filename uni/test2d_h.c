
//  gcc -Wall -Wextra test2d_h.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_h

//  Use MARCHE to set "gradient-limits" on mesh-spacing data
//
//  ensures: |GRAD(h)| <= slope_limit(x),
//
//  via a "fast-marching" solver for the Hamilton-Jacobi eq.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int main ()
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _hfun ;
        jigsaw_init_msh_t(&_hfun) ;

    /*
    --------------------------------------------------------
     * JIGSAW's "grid" uses a column-major numbering:
    --------------------------------------------------------
     *
     *                 v:5
     *      v:2 o-------o-------o v:8
     *          |       |       |
     *          |       |       |
     *          |      v:4      |
     *      v:1 o-------o-------o v:7
     *          |       |       |
     *          |       |       |
     *          |       |       |
     *      v:0 o-------o-------o v:6
     *                 v:3
     *
    --------------------------------------------------------
     */

        real_t         _hfun_xgrid[3] = {   // setup hfun.
            0., .5, 1.
            } ;

        real_t         _hfun_ygrid[3] = {
            0., .5, 1.
            } ;

        fp32_t         _hfun_value[9] = {
            2.f, 2.f, 2.f, 2.f, 1.f, 2.f, 2.f,
            2.f, 2.f
            } ;

        fp32_t         _hfun_slope[9] = {
            .4f, .4f, .4f, .4f, .1f, .4f, .4f,
            .4f, .4f
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_GRID;

        _hfun._xgrid._data = &_hfun_xgrid[0] ;
        _hfun._xgrid._size = +3 ;

        _hfun._ygrid._data = &_hfun_ygrid[0] ;
        _hfun._ygrid._size = +3 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +9 ;

        _hfun._slope._data = &_hfun_slope[0] ;
        _hfun._slope._size = +9 ;

    /*-------------------------------- build MARCHE hfun. */

        _jjig._verbosity =   +1 ;

        _retv = marche (
            &_jjig ,                // the config. opts
            &_hfun ) ;              // the spacing h(x)

    /*-------------------------------- print MARCHE hfun. */

        printf("\n VALUE: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _hfun._value._size ;
                   ++_ipos )
        {
            printf("%1.4f\n",
                _hfun._value._data[_ipos]
                ) ;
        }

        printf (
       "MARCHE returned code : %d \n",_retv);


        return _retv ;
    }



