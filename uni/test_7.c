
//  gcc -Wall test_7.c
//  -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test_7

//  Use MARCHE to set "gradient-limits" on mesh-spacing data

#   include "../inc/lib_jigsaw.h"

#   include "stdio.h"

    int main (
        int          _argc ,
        char       **_argv
        )
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _hfun ;
        jigsaw_init_msh_t(&_hfun) ;

    /*
    --------------------------------------------------------
     * JIGSAW's "mesh" is a piecewise linear complex:
    --------------------------------------------------------
     *
     *      v:3 o---------------o v:2
     *          | \    t:2    / |
     *          |   \       /   |
     *          |     \ 4 /     |
     *          | t:3   o   t:1 |
     *          |     /   \     |
     *          |   /       \   |
     *          | /    t:0    \ |
     *      v:0 o---------------o v:1
     *
    --------------------------------------------------------
     */

        jigsaw_VERT2_t _hfun_vert2[5] = {   // setup hfun.
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .5}, +0 }
            } ;

        jigsaw_TRIA3_t _hfun_tria3[4] = {
            { {+0, +1, +4}, +0 } ,
            { {+1, +2, +4}, +0 } ,
            { {+2, +3, +4}, +4 } ,
            { {+3, +0, +4}, +0 }
            } ;

        real_t _hfun_value[5] = {
            2., 2., 2., 2., 1.
            } ;

        real_t _hfun_slope[1] = {
            .1
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _hfun._vert2._data = &_hfun_vert2[0] ;
        _hfun._vert2._size = +5 ;

        _hfun._tria3._data = &_hfun_tria3[0] ;
        _hfun._tria3._size = +4 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +5 ;

        _hfun._slope._data = &_hfun_slope[0] ;
        _hfun._slope._size = +1 ;

    /*-------------------------------- build MARCHE hfun. */

        _jjig._verbosity =   +1 ;

        _retv = marche (
            &_jjig ,    // the config. opts
            &_hfun ) ;  // the spacing h(x)

    /*-------------------------------- print MARCHE hfun. */

        printf("\n VALUE: \n\n") ;

        for (indx_t _ipos = +0;
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



