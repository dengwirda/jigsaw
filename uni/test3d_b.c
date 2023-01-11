
//  gcc -Wall -Wextra test3d_b.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_b

//  Use JIGSAW to mesh a simple geometry with user-defined
//  mesh-spacing data defined on a "grid" in E^3.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test3d_b (int _verb)
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _geom ;
        jigsaw_init_msh_t(&_geom) ;

        jigsaw_msh_t _hfun ;
        jigsaw_init_msh_t(&_hfun) ;

        jigsaw_msh_t _mesh ;
        jigsaw_init_msh_t(&_mesh) ;

    /*-------------------------------- setup JIGSAW geom. */

        jigsaw_VERT3_t _vert3[8] = {       // setup geom.
            { {0., 0., 0.}, +0 } ,
            { {1., 0., 0.}, +0 } ,
            { {1., 1., 0.}, +0 } ,
            { {0., 1., 0.}, +0 } ,
            { {0., 0., 1.}, +0 } ,
            { {1., 0., 1.}, +0 } ,
            { {1., 1., 1.}, +0 } ,
            { {0., 1., 1.}, +0 }
            } ;

        jigsaw_TRIA3_t _tria3[12] = {
            { {+0, +1, +2}, +0 } ,
            { {+0, +2, +3}, +0 } ,
            { {+4, +5, +6}, +0 } ,
            { {+4, +6, +7}, +0 } ,
            { {+0, +1, +5}, +0 } ,
            { {+0, +5, +4}, +0 } ,
            { {+1, +2, +6}, +0 } ,
            { {+1, +6, +5}, +0 } ,
            { {+2, +3, +7}, +0 } ,
            { {+2, +7, +6}, +0 } ,
            { {+3, +7, +4}, +0 } ,
            { {+3, +4, +0}, +0 }
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _geom._vert3._data = &_vert3[0] ;
        _geom._vert3._size = +8 ;

        _geom._tria3._data = &_tria3[0] ;
        _geom._tria3._size = +12;

    /*-------------------------------- setup JIGSAW hfun. */

        real_t         _hfun_xgrid[2] = {  // setup hfun.
            0., 1.
            } ;

        real_t         _hfun_ygrid[2] = {
            0., 1.
            } ;

        real_t         _hfun_zgrid[2] = {
            0., 1.
            } ;

        fp32_t         _hfun_value[8] = {
            .6f, .6f, .6f, .6f,
            .3f, .3f, .3f, .3f
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_GRID;

        _hfun._xgrid._data = &_hfun_xgrid[0] ;
        _hfun._xgrid._size = +2 ;

        _hfun._ygrid._data = &_hfun_ygrid[0] ;
        _hfun._ygrid._size = +2 ;

        _hfun._zgrid._data = &_hfun_zgrid[0] ;
        _hfun._zgrid._size = +2 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +8 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity = _verb;

        _jjig._hfun_scal =
            JIGSAW_HFUN_ABSOLUTE;

        _jjig._hfun_hmax = 1. ;
        _jjig._hfun_hmin = 0. ;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;

        _jjig._mesh_dims =   +3 ;   // make 3-dim cell

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
            & _hfun ,               // hfun. data
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        if (_verb > 0 )
        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "[3d_b] JIGSAW returned code : %d \n", _retv) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test3d_b(1) ; }
#   endif//__SKIP_MAIN__



