
//  gcc -Wall -Wextra test3d_c.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_c

//  Use JIGSAW to mesh a simple geometry with user-defined
//  mesh-spacing data defined on a "mesh" in E^3.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int main ()
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

        jigsaw_VERT3_t _geom_vert3[8] = {  // setup geom.
            { {0., 0., 0.}, +0 } ,
            { {1., 0., 0.}, +0 } ,
            { {1., 1., 0.}, +0 } ,
            { {0., 1., 0.}, +0 } ,
            { {0., 0., 1.}, +0 } ,
            { {1., 0., 1.}, +0 } ,
            { {1., 1., 1.}, +0 } ,
            { {0., 1., 1.}, +0 }
            } ;

        jigsaw_TRIA3_t _geom_tria3[12] = {
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

        _geom._vert3._data = &_geom_vert3[0] ;
        _geom._vert3._size = +8 ;

        _geom._tria3._data = &_geom_tria3[0] ;
        _geom._tria3._size = +12;

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

        _hfun._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _hfun._vert3._data = &_hfun_vert3[0] ;
        _hfun._vert3._size = +9 ;

        _hfun._tria4._data = &_hfun_tria4[0] ;
        _hfun._tria4._size = +10;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +9 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_scal =
            JIGSAW_HFUN_ABSOLUTE;

        _jjig._hfun_hmax = 1. ;
        _jjig._hfun_hmin = 0. ;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;

        _jjig._mesh_dims =   +3 ;   // make 3-dim cell

        _jjig._mesh_kern =
            JIGSAW_KERN_DELAUNAY;

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
            & _hfun ,               // hfun. data
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



