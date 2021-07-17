
//  gcc -Wall -Wextra test2s_c.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2s_c

//  Use JIGSAW to mesh a simple geometry with user-defined
//  mesh-spacing data defined on a "mesh" in S^2.

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

        real_t _radii[3] = {+1., +1., +1.
            } ;

        _geom._flags =
            JIGSAW_ELLIPSOID_MESH ;

        _geom._radii._data = &_radii[0] ;
        _geom._radii._size = +3 ;

    /*-------------------------------- setup JIGSAW hfun. */

        real_t PI = 3.14159265358979323846 ; 

        jigsaw_VERT2_t _hfun_vert2[9] = {  // setup hfun.
            { {-1. * PI, -0.5 * PI}, +0 } ,
            { {+0. * PI, -0.5 * PI}, +0 } ,
            { {+1. * PI, -0.5 * PI}, +0 } ,
            { {-1. * PI, +0.0 * PI}, +0 } ,
            { {+0. * PI, +0.0 * PI}, +0 } ,
            { {+1. * PI, +0.0 * PI}, +0 } ,
            { {-1. * PI, +0.5 * PI}, +0 } ,
            { {+0. * PI, +0.5 * PI}, +0 } ,
            { {+1. * PI, +0.5 * PI}, +0 }
            } ;

        jigsaw_TRIA3_t _hfun_tria3[8] = {
            { {+0, +1, +3}, +0 } ,
            { {+3, +1, +4}, +0 } ,
            { {+1, +2, +4}, +0 } ,
            { {+4, +2, +5}, +0 } ,
            { {+3, +4, +6}, +0 } ,
            { {+6, +4, +7}, +0 } ,
            { {+4, +5, +7}, +0 } ,
            { {+7, +5, +8}, +0 }
            } ;

        fp32_t         _hfun_value[9] = {
            1., 1., 1., 1., 1., 1., .5,
            .5, .5
            } ;

        _hfun._flags
            = JIGSAW_ELLIPSOID_MESH;

        _hfun._vert2._data = &_hfun_vert2[0] ;
        _hfun._vert2._size = +9 ;

        _hfun._tria3._data = &_hfun_tria3[0] ;
        _hfun._tria3._size = +8 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +9 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_scal =
            JIGSAW_HFUN_ABSOLUTE;

        _jjig._hfun_hmax = 1. ;
        _jjig._hfun_hmin = 0. ;

        _jjig._mesh_dims =   +2 ;

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



