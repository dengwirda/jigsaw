
//  gcc -Wall -Wextra test2s_e.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2s_e

//  Use JIGSAW to mesh a simple domain in S^2, starting 
//  from user-defined initial-conditions.

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

        jigsaw_msh_t _init ;
        jigsaw_init_msh_t(&_init) ;

        jigsaw_msh_t _mesh ;
        jigsaw_init_msh_t(&_mesh) ;

    /*-------------------------------- setup JIGSAW geom. */

        real_t PI = 3.14159265358979323846 ;

        real_t _radii[3] = {+1., +1., +1.
            } ;

        _geom._flags =
            JIGSAW_ELLIPSOID_MESH ;

        _geom._radii._data = &_radii[0] ;
        _geom._radii._size = +3 ;

    /*-------------------------------- setup JIGSAW init. */

        jigsaw_VERT2_t _init_vert2[8] = {
            { {-1.0 * PI, +0.0 * PI}, +0 } ,
            { {-.33 * PI, +0.0 * PI}, +0 } ,
            { {+.33 * PI, +0.0 * PI}, +0 } ,
            { {+1.0 * PI, +0.0 * PI}, +0 } ,
            { {+0.0 * PI, +0.5 * PI}, +0 } ,
            { {+0.0 * PI, -0.5 * PI}, +0 } ,
            { {+0.0 * PI, +.25 * PI}, +0 } ,
            { {+0.0 * PI, -.25 * PI}, +0 } ,
            } ;

        jigsaw_EDGE2_t _init_edge2[1] = {
            { {+6, +7}, -1 } ,      // -1 => "un-refinable"
            } ;

        _init._flags
            = JIGSAW_ELLIPSOID_MESH ;

        _init._radii._data = &_radii[0] ;
        _init._radii._size = +3 ;

        _init._vert2._data = &_init_vert2[0] ;
        _init._vert2._size = +8 ;

        _init._edge2._data = &_init_edge2[0] ;
        _init._edge2._size = +1 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_hmax = 1./3.;
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._mesh_dims =   +2 ;

        _retv = jigsaw (
            &_jjig ,                // the config. opts
            &_geom ,                // geom. data
            &_init ,                // init. data
              NULL ,                // empty hfun. data
            &_mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



