
//  gcc -Wall -Wextra test2d_e.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_e

//  Use JIGSAW to mesh a simple domain in E^2, starting
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

    /*
    --------------------------------------------------------
     * JIGSAW's "mesh" is a piecewise linear complex:
    --------------------------------------------------------
     *
     *                 e:2
     *      v:3 o---------------o v:2
     *          |               |
     *          |               |
     *          |               |
     *      e:3 |               | e:1
     *          |               |
     *          |               |
     *          |               |
     *      v:0 o---------------o v:1
     *                 e:0
     *
    --------------------------------------------------------
     */

        jigsaw_VERT2_t _geom_vert2[4] = {  // setup geom.
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 }
            } ;

        jigsaw_EDGE2_t _geom_edge2[4] = {
            { {+0, +1}, +0 } ,
            { {+1, +2}, +0 } ,
            { {+2, +3}, +0 } ,
            { {+3, +0}, +0 }
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _geom._vert2._data = &_geom_vert2[0] ;
        _geom._vert2._size = +4 ;

        _geom._edge2._data = &_geom_edge2[0] ;
        _geom._edge2._size = +4 ;

    /*-------------------------------- form init. config. */

        jigsaw_VERT2_t _init_vert2[4] = {
            { {0., 0.}, +0 } ,
            { {0., .5}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .5}, +0 }
            } ;

        jigsaw_EDGE2_t _init_edge2[2] = {
            { {+0, +1}, -1 } ,      // -1 => "un-refinable"
            { {+1, +2}, -1 }
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _init._vert2._data = &_init_vert2[0] ;
        _init._vert2._size = +4 ;

        _init._edge2._data = &_init_edge2[0] ;
        _init._edge2._size = +2 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_hmax = 1./3.;
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._mesh_dims =   +2 ;
        _jjig._geom_feat =   +1 ;
        _jjig._mesh_top1 =   +1 ;

        _retv = jigsaw (
            &_jjig ,                // the config. opts
            &_geom ,                // geom. data
            &_init ,                // init. data
              NULL ,                // empty hfun. data
            &_mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_2(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



