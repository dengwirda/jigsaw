
//  gcc -Wall -Wextra test2d_f.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_f

//  Uses TRIPOD to build a "restricted" DT, given a set
//  points and a geometry in E^2. Compared to JIGSAW,
//  TRIPOD computes a "restricted" triangulation (given
//  a set of points a priori), rather than generating a
//  mesh dynamically.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2d_f (int _verb)
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _geom ;
        jigsaw_init_msh_t(&_geom) ;

        jigsaw_msh_t _init ;
        jigsaw_init_msh_t(&_init) ;

        jigsaw_msh_t _tria ;
        jigsaw_init_msh_t(&_tria) ;

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

        jigsaw_VERT2_t _geom_vert2[4] = {   // setup geom.
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 }
            } ;

        jigsaw_EDGE2_t _geom_edge2[4] = {
            { {+0, +1}, +1 } ,
            { {+1, +2}, +2 } ,
            { {+2, +3}, +3 } ,
            { {+3, +0}, +4 }
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _geom._vert2._data = &_geom_vert2[0] ;
        _geom._vert2._size = +4 ;

        _geom._edge2._data = &_geom_edge2[0] ;
        _geom._edge2._size = +4 ;

    /*-------------------------------- pts to triangulate */

        jigsaw_VERT2_t _init_vert2[9] = {
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .0}, +0 } ,
            { {1., .5}, +0 } ,
            { {.5, 1.}, +0 } ,
            { {.0, .5}, +0 } ,
            { {.3, .3}, +0 }
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _init._vert2._data = &_init_vert2[0] ;
        _init._vert2._size = +9 ;

    /*-------------------------------- build TRIPOD r-DT. */

        _jjig._verbosity = _verb;

        _jjig._mesh_dims =   +2 ;

        _retv = tripod (
            & _jjig ,               // the config. opts
            & _init ,               // init. data
            & _geom ,               // geom. data
            & _tria ) ;

    /*-------------------------------- print TRIPOD r-DT. */

        if (_verb > 0 )
        output_msh_data_2(&_tria);

        jigsaw_free_msh_t(&_tria);

        printf (
       "[2d_f] TRIPOD returned code : %d \n", _retv) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2d_f(1) ; }
#   endif//__SKIP_MAIN__



