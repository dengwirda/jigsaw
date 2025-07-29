
//  gcc -Wall -Wextra test2s_f.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2s_f

//  Uses TRIPOD to build a "restricted" DT, given a set
//  points and a geometry in S^2. Compared to JIGSAW,
//  TRIPOD computes a "restricted" triangulation (given
//  a set of points a priori), rather than generating a
//  mesh dynamically.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2s_f (int _verb)
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

    /*-------------------------------- setup JIGSAW geom. */

        real_t PI = 3.14159265358979323846 ;

        jigsaw_VERT2_t _point[6] = {
            { {0., -0.5 * PI}, +0 } ,
            { {0., +0.5 * PI}, +0 } ,
            { {-1.0 * PI, 0.}, +0 } ,
            { {-0.5 * PI, 0.}, +0 } ,
            { {+0.0 * PI, 0.}, +0 } ,
            { {+0.5 * PI, 0.}, +0 } ,
            } ;

        jigsaw_EDGE2_t _edges[4] = {
            { {2, 3}, +1} ,
            { {3, 4}, +2} ,
            { {4, 5}, +3} ,
            { {5, 2}, +4} ,
            } ;

        real_t _radii[3] = {+1., +1., +1.
            } ;

        _geom._flags =
            JIGSAW_ELLIPSOID_MESH ;

        _geom._radii._data = &_radii[0] ;
        _geom._radii._size = +3 ;

        _geom._vert2._data = &_point[0] ;
        _geom._vert2._size = +6 ;

        _geom._edge2._data = &_edges[0] ;
        _geom._edge2._size = +4 ;

    /*-------------------------------- form init. config. */

        _init._flags
            = JIGSAW_ELLIPSOID_MESH;

        _init._radii._data = &_radii[0] ;
        _init._radii._size = +3 ;

        _init._vert2._data = &_point[0] ;
        _init._vert2._size = +6 ;

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
        output_msh_data_3(&_tria);

        jigsaw_free_msh_t(&_tria);

        printf (
       "[2s_f] TRIPOD returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2s_f(1) ; }
#   endif//__SKIP_MAIN__



