
//  gcc -Wall -Wextra test3d_f.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_f

//  Uses TRIPOD to build a "restricted" DT, given a set
//  points and a geometry in E^3. Compared to JIGSAW,
//  TRIPOD computes a "restricted" triangulation (given
//  a set of points a priori), rather than generating a
//  mesh dynamically.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test3d_f (int _verb)
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

        jigsaw_VERT3_t _geom_vert3[8] = {  // setup geom.
            { {0., 0., 0.}, +0 } ,
            { {5., 0., 0.}, +0 } ,
            { {5., 5., 0.}, +0 } ,
            { {0., 5., 0.}, +0 } ,
            { {0., 0., 5.}, +0 } ,
            { {5., 0., 5.}, +0 } ,
            { {5., 5., 5.}, +0 } ,
            { {0., 5., 5.}, +0 } ,
            } ;

        jigsaw_TRIA3_t _geom_tria3[12] = {
            { { 0,  1,  2}, +0 } ,
            { { 0,  2,  3}, +0 } ,
            { { 4,  5,  6}, +0 } ,
            { { 4,  6,  7}, +0 } ,
            { { 0,  1,  5}, +0 } ,
            { { 0,  5,  4}, +0 } ,
            { { 1,  2,  6}, +0 } ,
            { { 1,  6,  5}, +0 } ,
            { { 2,  3,  7}, +0 } ,
            { { 2,  7,  6}, +0 } ,
            { { 3,  7,  4}, +0 } ,
            { { 3,  4,  0}, +0 } ,
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _geom._vert3._data = &_geom_vert3[0] ;
        _geom._vert3._size = +8 ;

        _geom._tria3._data = &_geom_tria3[0] ;
        _geom._tria3._size = +12;

    /*-------------------------------- form init. config. */

        jigsaw_VERT3_t _init_vert3[9] = {
            { {0., 0., 0.}, +0 } ,
            { {5., 0., 0.}, +0 } ,
            { {5., 5., 0.}, +0 } ,
            { {0., 5., 0.}, +0 } ,
            { {0., 0., 5.}, +0 } ,
            { {5., 0., 5.}, +0 } ,
            { {5., 5., 5.}, +0 } ,
            { {0., 5., 5.}, +0 } ,
            { {0., 0., 3.}, +0 } ,
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _init._vert3._data = &_init_vert3[0] ;
        _init._vert3._size = +9 ;

    /*-------------------------------- build TRIPOD r-DT. */

        _jjig._verbosity = _verb;

        _jjig._mesh_dims =   +3 ;

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
       "[3d_f] TRIPOD returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test3d_f(1) ; }
#   endif//__SKIP_MAIN__



