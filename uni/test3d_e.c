
//  gcc -Wall -Wextra test3d_e.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_e

//  Use JIGSAW to mesh a simple domain in E^3, starting 
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

        jigsaw_VERT3_t _geom_vert3[8] = {  // setup geom.
            { {0., 0., 0.}, +0 } ,
            { {3., 0., 0.}, +0 } ,
            { {3., 3., 0.}, +0 } ,
            { {0., 3., 0.}, +0 } ,
            { {0., 0., 3.}, +0 } ,
            { {3., 0., 3.}, +0 } ,
            { {3., 3., 3.}, +0 } ,
            { {0., 3., 3.}, +0 } ,
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

        jigsaw_VERT3_t _init_vert3[2] = {
            { {0., 1.5, 1.25}, +0 } ,
            { {0., 1.5, 1.75}, +0 } ,
            } ;

        jigsaw_EDGE2_t _init_edge2[1] = {
            { {+0, +1}, -1 } ,      // -1 => "un-refinable"
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _init._vert3._data = &_init_vert3[0] ;
        _init._vert3._size = +2 ;

        _init._edge2._data = &_init_edge2[0] ;
        _init._edge2._size = +1 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_hmax = 1./1.;
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._mesh_rad2 = 1.50 ;
        _jjig._mesh_rad3 = 2.50 ;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;
        
        _jjig._mesh_dims =   +3 ;   // make 3-dim cell

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
            & _init ,               // init. data
               NULL ,               // empty hfun. obj.
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



