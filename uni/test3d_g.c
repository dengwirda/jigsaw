
//  gcc -Wall -Wextra test3d_g.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_g

//  Uses TRIPOD to build an (unconstrained) DT in E^3,
//  comparable to more conventional Delaunay tessellation
//  approaches.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int main ()
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _init ;
        jigsaw_init_msh_t(&_init) ;

        jigsaw_msh_t _tria ;
        jigsaw_init_msh_t(&_tria) ;

    /*-------------------------------- pts to triangulate */

        jigsaw_VERT3_t _init_vert3[9] = {
            { {0., 0., 0.}, +0 } ,
            { {5., 1., 0.}, +0 } ,
            { {5., 5., 1.}, +0 } ,
            { {2., 5., 0.}, +0 } ,
            { {0., 2., 5.}, +0 } ,
            { {5., 1., 5.}, +0 } ,
            { {5., 6., 5.}, +0 } ,
            { {0., 5., 5.}, +0 } ,
            { {1., 0., 3.}, +0 } ,
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _init._vert3._data = &_init_vert3[0] ;
        _init._vert3._size = +9 ;

    /*-------------------------------- build TRIPOD r-DT. */

        _jjig._verbosity =   +2 ;

        _jjig._mesh_dims =   +3 ;

        _retv = tripod (
            &_jjig ,                // the config. opts
            &_init ,                // init. data
              NULL ,                // empty geom. data
            &_tria ) ;

    /*-------------------------------- print TRIPOD r-DT. */

        output_msh_data_3(&_tria);

        jigsaw_free_msh_t(&_tria);

        printf (
       "TRIPOD returned code : %d \n",_retv);


        return _retv ;
    }



