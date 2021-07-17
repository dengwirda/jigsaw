
//  gcc -Wall -Wextra test3d_a.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_a

//  A simple example to start: use JIGSAW to mesh a plain
//  domain in E^3.

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

        jigsaw_msh_t _mesh ;
        jigsaw_init_msh_t(&_mesh) ;

    /*-------------------------------- setup JIGSAW geom. */

        jigsaw_VERT3_t _vert3[8] = {    // setup geom.
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

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +2 ;

        _jjig._hfun_hmax = .675 ;   // uniform spacing
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;

        _jjig._mesh_dims =   +3 ;   // make 3-dim cell

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
               NULL ,               // empty hfun. obj.
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



