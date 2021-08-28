
//  gcc -Wall -Wextra test2d_c.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_c

//  Use JIGSAW to mesh a simple geometry with user-defined
//  mesh-spacing data defined on a "mesh" in E^2.

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

    /*
    --------------------------------------------------------
     * JIGSAW's "mesh" is a piecewise linear complex:
    --------------------------------------------------------
     *
     *      v:3 o---------------o v:2
     *          | \    t:2    / |
     *          |   \       /   |
     *          |     \ 4 /     |
     *          | t:3   o   t:1 |
     *          |     /   \     |
     *          |   /       \   |
     *          | /    t:0    \ |
     *      v:0 o---------------o v:1
     *
    --------------------------------------------------------
     */

        jigsaw_VERT2_t _hfun_vert2[5] = {  // setup hfun.
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .5}, +0 }
            } ;

        jigsaw_TRIA3_t _hfun_tria3[4] = {
            { {+0, +1, +4}, +0 } ,
            { {+1, +2, +4}, +0 } ,
            { {+2, +3, +4}, +0 } ,
            { {+3, +0, +4}, +0 }
            } ;

        fp32_t _hfun_value[5] = {
            .2f, .2f, .2f, .2f, .5f
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _hfun._vert2._data = &_hfun_vert2[0] ;
        _hfun._vert2._size = +5 ;

        _hfun._tria3._data = &_hfun_tria3[0] ;
        _hfun._tria3._size = +4 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +5 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_scal =
            JIGSAW_HFUN_ABSOLUTE;

        _jjig._hfun_hmax = 1. ;
        _jjig._hfun_hmin = 0. ;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;

        _jjig._mesh_dims =   +2 ;

        _jjig._mesh_kern =
            JIGSAW_KERN_DELAUNAY;

        _retv = jigsaw (
            & _jjig ,               // the config. opts
            & _geom ,               // geom. data
               NULL ,               // empty init. data
            & _hfun ,               // hfun. data
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_2(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



