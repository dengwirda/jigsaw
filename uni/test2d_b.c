
//  gcc -Wall -Wextra test2d_b.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_b

//  Use JIGSAW to mesh a simple geometry with user-defined
//  mesh-spacing data defined on a "grid" in E^2.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2d_b (int _verb)
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

        jigsaw_VERT2_t _geom_vert2[4] = {   // setup geom.
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
     * JIGSAW's "grid" uses a column-major numbering:
    --------------------------------------------------------
     *
     *                 v:5
     *      v:2 o-------o-------o v:8
     *          |       |       |
     *          |       |       |
     *          |      v:4      |
     *      v:1 o-------o-------o v:7
     *          |       |       |
     *          |       |       |
     *          |       |       |
     *      v:0 o-------o-------o v:6
     *                 v:3
     *
    --------------------------------------------------------
     */

        real_t         _hfun_xgrid[3] = {  // setup hfun.
            0., .5, 1.
            } ;

        real_t         _hfun_ygrid[3] = {
            0., .5, 1.
            } ;

        fp32_t         _hfun_value[9] = {
            .3f, .2f, .3f, .2f, .1f, .2f, .3f,
            .2f, .3f
            } ;

        _hfun._flags
            = JIGSAW_EUCLIDEAN_GRID;

        _hfun._xgrid._data = &_hfun_xgrid[0] ;
        _hfun._xgrid._size = +3 ;

        _hfun._ygrid._data = &_hfun_ygrid[0] ;
        _hfun._ygrid._size = +3 ;

        _hfun._value._data = &_hfun_value[0] ;
        _hfun._value._size = +9 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity = _verb;

        _jjig._hfun_scal =
            JIGSAW_HFUN_ABSOLUTE;

        _jjig._hfun_hmax = 1. ;
        _jjig._hfun_hmin = 0. ;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;

        _jjig._mesh_dims =   +2 ;

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
            & _hfun ,               // hfun. data
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        if (_verb > 0 )
        output_msh_data_2(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "[2d_b] JIGSAW returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2d_b(1) ; }
#   endif//__SKIP_MAIN__



