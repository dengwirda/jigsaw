
//  gcc -Wall -Wextra test2s_b.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2s_b

//  Use JIGSAW to mesh a simple geometry with user-defined
//  mesh-spacing data defined on a "grid" in S^2.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2s_b (int _verb)
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

    /*-------------------------------- setup JIGSAW geom. */

        real_t _radii[3] = {+1., +1., +1.
            } ;

        _geom._flags =
            JIGSAW_ELLIPSOID_MESH ;

        _geom._radii._data = &_radii[0] ;
        _geom._radii._size = +3 ;

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

        real_t PI = 3.14159265358979323846 ;

        real_t         _hfun_xgrid[3] = {  // setup hfun.
            -1.0 * PI, 0.0, +1.0 * PI
            } ;

        real_t         _hfun_ygrid[3] = {
            -0.5 * PI, 0.0, +0.5 * PI
            } ;

        fp32_t         _hfun_value[9] = {
            .5f, 1.f, .5f, .5f, 1.f, .5f, .5f,
            1.f, .5f
            } ;

        _hfun._flags
            = JIGSAW_ELLIPSOID_GRID;

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

        _jjig._mesh_dims =   +2 ;

        _jjig._mesh_kern =          // non-default opts
            JIGSAW_KERN_DELAUNAY;
        _jjig._optm_kern =
            JIGSAW_KERN_CVT_DQDX;
        _jjig._optm_cost =
            JIGSAW_KERN_SKEW_COS;

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
            & _hfun ,               // hfun. data
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        if (_verb > 0 )
        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "[2s_b] JIGSAW returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2s_b(1) ; }
#   endif//__SKIP_MAIN__



