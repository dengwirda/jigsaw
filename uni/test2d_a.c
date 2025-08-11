
//  gcc -Wall -Wextra test2d_a.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_a

//  A simple example to start: use JIGSAW to mesh a plain
//  domain in E^2.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2d_a (int _verb)
    {
        int _retv = 0;

    /*-------------------------------- setup JIGSAW types */
        jigsaw_jig_t _jjig ;
        jigsaw_init_jig_t(&_jjig) ;

        jigsaw_msh_t _geom ;
        jigsaw_init_msh_t(&_geom) ;

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

        jigsaw_VERT2_t _vert2[4] = {    // setup geom.
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 }
            } ;

        jigsaw_EDGE2_t _edge2[4] = {
            { {+0, +1}, +0 } ,
            { {+1, +2}, +0 } ,
            { {+2, +3}, +0 } ,
            { {+3, +0}, +0 }
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _geom._vert2._data = &_vert2[0] ;
        _geom._vert2._size = +4 ;

        _geom._edge2._data = &_edge2[0] ;
        _geom._edge2._size = +4 ;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity = _verb;

        _jjig._hfun_hmax = 0.25 ;
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._mesh_dims =   +2 ;

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
               NULL ,               // empty hfun. obj.
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        if (_verb > 0 )
        output_msh_data_2(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "[2d_a] JIGSAW returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2d_a(1) ; }
#   endif//__SKIP_MAIN__



