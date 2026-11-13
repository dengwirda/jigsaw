
//  gcc -Wall -Wextra test2s_g.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2s_g

//  Uses TRIPOD to build an (unconstrained) DT in S^2,
//  comparable to more conventional Delaunay tessellation
//  approaches.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test2s_g (int _verb)
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

        real_t _radii[3] = {+1., +1., +1.
            } ;

        _geom._flags =
            JIGSAW_ELLIPSOID_MESH ;

        _geom._radii._data = &_radii[0] ;
        _geom._radii._size = +3 ;

    /*-------------------------------- pts to triangulate */

        real_t PI = 3.14159265358979323846 ;
        indx_t _next = 0;

        jigsaw_VERT2_t _point[45] ;

        for (size_t _ipos = 0; _ipos < 9; ++_ipos)
        for (size_t _jpos = 0; _jpos < 5; ++_jpos)
        {
            _point[_next]._ppos[0] =
                -0.80 * PI + _ipos * 0.20 * PI ;
            _point[_next]._ppos[1] =
                -0.40 * PI + _jpos * 0.20 * PI ;

            _point[_next]._itag = +0 ;

            _next +=  +1;
        }

        _init._flags
            = JIGSAW_ELLIPSOID_MESH;

        _init._radii._data = &_radii[0] ;
        _init._radii._size = +3 ;

        _init._vert2._data = &_point[0] ;
        _init._vert2._size = +45;

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
       "[2s_g] TRIPOD returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test2s_g(1) ; }
#   endif//__SKIP_MAIN__



