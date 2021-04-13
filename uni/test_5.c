
//  gcc -Wall -Wextra test_5.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test_5

//  An example that uses TRIPOD to build a "restricted" DT.

#   include "../inc/lib_jigsaw.h"

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

        jigsaw_msh_t _tria ;
        jigsaw_init_msh_t(&_tria) ;

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

    /*-------------------------------- pts to triangulate */

        jigsaw_VERT2_t _point[9] = {
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .0}, +0 } ,
            { {1., .5}, +0 } ,
            { {.5, 1.}, +0 } ,
            { {.0, .5}, +0 } ,
            { {.3, .3}, +0 }
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _init._vert2._data = &_point[0] ;
        _init._vert2._size = +9 ;

    /*-------------------------------- build TRIPOD r-DT. */

        _jjig._verbosity =   +1 ;

        _jjig._mesh_dims =   +2 ;

        _retv = tripod (
            &_jjig ,    // the config. opts
            &_init ,    // init. data
            &_geom ,    // geom. data
            &_tria ) ;

    /*-------------------------------- print TRIPOD r-DT. */

        printf("\n VERT2: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _tria._vert2._size ;
                   ++_ipos )
        {
            printf("%1.4f, %1.4f\n",
            _tria._vert2.
                _data[_ipos]._ppos[0],
            _tria._vert2.
                _data[_ipos]._ppos[1]
                ) ;
        }

        printf("\n TRIA3: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _tria._tria3._size ;
                   ++_ipos )
        {
            printf("%d, %d, %d\n",
            _tria._tria3.
                _data[_ipos]._node[0],
            _tria._tria3.
                _data[_ipos]._node[1],
            _tria._tria3.
                _data[_ipos]._node[2]
                ) ;
        }

        jigsaw_free_msh_t(&_tria);

        printf (
       "TRIPOD returned code : %d \n",_retv);


        return _retv ;
    }



