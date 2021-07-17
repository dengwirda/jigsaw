
//  gcc -Wall -Wextra test2d_d.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test2d_d

//  Use JIGSAW to mesh "multi-part" geometry defined in
//  E^2. Cells inherit part information via ID-tags.

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

    /*
    --------------------------------------------------------
     * A domain with "interior" constraints
    --------------------------------------------------------
     *
     *                 e:2
     *      v:3 o---------------o v:2
     *          |               |
     *          |      e:6      |
     *          |     o---o     |
     *      e:3 | e:7 |   | e:5 | e:1
     *          |     o---o     |
     *          |      e:4      |
     *          |               |
     *      v:0 o---------------o v:1
     *                 e:0
     *
    --------------------------------------------------------
     */

        indx_t _ITAG
            = JIGSAW_EDGE2_TAG ;

        jigsaw_VERT2_t _vert2[8] = {    // setup geom.
            { {0., 0.}, +0 } ,
            { {3., 0.}, +0 } ,
            { {3., 3.}, +0 } ,
            { {0., 3.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {2., 1.}, +0 } ,
            { {2., 2.}, +0 } ,
            { {1., 2.}, +0 }
            } ;

        jigsaw_EDGE2_t _edge2[8] = {
            { {+0, +1}, +0 } ,          // outer geom.
            { {+1, +2}, +0 } ,
            { {+2, +3}, +0 } ,
            { {+3, +0}, +0 } ,
            { {+4, +5}, +0 } ,          // inner geom.
            { {+5, +6}, +0 } ,
            { {+6, +7}, +0 } ,
            { {+7, +4}, +0 }
            } ;

        jigsaw_BOUND_t _bound[12]= {
            {  +1, +0, _ITAG},
            {  +1, +1, _ITAG},
            {  +1, +2, _ITAG},
            {  +1, +3, _ITAG},
            {  +1, +4, _ITAG},
            {  +1, +5, _ITAG},
            {  +1, +6, _ITAG},
            {  +1, +7, _ITAG},
            {  +2, +4, _ITAG},          // inner geom.
            {  +2, +5, _ITAG},
            {  +2, +6, _ITAG},
            {  +2, +7, _ITAG},
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH ;

        _geom._vert2._data = &_vert2[0] ;
        _geom._vert2._size = +8 ;

        _geom._edge2._data = &_edge2[0] ;
        _geom._edge2._size = +8 ;

        _geom._bound._data = &_bound[0] ;
        _geom._bound._size = +12;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity =   +1 ;

        _jjig._hfun_hmax = 0.25 ;
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._mesh_dims =   +2 ;

        _retv = jigsaw (
            &_jjig ,                // the config. opts
            &_geom ,                // geom. data
              NULL ,                // empty init. data
              NULL ,                // empty hfun. data
            &_mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        output_msh_data_2(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "JIGSAW returned code : %d \n",_retv);


        return _retv ;
    }



