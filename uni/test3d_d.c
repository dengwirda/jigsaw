
//  gcc -Wall -Wextra test3d_d.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test3d_d

//  Use JIGSAW to mesh "multi-part" geometry defined in
//  E^3. Cells inherit part information via ID-tags.

#   include "../inc/lib_jigsaw.h"

#   include "print.h"
#   include "stdio.h"

    int test3d_d (int _verb)
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

        indx_t _ITAG
            = JIGSAW_TRIA3_TAG ;

        jigsaw_VERT3_t _geom_vert3[16] = { // setup geom.
            { {0., 0., 0.}, +0 } ,      // outer cube
            { {3., 0., 0.}, +0 } ,
            { {3., 3., 0.}, +0 } ,
            { {0., 3., 0.}, +0 } ,
            { {0., 0., 3.}, +0 } ,
            { {3., 0., 3.}, +0 } ,
            { {3., 3., 3.}, +0 } ,
            { {0., 3., 3.}, +0 } ,
            { {1., 1., 1.}, +0 } ,      // inner cube
            { {2., 1., 1.}, +0 } ,
            { {2., 2., 1.}, +0 } ,
            { {1., 2., 1.}, +0 } ,
            { {1., 1., 2.}, +0 } ,
            { {2., 1., 2.}, +0 } ,
            { {2., 2., 2.}, +0 } ,
            { {1., 2., 2.}, +0 } ,
            } ;

        jigsaw_TRIA3_t _geom_tria3[24] = {
            { { 0,  1,  2}, +0 } ,      // outer cube
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
            { { 8,  9, 10}, +0 } ,      // inner cube
            { { 8, 10, 11}, +0 } ,
            { {12, 13, 14}, +0 } ,
            { {12, 14, 15}, +0 } ,
            { { 8,  9, 13}, +0 } ,
            { { 8, 13, 12}, +0 } ,
            { { 9, 10, 14}, +0 } ,
            { { 9, 14, 13}, +0 } ,
            { {10, 11, 15}, +0 } ,
            { {10, 15, 14}, +0 } ,
            { {11, 15, 12}, +0 } ,
            { {11, 12,  8}, +0 } ,
            } ;

        jigsaw_BOUND_t _geom_bound[36]= {
            {  +1,  0, _ITAG},
            {  +1,  1, _ITAG},
            {  +1,  2, _ITAG},
            {  +1,  3, _ITAG},
            {  +1,  4, _ITAG},
            {  +1,  5, _ITAG},
            {  +1,  6, _ITAG},
            {  +1,  7, _ITAG},
            {  +1,  8, _ITAG},
            {  +1,  9, _ITAG},
            {  +1, 10, _ITAG},
            {  +1, 11, _ITAG},
            {  +1, 12, _ITAG},
            {  +1, 13, _ITAG},
            {  +1, 14, _ITAG},
            {  +1, 15, _ITAG},
            {  +1, 16, _ITAG},
            {  +1, 17, _ITAG},
            {  +1, 18, _ITAG},
            {  +1, 19, _ITAG},
            {  +1, 20, _ITAG},
            {  +1, 21, _ITAG},
            {  +1, 22, _ITAG},
            {  +1, 23, _ITAG},
            {  +2, 12, _ITAG},          // inner geom.
            {  +2, 13, _ITAG},
            {  +2, 14, _ITAG},
            {  +2, 15, _ITAG},
            {  +2, 16, _ITAG},
            {  +2, 17, _ITAG},
            {  +2, 18, _ITAG},
            {  +2, 19, _ITAG},
            {  +2, 20, _ITAG},
            {  +2, 21, _ITAG},
            {  +2, 22, _ITAG},
            {  +2, 23, _ITAG},
            } ;

        _geom._flags
            = JIGSAW_EUCLIDEAN_MESH;

        _geom._vert3._data = &_geom_vert3[0] ;
        _geom._vert3._size = +16;

        _geom._tria3._data = &_geom_tria3[0] ;
        _geom._tria3._size = +24;

        _geom._bound._data = &_geom_bound[0] ;
        _geom._bound._size = +36;

    /*-------------------------------- build JIGSAW tria. */

        _jjig._verbosity = _verb;

        _jjig._hfun_hmax = 3./3.;
        _jjig._hfun_scal =
            JIGSAW_HFUN_RELATIVE;

        _jjig._mesh_rad2 = 1.50 ;
        _jjig._mesh_rad3 = 2.50 ;

        _jjig._geom_feat =   +1 ;   // do "sharp" geom.
        _jjig._mesh_top1 =   +1 ;

        _jjig._mesh_dims =   +3 ;   // make 3-dim cell

        _jjig._mesh_kern =
            JIGSAW_KERN_DELAUNAY;

        _retv = jigsaw (
            & _jjig ,               // the config. opt.
            & _geom ,               // geom. data
               NULL ,               // empty init. obj.
               NULL ,               // empty hfun. obj.
            & _mesh ) ;

    /*-------------------------------- print JIGSAW tria. */

        if (_verb > 0 )
        output_msh_data_3(&_mesh);

        jigsaw_free_msh_t(&_mesh);

        printf (
       "[3d_d] JIGSAW returned code : %d \n", _retv) ;

        fflush (stdout) ;


        return _retv ;
    }

#   ifndef __SKIP_MAIN__
    int main () { return test3d_d(1) ; }
#   endif//__SKIP_MAIN__



