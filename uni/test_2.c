
//  gcc -Wall test_2.c 
//  -Xlinker -rpath=../lib/LNX-64 
//  -L ../lib/LNX-64 -ljigsaw64r -o test_2

#   include "../inc/lib_jigsaw.h"
    
#   include "stdio.h"
    
    int main (
        int          _argc , 
        char       **_argv
        )
    {
        int _retv = 0;
        
    /*-------------------------------- setup JIGSAW types */      
        jigsaw_jig_t _jjig ;
        jigsaw_msh_t _geom ;
        jigsaw_msh_t _hfun ;
        jigsaw_msh_t _mesh ;
  
        jigsaw_init_jig_t(&_jjig) ;
        jigsaw_init_msh_t(&_geom) ;
        jigsaw_init_msh_t(&_hfun) ;
        jigsaw_init_msh_t(&_mesh) ;
 
    /* 
    --------------------------------------------------------
     * JIGSAW's "mesh" is a piecewise linear complex:
    --------------------------------------------------------
     *
     *             e:2
     *  v:3 o---------------o v:2
     *      |               |
     *      |               |
     *      |               |
     *  e:3 |               | e:1
     *      |               |
     *      |               |
     *      |               |
     *  v:0 o---------------o v:1
     *             e:0
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
     * JIGSAW's "mesh" is a piecewise linear complex:
    --------------------------------------------------------
     *
     *  v:3 o---------------o v:2
     *      | \    t:2    / |
     *      |   \       /   |
     *      |     \ 4 /     |
     *      | t:3   o   t:1 |
     *      |     /   \     |
     *      |   /       \   |
     *      | /    t:0    \ |
     *  v:0 o---------------o v:1
     *
    --------------------------------------------------------
     */
     
        jigsaw_VERT2_t _hfun_vert2[5] = {   // setup hfun.
            { {0., 0.}, +0 } ,
            { {1., 0.}, +0 } ,
            { {1., 1.}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .5}, +0 }
            } ;
            
        jigsaw_TRIA3_t _hfun_tria3[4] = {
            { {+0, +1, +4}, +0 } ,
            { {+1, +2, +4}, +0 } ,
            { {+2, +3, +4}, +4 } ,
            { {+3, +0, +4}, +0 }
            } ;
            
        real_t         _hfun_value[5] = {
             .1 ,
             .1 ,
             .1 ,
             .1 ,
             .2
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
        
        _jjig._verbosity = +1 ;
        
        _jjig._hfun_scal = 
            JIGSAW_HFUN_ABSOLUTE ;
            
        _jjig._hfun_hmax = 1. ;
        _jjig._hfun_hmin = 0. ;
 
 
        _retv = jigsaw_make_mesh (
            &_jjig, &_geom,
                &_hfun, &_mesh
             ) ;
 
    /*-------------------------------- print JIGSAW tria. */

        printf("\n VERT2: \n\n") ;

        for (indx_t _ipos = +0; 
                _ipos != _mesh._vert2._size ; 
                   ++_ipos )
        {
            printf("%1.4f, %1.4f\n",
            _mesh._vert2.
                _data[_ipos]._ppos[0],
            _mesh._vert2.
                _data[_ipos]._ppos[1]
                ) ;
        }
        
        printf("\n TRIA3: \n\n") ;
 
        for (indx_t _ipos = +0; 
                _ipos != _mesh._tria3._size ; 
                   ++_ipos )
        {
            printf("%i, %i, %i\n",
            _mesh._tria3.
                _data[_ipos]._node[0],
            _mesh._tria3.
                _data[_ipos]._node[1],
            _mesh._tria3.
                _data[_ipos]._node[2]
                ) ;
        }
 
        jigsaw_free_msh_t(&_mesh);
 
        printf (
    "JIGSAW returned code: %i \n", _retv) ;
 
 
        return _retv ;
    }



