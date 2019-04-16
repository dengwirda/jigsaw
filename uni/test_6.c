
//  gcc -Wall test_6.c 
//  -Xlinker -rpath=../lib 
//  -L ../lib -ljigsaw -o test_6

//  Use JIGSAW to mesh a simple domain, but starting from
//  user-defined initial-conditions. 

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
        jigsaw_msh_t _init ;
        jigsaw_msh_t _mesh ;
  
        jigsaw_init_jig_t(&_jjig) ;
        jigsaw_init_msh_t(&_geom) ;
        jigsaw_init_msh_t(&_init) ;
        jigsaw_init_msh_t(&_mesh) ;
 
    /* 
    --------------------------------------------------------
     * A simple square domain
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
    
    /*-------------------------------- form init. config. */
  
        jigsaw_VERT2_t _point[4] = {
            { {0., 0.}, +0 } ,
            { {0., .5}, +0 } ,
            { {0., 1.}, +0 } ,
            { {.5, .5}, +0 }
            } ;
        
        jigsaw_EDGE2_t _edges[3] = {
            { {+0, +1}, -1 } ,      // -1 => "un-refinable"
            { {+1, +2}, -1 } ,
            { {+2, +3}, -1 }
            } ;

        _init._flags
            = JIGSAW_EUCLIDEAN_MESH;
        
        _init._vert2._data = &_point[0] ;
        _init._vert2._size = +4 ;

        _init._edge2._data = &_edges[0] ;
        _init._edge2._size = +3 ;  

    /*-------------------------------- build JIGSAW tria. */
        
        _jjig._verbosity =   +1 ;
        
        _jjig._hfun_hmax = 0.33 ;
        _jjig._hfun_scal = 
            JIGSAW_HFUN_RELATIVE;
            
        _jjig._mesh_dims =   +2 ;
        _jjig._geom_feat =   +1 ;
        _jjig._mesh_top1 =   +1 ;

        _jjig._optm_iter =   +0 ;
        
        _retv = jigsaw (
            &_jjig ,    // the config. opts
            &_geom ,    // geom. data
            &_init ,    // init. data 
              NULL ,    // empty hfun. data
            &_mesh ) ;
 
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
            printf("%d, %d, %d\n",
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
       "JIGSAW returned code : %d \n",_retv);
 
 
        return _retv ;
    }



