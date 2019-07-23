
//  gcc -Wall test_9.c 
//  -Xlinker -rpath=../lib 
//  -L ../lib -ljigsaw -o test_9

//  Use JIGSAW to generate a uniform mesh for a spheroidal 
//  domain. 

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
        jigsaw_init_jig_t(&_jjig) ;        

        jigsaw_msh_t _geom ;
        jigsaw_init_msh_t(&_geom) ;

        jigsaw_msh_t _mesh ;
        jigsaw_init_msh_t(&_mesh) ;        
 
    /*-------------------------------- setup JIGSAW geom. */

        real_t _radii[3] = {+1., +1., +1.
            } ;

        _geom._flags = 
            JIGSAW_ELLIPSOID_MESH ;

        _geom._radii._data = &_radii[0] ;
        _geom._radii._size = +3 ;

    /*-------------------------------- build JIGSAW tria. */
        
        _jjig._verbosity =   +1 ;
        
        _jjig._hfun_hmax = 0.80 ;
        _jjig._hfun_scal = 
            JIGSAW_HFUN_ABSOLUTE;
            
        _jjig._mesh_dims =   +2 ;
        
        _retv = jigsaw (
            &_jjig ,    // the config. opts
            &_geom ,    // geom. data
              NULL ,    // empty init. data 
              NULL ,    // empty hfun. data
            &_mesh ) ;
 
    /*-------------------------------- print JIGSAW tria. */

        printf("\n VERT3: \n\n") ;

        for (indx_t _ipos = +0; 
                _ipos != _mesh._vert3._size ; 
                   ++_ipos )
        {
            printf(
             "%1.4f, %1.4f, %1.4f\n",
            _mesh._vert3.
                _data[_ipos]._ppos[0],
            _mesh._vert3.
                _data[_ipos]._ppos[1],
            _mesh._vert3.
                _data[_ipos]._ppos[2]
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



