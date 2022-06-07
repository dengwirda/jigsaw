
#   include "../inc/lib_jigsaw.h"

#   include "stdio.h"

    void output_msh_data_2 (
        jigsaw_msh_t *_mesh
        )
    {
    /*-------------------------------- helper: print msh. */
        if (_mesh->_vert2._size > 0)
            printf("\n VERT2: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_vert2._size ;
                   ++_ipos )
        {
            printf("%1.4f, %1.4f\n" ,
            _mesh->_vert2.
                _data[_ipos]._ppos[0],
            _mesh->_vert2.
                _data[_ipos]._ppos[1]
                ) ;
        }

        if (_mesh->_edge2._size > 0)
            printf("\n EDGE2: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_edge2._size ;
                   ++_ipos )
        {
            printf("%d, %d, %d\n",
            _mesh->_edge2.
                _data[_ipos]._node[0],
            _mesh->_edge2.
                _data[_ipos]._node[1],
            _mesh->_edge2.
                _data[_ipos]._itag
                ) ;
        }

        if (_mesh->_tria3._size > 0)
            printf("\n TRIA3: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_tria3._size ;
                   ++_ipos )
        {
            printf("%d, %d, %d, %d\n",
            _mesh->_tria3.
                _data[_ipos]._node[0],
            _mesh->_tria3.
                _data[_ipos]._node[1],
            _mesh->_tria3.
                _data[_ipos]._node[2],
            _mesh->_tria3.
                _data[_ipos]._itag
                ) ;
        }
    }

    void output_msh_data_3 (
        jigsaw_msh_t *_mesh
        )
    {
    /*-------------------------------- helper: print msh. */
        if (_mesh->_vert3._size > 0)
            printf("\n VERT3: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_vert3._size ;
                   ++_ipos )
        {
            printf("%1.4f, %1.4f, %1.4f\n" ,
            _mesh->_vert3.
                _data[_ipos]._ppos[0],
            _mesh->_vert3.
                _data[_ipos]._ppos[1],
            _mesh->_vert3.
                _data[_ipos]._ppos[2]
                ) ;
        }

        if (_mesh->_edge2._size > 0)
            printf("\n EDGE2: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_edge2._size ;
                   ++_ipos )
        {
            printf("%d, %d, %d\n",
            _mesh->_edge2.
                _data[_ipos]._node[0],
            _mesh->_edge2.
                _data[_ipos]._node[1],
            _mesh->_edge2.
                _data[_ipos]._itag
                ) ;
        }

        if (_mesh->_tria3._size > 0)
            printf("\n TRIA3: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_tria3._size ;
                   ++_ipos )
        {
            printf("%d, %d, %d, %d\n",
            _mesh->_tria3.
                _data[_ipos]._node[0],
            _mesh->_tria3.
                _data[_ipos]._node[1],
            _mesh->_tria3.
                _data[_ipos]._node[2],
            _mesh->_tria3.
                _data[_ipos]._itag
                ) ;
        }

        if (_mesh->_tria4._size > 0)
            printf("\n TRIA4: \n\n") ;

        for (size_t _ipos = +0;
                _ipos != _mesh->_tria4._size ;
                   ++_ipos )
        {
            printf("%d, %d, %d, %d, %d \n" ,
            _mesh->_tria4.
                _data[_ipos]._node[0],
            _mesh->_tria4.
                _data[_ipos]._node[1],
            _mesh->_tria4.
                _data[_ipos]._node[2],
            _mesh->_tria4.
                _data[_ipos]._node[3],
            _mesh->_tria4.
                _data[_ipos]._itag
                ) ;
        }
    }



