
    /*
    --------------------------------------------------------
     * INIT-MSH_T: init. msh_t for lib_jigsaw.
    --------------------------------------------------------
     *
     * This program may be freely redistributed under the
     * condition that the copyright notices (including this
     * entire header) are not removed, and no compensation
     * is received through use of the software.  Private,
     * research, and institutional use is free.  You may
     * distribute modified versions of this code UNDER THE
     * CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE
     * TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE
     * ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE
     * MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR
     * NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution
     * of this code as part of a commercial system is
     * permissible ONLY BY DIRECT ARRANGEMENT WITH THE
     * AUTHOR.  (If you are not directly supplying this
     * code to a customer, and you are instead telling them
     * how they can obtain it for free, then you are not
     * required to make any arrangement with me.)
     *
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 30 May, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __INIT_MSH_T__
#   define __INIT_MSH_T__

    /*
    --------------------------------------------------------
     * INIT-MSH_T: init. msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    void  jigsaw_init_msh_t (               // init. msh_t
        jigsaw_msh_t *_mesh
        )
    {
        if (_mesh != nullptr)
        {
        _mesh->_flags = JIGSAW_NULL_FLAG ;

        _mesh->_vert2._size = +0 ;
        _mesh->_vert2._data = nullptr ;

        _mesh->_vert3._size = +0 ;
        _mesh->_vert3._data = nullptr ;

        _mesh->_radii._size = +0 ;
        _mesh->_radii._data = nullptr ;

        _mesh->_power._size = +0 ;
        _mesh->_power._data = nullptr ;

        _mesh->_edge2._size = +0 ;
        _mesh->_edge2._data = nullptr ;

        _mesh->_tria3._size = +0 ;
        _mesh->_tria3._data = nullptr ;

        _mesh->_quad4._size = +0 ;
        _mesh->_quad4._data = nullptr ;

        _mesh->_tria4._size = +0 ;
        _mesh->_tria4._data = nullptr ;

        _mesh->_hexa8._size = +0 ;
        _mesh->_hexa8._data = nullptr ;

        _mesh->_wedg6._size = +0 ;
        _mesh->_wedg6._data = nullptr ;

        _mesh->_pyra5._size = +0 ;
        _mesh->_pyra5._data = nullptr ;

        _mesh->_bound._size = +0 ;
        _mesh->_bound._data = nullptr ;

        _mesh->_seed2._size = +0 ;
        _mesh->_seed2._data = nullptr ;

        _mesh->_seed3._size = +0 ;
        _mesh->_seed3._data = nullptr ;

        _mesh->_value._size = +0 ;
        _mesh->_value._data = nullptr ;

        _mesh->_slope._size = +0 ;
        _mesh->_slope._data = nullptr ;

        _mesh->_xgrid._size = +0 ;
        _mesh->_xgrid._data = nullptr ;

        _mesh->_ygrid._size = +0 ;
        _mesh->_ygrid._data = nullptr ;

        _mesh->_zgrid._size = +0 ;
        _mesh->_zgrid._data = nullptr ;
        }
    }

#   endif//__lib_jigsaw


    /*
    --------------------------------------------------------
     * FREE-MSH_T: _free msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    void  jigsaw_free_msh_t (                // free-alloc
        jigsaw_msh_t *_mesh
        )
    {
    /*-------------------------- dealloc. msh_t data type */
        if (_mesh != nullptr)
        {
        _mesh->_flags = JIGSAW_NULL_FLAG ;

        jigsaw_free_vert2(&_mesh->_vert2);
        jigsaw_free_vert3(&_mesh->_vert3);
        jigsaw_free_reals(&_mesh->_power);

        jigsaw_free_reals(&_mesh->_radii);

        jigsaw_free_edge2(&_mesh->_edge2);
        jigsaw_free_tria3(&_mesh->_tria3);
        jigsaw_free_quad4(&_mesh->_quad4);
        jigsaw_free_tria4(&_mesh->_tria4);
        jigsaw_free_hexa8(&_mesh->_hexa8);
        jigsaw_free_wedg6(&_mesh->_wedg6);
        jigsaw_free_pyra5(&_mesh->_pyra5);

        jigsaw_free_bound(&_mesh->_bound);

        jigsaw_free_reals(&_mesh->_xgrid);
        jigsaw_free_reals(&_mesh->_ygrid);
        jigsaw_free_reals(&_mesh->_zgrid);

        jigsaw_free_reals(&_mesh->_value);
        jigsaw_free_reals(&_mesh->_slope);

        }
    }

#   endif//__lib_jigsaw

    /*
    --------------------------------------------------------
     * ALLOC-MSH_T: alloc. msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    void jigsaw_alloc_vert2 (               // _new-alloc
        jigsaw_VERT2_array_t *_vert2 ,
        size_t _size
        )
    {
        _vert2->_size = _size ;
        _vert2->_data =(jigsaw_VERT2_t*)
            std::malloc(
        _size * sizeof (jigsaw_VERT2_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _vert2->_data[_item].
                _ppos[ 0 ] = (real_t) +0.;
            _vert2->_data[_item].
                _ppos[ 1 ] = (real_t) +0.;

            _vert2->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_vert3 (               // _new-alloc
        jigsaw_VERT3_array_t *_vert3 ,
        size_t _size
        )
    {
        _vert3->_size = _size ;
        _vert3->_data =(jigsaw_VERT3_t*)
            std::malloc(
        _size * sizeof (jigsaw_VERT3_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _vert3->_data[_item].
                _ppos[ 0 ] = (real_t) +0.;
            _vert3->_data[_item].
                _ppos[ 1 ] = (real_t) +0.;
            _vert3->_data[_item].
                _ppos[ 2 ] = (real_t) +0.;

            _vert3->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_edge2 (               // _new-alloc
        jigsaw_EDGE2_array_t *_edge2 ,
        size_t _size
        )
    {
        _edge2->_size = _size ;
        _edge2->_data =(jigsaw_EDGE2_t*)
            std::malloc(
        _size * sizeof (jigsaw_EDGE2_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _edge2->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _edge2->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;

            _edge2->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_tria3 (               // _new-alloc
        jigsaw_TRIA3_array_t *_tria3 ,
        size_t _size
        )
    {
        _tria3->_size = _size ;
        _tria3->_data =(jigsaw_TRIA3_t*)
            std::malloc(
        _size * sizeof (jigsaw_TRIA3_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _tria3->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _tria3->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;
            _tria3->_data[_item].
                _node[ 2 ] = (indx_t) +0 ;

            _tria3->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_quad4 (               // _new-alloc
        jigsaw_QUAD4_array_t *_quad4 ,
        size_t _size
        )
    {
        _quad4->_size = _size ;
        _quad4->_data =(jigsaw_QUAD4_t*)
            std::malloc(
        _size * sizeof (jigsaw_QUAD4_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _quad4->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _quad4->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;
            _quad4->_data[_item].
                _node[ 2 ] = (indx_t) +0 ;
            _quad4->_data[_item].
                _node[ 3 ] = (indx_t) +0 ;

            _quad4->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_tria4 (               // _new-alloc
        jigsaw_TRIA4_array_t *_tria4 ,
        size_t _size
        )
    {
        _tria4->_size = _size ;
        _tria4->_data =(jigsaw_TRIA4_t*)
            std::malloc(
        _size * sizeof (jigsaw_TRIA4_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _tria4->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _tria4->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;
            _tria4->_data[_item].
                _node[ 2 ] = (indx_t) +0 ;
            _tria4->_data[_item].
                _node[ 3 ] = (indx_t) +0 ;

            _tria4->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_hexa8 (               // _new-alloc
        jigsaw_HEXA8_array_t *_hexa8 ,
        size_t _size
        )
    {
        _hexa8->_size = _size ;
        _hexa8->_data =(jigsaw_HEXA8_t*)
            std::malloc(
        _size * sizeof (jigsaw_HEXA8_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _hexa8->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 2 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 3 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 4 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 5 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 6 ] = (indx_t) +0 ;
            _hexa8->_data[_item].
                _node[ 7 ] = (indx_t) +0 ;

            _hexa8->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_wedg6 (               // _new-alloc
        jigsaw_WEDG6_array_t *_wedg6 ,
        size_t _size
        )
    {
        _wedg6->_size = _size ;
        _wedg6->_data =(jigsaw_WEDG6_t*)
            std::malloc(
        _size * sizeof (jigsaw_WEDG6_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _wedg6->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _wedg6->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;
            _wedg6->_data[_item].
                _node[ 2 ] = (indx_t) +0 ;
            _wedg6->_data[_item].
                _node[ 3 ] = (indx_t) +0 ;
            _wedg6->_data[_item].
                _node[ 4 ] = (indx_t) +0 ;
            _wedg6->_data[_item].
                _node[ 5 ] = (indx_t) +0 ;

            _wedg6->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_pyra5 (               // _new-alloc
        jigsaw_PYRA5_array_t *_pyra5 ,
        size_t _size
        )
    {
        _pyra5->_size = _size ;
        _pyra5->_data =(jigsaw_PYRA5_t*)
            std::malloc(
        _size * sizeof (jigsaw_PYRA5_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _pyra5->_data[_item].
                _node[ 0 ] = (indx_t) +0 ;
            _pyra5->_data[_item].
                _node[ 1 ] = (indx_t) +0 ;
            _pyra5->_data[_item].
                _node[ 2 ] = (indx_t) +0 ;
            _pyra5->_data[_item].
                _node[ 3 ] = (indx_t) +0 ;
            _pyra5->_data[_item].
                _node[ 4 ] = (indx_t) +0 ;

            _pyra5->_data[_item].
                _itag      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_bound (               // _new-alloc
        jigsaw_BOUND_array_t *_bound ,
        size_t _size
        )
    {
        _bound->_size = _size ;
        _bound->_data =(jigsaw_BOUND_t*)
            std::malloc(
        _size * sizeof (jigsaw_BOUND_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _bound->_data[_item].
                _itag      = (indx_t) +0 ;
            _bound->_data[_item].
                _indx      = (indx_t) +0 ;
            _bound->_data[_item].
                _kind      = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_index (               // _new-alloc
        jigsaw_INDEX_array_t *_index ,
        size_t _size
        )
    {
        _index->_size = _size ;
        _index->_data = (indx_t *)
            std::malloc (
                _size * sizeof (indx_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _index->
              _data[_item] = (indx_t) +0 ;
        }
    }

    void jigsaw_alloc_reals (               // _new-alloc
        jigsaw_REALS_array_t *_reals ,
        size_t _size
        )
    {
        _reals->_size = _size ;
        _reals->_data = (real_t *)
            std::malloc (
                _size * sizeof (real_t)) ;

        for (size_t _item  = (size_t) +0 ;
                    _item != _size ;
                  ++_item  )
        {
            _reals->
              _data[_item] = (real_t) +0.;
        }
    }

#   endif//__lib_jigsaw

    /*
    --------------------------------------------------------
     * FREE-MSH_T: _free msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    void jigsaw_free_vert2 (               // free-alloc
        jigsaw_VERT2_array_t *_vert2
        )
    {
        if (_vert2 != nullptr &&
            _vert2->_size > 0 )
        {
            _vert2->_size = 0 ;
            std::free(_vert2->_data) ;
        }
    }

    void jigsaw_free_vert3 (               // free-alloc
        jigsaw_VERT3_array_t *_vert3
        )
    {
        if (_vert3 != nullptr &&
            _vert3->_size > 0 )
        {
            _vert3->_size = 0 ;
            std::free(_vert3->_data) ;
        }
    }

    void jigsaw_free_edge2 (               // free-alloc
        jigsaw_EDGE2_array_t *_edge2
        )
    {
        if (_edge2 != nullptr &&
            _edge2->_size > 0 )
        {
            _edge2->_size = 0 ;
            std::free(_edge2->_data) ;
        }
    }

    void jigsaw_free_tria3 (               // free-alloc
        jigsaw_TRIA3_array_t *_tria3
        )
    {
        if (_tria3 != nullptr &&
            _tria3->_size > 0 )
        {
            _tria3->_size = 0 ;
            std::free(_tria3->_data) ;
        }
    }

    void jigsaw_free_quad4 (               // free-alloc
        jigsaw_QUAD4_array_t *_quad4
        )
    {
        if (_quad4 != nullptr &&
            _quad4->_size > 0 )
        {
            _quad4->_size = 0 ;
            std::free(_quad4->_data) ;
        }
    }

    void jigsaw_free_tria4 (               // free-alloc
        jigsaw_TRIA4_array_t *_tria4
        )
    {
        if (_tria4 != nullptr &&
            _tria4->_size > 0 )
        {
            _tria4->_size = 0 ;
            std::free(_tria4->_data) ;
        }
    }

    void jigsaw_free_hexa8 (               // free-alloc
        jigsaw_HEXA8_array_t *_hexa8
        )
    {
        if (_hexa8 != nullptr &&
            _hexa8->_size > 0 )
        {
            _hexa8->_size = 0 ;
            std::free(_hexa8->_data) ;
        }
    }

    void jigsaw_free_wedg6 (               // free-alloc
        jigsaw_WEDG6_array_t *_wedg6
        )
    {
        if (_wedg6 != nullptr &&
            _wedg6->_size > 0 )
        {
            _wedg6->_size = 0 ;
            std::free(_wedg6->_data) ;
        }
    }

    void jigsaw_free_pyra5 (               // free-alloc
        jigsaw_PYRA5_array_t *_pyra5
        )
    {
        if (_pyra5 != nullptr &&
            _pyra5->_size > 0 )
        {
            _pyra5->_size = 0 ;
            std::free(_pyra5->_data) ;
        }
    }

    void jigsaw_free_bound (               // free-alloc
        jigsaw_BOUND_array_t *_bound
        )
    {
        if (_bound != nullptr &&
            _bound->_size > 0 )
        {
            _bound->_size = 0 ;
            std::free(_bound->_data) ;
        }
    }

    void jigsaw_free_index (               // free-alloc
        jigsaw_INDEX_array_t *_index
        )
    {
        if (_index != nullptr &&
            _index->_size > 0 )
        {
            _index->_size = 0 ;
            std::free(_index->_data) ;
        }
    }

    void jigsaw_free_reals (               // free-alloc
        jigsaw_REALS_array_t *_reals
        )
    {
        if (_reals != nullptr &&
            _reals->_size > 0 )
        {
            _reals->_size = 0 ;
            std::free(_reals->_data) ;
        }
    }

#   endif//__lib_jigsaw

#   endif



