
    /*
    --------------------------------------------------------
     * MSH-READ-API: parse MSH_t data into mesh data.
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
     * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant
     * this code in any way whatsoever.  This code is
     * provided "as-is" to be used at your own risk.
     *
     * THE CONTRIBUTORS include:
     * (a) The University of Sydney
     * (b) The Massachusetts Institute of Technology
     * (c) Columbia University
     * (d) The National Aeronautics & Space Administration
     * (e) Los Alamos National Laboratory
     *
    --------------------------------------------------------
     *
     * Last updated: 09 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

     // from msh_read.hpp


    /*
    --------------------------------------------------------
     * MSH-READER-API: read MSH_t mesh data
    --------------------------------------------------------
     */

    class msh_reader_api
    {
    public  :

    typedef containers::array <
            std::string > string_tokens;

    string_tokens           _errs ;

    std::int32_t            _ftag = -1 ;
    std:: size_t            _ndim = +0 ;

    jmsh_kind::enum_data    _kind = 
               jmsh_kind::euclidean_mesh ;

    public  :

    template <
        typename  dest_type
             >
    __normal_call void_type read_data (
        jigsaw_msh_t const&_msht,
        dest_type        &&_dest
        )
    {
        switch (_msht._flags)
        {
        case JIGSAW_EUCLIDEAN_MESH:
            this->_kind = jmsh_kind::euclidean_mesh;
            break;

        case JIGSAW_ELLIPSOID_MESH:
            this->_kind = jmsh_kind::ellipsoid_mesh;
            break;

        case JIGSAW_EUCLIDEAN_GRID:
            this->_kind = jmsh_kind::euclidean_grid;
            break;

        case JIGSAW_ELLIPSOID_GRID:
            this->_kind = jmsh_kind::ellipsoid_grid;
            break;

        default:
            this->_errs.push_tail(
                "Unsupported msh kind") ;
            return;
        }

        if (std::max({_msht._vert2._size,
                      _msht._vert3._size}) > 0 &&
            std::max({_msht._xgrid._size,
                      _msht._ygrid._size,
                      _msht._zgrid._size}) > 0 )
        {
            _errs.push_tail(
                "POINT/COORD mismatch") ;
        }

        if (_msht._xgrid._size > 0) this->_ndim =  +1;
        if (_msht._ygrid._size > 0) this->_ndim =  +2;
        if (_msht._zgrid._size > 0) this->_ndim =  +3;

        if (_msht._seed2._size > 0) this->_ndim =  +2;
        if (_msht._vert2._size > 0) this->_ndim =  +2;
        if (_msht._seed3._size > 0) this->_ndim =  +3;
        if (_msht._vert3._size > 0) this->_ndim =  +3;

        _dest.initialise() ;

        _dest.push_mshid(this->_ftag, 
                         this->_kind) ;
        _dest.push_ndims(this->_ndim) ;

        if (this->_kind == jmsh_kind::ellipsoid_mesh ||
            this->_kind == jmsh_kind::ellipsoid_grid )
        {
            double _rrad[3] = {1., 1., 1.} ;
            if (_msht._radii._size == 3)
            {
                _rrad[0] = _msht._radii._data[0];
                _rrad[1] = _msht._radii._data[1];
                _rrad[2] = _msht._radii._data[2];
            }
            else
            if (_msht._radii._size == 1)
            {
                _rrad[0] = _msht._radii._data[0];
                _rrad[1] = _msht._radii._data[0];
                _rrad[2] = _msht._radii._data[0];
            }
            _dest.push_radii(_rrad);
        }

        if (this->_ndim != 2 && _msht._vert2._size > 0)
        {
            _errs.push_tail(
                "NDIMS/POINT mismatch") ;
        }
        else
        {
        _dest.open_point( _msht._vert2._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._vert2._size ; 
                ++_ipos )
        {
            _dest.push_point( _ipos, 
                 _msht._vert2._data[_ipos]._ppos,
                 _msht._vert2._data[_ipos]._itag
            ) ;
        }
        }

        if (this->_ndim != 2 && _msht._seed2._size > 0)
        {
            _errs.push_tail(
                "NDIMS/SEEDS mismatch") ;
        }
        else
        {
        _dest.open_seeds( _msht._seed2._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._seed2._size ; 
                ++_ipos )
        {
            _dest.push_seeds( _ipos, 
                 _msht._seed2._data[_ipos]._ppos,
                 _msht._seed2._data[_ipos]._itag
            ) ;
        }
        }

        if (this->_ndim != 3 && _msht._vert3._size > 0)
        {
            _errs.push_tail(
                "NDIMS/POINT mismatch") ;
        }
        else
        {
        _dest.open_point( _msht._vert3._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._vert3._size ; 
                ++_ipos )
        {
            _dest.push_point( _ipos, 
                 _msht._vert3._data[_ipos]._ppos,
                 _msht._vert3._data[_ipos]._itag
            ) ;
        }
        }

        if (this->_ndim != 3 && _msht._seed3._size > 0)
        {
            _errs.push_tail(
                "NDIMS/SEEDS mismatch") ;
        }
        else
        {
        _dest.open_seeds( _msht._seed3._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._seed3._size ; 
                ++_ipos )
        {
            _dest.push_seeds( _ipos, 
                 _msht._seed3._data[_ipos]._ppos,
                 _msht._seed3._data[_ipos]._itag
            ) ;
        }
        }

        _dest.open_edge2( _msht._edge2._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._edge2._size ; 
                ++_ipos )
        {
            _dest.push_edge2( _ipos, 
                 _msht._edge2._data[_ipos]._node,
                 _msht._edge2._data[_ipos]._itag
            ) ;
        }

        _dest.open_tria3( _msht._tria3._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._tria3._size ; 
                ++_ipos )
        {
            _dest.push_tria3( _ipos, 
                 _msht._tria3._data[_ipos]._node,
                 _msht._tria3._data[_ipos]._itag
            ) ;
        }

        _dest.open_quad4( _msht._quad4._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._quad4._size ; 
                ++_ipos )
        {
            _dest.push_quad4( _ipos, 
                 _msht._quad4._data[_ipos]._node,
                 _msht._quad4._data[_ipos]._itag
            ) ;
        }

        _dest.open_tria4( _msht._tria4._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._tria4._size ; 
                ++_ipos )
        {
            _dest.push_tria4( _ipos, 
                 _msht._tria4._data[_ipos]._node,
                 _msht._tria4._data[_ipos]._itag
            ) ;
        }

        _dest.open_hexa8( _msht._hexa8._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._hexa8._size ; 
                ++_ipos )
        {
            _dest.push_hexa8( _ipos, 
                 _msht._hexa8._data[_ipos]._node,
                 _msht._hexa8._data[_ipos]._itag
            ) ;
        }

        _dest.open_wedg6( _msht._wedg6._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._wedg6._size ; 
                ++_ipos )
        {
            _dest.push_wedg6( _ipos, 
                 _msht._wedg6._data[_ipos]._node,
                 _msht._wedg6._data[_ipos]._itag
            ) ;
        }

        _dest.open_pyra5( _msht._pyra5._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._pyra5._size ; 
                ++_ipos )
        {
            _dest.push_pyra5( _ipos, 
                 _msht._pyra5._data[_ipos]._node,
                 _msht._pyra5._data[_ipos]._itag
            ) ;
        }

        _dest.open_coord( _msht._xgrid._size, 1);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._xgrid._size ; 
                ++_ipos )
        {
            _dest.push_coord( 1, _ipos, 
                 _msht._xgrid._data[_ipos]
            ) ;
        }

        _dest.open_coord( _msht._ygrid._size, 2);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._ygrid._size ; 
                ++_ipos )
        {
            _dest.push_coord( 2, _ipos, 
                 _msht._ygrid._data[_ipos]
            ) ;
        }

        _dest.open_coord( _msht._zgrid._size, 3);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._zgrid._size ; 
                ++_ipos )
        {
            _dest.push_coord( 3, _ipos, 
                 _msht._zgrid._data[_ipos]
            ) ;
        }

        _dest.open_power( _msht._power._size, 1);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._power._size ; 
                ++_ipos )
        {
            _dest.push_power( _ipos, 
                &_msht._power._data[_ipos]
            ) ;
        }

        if (_msht._power._size  >  +0 &&
            _msht._power._size != 
            _msht._vert2._size &&
            _msht._power._size != 
            _msht._vert3._size  )
        {
            _errs.push_tail(
                "POWER/POINT mismatch") ;
        }

        _dest.open_value( _msht._value._size, 1);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._value._size ; 
                ++_ipos )
        {
            _dest.push_value( _ipos, 
                &_msht._value._data[_ipos]
            ) ;
        }

        if (_msht._value._size  >  +1 &&
            _msht._value._size != 
            _msht._vert2._size &&
            _msht._value._size != 
            _msht._vert3._size &&
            _msht._value._size !=
        std::max<size_t>(+1u, _msht._xgrid._size) *
        std::max<size_t>(+1u, _msht._ygrid._size) *
        std::max<size_t>(+1u, _msht._zgrid._size) )
        {
            _errs.push_tail(
                "VALUE/POINT mismatch") ;
        }

        _dest.open_slope( _msht._slope._size, 1);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._slope._size ; 
                ++_ipos )
        {
            _dest.push_slope( _ipos, 
                &_msht._slope._data[_ipos]
            ) ;
        }

        if (_msht._slope._size  >  +1 &&
            _msht._slope._size != 
            _msht._vert2._size &&
            _msht._slope._size != 
            _msht._vert3._size &&
            _msht._slope._size !=
        std::max<size_t>(+1u, _msht._xgrid._size) *
        std::max<size_t>(+1u, _msht._ygrid._size) *
        std::max<size_t>(+1u, _msht._zgrid._size) )
        {
            _errs.push_tail(
                "SLOPE/POINT mismatch") ;
        }

        _dest.open_bound( _msht._bound._size);
        for (auto _ipos = (size_t) +0 ;
                  _ipos!= _msht._bound._size ; 
                ++_ipos )
        {
            _dest.push_bound( _ipos, 
                 _msht._bound._data[_ipos]._itag,
                 _msht._bound._data[_ipos]._indx,
                 _msht._bound._data[_ipos]._kind
            ) ;
        }

        {
           _dest.postscript();
        }
    }

    } ;



