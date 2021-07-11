
    /*
    --------------------------------------------------------
     * HFN-LOAD: parse *.MSH file into HFUN data.
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
     * Last updated: 07 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   include "msh_read.hpp"

#   ifndef __HFN_LOAD__
#   define __HFN_LOAD__

    /*
    --------------------------------------------------------
     * HFUN-FROM-JMSH: read *.JMSH file into HFUN data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type hfun_from_jmsh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_hfun
        )
    {
        class hfun_reader:
            public jmsh_reader_base
        {
        public  :
            hfun_data           *_hfun ;

            std::int32_t         _ftag ;
            jmsh_kind::
            enum_data            _kind ;
            std:: size_t         _ndim ;

        public  :
        __normal_call hfun_reader (
            hfun_data*_hsrc =  nullptr
            ) : _hfun(_hsrc) {}
    /*-------------------------------- read MSHID section */
        __normal_call void_type push_mshid (
            std::int32_t  _FTAG ,
            jmsh_kind::enum_data _KIND
            )
        {
            this->_ftag = _FTAG ;
            this->_kind = _KIND ;
            this->
           _hfun->_kind = _KIND ;
        }
    /*-------------------------------- read NDIMS section */
        __normal_call void_type push_ndims (
            std:: size_t  _NDIM
            )
        {
            this->_ndim = _NDIM ;
            this->
           _hfun->_ndim = _NDIM ;
        }
    /*-------------------------------- read POINT section */
        __normal_call void_type push_point (
            std:: size_t/*_ipos*/ ,
            double       *_pval ,
            std::int32_t/*_itag*/
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                hfun_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];

                this->_hfun->
                    _euclidean_mesh_2d.
                _mesh.push_node(_ndat, false) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                hfun_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) = _pval[2];

                this->_hfun->
                    _euclidean_mesh_3d.
                _mesh.push_node(_ndat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                typename
                hfun_data::ellipsoid_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) = +0.0 ;

                this->_hfun->
                    _ellipsoid_mesh_3d.
                _mesh.push_node(_ndat, false) ;
            }
        }
    /*-------------------------------- read TRIA3 section */
        __normal_call void_type push_tria3 (
            std:: size_t/*_ipos*/ ,
            std::int32_t *_node ,
            std::int32_t/*_itag*/
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                hfun_data::euclidean_mesh_2d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];

                this->_hfun->
                    _euclidean_mesh_2d.
               _mesh.push_tri3(_tdat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                typename
                hfun_data::ellipsoid_mesh_3d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];

                this->_hfun->
                    _ellipsoid_mesh_3d.
               _mesh.push_tri3(_tdat, false) ;
            }
        }
    /*-------------------------------- read TRIA4 section */
        __normal_call void_type push_tria4 (
            std:: size_t/*_ipos*/ ,
            std::int32_t *_node ,
            std::int32_t/*_itag*/
            )
        {
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                hfun_data::euclidean_mesh_3d
                    ::tri4_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                _tdat.node(3) = _node[3];

                this->_hfun->
                    _euclidean_mesh_3d.
               _mesh.push_tri4(_tdat, false) ;
            }
        }
    /*-------------------------------- read COORD section */
        __normal_call void_type push_coord (
            std:: size_t  _idim ,
            std:: size_t/*_irow*/ ,
            double        _ppos
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                if (_idim == +1)
                {
                    this->_hfun->
                    _euclidean_grid_2d.
                        _xpos.push_tail(_ppos) ;
                }
                else
                if (_idim == +2)
                {
                    this->_hfun->
                    _euclidean_grid_2d.
                        _ypos.push_tail(_ppos) ;
                }
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                if (_idim == +1)
                {
                    this->_hfun->
                    _euclidean_grid_3d.
                        _xpos.push_tail(_ppos) ;
                }
                else
                if (_idim == +2)
                {
                    this->_hfun->
                    _euclidean_grid_3d.
                        _ypos.push_tail(_ppos) ;
                }
                else
                if (_idim == +3)
                {
                    this->_hfun->
                    _euclidean_grid_3d.
                        _zpos.push_tail(_ppos) ;
                }
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_grid)
            {
                if (_idim == +1)
                {
                    this->_hfun->
                    _ellipsoid_grid_3d.
                        _xpos.push_tail(_ppos) ;
                }
                else
                if (_idim == +2)
                {
                     this->_hfun->
                    _ellipsoid_grid_3d.
                        _ypos.push_tail(_ppos) ;
                }
            }
        }
    /*-------------------------------- open VALUE section */
        __normal_call void_type open_value (
            std:: size_t  _nrow ,
            std:: size_t/*_nval*/
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_2d.
                        _hval.set_count(_nrow) ;
            }
            else
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_2d.
                        _hmat.set_count(_nrow) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_3d.
                        _hval.set_count(_nrow) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_3d.
                        _hmat.set_count(_nrow) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                this->_hfun->
                    _ellipsoid_mesh_3d.
                        _hval.set_count(_nrow) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_grid)
            {
                this->_hfun->
                    _ellipsoid_grid_3d.
                        _hmat.set_count(_nrow) ;
            }
        }
    /*-------------------------------- read VALUE section */
        __normal_call void_type push_value (
            std:: size_t  _ipos ,
            float        *_vval
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_2d.
                        _hval[_ipos] = *_vval;
            }
            else
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_2d.
                        _hmat[_ipos] = *_vval;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_3d.
                        _hval[_ipos] = *_vval;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_3d.
                        _hmat[_ipos] = *_vval;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                this->_hfun->
                    _ellipsoid_mesh_3d.
                        _hval[_ipos] = *_vval;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_grid)
            {
                this->_hfun->
                    _ellipsoid_grid_3d.
                        _hmat[_ipos] = *_vval;
            }
        }
    /*-------------------------------- open SLOPE section */
        __normal_call void_type open_slope (
            std:: size_t  _nrow ,
            std:: size_t/*_nval*/
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_2d.
                        _dhdx.set_count(_nrow) ;
            }
            else
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_2d.
                        _dhdx.set_count(_nrow) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_3d.
                        _dhdx.set_count(_nrow) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_3d.
                        _dhdx.set_count(_nrow) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                this->_hfun->
                    _ellipsoid_mesh_3d.
                        _dhdx.set_count(_nrow) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_grid)
            {
                this->_hfun->
                    _ellipsoid_grid_3d.
                        _dhdx.set_count(_nrow) ;
            }
        }
    /*-------------------------------- read SLOPE section */
        __normal_call void_type push_slope (
            std:: size_t  _ipos ,
            float        *_vval
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_2d.
                        _dhdx[_ipos] = *_vval;
            }
            else
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_2d.
                        _dhdx[_ipos] = *_vval;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                this->_hfun->
                    _euclidean_mesh_3d.
                        _dhdx[_ipos] = *_vval;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_3d.
                        _dhdx[_ipos] = *_vval;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                this->_hfun->
                    _ellipsoid_mesh_3d.
                        _dhdx[_ipos] = *_vval;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_grid)
            {
                this->_hfun->
                    _ellipsoid_grid_3d.
                        _dhdx[_ipos] = *_vval;
            }
        }
    /*---------------------------------- parse RADII data */
        __normal_call void_type push_radii (
            double       *_erad
            )
        {
            this->_hfun->_ellipsoid_mesh_3d.
                _radA = _erad[ 0] ;
            this->_hfun->_ellipsoid_mesh_3d.
                _radB = _erad[ 1] ;
            this->_hfun->_ellipsoid_mesh_3d.
                _radC = _erad[ 2] ;

            this->_hfun->_ellipsoid_grid_3d.
                _radA = _erad[ 0] ;
            this->_hfun->_ellipsoid_grid_3d.
                _radB = _erad[ 1] ;
            this->_hfun->_ellipsoid_grid_3d.
                _radC = _erad[ 2] ;
        }
        } ;

    /*---------------------------------- parse HFUN. file */
        iptr_type _errv  = __no_error ;

        try
        {
            jmsh_reader   _jmsh ;
            std::ifstream _file ;
            _file. open(
            _jcfg._hfun_file, std::ifstream::in) ;

            if (_file.is_open() )
            {
                _jmsh.read_file (
                _file, hfun_reader(&_hfun));
            }
            else
            {
                _jlog.push(
            "**parse error: file not found!\n" ) ;

                _errv = __file_not_located ;
            }
            _file.close ();

            for (auto _iter  =
                      _jmsh._errs.head();
                      _iter !=
                      _jmsh._errs.tend();
                    ++_iter  )
            {
                _jlog.push(
            "**parse error: " + * _iter + "\n" ) ;

                _errv = __invalid_argument ;
            }
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * HFUN-FROM-MSHT: read MSH_t data into HFUN data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type hfun_from_msht (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_hfun ,
        jigsaw_msh_t const&_hmsh
        )
    {
        iptr_type _errv  = __no_error ;

        try
        {

        __unreferenced (_jlog) ;
        __unreferenced (_jcfg) ;

        if (_hmsh._flags ==
                JIGSAW_EUCLIDEAN_MESH )
        {
            if (_hmsh._vert2._size > 0)
            {
    /*--------------------------------- euclidean-mesh-2d */
            _hfun._kind
                = jmsh_kind::euclidean_mesh ;
            _hfun._ndim = +2;

            if (_hmsh._vert2._size !=
                _hmsh._value._size )
            return __invalid_argument ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._vert2._size ;
                    ++_ipos )
            {
                typename
                hfun_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _hmsh.
                _vert2._data[_ipos]._ppos[0];
                _ndat.pval(1) = _hmsh.
                _vert2._data[_ipos]._ppos[1];

                _hfun._euclidean_mesh_2d.
                _mesh.push_node(_ndat, false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._tria3._size ;
                    ++_ipos )
            {
                typename
                hfun_data::euclidean_mesh_2d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _hmsh.
                _tria3._data[_ipos]._node[0];
                _tdat.node(1) = _hmsh.
                _tria3._data[_ipos]._node[1];
                _tdat.node(2) = _hmsh.
                _tria3._data[_ipos]._node[2];

                _hfun._euclidean_mesh_2d.
                _mesh.push_tri3(_tdat, false) ;
            }

            _hfun._euclidean_mesh_2d._hval.
                set_count(_hmsh._value._size) ;

            _hfun._euclidean_mesh_2d._dhdx.
                set_count(_hmsh._slope._size) ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._value._size ;
                    ++_ipos )
            {
                _hfun._euclidean_mesh_2d.
                    _hval[_ipos] =
                    _hmsh._value._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._slope._size ;
                    ++_ipos )
            {
                _hfun._euclidean_mesh_2d.
                    _dhdx[_ipos] =
                    _hmsh._slope._data[_ipos] ;
            }

            }
            else
            if (_hmsh._vert3._size > 0)
            {
    /*--------------------------------- euclidean-mesh-3d */
            _hfun._kind
                = jmsh_kind::euclidean_mesh ;
            _hfun._ndim = +3;

            if (_hmsh._vert3._size !=
                _hmsh._value._size )
            return __invalid_argument ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._vert3._size ;
                    ++_ipos )
            {
                typename
                hfun_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _hmsh.
                _vert3._data[_ipos]._ppos[0];
                _ndat.pval(1) = _hmsh.
                _vert3._data[_ipos]._ppos[1];
                _ndat.pval(2) = _hmsh.
                _vert3._data[_ipos]._ppos[2];

                _hfun._euclidean_mesh_3d.
                _mesh.push_node(_ndat, false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._tria4._size ;
                    ++_ipos )
            {
                typename
                hfun_data::euclidean_mesh_3d
                    ::tri4_type _tdat ;
                _tdat.node(0) = _hmsh.
                _tria4._data[_ipos]._node[0];
                _tdat.node(1) = _hmsh.
                _tria4._data[_ipos]._node[1];
                _tdat.node(2) = _hmsh.
                _tria4._data[_ipos]._node[2];
                _tdat.node(3) = _hmsh.
                _tria4._data[_ipos]._node[3];

                _hfun._euclidean_mesh_3d.
                _mesh.push_tri4(_tdat, false) ;
            }

            _hfun._euclidean_mesh_3d._hval.
                set_count(_hmsh._value._size) ;

            _hfun._euclidean_mesh_3d._dhdx.
                set_count(_hmsh._slope._size) ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._value._size ;
                    ++_ipos )
            {
                _hfun._euclidean_mesh_3d.
                    _hval[_ipos] =
                    _hmsh._value._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._slope._size ;
                    ++_ipos )
            {
                _hfun._euclidean_mesh_3d.
                    _dhdx[_ipos] =
                    _hmsh._slope._data[_ipos] ;
            }

            }
        }
        else
        if (_hmsh._flags ==
                JIGSAW_ELLIPSOID_MESH )
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            _hfun._kind
                = jmsh_kind::ellipsoid_mesh ;
            _hfun._ndim = +2;

            if (_hmsh._vert2._size !=
                _hmsh._value._size )
            return __invalid_argument ;

            if (_hmsh._radii._size==+3)
            {
            _hfun._ellipsoid_mesh_3d.
                _radA = _hmsh._radii._data[0] ;
            _hfun._ellipsoid_mesh_3d.
                _radB = _hmsh._radii._data[1] ;
            _hfun._ellipsoid_mesh_3d.
                _radC = _hmsh._radii._data[2] ;
            }
            else
            if (_hmsh._radii._size==+1)
            {
            _hfun._ellipsoid_mesh_3d.
                _radA = _hmsh._radii._data[0] ;
            _hfun._ellipsoid_mesh_3d.
                _radB = _hmsh._radii._data[0] ;
            _hfun._ellipsoid_mesh_3d.
                _radC = _hmsh._radii._data[0] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._vert2._size ;
                    ++_ipos )
            {
                typename
                hfun_data::ellipsoid_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _hmsh.
                _vert2._data[_ipos]._ppos[0];
                _ndat.pval(1) = _hmsh.
                _vert2._data[_ipos]._ppos[1];
                _ndat.pval(2) = + 0.0 ;

                _hfun._ellipsoid_mesh_3d.
                _mesh.push_node(_ndat, false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._tria3._size ;
                    ++_ipos )
            {
                typename
                hfun_data::ellipsoid_mesh_3d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _hmsh.
                _tria3._data[_ipos]._node[0];
                _tdat.node(1) = _hmsh.
                _tria3._data[_ipos]._node[1];
                _tdat.node(2) = _hmsh.
                _tria3._data[_ipos]._node[2];

                _hfun._ellipsoid_mesh_3d.
                _mesh.push_tri3(_tdat, false) ;
            }

            _hfun._ellipsoid_mesh_3d._hval.
                set_count(_hmsh._value._size) ;

            _hfun._ellipsoid_mesh_3d._dhdx.
                set_count(_hmsh._slope._size) ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._value._size ;
                    ++_ipos )
            {
                _hfun._ellipsoid_mesh_3d.
                    _hval[_ipos] =
                    _hmsh._value._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._slope._size ;
                    ++_ipos )
            {
                _hfun._ellipsoid_mesh_3d.
                    _dhdx[_ipos] =
                    _hmsh._slope._data[_ipos] ;
            }

        }
        else
        if (_hmsh._flags ==
                JIGSAW_EUCLIDEAN_GRID )
        {
            if (_hmsh._zgrid._size== 0)
            {
    /*--------------------------------- euclidean-grid-2d */
            _hfun._kind
                = jmsh_kind::euclidean_grid ;
            _hfun._ndim = +2;

            _hfun._euclidean_grid_2d._xpos.
                set_count(_hmsh._xgrid._size) ;
            _hfun._euclidean_grid_2d._ypos.
                set_count(_hmsh._ygrid._size) ;

            _hfun._euclidean_grid_2d._hmat.
                set_count(_hmsh._value._size) ;

            _hfun._euclidean_grid_2d._dhdx.
                set_count(_hmsh._slope._size) ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._xgrid._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d.
                    _xpos[_ipos] =
                    _hmsh._xgrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._ygrid._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d.
                    _ypos[_ipos] =
                    _hmsh._ygrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._value._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d.
                    _hmat[_ipos] =
                    _hmsh._value._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._slope._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d.
                    _dhdx[_ipos] =
                    _hmsh._slope._data[_ipos] ;
            }

            }
            else
            if (_hmsh._zgrid._size!= 0)
            {
    /*--------------------------------- euclidean-grid-3d */
            _hfun._kind
                = jmsh_kind::euclidean_grid ;
            _hfun._ndim = +3;

            _hfun._euclidean_grid_3d._xpos.
                set_count(_hmsh._xgrid._size) ;
            _hfun._euclidean_grid_3d._ypos.
                set_count(_hmsh._ygrid._size) ;
            _hfun._euclidean_grid_3d._zpos.
                set_count(_hmsh._zgrid._size) ;

            _hfun._euclidean_grid_3d._hmat.
                set_count(_hmsh._value._size) ;

            _hfun._euclidean_grid_3d._dhdx.
                set_count(_hmsh._slope._size) ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._xgrid._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d.
                    _xpos[_ipos] =
                    _hmsh._xgrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._ygrid._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d.
                    _ypos[_ipos] =
                    _hmsh._ygrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._zgrid._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d.
                    _zpos[_ipos] =
                    _hmsh._zgrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._value._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d.
                    _hmat[_ipos] =
                    _hmsh._value._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._slope._size ;
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d.
                    _dhdx[_ipos] =
                    _hmsh._slope._data[_ipos] ;
            }

            }
        }
        else
        if (_hmsh._flags ==
                JIGSAW_ELLIPSOID_GRID )
        {
            if (_hmsh._zgrid._size== 0)
            {
    /*--------------------------------- ellipsoid-grid-3d */
            _hfun._kind
                = jmsh_kind::ellipsoid_grid ;
            _hfun._ndim = +2;

            if (_hmsh._radii._size==+3)
            {
            _hfun._ellipsoid_grid_3d.
                _radA = _hmsh._radii._data[0] ;
            _hfun._ellipsoid_grid_3d.
                _radB = _hmsh._radii._data[1] ;
            _hfun._ellipsoid_grid_3d.
                _radC = _hmsh._radii._data[2] ;
            }
            else
            if (_hmsh._radii._size==+1)
            {
            _hfun._ellipsoid_grid_3d.
                _radA = _hmsh._radii._data[0] ;
            _hfun._ellipsoid_grid_3d.
                _radB = _hmsh._radii._data[0] ;
            _hfun._ellipsoid_grid_3d.
                _radC = _hmsh._radii._data[0] ;
            }

            _hfun._ellipsoid_grid_3d._xpos.
                set_count(_hmsh._xgrid._size) ;
            _hfun._ellipsoid_grid_3d._ypos.
                set_count(_hmsh._ygrid._size) ;

            _hfun._ellipsoid_grid_3d._hmat.
                set_count(_hmsh._value._size) ;

            _hfun._ellipsoid_grid_3d._dhdx.
                set_count(_hmsh._slope._size) ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._xgrid._size ;
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d.
                    _xpos[_ipos] =
                    _hmsh._xgrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._ygrid._size ;
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d.
                    _ypos[_ipos] =
                    _hmsh._ygrid._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._value._size ;
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d.
                    _hmat[_ipos] =
                    _hmsh._value._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _hmsh._slope._size ;
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d.
                    _dhdx[_ipos] =
                    _hmsh._slope._data[_ipos] ;
            }

            }
        }

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * READ-HFUN: read HFUN input file.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type read_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_hfun
        )
    {
        return hfun_from_jmsh (
                _jcfg, _jlog, _hfun ) ;
    }

    /*
    --------------------------------------------------------
     * COPY-HFUN: copy HFUN input data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type copy_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_hfun ,
        jigsaw_msh_t const&_hmsh
        )
    {
        return hfun_from_msht (
           _jcfg, _jlog, _hfun, _hmsh);
    }

    /*
    --------------------------------------------------------
     * TEST-HFUN: test HFUN data validity.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type test_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_hfun
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;

        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            fp32_type _hmin =
            std::numeric_limits<fp32_type>::max() ;
            fp32_type _smin =
            std::numeric_limits<fp32_type>::max() ;

            iptr_type _imin =
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min() ;

            iptr_type _nmax = +0 ;

            for (auto _iter  = _hfun.
                _euclidean_mesh_2d._hval.head() ;
                      _iter != _hfun.
                _euclidean_mesh_2d._hval.tend() ;
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
            }

            for (auto _iter  = _hfun.
                _euclidean_mesh_2d._dhdx.head() ;
                      _iter != _hfun.
                _euclidean_mesh_2d._dhdx.tend() ;
                    ++_iter  )
            {
                _smin = std::min(_smin, *_iter) ;
            }

            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh.node().head();
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh.node().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _nmax += +1  ;
            }

            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh.tri3().head();
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh.tri3().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imin = std::min(
                    _imin, _iter->node(2)) ;

                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(2)) ;
            }

            auto _hnum = _hfun.
                _euclidean_mesh_2d._hval.count();

            auto _gnum = _hfun.
                _euclidean_mesh_2d._dhdx.count();

            if (_gnum > +0 &&
                    _gnum != 1 &&
                        _gnum != _hnum)
            {
                _jlog.push (
    "**input error: DHDX. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_hmin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: HFUN. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_smin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: DHDX. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_imin < +0 || _imax>=_nmax)
            {
                _jlog.push (
    "**input error: HFUN. tria. indexing is incorrect.\n") ;

                _errv = __invalid_indexing ;
            }
        }
        else
        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-2d */
            fp32_type _hmin =
            std::numeric_limits<fp32_type>::infinity();
            fp32_type _smin =
            std::numeric_limits<fp32_type>::infinity();

            for (auto _iter  = _hfun.
            _euclidean_grid_2d._hmat.head();
                      _iter != _hfun.
            _euclidean_grid_2d._hmat.tend();
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
            }

            for (auto _iter  = _hfun.
            _euclidean_grid_2d._dhdx.head();
                      _iter != _hfun.
            _euclidean_grid_2d._dhdx.tend();
                    ++_iter  )
            {
                _smin = std::min(_smin, *_iter) ;
            }

            bool_type _mono = true;

            if(!_hfun._euclidean_grid_2d._xpos.empty())
            for (auto _iter  = _hfun.
            _euclidean_grid_2d._xpos.head();
                      _iter != _hfun.
            _euclidean_grid_2d._xpos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }
            }

            if(!_hfun._euclidean_grid_2d._ypos.empty())
            for (auto _iter  = _hfun.
            _euclidean_grid_2d._ypos.head();
                      _iter != _hfun.
            _euclidean_grid_2d._ypos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }
            }

            auto _xnum = _hfun.
                _euclidean_grid_2d._xpos.count();

            auto _ynum = _hfun.
                _euclidean_grid_2d._ypos.count();

            auto _hnum = _hfun.
                _euclidean_grid_2d._hmat.count();

            auto _gnum = _hfun.
                _euclidean_grid_2d._dhdx.count();

            if (_hnum != _xnum * _ynum)
            {
                _jlog.push (
    "**input error: HFUN. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_gnum > +0 &&
                    _gnum != 1 &&
                        _gnum != _hnum)
            {
                _jlog.push (
    "**input error: DHDX. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_hmin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: HFUN. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_smin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: DHDX. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (!_mono)
            {
                _jlog.push (
    "**input error: grid must be monotonic increasing.\n") ;

                _errv = __invalid_argument ;
            }
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            fp32_type _hmin =
            std::numeric_limits<fp32_type>::max() ;
            fp32_type _smin =
            std::numeric_limits<fp32_type>::max() ;

            iptr_type _imin =
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min() ;

            iptr_type _nmax = +0 ;

            for (auto _iter  = _hfun.
                _euclidean_mesh_3d._hval.head() ;
                      _iter != _hfun.
                _euclidean_mesh_3d._hval.tend() ;
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
            }

            for (auto _iter  = _hfun.
                _euclidean_mesh_3d._dhdx.head() ;
                      _iter != _hfun.
                _euclidean_mesh_3d._dhdx.tend() ;
                    ++_iter  )
            {
                _smin = std::min(_smin, *_iter) ;
            }

            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh.node().head();
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh.node().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _nmax += +1  ;
            }

            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh.tri4().head();
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh.tri4().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imin = std::min(
                    _imin, _iter->node(2)) ;
                _imin = std::min(
                    _imin, _iter->node(3)) ;

                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(2)) ;
                _imax = std::max(
                    _imax, _iter->node(3)) ;
            }

            auto _hnum = _hfun.
                _euclidean_mesh_3d._hval.count();

            auto _gnum = _hfun.
                _euclidean_mesh_3d._dhdx.count();

            if (_gnum > +0 &&
                    _gnum != 1 &&
                        _gnum != _hnum)
            {
                _jlog.push (
    "**input error: DHDX. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_hmin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: HFUN. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_smin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: DHDX. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_imin < +0 || _imax>=_nmax)
            {
                _jlog.push (
    "**input error: HFUN. tria. indexing is incorrect.\n") ;

                _errv = __invalid_indexing ;
            }
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-3d */
            fp32_type _hmin =
            std::numeric_limits<fp32_type>::infinity();
            fp32_type _smin =
            std::numeric_limits<fp32_type>::infinity();

            for (auto _iter  = _hfun.
            _euclidean_grid_3d._hmat.head();
                      _iter != _hfun.
            _euclidean_grid_3d._hmat.tend();
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
            }

            for (auto _iter  = _hfun.
            _euclidean_grid_3d._dhdx.head();
                      _iter != _hfun.
            _euclidean_grid_3d._dhdx.tend();
                    ++_iter  )
            {
                _smin = std::min(_smin, *_iter) ;
            }

            bool_type _mono = true;

            if(!_hfun._euclidean_grid_3d._xpos.empty())
            for (auto _iter  = _hfun.
            _euclidean_grid_3d._xpos.head();
                      _iter != _hfun.
            _euclidean_grid_3d._xpos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }
            }

            if(!_hfun._euclidean_grid_3d._ypos.empty())
            for (auto _iter  = _hfun.
            _euclidean_grid_3d._ypos.head();
                      _iter != _hfun.
            _euclidean_grid_3d._ypos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }
            }

            if(!_hfun._euclidean_grid_3d._zpos.empty())
            for (auto _iter  = _hfun.
            _euclidean_grid_3d._zpos.head();
                      _iter != _hfun.
            _euclidean_grid_3d._zpos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }
            }

            auto _xnum = _hfun.
                _euclidean_grid_3d._xpos.count();

            auto _ynum = _hfun.
                _euclidean_grid_3d._ypos.count();

            auto _znum = _hfun.
                _euclidean_grid_3d._zpos.count();

            auto _hnum = _hfun.
                _euclidean_grid_3d._hmat.count();

            auto _gnum = _hfun.
                _euclidean_grid_3d._dhdx.count();

            if (_hnum!=_xnum*_ynum*_znum)
            {
                _jlog.push (
    "**input error: HFUN. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_gnum > +0 &&
                    _gnum != 1 &&
                        _gnum != _hnum)
            {
                _jlog.push (
    "**input error: DHDX. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_hmin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: HFUN. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_smin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: DHDX. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (!_mono)
            {
                _jlog.push (
    "**input error: grid must be monotonic increasing.\n") ;

                _errv = __invalid_argument ;
            }
        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            if (_hfun._ellipsoid_mesh_3d.
                _radA <= (real_type)  +0. ||
                _hfun._ellipsoid_mesh_3d.
                _radB <= (real_type)  +0. ||
                _hfun._ellipsoid_mesh_3d.
                _radC <= (real_type)  +0. )
            {
                _jlog.push (
    "**input error: HFUN. RADII entries are incorrect.\n") ;

                _errv = __invalid_argument ;
            }

            fp32_type _hmin =
            std::numeric_limits<fp32_type>::max() ;
            fp32_type _smin =
            std::numeric_limits<fp32_type>::max() ;

            real_type static const _PI   =
           (real_type)std::atan(+1.0) * 4. ;

        //  careful with the way PI truncations onto float

        //  expanded range so that we don't throw warnings
        //  due to rounding issues...

            real_type static const _XMIN =
                    (real_type) -2.1 * _PI ;
            real_type static const _XMAX =
                    (real_type) +2.1 * _PI ;
            real_type static const _YMIN =
                    (real_type) -1.1 * _PI ;
            real_type static const _YMAX =
                    (real_type) +1.1 * _PI ;

            real_type _xmin =  _XMAX;
            real_type _xmax =  _XMIN;
            real_type _ymin =  _YMAX;
            real_type _ymax =  _YMIN;

            iptr_type _imin =
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min() ;

            iptr_type _nmax = +0 ;

            for (auto _iter  = _hfun.
                _ellipsoid_mesh_3d._hval.head() ;
                      _iter != _hfun.
                _ellipsoid_mesh_3d._hval.tend() ;
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
            }

            for (auto _iter  = _hfun.
                _ellipsoid_mesh_3d._dhdx.head() ;
                      _iter != _hfun.
                _ellipsoid_mesh_3d._dhdx.tend() ;
                    ++_iter  )
            {
                _smin = std::min(_smin, *_iter) ;
            }

            for (auto _iter  = _hfun.
            _ellipsoid_mesh_3d._mesh.node().head();
                      _iter != _hfun.
            _ellipsoid_mesh_3d._mesh.node().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _xmin = std::min(
                    _xmin, _iter->pval(0)) ;
                _xmax = std::max(
                    _xmax, _iter->pval(0)) ;

                _ymin = std::min(
                    _ymin, _iter->pval(1)) ;
                _ymax = std::max(
                    _ymax, _iter->pval(1)) ;

                _nmax += +1  ;
            }

            for (auto _iter  = _hfun.
            _ellipsoid_mesh_3d._mesh.tri3().head();
                      _iter != _hfun.
            _ellipsoid_mesh_3d._mesh.tri3().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imin = std::min(
                    _imin, _iter->node(2)) ;

                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(2)) ;
            }

            auto _hnum = _hfun.
                _ellipsoid_mesh_3d._hval.count();

            auto _gnum = _hfun.
                _ellipsoid_mesh_3d._dhdx.count();

            if (_xmin <  _XMIN ||
                _xmax >  _XMAX )
            {
                _jlog.push (
    "**input error: XPOS. must be in [-1.*pi, +1.*pi].\n") ;

                _errv = __invalid_argument ;
            }

            if (_ymin <  _YMIN ||
                _ymax >  _YMAX )
            {
                _jlog.push (
    "**input error: YPOS. must be in [-.5*pi, +.5*pi].\n") ;

                _errv = __invalid_argument ;
            }

            if (_gnum > +0 &&
                    _gnum != 1 &&
                        _gnum != _hnum)
            {
                _jlog.push (
    "**input error: DHDX. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_hmin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: HFUN. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_smin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: DHDX. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_imin < +0 || _imax>=_nmax)
            {
                _jlog.push (
    "**input error: HFUN. tria. indexing is incorrect.\n") ;

                _errv = __invalid_indexing ;
            }
        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_grid)
        {
    /*--------------------------------- ellipsoid-grid-3d */
            if (_hfun._ellipsoid_grid_3d.
                _radA <= (real_type)  +0. ||
                _hfun._ellipsoid_grid_3d.
                _radB <= (real_type)  +0. ||
                _hfun._ellipsoid_grid_3d.
                _radC <= (real_type)  +0. )
            {
                _jlog.push (
    "**input error: HFUN. RADII entries are incorrect.\n") ;

                _errv = __invalid_argument ;
            }

            fp32_type _hmin =
            std::numeric_limits<fp32_type>::infinity();
            fp32_type _smin =
            std::numeric_limits<fp32_type>::infinity();

            real_type static const _PI   =
           (real_type)std::atan(+1.0) * 4. ;

        //  careful with the way PI truncations onto float
        //  expanded range so that we don't throw warnings
        //  due to rounding issues...

            real_type static const _XMIN =
                    (real_type) -2.1 * _PI ;
            real_type static const _XMAX =
                    (real_type) +2.1 * _PI ;
            real_type static const _YMIN =
                    (real_type) -1.1 * _PI ;
            real_type static const _YMAX =
                    (real_type) +1.1 * _PI ;

            real_type _xmin =  _XMAX;
            real_type _xmax =  _XMIN;
            real_type _ymin =  _YMAX;
            real_type _ymax =  _YMIN;

            for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._hmat.head();
                      _iter != _hfun.
            _ellipsoid_grid_3d._hmat.tend();
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
            }

            for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._dhdx.head();
                      _iter != _hfun.
            _ellipsoid_grid_3d._dhdx.tend();
                    ++_iter  )
            {
                _smin = std::min(_smin, *_iter) ;
            }

            bool_type _mono = true;

            if(!_hfun._ellipsoid_grid_3d._xpos.empty())
            for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._xpos.head();
                      _iter != _hfun.
            _ellipsoid_grid_3d._xpos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }

                _xmin = std::min(
                        _xmin , *_iter) ;
                _xmax = std::max(
                        _xmax , *_iter) ;
            }

            if(!_hfun._ellipsoid_grid_3d._ypos.empty())
            for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._ypos.head();
                      _iter != _hfun.
            _ellipsoid_grid_3d._ypos.tail();
                    ++_iter  )
            {
                if (*(_iter+1) < *(_iter+0))
                {
                    _mono = false; break;
                }

                _ymin = std::min(
                        _ymin , *_iter) ;
                _ymax = std::max(
                        _ymax , *_iter) ;
            }

            auto _xnum = _hfun.
                _ellipsoid_grid_3d._xpos.count();

            auto _ynum = _hfun.
                _ellipsoid_grid_3d._ypos.count();

            auto _hnum = _hfun.
                _ellipsoid_grid_3d._hmat.count();

            auto _gnum = _hfun.
                _ellipsoid_grid_3d._dhdx.count();

            if (_xmin <  _XMIN ||
                _xmax >  _XMAX )
            {
                _jlog.push (
    "**input error: XPOS. must be in [-1.*pi, +1.*pi].\n") ;

                _errv = __invalid_argument ;
            }

            if (_ymin <  _YMIN ||
                _ymax >  _YMAX )
            {
                _jlog.push (
    "**input error: YPOS. must be in [-.5*pi, +.5*pi].\n") ;

                _errv = __invalid_argument ;
            }

            if (_hnum != _xnum * _ynum)
            {
                _jlog.push (
    "**input error: HFUN. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_gnum > +0 &&
                    _gnum != 1 &&
                        _gnum != _hnum)
            {
                _jlog.push (
    "**input error: DHDX. matrix incorrect dimensions.\n") ;

                _errv = __invalid_arraydim ;
            }

            if (_hmin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: HFUN. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (_smin <  (fp32_type) +0.)
            {
                _jlog.push (
    "**input error: DHDX. values must be non-negative.\n") ;

                _errv = __invalid_argument ;
            }

            if (!_mono)
            {
                _jlog.push (
    "**input error: grid must be monotonic increasing.\n") ;

                _errv = __invalid_argument ;
            }
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * ECHO-HFUN: print summary of HFUN data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type echo_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_hfun
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;

        fp32_type _hmin =
            +std::numeric_limits
                <real_type>::infinity() ;

        fp32_type _hmax =
            -std::numeric_limits
                <real_type>::infinity() ;

        if (_hfun._ndim == +0)
        {
    /*--------------------------------- constant-value-kd */
            _jlog.push(
                "  CONSTANT-VALUE\n\n") ;

            _jlog.push(
                "  .VAL(H). = " + 
                to_string_prec (
                _hfun._constant_value_kd._hval, 2) + "\n") ;
        }
        else
        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(2) + "\n");

            _jlog.push("\n") ;

            iptr_type _num1 = +0 ;
            iptr_type _num3 = +0 ;

            for (auto _iter  = _hfun.
                _euclidean_mesh_2d._hval.head() ;
                      _iter != _hfun.
                _euclidean_mesh_2d._hval.tend() ;
                    ++_iter  )
            {
                _hmin = std::min(
                            _hmin,*_iter) ;
                _hmax = std::max(
                            _hmax,*_iter) ;
            }

            _jlog.push(
                "  .MIN(H). = " +
                to_string_prec(_hmin, 2) + "\n");
            _jlog.push(
                "  .MAX(H). = " +
                to_string_prec(_hmax, 2) + "\n");

            _jlog.push("  \n") ;

             auto _mptr =
           &_hfun._euclidean_mesh_2d._mesh;

            for (auto _iter  =
                      _mptr->node().head();
                      _iter !=
                      _mptr->node().tend();
                    ++_iter  )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }

            _jlog.push(
                "  |COORD.| = " 
                + std::to_string(_num1) + "\n");

            for (auto _iter  =
                      _mptr->tri3().head();
                      _iter !=
                      _mptr->tri3().tend();
                    ++_iter  )
            {
            if (_iter->mark()>=+0) _num3 += +1 ;
            }

            _jlog.push(
                "  |TRIA-3| = " 
                + std::to_string(_num3) + "\n");

        }
        else
        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-2d */
            _jlog.push(
                "  EUCLIDEAN-GRID\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(2) + "\n");

            _jlog.push("\n") ;

            for (auto _iter  = _hfun.
            _euclidean_grid_2d._hmat.head() ;
                      _iter != _hfun.
            _euclidean_grid_2d._hmat.tend() ;
                    ++_iter )
            {
                _hmin =
                    std::min(_hmin, *_iter) ;
                _hmax =
                    std::max(_hmax, *_iter) ;
            }

            _jlog.push(
                "  .MIN(H). = " +
                to_string_prec(_hmin, 2) + "\n");
            _jlog.push(
                "  .MAX(H). = " +
                to_string_prec(_hmax, 2) + "\n");

            _jlog.push("  \n") ;

            auto _xnum = _hfun.
                _euclidean_grid_2d._xpos.count();

            auto _ynum = _hfun.
                _euclidean_grid_2d._ypos.count();

            _jlog.push(
                "  |XGRID.| = " 
                + std::to_string(_xnum) + "\n") ;
            _jlog.push(
                "  |YGRID.| = " 
                + std::to_string(_ynum) + "\n") ;
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(3) + "\n");

            _jlog.push("\n") ;

            iptr_type _num1 = +0 ;
            iptr_type _num4 = +0 ;

            for (auto _iter  = _hfun.
                _euclidean_mesh_3d._hval.head() ;
                      _iter != _hfun.
                _euclidean_mesh_3d._hval.tend() ;
                    ++_iter  )
            {
                _hmin = std::min(
                            _hmin,*_iter) ;
                _hmax = std::max(
                            _hmax,*_iter) ;
            }

            _jlog.push(
                "  .MIN(H). = " +
                to_string_prec(_hmin, 2) + "\n");
            _jlog.push(
                "  .MAX(H). = " +
                to_string_prec(_hmax, 2) + "\n");

            _jlog.push("  \n") ;

             auto _mptr =
           &_hfun._euclidean_mesh_3d._mesh;

            for (auto _iter  =
                      _mptr->node().head();
                      _iter !=
                      _mptr->node().tend();
                    ++_iter  )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }

            _jlog.push(
                "  |COORD.| = " 
                + std::to_string(_num1) + "\n");

            for (auto _iter  =
                      _mptr->tri4().head();
                      _iter !=
                      _mptr->tri4().tend();
                    ++_iter  )
            {
            if (_iter->mark()>=+0) _num4 += +1 ;
            }

            _jlog.push(
                "  |TRIA-4| = " 
                + std::to_string(_num4) + "\n");

        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-3d */
            _jlog.push(
                "  EUCLIDEAN-GRID\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(3) + "\n");

            _jlog.push("\n") ;

            for (auto _iter  = _hfun.
            _euclidean_grid_3d._hmat.head() ;
                      _iter != _hfun.
            _euclidean_grid_3d._hmat.tend() ;
                    ++_iter )
            {
                _hmin =
                    std::min(_hmin, *_iter) ;
                _hmax =
                    std::max(_hmax, *_iter) ;
            }

            _jlog.push(
                "  .MIN(H). = " +
                to_string_prec(_hmin, 2) + "\n");
            _jlog.push(
                "  .MAX(H). = " +
                to_string_prec(_hmax, 2) + "\n");

            _jlog.push("  \n") ;

            auto _xnum = _hfun.
                _euclidean_grid_3d._xpos.count();

            auto _ynum = _hfun.
                _euclidean_grid_3d._ypos.count();

            auto _znum = _hfun.
                _euclidean_grid_3d._zpos.count();

            _jlog.push(
                "  |XGRID.| = " 
                + std::to_string(_xnum) + "\n") ;
            _jlog.push(
                "  |YGRID.| = " 
                + std::to_string(_ynum) + "\n") ;
            _jlog.push(
                "  |ZGRID.| = " 
                + std::to_string(_znum) + "\n") ;

        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            _jlog.push(
                "  ELLIPSOID-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(2) + "\n");

            iptr_type _num1 = +0 ;
            iptr_type _num3 = +0 ;

            for (auto _iter  = _hfun.
                _ellipsoid_mesh_3d._hval.head() ;
                      _iter != _hfun.
                _ellipsoid_mesh_3d._hval.tend() ;
                    ++_iter  )
            {
                _hmin = std::min(
                            _hmin,*_iter) ;
                _hmax = std::max(
                            _hmax,*_iter) ;
            }

            _jlog.push(
                "  .MIN(H). = " +
                to_string_prec(_hmin, 2) + "\n");
            _jlog.push(
                "  .MAX(H). = " +
                to_string_prec(_hmax, 2) + "\n");

            _jlog.push("  \n") ;

             auto _mptr =
           &_hfun._ellipsoid_mesh_3d._mesh;

            for (auto _iter  =
                      _mptr->node().head();
                      _iter !=
                      _mptr->node().tend();
                    ++_iter  )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }

            _jlog.push(
                "  |COORD.| = " 
                + std::to_string(_num1) + "\n");

            for (auto _iter  =
                      _mptr->tri3().head();
                      _iter !=
                      _mptr->tri3().tend();
                    ++_iter  )
            {
            if (_iter->mark()>=+0) _num3 += +1 ;
            }

            _jlog.push(
                "  |TRIA-3| = " 
                + std::to_string(_num3) + "\n");

        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_grid)
        {
    /*--------------------------------- ellipsoid-grid-3d */
            _jlog.push(
                "  ELLIPSOID-GRID\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(2) + "\n");

            for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._hmat.head() ;
                      _iter != _hfun.
            _ellipsoid_grid_3d._hmat.tend() ;
                    ++_iter )
            {
                _hmin =
                    std::min(_hmin, *_iter) ;
                _hmax =
                    std::max(_hmax, *_iter) ;
            }

            _jlog.push(
                "  .MIN(H). = " +
                to_string_prec(_hmin, 2) + "\n");
            _jlog.push(
                "  .MAX(H). = " +
                to_string_prec(_hmax, 2) + "\n");

            _jlog.push("  \n") ;

            auto _xnum = _hfun.
                _ellipsoid_grid_3d._xpos.count();

            auto _ynum = _hfun.
                _ellipsoid_grid_3d._ypos.count();

            _jlog.push(
                "  |XGRID.| = " 
                + std::to_string(_xnum) + "\n") ;
            _jlog.push(
                "  |YGRID.| = " 
                + std::to_string(_ynum) + "\n") ;

            _jlog.push("  \n") ;

            if (_hfun._ellipsoid_grid_3d._wrap)
            _jlog.push("  PERIODIC = TRUE\n"  ) ;

        }

        _jlog.push("\n") ;

        return (  _errv) ;
    }


#   endif   //__HFN_LOAD__



