
    /*
    --------------------------------------------------------
     * INI-LOAD: parse *.MSH file into INIT data.
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

#   ifndef __INI_LOAD__
#   define __INI_LOAD__

    /*
    --------------------------------------------------------
     * INIT-FROM-JMSH: read *.JMSH file into INIT data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type init_from_jmsh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_init
        )
    {
        class init_reader:
            public jmsh_reader_base
        {
        public  :
            mesh_data           *_init ;

            std::int32_t         _ftag ;
            jmsh_kind::
            enum_data            _kind ;
            std:: size_t         _ndim ;

            real_type            _rrad[3] = {1.} ;

            containers::
                array<real_type> _vpwr ;
        public  :
    /*---------------------------------- construct reader */
        __normal_call init_reader (
            mesh_data*_isrc =  nullptr
            ) : _init(_isrc) {}
    /*---------------------------------- parse MSHID data */
        __normal_call void_type push_mshid (
            std::int32_t  _FTAG ,
            jmsh_kind::enum_data _KIND
            )
        {
            this->_ftag = _FTAG ;
            this->_kind = _KIND ;
            this->
           _init->_kind = _KIND ;
        }
    /*---------------------------------- parse NDIMS data */
        __normal_call void_type push_ndims (
            std:: size_t  _NDIM
            )
        {
            this->_ndim = _NDIM ;
            this->
           _init->_ndim = _NDIM ;
        }
    /*---------------------------------- parse RADII data */
        __normal_call void_type push_radii (
            double       *_erad
            )
        {
            this->_rrad[0] = _erad[ 0] ;
            this->_rrad[1] = _erad[ 1] ;
            this->_rrad[2] = _erad[ 2] ;
        }
    /*---------------------------------- parse POINT data */
        __normal_call void_type push_point (
            std:: size_t/*_ipos*/ ,
            double       *_pval ,
            std::int32_t  _itag
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.itag () = _itag ;

                _ndat.pval(2) =
                        (real_type) +0. ;

                _ndat.fdim () =     +0  ;

                if (_itag < +0)
                    _ndat.feat () =
                        mesh::user_feat ;
                else
                    _ndat.feat () =
                        mesh::null_feat ;

                this->_init->
                   _euclidean_mesh_2d.
               _mesh.push_node(_ndat, false) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) = _pval[2];
                _ndat.itag () = _itag ;

                _ndat.pval(3) =
                        (real_type) +0. ;

                _ndat.fdim () =     +0  ;

                if (_itag < +0)
                    _ndat.feat () =
                        mesh::user_feat ;
                else
                    _ndat.feat () =
                        mesh::null_feat ;

                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_node(_ndat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) =
                        (real_type) +0. ;
                _ndat.itag () = _itag ;

                _ndat.pval(3) =
                        (real_type) +0. ;

                _ndat.fdim () =     +0  ;

                if (_itag < +0)
                    _ndat.feat () =
                        mesh::user_feat ;
                else
                    _ndat.feat () =
                        mesh::null_feat ;

                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_node(_ndat, false) ;
            }
        }
    /*---------------------------------- store POINT data */
        __normal_call void_type save_point (
            )
        {
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
            iptr_type _ipos = + 0 ;

            for (auto
                _iter  = this->_init->
            _euclidean_mesh_3d._mesh.node().head();
                _iter != this->_init->
            _euclidean_mesh_3d._mesh.node().tend();
              ++_iter, ++_ipos )
            {
                real_type _apos[2];
                real_type _ppos[3];

                _apos[0] = _iter->pval (0) ;
                _apos[1] = _iter->pval (1) ;

                math::toR3(
                    _rrad, _apos, _ppos) ;

                _iter->pval(0) = _ppos [0] ;
                _iter->pval(1) = _ppos [1] ;
                _iter->pval(2) = _ppos [2] ;
            }
            }
        }
    /*---------------------------------- parse POWER data */
        __normal_call void_type push_power (
            std:: size_t/*_ipos*/ ,
            double       *_ppwr
            )
        {
            this->_vpwr.push_tail(_ppwr[0]);
        }
    /*---------------------------------- store POWER data */
        __normal_call void_type save_power (
            )
        {
            if (this->_vpwr.count() == +0) return ;

            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
            if (this->_vpwr.count() !=
                this->_init->
            _euclidean_mesh_2d._mesh.node().count())
                return ;

            iptr_type _ipos = +0 ;
            for (auto
                _iter  = this->_init->
            _euclidean_mesh_2d._mesh.node().head();
                _iter != this->_init->
            _euclidean_mesh_2d._mesh.node().tend();
              ++_iter, ++_ipos )
            {
                _iter->pval(2) =
                    this->_vpwr[_ipos] ;
            }
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh ||
               (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh) )
            {
            if (this->_vpwr.count() !=
                this->_init->
            _euclidean_mesh_3d._mesh.node().count())
                return ;

            iptr_type _ipos = +0 ;
            for (auto
                _iter  = this->_init->
            _euclidean_mesh_3d._mesh.node().head();
                _iter != this->_init->
            _euclidean_mesh_3d._mesh.node().tend();
              ++_iter, ++_ipos )
            {
                _iter->pval(3) =
                    this->_vpwr[_ipos] ;
            }
            }
        }
    /*---------------------------------- parse EDGE2 data */
        __normal_call void_type push_edge2 (
            std:: size_t/*_ipos*/ ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                _edat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_2d.
               _mesh.push_edge(_edat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh ||
               (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh) )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                _edat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_edge(_edat, false) ;
            }
        }
    /*---------------------------------- parse TRIA3 data */
        __normal_call void_type push_tria3 (
            std:: size_t/*_ipos*/ ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                _tdat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_2d.
               _mesh.push_tri3(_tdat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh ||
               (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh) )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                _tdat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_tri3(_tdat, false) ;
            }
        }
    /*---------------------------------- parse QUAD4 data */
        __normal_call void_type push_quad4 (
            std:: size_t/*_ipos*/ ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::quad_type _qdat ;
                _qdat.node(0) = _node[0];
                _qdat.node(1) = _node[1];
                _qdat.node(2) = _node[2];
                _qdat.node(3) = _node[3];
                _qdat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_2d.
               _mesh.push_quad(_qdat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh ||
               (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh) )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::quad_type _qdat ;
                _qdat.node(0) = _node[0];
                _qdat.node(1) = _node[1];
                _qdat.node(2) = _node[2];
                _qdat.node(3) = _node[3];
                _qdat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_quad(_qdat, false) ;
            }
        }
    /*---------------------------------- parse TRIA4 data */
        __normal_call void_type push_tria4 (
            std:: size_t/*_ipos*/ ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                // do nothing...
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh ||
               (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh) )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::tri4_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                _tdat.node(3) = _node[3];
                _tdat.itag () = _itag ;

                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_tri4(_tdat, false) ;
            }
        }
    /*---------------------------------- finalise objects */
        __normal_call void_type postscript (
            )
        {
            save_point () ;           // lon-lat to R^3
            save_power () ;           // match to point
        }
        } ;

    /*---------------------------------- parse INIT. file */
        iptr_type _errv  = __no_error ;

        try
        {
            jmsh_reader   _jmsh ;
            std::ifstream _file ;
            _file. open(
            _jcfg._init_file, std::ifstream::in) ;

            if (_file.is_open() )
            {
                _jmsh.read_file (
                _file, init_reader(&_init));
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
     * INIT-FROM-MSHT: read MSH_t data into INIT data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type init_from_msht (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_init ,
        jigsaw_msh_t const&_imsh
        )
    {
        iptr_type _errv  = __no_error ;

        try
        {

        __unreferenced (_jlog) ;
        __unreferenced (_jcfg) ;

        if (_imsh._vert2._size > 0 &&
            _imsh._flags ==
                JIGSAW_EUCLIDEAN_MESH )
        {
    /*--------------------------------- euclidean-mesh-2d */
            _init._kind
                = jmsh_kind::euclidean_mesh ;
            _init._ndim = +2;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._vert2._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _imsh.
                    _vert2._data[_ipos]._ppos[0];
                _ndat.pval(1) = _imsh.
                    _vert2._data[_ipos]._ppos[1];
                _ndat.itag () = _imsh.
                    _vert2._data[_ipos]._itag ;

                _ndat.pval(2) =
                        (real_type) +0. ;

                _ndat.fdim () =     +0  ;

                if (_ndat.itag () < +0)
                    _ndat.feat () =
                        mesh::user_feat ;
                else
                    _ndat.feat () =
                        mesh::null_feat ;

                _init._euclidean_mesh_2d.
                    _mesh.push_node(_ndat,false);
            }

            if (_imsh._power._size != +0 &&
                _imsh._power._size !=
                _imsh._vert2._size  )
                return __invalid_argument ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._power._size ;
                    ++_ipos )
            {
                _init._euclidean_mesh_2d.
                    _mesh.node(_ipos).pval(2) =
                        _imsh._power._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._edge2._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::edge_type _edat ;
                _edat.node(0) = _imsh.
                    _edge2._data[_ipos]._node[0];
                _edat.node(1) = _imsh.
                    _edge2._data[_ipos]._node[1];
                _edat.itag () = _imsh.
                    _edge2._data[_ipos]._itag ;

                _init._euclidean_mesh_2d.
                    _mesh.push_edge(_edat,false);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._tria3._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _imsh.
                    _tria3._data[_ipos]._node[0];
                _tdat.node(1) = _imsh.
                    _tria3._data[_ipos]._node[1];
                _tdat.node(2) = _imsh.
                    _tria3._data[_ipos]._node[2];
                _tdat.itag () = _imsh.
                    _tria3._data[_ipos]._itag ;

                _init._euclidean_mesh_2d.
                    _mesh.push_tri3(_tdat,false);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._quad4._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_2d
                    ::quad_type _qdat ;
                _qdat.node(0) = _imsh.
                    _quad4._data[_ipos]._node[0];
                _qdat.node(1) = _imsh.
                    _quad4._data[_ipos]._node[1];
                _qdat.node(2) = _imsh.
                    _quad4._data[_ipos]._node[2];
                _qdat.node(3) = _imsh.
                    _quad4._data[_ipos]._node[3];
                _qdat.itag () = _imsh.
                    _quad4._data[_ipos]._itag ;

                _init._euclidean_mesh_2d.
                    _mesh.push_quad(_qdat,false);
            }

        }
        else
        if (_imsh._flags ==
                JIGSAW_ELLIPSOID_MESH ||
           (_imsh._vert3._size > 0 &&
            _imsh._flags ==
                JIGSAW_EUCLIDEAN_MESH ) )
        {
    /*--------------------------------- euclidean-mesh-3d */
            _init._kind
                = jmsh_kind::euclidean_mesh ;
            _init._ndim = + 3;

            real_type _rrad[3] = { 1. } ;

            if (_imsh._radii._size == 3 )
            {
                _rrad[0] = _imsh._radii._data[0];
                _rrad[1] = _imsh._radii._data[1];
                _rrad[2] = _imsh._radii._data[2];
            }
            else
            if (_imsh._radii._size == 1 )
            {
                _rrad[0] = _imsh._radii._data[0];
                _rrad[1] = _imsh._radii._data[0];
                _rrad[2] = _imsh._radii._data[0];
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._vert3._size ;
                    ++_ipos )   // from euclidean...
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _imsh.
                    _vert3._data[_ipos]._ppos[0];
                _ndat.pval(1) = _imsh.
                    _vert3._data[_ipos]._ppos[1];
                _ndat.pval(2) = _imsh.
                    _vert3._data[_ipos]._ppos[2];
                _ndat.itag () = _imsh.
                    _vert3._data[_ipos]._itag ;

                _ndat.pval(3) =
                        (real_type) +0. ;

                _ndat.fdim () =     +0  ;

                if (_ndat.itag () < +0)
                    _ndat.feat () =
                        mesh::user_feat ;
                else
                    _ndat.feat () =
                        mesh::null_feat ;

                _init._euclidean_mesh_3d.
                    _mesh.push_node(_ndat,false);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._vert2._size ;
                    ++_ipos )   // from ellipsoid...
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _imsh.
                    _vert2._data[_ipos]._ppos[0];
                _ndat.pval(1) = _imsh.
                    _vert2._data[_ipos]._ppos[1];
                _ndat.itag () = _imsh.
                    _vert2._data[_ipos]._itag ;

                _ndat.pval(3) =
                        (real_type) +0. ;

                _ndat.fdim () =     +0  ;

                real_type _apos[2];
                real_type _ppos[3];

                _apos[0] = _ndat. pval (0) ;
                _apos[1] = _ndat. pval (1) ;

                math::toR3(
                    _rrad, _apos, _ppos) ;

                _ndat. pval(0) = _ppos [0] ;
                _ndat. pval(1) = _ppos [1] ;
                _ndat. pval(2) = _ppos [2] ;

                if (_ndat.itag () < +0)
                    _ndat.feat () =
                        mesh::user_feat ;
                else
                    _ndat.feat () =
                        mesh::null_feat ;

                _init._euclidean_mesh_3d.
                    _mesh.push_node(_ndat,false);
            }

            if (_imsh._power._size != +0 &&
                _imsh._power._size !=
                _imsh._vert3._size &&
                _imsh._power._size !=
                _imsh._vert2._size  )
                return __invalid_argument ;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._power._size ;
                    ++_ipos )
            {
                _init._euclidean_mesh_3d.
                    _mesh.node(_ipos).pval(3) =
                        _imsh._power._data[_ipos] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._edge2._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _imsh.
                    _edge2._data[_ipos]._node[0];
                _edat.node(1) = _imsh.
                    _edge2._data[_ipos]._node[1];
                _edat.itag () = _imsh.
                    _edge2._data[_ipos]._itag ;

                _init._euclidean_mesh_3d.
                    _mesh.push_edge(_edat,false);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._tria3._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _imsh.
                    _tria3._data[_ipos]._node[0];
                _tdat.node(1) = _imsh.
                    _tria3._data[_ipos]._node[1];
                _tdat.node(2) = _imsh.
                    _tria3._data[_ipos]._node[2];
                _tdat.itag () = _imsh.
                    _tria3._data[_ipos]._itag ;

                _init._euclidean_mesh_3d.
                    _mesh.push_tri3(_tdat,false);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._quad4._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::quad_type _qdat ;
                _qdat.node(0) = _imsh.
                    _quad4._data[_ipos]._node[0];
                _qdat.node(1) = _imsh.
                    _quad4._data[_ipos]._node[1];
                _qdat.node(2) = _imsh.
                    _quad4._data[_ipos]._node[2];
                _qdat.node(3) = _imsh.
                    _quad4._data[_ipos]._node[3];
                _qdat.itag () = _imsh.
                    _quad4._data[_ipos]._itag ;

                _init._euclidean_mesh_3d.
                    _mesh.push_quad(_qdat,false);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _imsh._tria4._size ;
                    ++_ipos )
            {
                typename
                mesh_data::euclidean_mesh_3d
                    ::tri4_type _tdat ;
                _tdat.node(0) = _imsh.
                    _tria4._data[_ipos]._node[0];
                _tdat.node(1) = _imsh.
                    _tria4._data[_ipos]._node[1];
                _tdat.node(2) = _imsh.
                    _tria4._data[_ipos]._node[2];
                _tdat.node(3) = _imsh.
                    _tria4._data[_ipos]._node[3];
                _tdat.itag () = _imsh.
                    _tria4._data[_ipos]._itag ;

                _init._euclidean_mesh_3d.
                    _mesh.push_tri4(_tdat,false);
            }
        }
        else
        if (_imsh._flags ==
                JIGSAW_EUCLIDEAN_GRID)
        {
            _errv =
                 __invalid_argument ;
        }
        else
        if (_imsh._flags ==
                JIGSAW_ELLIPSOID_GRID)
        {
            _errv =
                 __invalid_argument ;
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
     * READ-INIT: read init. distribution file.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type read_init (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_init
        )
    {
        return init_from_jmsh (
                _jcfg, _jlog, _init ) ;
    }

    /*
    --------------------------------------------------------
     * COPY-INIT: copy INIT input data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type copy_init (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_init ,
        jigsaw_msh_t const&_imsh
        )
    {
        return init_from_msht (
           _jcfg, _jlog, _init, _imsh);
    }

    /*
    --------------------------------------------------------
     * TEST-INIT: test INIT data validity.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type test_init (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_init
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;

        if (_init._ndim == +2 &&
            _init._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            iptr_type _imin =
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min() ;

            iptr_type _nmax = +0 ;

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.node().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.node().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _nmax += + 1 ;
            }

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.edge().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.edge().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;
            }

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.tri3().head();
                      _iter != _init.
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

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.quad().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.quad().tend();
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

            if (_imin < +0 || _imax>=_nmax)
            {
                _jlog.push (
    "**input error: INIT. cell. indexing is incorrect.\n") ;

                _errv = __invalid_indexing ;
            }
        }
        else
        if (_init._ndim == +3 &&
            _init._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            iptr_type _imin =
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min() ;

            iptr_type _nmax = +0 ;

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.node().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.node().tend();
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;

                _nmax += + 1 ;
            }

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.edge().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.edge().tend();
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;
            }

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.tri3().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.tri3().tend();
                    ++_iter )
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

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.quad().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.quad().tend();
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

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.tri4().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.tri4().tend();
                    ++_iter )
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

            if (_imin < +0 || _imax>=_nmax)
            {
                _jlog.push (
    "**input error: INIT. cell. indexing is incorrect.\n") ;

                _errv = __invalid_indexing ;
            }
        }
        else
        {
            _jlog.  push (
    "**input error: couldn't determine format of INIT.\n") ;

            return __invalid_argument ;
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * ECHO-INIT: print summary of INIT data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type echo_init (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_init
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;

        if (_init._ndim == +2 &&
            _init._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(2) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnN1 = +0 ;
            iptr_type _nnE2 = +0 ;
            iptr_type _nnT3 = +0 , _nnQ4 =  +0 ;

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.node().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.node().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnN1 += +1 ;
            }

            _jlog.push(
                "  |COORD.| = "
                + std::to_string(_nnN1) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.edge().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.edge().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnE2 += +1 ;
            }

            _jlog.push(
                "  |EDGE-2| = "
                + std::to_string(_nnE2) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.tri3().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.tri3().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnT3 += +1 ;
            }

            _jlog.push(
                "  |TRIA-3| = "
                + std::to_string(_nnT3) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh.quad().head();
                      _iter != _init.
            _euclidean_mesh_2d._mesh.quad().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnQ4 += +1 ;
            }

            _jlog.push(
                "  |QUAD-4| = "
                + std::to_string(_nnQ4) + "\n");
        }
        else
        if (_init._ndim == +3 &&
            _init._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(3) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnN1 = +0 ;
            iptr_type _nnE2 = +0 ;
            iptr_type _nnT3 = +0 , _nnQ4 =  +0 ;
            iptr_type _nnT4 = +0 ;

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.node().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.node().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnN1 += +1 ;
            }

            _jlog.push(
                "  |COORD.| = "
                + std::to_string(_nnN1) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.edge().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.edge().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnE2 += +1 ;
            }

            _jlog.push(
                "  |EDGE-2| = "
                + std::to_string(_nnE2) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.tri3().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.tri3().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnT3 += +1 ;
            }

            _jlog.push(
                "  |TRIA-3| = "
                + std::to_string(_nnT3) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.quad().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.quad().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnQ4 += +1 ;
            }

            _jlog.push(
                "  |QUAD-4| = "
                + std::to_string(_nnQ4) + "\n");

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh.tri4().head();
                      _iter != _init.
            _euclidean_mesh_3d._mesh.tri4().tend();
                    ++_iter )
            {
            if (_iter->mark()>=+0) _nnT4 += +1 ;
            }

            _jlog.push(
                "  |TRIA-4| = "
                + std::to_string(_nnT4) + "\n");
        }

        _jlog.push("\n") ;

        return (  _errv) ;
    }


#   endif   //__INI_LOAD__



