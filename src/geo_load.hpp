
    /*
    --------------------------------------------------------
     * GEO-LOAD: parse *.MSH file into GEOM data.
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
     * Last updated: 10 Jul., 2021
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

#   ifndef __GEO_LOAD__
#   define __GEO_LOAD__

    /*
    --------------------------------------------------------
     * GEOM-FROM-JMSH: read *.JMSH file into GEOM data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type geom_from_jmsh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom
        )
    {
        class geom_reader:
            public jmsh_reader_base
        {
        public  :
            geom_data           *_geom ;
            std::int32_t         _ftag ;
            jmsh_kind::
            enum_data            _kind ;
            std:: size_t         _ndim ;
        public  :
    /*---------------------------------- construct reader */
        __normal_call geom_reader (
            geom_data*_gsrc =  nullptr
            ) : _geom(_gsrc) {}
    /*---------------------------------- parse MSHID data */
        __normal_call void_type push_mshid (
            std::int32_t  _FTAG ,
            jmsh_kind::enum_data _KIND
            )
        {
            this->_ftag = _FTAG ;
            this->_kind = _KIND ;
            this->
           _geom->_kind = _KIND ;
        }
    /*---------------------------------- parse RADII data */
        __normal_call void_type push_radii (
            double       *_erad
            )
        {
            this->_geom->_ellipsoid_mesh_3d.
                _radA = _erad[ 0] ;
            this->_geom->_ellipsoid_mesh_3d.
                _radB = _erad[ 1] ;
            this->_geom->_ellipsoid_mesh_3d.
                _radC = _erad[ 2] ;
        }
    /*---------------------------------- parse NDIMS data */
        __normal_call void_type push_ndims (
            std:: size_t  _NDIM
            )
        {
            this->_ndim = _NDIM ;
            this->
           _geom->_ndim = _NDIM ;
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
                geom_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.itag () = _itag ;

                this->_geom->
                   _euclidean_mesh_2d.
               _tria.push_node(_ndat, false) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) = _pval[2];
                _ndat.itag () = _itag ;

                this->_geom->
                   _euclidean_mesh_3d.
               _tria.push_node(_ndat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                typename
                geom_data::ellipsoid_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) = + 0.0 ;
                _ndat.itag () = _itag ;

                this->_geom->
                   _ellipsoid_mesh_3d.
               _mesh.push_node(_ndat, false) ;
            }
        }
    /*---------------------------------- parse SEEDS data */
        __normal_call void_type push_seeds (
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
                geom_data::euclidean_mesh_2d
                    ::seed_type _sdat ;
                _sdat.mark () = +   0 ;
                _sdat.self () = +   1 ;
                _sdat.pval(0) = _pval[0];
                _sdat.pval(1) = _pval[1];
                _sdat.itag () = _itag ;

                this->_geom->_euclidean_mesh_2d.
                   _seed.push_tail(_sdat) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::seed_type _sdat ;
                _sdat.mark () = +   0 ;
                _sdat.self () = +   1 ;
                _sdat.pval(0) = _pval[0];
                _sdat.pval(1) = _pval[1];
                _sdat.pval(2) = _pval[2];
                _sdat.itag () = _itag ;

                this->_geom->_euclidean_mesh_3d.
                   _seed.push_tail(_sdat) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                typename
                geom_data::ellipsoid_mesh_3d
                    ::seed_type _sdat ;
                _sdat.mark () = +   0 ;
                _sdat.self () = +   1 ;
                _sdat.pval(0) = _pval[0];
                _sdat.pval(1) = _pval[1];
                _sdat.pval(2) = + 0.0 ;
                _sdat.itag () = _itag ;

                this->_geom->_ellipsoid_mesh_3d.
                   _seed.push_tail(_sdat) ;
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
                geom_data::euclidean_mesh_2d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                _edat.itag () = _itag ;

                this->_geom->
                   _euclidean_mesh_2d.
               _tria.push_edge(_edat, false) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                _edat.itag () = _itag ;

                this->_geom->
                   _euclidean_mesh_3d.
               _tria.push_edge(_edat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                typename
                geom_data::ellipsoid_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                _edat.itag () = _itag ;

                this->_geom->
                   _ellipsoid_mesh_3d.
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
                // do nothing...
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                _tdat.itag () = _itag ;

                this->_geom->
                   _euclidean_mesh_3d.
               _tria.push_tri3(_tdat, false) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                // do nothing...
            }
        }
    /*---------------------------------- parse BOUND data */
        __normal_call void_type push_bound (
            std:: size_t/*_ipos*/ ,
            std::int32_t  _itag ,
            std::int32_t  _inum ,
            std::int32_t  _KIND
            )
        {
            if (this->_ndim == +2 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                geom_data::euclidean_mesh_2d
                    ::part_data _pdat ;
                _pdat.itag () = _itag ;
                _pdat.indx () = _inum ;
                _pdat.kind () = _KIND ;

                this->_geom->_euclidean_mesh_2d.
                   _part.push  (_pdat) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind ==
                    jmsh_kind::euclidean_mesh)
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::part_data _pdat ;
                _pdat.itag () = _itag ;
                _pdat.indx () = _inum ;
                _pdat.kind () = _KIND ;

                this->_geom->_euclidean_mesh_3d.
                   _part.push  (_pdat) ;
            }
            else
            if (this->_kind ==
                    jmsh_kind::ellipsoid_mesh)
            {
                //!! do things here...
            }
        }
        } ;

    /*---------------------------------- parse GEOM. file */
        iptr_type _errv  = __no_error ;

        try
        {
            jmsh_reader   _jmsh ;
            std::ifstream _file ;
            _file. open(
            _jcfg._geom_file, std::ifstream::in) ;

            if (_file.is_open() )
            {
                _jmsh.read_file (
                _file, geom_reader(&_geom));
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
     * GEOM-FROM-MSHT: read MSH_t data into GEOM data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type geom_from_msht (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom ,
        jigsaw_msh_t const&_gmsh
        )
    {
        iptr_type _errv  = __no_error ;

        try
        {

        __unreferenced (_jlog) ;
        __unreferenced (_jcfg) ;

        if (_gmsh._flags ==
                JIGSAW_EUCLIDEAN_MESH )
        {
            if (_gmsh._vert2._size > 0)
            {
    /*--------------------------------- euclidean-mesh-2d */
            _geom._kind
                = jmsh_kind::euclidean_mesh ;
            _geom._ndim = +2;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._vert2._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _gmsh.
                    _vert2._data[_ipos]._ppos[0];
                _ndat.pval(1) = _gmsh.
                    _vert2._data[_ipos]._ppos[1];
                _ndat.itag () = _gmsh.
                    _vert2._data[_ipos]._itag ;

                _geom._euclidean_mesh_2d.
                _tria.push_node(_ndat , false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._seed2._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_2d
                    ::seed_type _sdat ;
                _sdat.pval(0) = _gmsh.
                    _seed2._data[_ipos]._ppos[0];
                _sdat.pval(1) = _gmsh.
                    _seed2._data[_ipos]._ppos[1];
                _sdat.itag () = _gmsh.
                    _seed2._data[_ipos]._itag ;

                _sdat.mark () = +   0 ;
                _sdat.self () = +   1 ;

                _geom._euclidean_mesh_2d.
                      _seed.push_tail(_sdat);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._edge2._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_2d
                    ::edge_type _edat ;
                _edat.node(0) = _gmsh.
                    _edge2._data[_ipos]._node[0];
                _edat.node(1) = _gmsh.
                    _edge2._data[_ipos]._node[1];
                _edat.itag () = _gmsh.
                    _edge2._data[_ipos]._itag ;

                _geom._euclidean_mesh_2d.
                _tria.push_edge(_edat , false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._bound._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_2d
                    ::part_data _pdat ;
                _pdat.itag () = _gmsh.
                    _bound._data[_ipos]._itag ;
                _pdat.indx () = _gmsh.
                    _bound._data[_ipos]._indx ;
                _pdat.kind () = _gmsh.
                    _bound._data[_ipos]._kind ;

                _geom._euclidean_mesh_2d.
                      _part.push(_pdat) ;
            }

            }
            else
            if (_gmsh._vert3._size > 0)
            {
    /*--------------------------------- euclidean-mesh-3d */
            _geom._kind
                = jmsh_kind::euclidean_mesh ;
            _geom._ndim = +3;

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._vert3._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _gmsh.
                    _vert3._data[_ipos]._ppos[0];
                _ndat.pval(1) = _gmsh.
                    _vert3._data[_ipos]._ppos[1];
                _ndat.pval(2) = _gmsh.
                    _vert3._data[_ipos]._ppos[2];
                _ndat.itag () = _gmsh.
                    _vert3._data[_ipos]._itag ;

                _geom._euclidean_mesh_3d.
                _tria.push_node(_ndat , false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._seed3._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::seed_type _sdat ;
                _sdat.pval(0) = _gmsh.
                    _seed3._data[_ipos]._ppos[0];
                _sdat.pval(1) = _gmsh.
                    _seed3._data[_ipos]._ppos[1];
                _sdat.pval(2) = _gmsh.
                    _seed3._data[_ipos]._ppos[2];
                _sdat.itag () = _gmsh.
                    _seed3._data[_ipos]._itag ;

                _sdat.mark () = +   0 ;
                _sdat.self () = +   1 ;

                _geom._euclidean_mesh_3d.
                      _seed.push_tail(_sdat);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._edge2._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _gmsh.
                    _edge2._data[_ipos]._node[0];
                _edat.node(1) = _gmsh.
                    _edge2._data[_ipos]._node[1];
                _edat.itag () = _gmsh.
                    _edge2._data[_ipos]._itag ;

                _geom._euclidean_mesh_3d.
                _tria.push_edge(_edat , false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._tria3._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _gmsh.
                    _tria3._data[_ipos]._node[0];
                _tdat.node(1) = _gmsh.
                    _tria3._data[_ipos]._node[1];
                _tdat.node(2) = _gmsh.
                    _tria3._data[_ipos]._node[2];
                _tdat.itag () = _gmsh.
                    _tria3._data[_ipos]._itag ;

                _geom._euclidean_mesh_3d.
                _tria.push_tri3(_tdat , false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._bound._size ;
                    ++_ipos )
            {
                typename
                geom_data::euclidean_mesh_3d
                    ::part_data _pdat ;
                _pdat.itag () = _gmsh.
                    _bound._data[_ipos]._itag ;
                _pdat.indx () = _gmsh.
                    _bound._data[_ipos]._indx ;
                _pdat.kind () = _gmsh.
                    _bound._data[_ipos]._kind ;

                _geom._euclidean_mesh_3d.
                      _part.push(_pdat) ;
            }

            }
        }
        else
        if (_gmsh._flags ==
                JIGSAW_ELLIPSOID_MESH )
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            _geom._kind
                = jmsh_kind::ellipsoid_mesh ;
            _geom._ndim = +2;

            if (_gmsh._radii._size==+3)
            {
            _geom._ellipsoid_mesh_3d.
                _radA = _gmsh._radii._data[0] ;
            _geom._ellipsoid_mesh_3d.
                _radB = _gmsh._radii._data[1] ;
            _geom._ellipsoid_mesh_3d.
                _radC = _gmsh._radii._data[2] ;
            }
            else
            if (_gmsh._radii._size==+1)
            {
            _geom._ellipsoid_mesh_3d.
                _radA = _gmsh._radii._data[0] ;
            _geom._ellipsoid_mesh_3d.
                _radB = _gmsh._radii._data[0] ;
            _geom._ellipsoid_mesh_3d.
                _radC = _gmsh._radii._data[0] ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._vert2._size ;
                    ++_ipos )
            {
                typename
                geom_data::ellipsoid_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _gmsh.
                    _vert2._data[_ipos]._ppos[0];
                _ndat.pval(1) = _gmsh.
                    _vert2._data[_ipos]._ppos[1];
                _ndat.pval(2) = + 0.0 ;
                _ndat.itag () = _gmsh.
                    _vert2._data[_ipos]._itag ;

                _geom._ellipsoid_mesh_3d.
                _mesh.push_node(_ndat , false) ;
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._seed2._size ;
                    ++_ipos )
            {
                typename
                geom_data::ellipsoid_mesh_3d
                    ::seed_type _sdat ;
                _sdat.pval(0) = _gmsh.
                    _seed2._data[_ipos]._ppos[0];
                _sdat.pval(1) = _gmsh.
                    _seed2._data[_ipos]._ppos[1];
                _sdat.pval(2) = + 0.0 ;
                _sdat.itag () = _gmsh.
                    _seed2._data[_ipos]._itag ;

                _sdat.mark () = +   0 ;
                _sdat.self () = +   1 ;

                _geom._ellipsoid_mesh_3d.
                      _seed.push_tail(_sdat);
            }

            for (auto _ipos = (size_t) +0 ;
                _ipos != _gmsh._edge2._size ;
                    ++_ipos )
            {
                typename
                geom_data::ellipsoid_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _gmsh.
                    _edge2._data[_ipos]._node[0];
                _edat.node(1) = _gmsh.
                    _edge2._data[_ipos]._node[1];
                _edat.itag () = _gmsh.
                    _edge2._data[_ipos]._itag ;

                _geom._ellipsoid_mesh_3d.
                _mesh.push_edge(_edat , false) ;
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
     * READ-GEOM: read geometry input file.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type read_geom (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom
        )
    {
        return geom_from_jmsh (
                _jcfg, _jlog, _geom ) ;
    }

    /*
    --------------------------------------------------------
     * COPY-GEOM: read geometry input data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type copy_geom (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom ,
        jigsaw_msh_t const&_gmsh
        )
    {
        return geom_from_msht (
           _jcfg, _jlog, _geom, _gmsh);
    }

    /*
    --------------------------------------------------------
     * TEST-GEOM: test GEOM data validity.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type test_geom (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;

        if (_geom._ndim == +2 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            iptr_type _imin =
            std::numeric_limits<iptr_type>::max () ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min () ;

            iptr_type _nnPT = +0 ;
            iptr_type _nnE2 = +0 ;

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria.node().head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria.node().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

                _nnPT += +1  ;
            }

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria.edge().head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria.edge().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;

                _nnE2 +=    +1 ;

                if (_imin < +0 ||
                        _imax >= _nnPT)
                {
                    _errv = __invalid_indexing ;
                }
            }

            if (_errv != __no_error)
            {
                _jlog.  push (
    "**input error: GEOM. EDGE2 indexing is incorrect.\n") ;

                return _errv ;
            }

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._part._lptr.head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._part._lptr.tend() ;
                    ++_iter  )
            {
                if (*_iter == nullptr) continue;

                for (auto _pptr  = *_iter ;
                    _pptr != nullptr;
                    _pptr  = _pptr->_next )
                {

                if (_pptr->_data.kind()
                         == JIGSAW_EDGE2_TAG )
                {
                if (_pptr->_data.indx() <  +0 ||
                    _pptr->_data.indx() >= _nnE2
                    )
                {
                    _errv = __invalid_indexing ;
                }
                }
                else
                {
                    _errv = __invalid_indexing ;
                }

                }
            }

            if (_errv != __no_error)
            {
                _jlog.  push (
    "**input error: GEOM. BOUND indexing is incorrect.\n") ;

                return _errv ;
            }
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            iptr_type _imin =
            std::numeric_limits<iptr_type>::max () ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min () ;

            iptr_type _nnPT = +0 ;
            iptr_type _nnE2 = +0 ;
            iptr_type _nnT3 = +0 ;

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.node().head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria.node().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

                _nnPT += +1  ;
            }

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.edge().head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria.edge().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;

                _nnE2 +=    +1 ;

                if (_imin < +0 ||
                        _imax >= _nnPT)
                {
                    _errv = __invalid_indexing ;
                }
            }

            if (_errv != __no_error)
            {
                _jlog.  push (
    "**input error: GEOM. EDGE2 indexing is incorrect.\n") ;

                return _errv ;
            }

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.tri3().head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria.tri3().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

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

                _nnT3 +=    +1 ;

                if (_imin < +0 ||
                        _imax >= _nnPT)
                {
                    _errv = __invalid_indexing ;
                }
            }

            if (_errv != __no_error)
            {
                _jlog.  push (
    "**input error: GEOM. TRIA3 indexing is incorrect.\n") ;

                return _errv ;
            }

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._part._lptr.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._part._lptr.tend() ;
                    ++_iter  )
            {
                if (*_iter == nullptr) continue;

                for (auto _pptr  = *_iter ;
                    _pptr != nullptr;
                    _pptr  = _pptr->_next )
                {

                if (_pptr->_data.kind()
                         == JIGSAW_TRIA3_TAG )
                {
                if (_pptr->_data.indx() <  +0 ||
                    _pptr->_data.indx() >= _nnT3
                    )
                {
                    _errv = __invalid_indexing ;
                }
                }
                else
                {
                    _errv = __invalid_indexing ;
                }

                }
            }

            if (_errv != __no_error)
            {
                _jlog.  push (
    "**input error: GEOM. BOUND indexing is incorrect.\n") ;

                return _errv ;
            }
        }
        else
        if (_geom._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            if (_geom._ellipsoid_mesh_3d.
                _radA <= (real_type)  +0. ||
                _geom._ellipsoid_mesh_3d.
                _radB <= (real_type)  +0. ||
                _geom._ellipsoid_mesh_3d.
                _radC <= (real_type)  +0. )
            {
                _jlog.push (
    "**input error: GEOM. RADII entries are incorrect.\n") ;

                _errv = __invalid_argument ;
            }

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
            std::numeric_limits<iptr_type>::max () ;
            iptr_type _imax =
            std::numeric_limits<iptr_type>::min () ;

            iptr_type _nnPT = +0 ;
            iptr_type _nnE2 = +0 ;

            for (auto _iter  = _geom.
            _ellipsoid_mesh_3d._mesh.node().head() ;
                      _iter != _geom.
            _ellipsoid_mesh_3d._mesh.node().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

                _xmin = std::min(
                    _xmin, _iter->pval(0)) ;
                _xmax = std::max(
                    _xmax, _iter->pval(0)) ;

                _ymin = std::min(
                    _ymin, _iter->pval(1)) ;
                _ymax = std::max(
                    _ymax, _iter->pval(1)) ;

                _nnPT += +1  ;
            }

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

            for (auto _iter  = _geom.
            _ellipsoid_mesh_3d._mesh.edge().head() ;
                      _iter != _geom.
            _ellipsoid_mesh_3d._mesh.edge().tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue;

                _imin = std::min(
                    _imin, _iter->node(0)) ;
                _imin = std::min(
                    _imin, _iter->node(1)) ;
                _imax = std::max(
                    _imax, _iter->node(0)) ;
                _imax = std::max(
                    _imax, _iter->node(1)) ;

                _nnE2 +=    +1 ;

                if (_imin < +0 ||
                        _imax >= _nnPT)
                {
                    _errv = __invalid_indexing ;
                }
            }

            if (_errv != __no_error)
            {
                _jlog.  push (
    "**input error: GEOM. EDGE2 indexing is incorrect.\n") ;

                return _errv ;
            }
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * ECHO-GEOM: print summary of GEOM data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type echo_geom (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;

        if (_geom._ndim == +2 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(2) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnPT = +0 ;
            iptr_type _nnE2 = +0 ;

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria.node().head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria.node().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnPT += +1 ;
            }

            _jlog.push(
                "  |COORD.| = " 
                + std::to_string(_nnPT) + "\n");

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria.edge().head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria.edge().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnE2 += +1 ;
            }

            _jlog.push(
                "  |EDGE-2| = " 
                + std::to_string(_nnE2) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnSD = +0 ;

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._seed.head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._seed.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnSD += +1 ;
            }

            _jlog.push(
                "  |SEEDS.| = " 
                + std::to_string(_nnSD) + "\n");

            _jlog.push("\n") ;

            _jlog.push("  |BOUND.| = " +
                std::to_string(_geom.
                _euclidean_mesh_2d._ptag.count()) + 
                " (" + 
                std::to_string(_geom.
                _euclidean_mesh_2d._part.count()) +
                ") " + "\n") ;
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(3) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnPT = +0 ;
            iptr_type _nnE2 = +0 ;
            iptr_type _nnT3 = +0 ;

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.node().head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria.node().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnPT += +1 ;
            }

            _jlog.push(
                "  |COORD.| = " 
                + std::to_string(_nnPT) + "\n");

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.edge().head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria.edge().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnE2 += +1 ;
            }

            _jlog.push(
                "  |EDGE-2| = " 
                + std::to_string(_nnE2) + "\n");

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.tri3().head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria.tri3().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnT3 += +1 ;
            }

            _jlog.push(
                "  |TRIA-3| = " 
                + std::to_string(_nnT3) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnSD = +0 ;

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._seed.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._seed.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnSD += +1 ;
            }

            _jlog.push(
                "  |SEEDS.| = " 
                + std::to_string(_nnSD) + "\n");

            _jlog.push("\n") ;

            _jlog.push("  |BOUND.| = " +
                std::to_string(_geom.
                _euclidean_mesh_3d._ptag.count()) + 
                " (" + 
                std::to_string(_geom.
                _euclidean_mesh_3d._part.count()) +
                ") " + "\n") ;
        }
        else
        if (_geom._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            _jlog.push(
                "  ELLIPSOID-MESH\n\n") ;

            _jlog.push(
                "  |NDIMS.| = " + std::to_string(3) + "\n");

            _jlog.push("\n") ;

            _jlog.push("  |1-RAD.| = " +
                to_string_prec (
            _geom._ellipsoid_mesh_3d._radA, 2) + "\n");

            _jlog.push("  |2-RAD.| = " +
                to_string_prec (
            _geom._ellipsoid_mesh_3d._radB, 2) + "\n");

            _jlog.push("  |3-RAD.| = " +
                to_string_prec (
            _geom._ellipsoid_mesh_3d._radC, 2) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnPT = +0 ;
            iptr_type _nnE2 = +0 ;

            for (auto _iter  = _geom.
            _ellipsoid_mesh_3d._mesh.node().head() ;
                      _iter != _geom.
            _ellipsoid_mesh_3d._mesh.node().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnPT += +1 ;
            }

            _jlog.push(
                "  |COORD.| = " 
                + std::to_string(_nnPT) + "\n");

            for (auto _iter  = _geom.
            _ellipsoid_mesh_3d._mesh.edge().head() ;
                      _iter != _geom.
            _ellipsoid_mesh_3d._mesh.edge().tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnE2 += +1 ;
            }

            _jlog.push(
                "  |EDGE-2| = " 
                + std::to_string(_nnE2) + "\n");

            _jlog.push("\n") ;

            iptr_type _nnSD = +0 ;

            for (auto _iter  = _geom.
            _ellipsoid_mesh_3d._seed.head() ;
                      _iter != _geom.
            _ellipsoid_mesh_3d._seed.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=0 ) _nnSD += +1 ;
            }

            _jlog.push(
                "  |SEEDS.| = " 
                + std::to_string(_nnSD) + "\n");
        }

        _jlog.push("\n") ;

        #undef  __dumpINTS
        #undef  __dumpREAL

        return (  _errv) ;
    }


#   endif   //__GEO_LOAD__



