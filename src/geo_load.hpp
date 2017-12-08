
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
     * Last updated: 16 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

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
            std::int32_t         _ndim ;
        public  :
    /*---------------------------------- construct reader */
        __normal_call geom_reader (
            geom_data*_gsrc =  nullptr
            ) : _geom(_gsrc) {}
    /*---------------------------------- parse MSHID data */
        __normal_call void_type push_mshid (
            std::int32_t  _ftag ,
            jmsh_kind::enum_data _kind
            )
        {   
            this->_ftag = _ftag ;
            this->_kind = _kind ;
            this->
           _geom->_kind = _kind ;
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
            std::int32_t  _ndim
            )
        {   
            this->_ndim = _ndim ;
            this->
           _geom->_ndim = _ndim ;
        }
    /*---------------------------------- parse POINT data */
        __normal_call void_type push_point (
            std::int32_t  _ipos ,
            double       *_pval ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;

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
                        _tria.push_node(_ndat) ;
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
                        _tria.push_node(_ndat) ;
            }
            else
            if (this->_kind == 
                    jmsh_kind::ellipsoid_mesh)
            {
                //!! do things here...
            }
        }
    /*---------------------------------- parse EDGE2 data */
        __normal_call void_type push_edge2 (
            std::int32_t  _ipos ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;

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
                        _tria.push_edge(_edat) ;
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
                        _tria.push_edge(_edat) ;
            }
            else
            if (this->_kind == 
                    jmsh_kind::ellipsoid_mesh)
            {
                //!! do things here...
            } 
        }
    /*---------------------------------- parse TRIA3 data */
        __normal_call void_type push_tria3 (
            std::int32_t  _ipos ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;

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
                        _tria.push_tri3(_tdat) ;
            }
            else
            if (this->_kind == 
                    jmsh_kind::ellipsoid_mesh)
            {
                // do nothing...
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
            _jcfg._geom_file, std::ifstream::in);

            if (_file.is_open() )
            {
                _jmsh.read_file (
                _file, geom_reader(&_geom));
            }
            else
            {           
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
            "  **parse error: " + *_iter + "\n");
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
        iptr_type _errv  = __no_error ;

        std::string _path ;
        std::string _name ;
        std::string _fext ;
        file_part(_jcfg._geom_file, 
            _path, _name, _fext);

        if (_fext.find("msh") == +0)
        {
        return geom_from_jmsh (
                _jcfg, _jlog, _geom)  ;
        }
        else
        {   
            _errv =__file_not_located ;
        }

        return ( _errv ) ;
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
        
        if (_gmsh._flags == 
                JIGSAW_EUCLIDEAN_MESH )
        {
            if (_gmsh._vert2._size > 0)
            {
    /*--------------------------------- euclidean-mesh-2d */
            _geom._kind 
                = jmsh_kind::euclidean_mesh ;
            _geom._ndim = +2;
    
            for (auto _ipos = (iptr_type)+0 ;
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
                        _tria.push_node(_ndat);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                        _tria.push_edge(_edat);
            }
            
            }
            else
            if (_gmsh._vert3._size > 0)
            {
    /*--------------------------------- euclidean-mesh-3d */
            _geom._kind 
                = jmsh_kind::euclidean_mesh ;
            _geom._ndim = +3;
            
            for (auto _ipos = (iptr_type)+0 ;
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
                        _tria.push_node(_ndat);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                        _tria.push_edge(_edat);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                        _tria.push_tri3(_tdat);
            }
     
            }        
        }
        else
        if (_gmsh._flags == 
                JIGSAW_ELLIPSOID_MESH )
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            if (_gmsh._radii._size==+3)
            {
            _geom._ellipsoid_mesh_3d.
                _radA = _gmsh._radii._data[0] ;
            _geom._ellipsoid_mesh_3d.
                _radB = _gmsh._radii._data[1] ;
            _geom._ellipsoid_mesh_3d.
                _radC = _gmsh._radii._data[2] ;
                
            //!! blah, also need coastlines...
                
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
            _jcfg, _jlog, _geom, _gmsh) ;
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
            #define __max2 (iptr_type) \
        _geom._euclidean_mesh_2d._tria._set1.count()
        
            iptr_type _imin = 
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax = 
            std::numeric_limits<iptr_type>::min() ;

            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria._set2.head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria._set2.tend() ;
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

            if (_imin < +0 || _imax >= __max2)
            {
                _jlog.push (
            "  **input error: invalid indexing\n");
        
                _errv = __invalid_argument ;
            }
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            #define __max3 (iptr_type) \
        _geom._euclidean_mesh_3d._tria._set1.count()
        
            iptr_type _imin = 
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax = 
            std::numeric_limits<iptr_type>::min() ;

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria._set2.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria._set2.tend() ;
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

            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria._set3.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria._set3.tend() ;
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

            if (_imin < +0 || _imax >= __max3)
            {
                _jlog.push (
            "  **input error: invalid indexing\n");
        
                _errv = __invalid_argument ;
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
            "  **input error: bad radius value\n");
        
                _errv = __invalid_argument ;
            }
        }

        #undef  __max2
        #undef  __max3

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
        std::stringstream  _sstr;

    /*---------------------------------- push "int_" data */
        #define __dumpINTS(__tag,__var)     \
            _sstr.str("");                  \
            _sstr.clear();                  \
            _sstr << "  " __tag " = "       \
                  << __var << "\n" ;        \
            _jlog.push(_sstr.str());
            
    /*---------------------------------- push "real" data */
        #define __dumpREAL(__tag,__var)     \
            _sstr.str("");                  \
            _sstr.clear();                  \
            _sstr << "  " __tag " = "       \
                  << std::scientific        \
                  << std::setprecision(2)   \
                  << __var                  \
                  << " \n" ;                \
            _jlog.push(_sstr.str()) ;       \

        iptr_type _errv  = __no_error ;
        
        __unreferenced(_jcfg) ;

        if (_geom._ndim == +2 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +2) ;
            
            _jlog.push("\n") ;
            
            iptr_type _num1 = +0 ;
            iptr_type _num2 = +0 ;
            
            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria._set1.head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria._set1.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }
            
            __dumpINTS("|COORD.|", _num1)
            
            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria._set2.head() ;
                      _iter != _geom.
            _euclidean_mesh_2d._tria._set2.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num2 += +1 ;
            }
            
            __dumpINTS("|EDGE-2|", _num2)   
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +3) ;
            
            _jlog.push("\n") ;
            
            iptr_type _num1 = +0 ;
            iptr_type _num2 = +0 ;
            iptr_type _num3 = +0 ;
            
            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria._set1.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria._set1.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }
            
            __dumpINTS("|COORD.|", _num1)
            
            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria._set2.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria._set2.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num2 += +1 ;
            }
            
            __dumpINTS("|EDGE-2|", _num2)
            
            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria._set3.head() ;
                      _iter != _geom.
            _euclidean_mesh_3d._tria._set3.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num3 += +1 ;
            }
            
            __dumpINTS("|TRIA-3|", _num3)           
        }
        else
        if (_geom._kind ==
             jmsh_kind::ellipsoid_mesh)
        {
    /*--------------------------------- ellipsoid-mesh-3d */
            _jlog.push(
                "  ELLIPSOID-MESH\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +3) ;
            
            _jlog.push("\n") ;
            
            __dumpREAL("|1-RAD.|", 
                _geom._ellipsoid_mesh_3d._radA) ;
                
            __dumpREAL("|2-RAD.|", 
                _geom._ellipsoid_mesh_3d._radB) ;
            
            __dumpREAL("|3-RAD.|", 
                _geom._ellipsoid_mesh_3d._radC) ;                       
        }

        _jlog.push("\n") ;

        #undef  __dumpINTS
        #undef  __dumpREAL

        return (  _errv) ;
    }


#   endif   //__GEO_LOAD__



