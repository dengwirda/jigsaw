
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
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The 
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any 
     * way whatsoever.  This code is provided "as-is" to be 
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 10 June, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
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
            std::int32_t         _ndim ;
        public  :
    /*---------------------------------- construct reader */
        __normal_call init_reader (
            mesh_data*_isrc =  nullptr
            ) : _init(_isrc) {}
    /*---------------------------------- parse MSHID data */
        __normal_call void_type push_mshid (
            std::int32_t  _ftag ,
            jmsh_kind::enum_data _kind
            )
        {   
            this->_ftag = _ftag ;
            this->_kind = _kind ;
            this->
           _init->_kind = _kind ;
        }
    /*---------------------------------- parse NDIMS data */
        __normal_call void_type push_ndims (
            std::int32_t  _ndim
            )
        {   
            this->_ndim = _ndim ;
            this->
           _init->_ndim = _ndim ;
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
                mesh_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                
                _ndat.pval(2) = 
                        (real_type) +0. ;
                
                _ndat.fdim () = +0 ;
                
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
                
                _ndat.pval(3) = 
                        (real_type) +0. ;
                
                _ndat.fdim () = +0 ;
                
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
    /*---------------------------------- parse EDGE2 data */
        __normal_call void_type push_edge2 (
            std::int32_t  _ipos ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;
            __unreferenced(_itag) ;

            if (this->_ndim == +2 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                typename 
                mesh_data::euclidean_mesh_2d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                
                this->_init->
                   _euclidean_mesh_2d.
               _mesh.push_edge(_edat, false) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                typename 
                mesh_data::euclidean_mesh_3d
                    ::edge_type _edat ;
                _edat.node(0) = _node[0];
                _edat.node(1) = _node[1];
                
                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_edge(_edat, false) ;
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
            __unreferenced(_itag) ;

            if (this->_ndim == +2 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                typename 
                mesh_data::euclidean_mesh_2d
                    ::tria_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                
                this->_init->
                   _euclidean_mesh_2d.
               _mesh.push_tri3(_tdat, false) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                typename 
                mesh_data::euclidean_mesh_3d
                    ::face_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                
                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_tri3(_tdat, false) ;
            }
        }
    /*---------------------------------- parse TRIA4 data */
        __normal_call void_type push_tria4 (
            std::int32_t  _ipos ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;
            __unreferenced(_itag) ;

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
                mesh_data::euclidean_mesh_3d
                    ::tria_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                _tdat.node(3) = _node[3];
                
                this->_init->
                   _euclidean_mesh_3d.
               _mesh.push_tri4(_tdat, false) ;
            }
        }      
        } ;
    
    /*---------------------------------- parse INIT. file */
        iptr_type _errv  = __no_error ;
        
        try
        {
            jmsh_reader   _jmsh ;
            std::ifstream _file ; 
            _file. open(
            _jcfg._init_file, std::ifstream::in);

            if (_file.is_open() )
            {
                _jmsh.read_file (
                _file, init_reader(&_init));
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
            "**parse error: " + * _iter + "\n") ;
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
        
        if (_imsh._flags == 
                JIGSAW_EUCLIDEAN_MESH )
        {
            if (_imsh._vert2._size > 0)
            {
    /*--------------------------------- euclidean-mesh-2d */
            _init._kind 
                = jmsh_kind::euclidean_mesh ;
            _init._ndim = +2;
    
            for (auto _ipos = (iptr_type)+0 ;
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
                
                _ndat.pval(2) = 
                        (real_type) +0. ;
                
                _ndat.fdim () = +0;
                
                if (_imsh._vert2.
                        _data[_ipos]._itag < +0)
                    _ndat.feat () = 
                        mesh::user_feat ;
                else
                    _ndat.feat () = 
                        mesh::null_feat ;
                
                _init._euclidean_mesh_2d.
                    _mesh.push_node(_ndat,false);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                
                _init._euclidean_mesh_2d.
                    _mesh.push_edge(_edat,false);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _imsh._tria3._size ; 
                    ++_ipos )
            {
                typename 
                mesh_data::euclidean_mesh_2d
                    ::tria_type _tdat ;
                _tdat.node(0) = _imsh.
                    _tria3._data[_ipos]._node[0];
                _tdat.node(1) = _imsh.
                    _tria3._data[_ipos]._node[1];
                _tdat.node(2) = _imsh.
                    _tria3._data[_ipos]._node[2];
                
                _init._euclidean_mesh_2d.
                    _mesh.push_tri3(_tdat,false);
            }
    
            }
            else
            if (_imsh._vert3._size > 0)
            {
    /*--------------------------------- euclidean-mesh-3d */
            _init._kind 
                = jmsh_kind::euclidean_mesh ;
            _init._ndim = +3;
    
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _imsh._vert3._size ; 
                    ++_ipos )
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
                
                _ndat.pval(3) = 
                        (real_type) +0. ;
                
                _ndat.fdim () = +0;
                
                if (_imsh._vert3.
                        _data[_ipos]._itag < +0)
                    _ndat.feat () = 
                        mesh::user_feat ;
                else
                    _ndat.feat () = 
                        mesh::null_feat ;
                
                _init._euclidean_mesh_3d.
                    _mesh.push_node(_ndat,false);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                
                _init._euclidean_mesh_3d.
                    _mesh.push_edge(_edat,false);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _imsh._tria3._size ; 
                    ++_ipos )
            {
                typename 
                mesh_data::euclidean_mesh_3d
                    ::face_type _tdat ;
                _tdat.node(0) = _imsh.
                    _tria4._data[_ipos]._node[0];
                _tdat.node(1) = _imsh.
                    _tria4._data[_ipos]._node[1];
                _tdat.node(2) = _imsh.
                    _tria4._data[_ipos]._node[2];
                
                _init._euclidean_mesh_3d.
                    _mesh.push_tri3(_tdat,false);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _imsh._tria4._size ; 
                    ++_ipos )
            {
                typename 
                mesh_data::euclidean_mesh_3d
                    ::tria_type _tdat ;
                _tdat.node(0) = _imsh.
                    _tria4._data[_ipos]._node[0];
                _tdat.node(1) = _imsh.
                    _tria4._data[_ipos]._node[1];
                _tdat.node(2) = _imsh.
                    _tria4._data[_ipos]._node[2];
                _tdat.node(3) = _imsh.
                    _tria4._data[_ipos]._node[3];
                
                _init._euclidean_mesh_3d.
                    _mesh.push_tri4(_tdat,false);
            }
    
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
        iptr_type _errv  = __no_error ;

        std::string _path ;
        std::string _name ;
        std::string _fext ;
        file_part (
            _jcfg._init_file, 
                _path, _name, _fext ) ;

        if (_fext.find("msh") == +0 )
        {
        return init_from_jmsh (
                _jcfg, _jlog, _init ) ;
        }
        else
        {   
            _errv =__file_not_located ;
        }

        return ( _errv ) ;
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
            _euclidean_mesh_2d._mesh._set1.head() ;
                      _iter != _init.
            _euclidean_mesh_2d._mesh._set1.tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;
            
                _nmax += + 1 ;
            }

            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh._set2.head() ;
                      _iter != _init.
            _euclidean_mesh_2d._mesh._set2.tend() ;
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
            _euclidean_mesh_2d._mesh._set3.head() ;
                      _iter != _init.
            _euclidean_mesh_2d._mesh._set3.tend() ;
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

            if (_imin < +0 || _imax>=_nmax)
            {
                _jlog.push (
    "**input error: GEOM. tria. indexing is incorrect.\n") ;
        
                _errv = __invalid_argument ;
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
            _euclidean_mesh_3d._mesh._set1.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set1.tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;
            
                _nmax += + 1 ;
            }

            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh._set2.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set2.tend() ;
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
            _euclidean_mesh_3d._mesh._set3.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set3.tend() ;
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
            _euclidean_mesh_3d._mesh._set4.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set4.tend() ;
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
    "**input error: GEOM. tria. indexing is incorrect.\n") ;
        
                _errv = __invalid_argument ;
            }
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

        if (_init._ndim == +2 &&
            _init._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +2) ;
            
            _jlog.push("\n") ;
            
            iptr_type _num1 = +0 ;
            iptr_type _num2 = +0 ;
            iptr_type _num3 = +0 ;
            
            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh._set1.head() ;
                      _iter != _init.
            _euclidean_mesh_2d._mesh._set1.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }
            
            __dumpINTS("|COORD.|", _num1)
            
            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh._set2.head() ;
                      _iter != _init.
            _euclidean_mesh_2d._mesh._set2.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num2 += +1 ;
            }
            
            __dumpINTS("|EDGE-2|", _num2)
            
            for (auto _iter  = _init.
            _euclidean_mesh_2d._mesh._set3.head() ;
                      _iter != _init.
            _euclidean_mesh_2d._mesh._set3.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num3 += +1 ;
            }
            
            __dumpINTS("|TRIA-3|", _num3)   
        }
        else
        if (_init._ndim == +3 &&
            _init._kind ==
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
            iptr_type _num4 = +0 ;
            
            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh._set1.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set1.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }
            
            __dumpINTS("|COORD.|", _num1)
            
            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh._set2.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set2.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num2 += +1 ;
            }
            
            __dumpINTS("|EDGE-2|", _num2)
            
            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh._set3.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set3.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num3 += +1 ;
            }
            
            __dumpINTS("|TRIA-3|", _num3)  
            
            for (auto _iter  = _init.
            _euclidean_mesh_3d._mesh._set4.head() ;
                      _iter != _init.
            _euclidean_mesh_3d._mesh._set4.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num4 += +1 ;
            }
            
            __dumpINTS("|TRIA-4|", _num4)
        }

        _jlog.push("\n") ;

        #undef  __dumpINTS
        #undef  __dumpREAL

        return (  _errv) ;
    }
    
    
#   endif   //__INI_LOAD__
    
    
    
