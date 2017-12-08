
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
            
            containers::array<
                std::int32_t >   _pmap ;
            
            std::int32_t         _ftag ;
            jmsh_kind::
            enum_data            _kind ;
            std::int32_t         _ndim ;
            std::int32_t         _idim ;
            
        public  :
        __normal_call hfun_reader (
            hfun_data*_hsrc =  nullptr
            ) : _hfun(_hsrc) {}            
    /*-------------------------------- read MSHID section */
        __normal_call void_type push_mshid (
            std::int32_t  _ftag ,
            jmsh_kind::enum_data _kind
            )
        {   
            this->_ftag = _ftag ;
            this->_kind = _kind ;
            this->
           _hfun->_kind = _kind ;
        }
    /*-------------------------------- read NDIMS section */
        __normal_call void_type push_ndims (
            std::int32_t  _ndim
            )
        {   
            this->_ndim = _ndim ;
            this->
           _hfun->_ndim = _ndim ;
        }
    /*-------------------------------- read POINT section */
        __normal_call void_type push_point (
            std::int32_t  _ipos ,
            double       *_pval ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;
            __unreferenced(_itag) ;

            if (this->_ndim == +2 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                std::int32_t _npos ;
                typename 
                hfun_data::euclidean_mesh_2d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                
                _npos = this->_hfun->
                   _euclidean_mesh_2d.
                        _mesh.push_node(_ndat) ;
                        
                this->_pmap.
                    push_tail(_npos) ; 
            }           
            else
            if (this->_ndim == +3 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                std::int32_t _npos ;
                typename 
                hfun_data::euclidean_mesh_3d
                    ::node_type _ndat ;
                _ndat.pval(0) = _pval[0];
                _ndat.pval(1) = _pval[1];
                _ndat.pval(2) = _pval[2];
                
                _npos = this->_hfun->
                   _euclidean_mesh_3d.
                        _mesh.push_node(_ndat) ;
                        
                this->_pmap.
                    push_tail(_npos) ;
            }
            else
            if (this->_kind == 
                    jmsh_kind::ellipsoid_mesh)
            {
                //!! do things
            }
        }      
    /*-------------------------------- read COORD section */
        __normal_call void_type open_coord (
            std::int32_t _idim,
            std::int32_t
            )
        {   this->_idim =_idim;
        }
        __normal_call void_type push_coord (
            std::int32_t _irow,
            double       _ppos
            )
        {
            __unreferenced(_irow) ;

            if (this->_ndim == +2 &&
                this->_kind == 
                    jmsh_kind::euclidean_grid)
            {
                if (this->_idim == +1)
                {
                    this->_hfun->
                    _euclidean_grid_2d.
                        _xpos.push_tail(_ppos) ;
                }
                else
                if (this->_idim == +2)
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
                if (this->_idim == +1)
                {
                    this->_hfun->
                    _euclidean_grid_3d.
                        _xpos.push_tail(_ppos) ;
                }
                else
                if (this->_idim == +2)
                {
                    this->_hfun->
                    _euclidean_grid_3d.
                        _ypos.push_tail(_ppos) ;
                }
                else
                if (this->_idim == +3)
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
                if (this->_idim == +1)
                {
                    this->_hfun->
                    _ellipsoid_grid_3d.
                        _xpos.push_tail(_ppos) ;
                }
                else
                if (this->_idim == +2)
                {
                    this->_hfun->
                    _ellipsoid_grid_3d.
                        _ypos.push_tail(_ppos) ;
                }            
            }
        }           
    /*-------------------------------- read VALUE section */
        __normal_call void_type push_value (
            std::int32_t  _ipos ,
            double       *_vval
            )
        {
            if (this->_ndim == +2 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                if (_ipos  < 
               (std::int32_t)this->_pmap.count())
                {
                std::int32_t  _nmap ;
                _nmap = this->_pmap [_ipos];
                
                this->_hfun->
                _euclidean_mesh_2d._mesh.
                    _set1[_nmap].hval() = *_vval;
                }
            }
            else
            if (this->_ndim == +2 &&
                this->_kind == 
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_2d.
                        _hmat.push_tail(*_vval) ;
            }
            else
            if (this->_ndim == +3 &&
                this->_kind == 
                    jmsh_kind::euclidean_mesh)
            {
                if (_ipos < 
               (std::int32_t)this->_pmap.count())
                {
                std::int32_t  _nmap ;
                _nmap = this->_pmap [_ipos];
                
                this->_hfun->
                _euclidean_mesh_3d._mesh.
                    _set1[_nmap].hval() = *_vval;
                }
            }
            else
            if (this->_ndim == +3 &&
                this->_kind == 
                    jmsh_kind::euclidean_grid)
            {
                this->_hfun->
                    _euclidean_grid_3d.
                        _hmat.push_tail(*_vval) ;
            }
            else
            if (this->_kind == 
                    jmsh_kind::ellipsoid_grid)
            {
                this->_hfun->
                    _ellipsoid_grid_3d.
                        _hmat.push_tail(*_vval) ;
            }
        }
    /*-------------------------------- read TRIA3 section */
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
                hfun_data::euclidean_mesh_2d
                    ::tri3_type _tdat ;
                _tdat.node(0) = _node[0];
                _tdat.node(1) = _node[1];
                _tdat.node(2) = _node[2];
                
                this->_hfun->
                    _euclidean_mesh_2d.
                        _mesh.push_tri3(_tdat) ;
            }
            else
            if (this->_ndim == +3)
            {
                // 3-dimensional!!          
            }  
        }
    /*-------------------------------- read TRIA4 section */
        __normal_call void_type push_tria4 (
            std::int32_t  _ipos ,
            std::int32_t *_node ,
            std::int32_t  _itag
            )
        {
            __unreferenced(_ipos) ;
            __unreferenced(_itag) ;

            if (this->_ndim == +2)
            {
                // 2-dimensional!!            
            }
            else
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
                        _mesh.push_tri4(_tdat) ;
            }
        }
        } ;
    
    /*---------------------------------- parse HFUN. file */
        iptr_type _errv  = __no_error ;
        
        try
        {
            jmsh_reader   _jmsh ;
            std::ifstream _file ; 
            _file. open(
            _jcfg._hfun_file, std::ifstream::in);

            if (_file.is_open() )
            {
                _jmsh.read_file (
                _file, hfun_reader(&_hfun));
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
        iptr_type _errv  = __no_error ;

        std::string _path ;
        std::string _name ;
        std::string _fext ;
        file_part(_jcfg._hfun_file , 
            _path, _name, _fext);

        if (_fext.find("msh") == +0)
        {
        return hfun_from_jmsh (
                _jcfg, _jlog, _hfun)  ;
        }
        else
        {   
            _errv =__file_not_located ;
        }
        
        return ( _errv ) ;
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
    
            for (auto _ipos = (iptr_type)+0 ;
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
                _ndat.hval () = _hmsh.
                    _value._data[_ipos] ;
                
                _hfun._euclidean_mesh_2d.
                        _mesh.push_node(_ndat);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                        _mesh.push_tri3(_tdat);
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
    
            for (auto _ipos = (iptr_type)+0 ;
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
                _ndat.hval () = _hmsh.
                    _value._data[_ipos] ;
                
                _hfun._euclidean_mesh_3d.
                        _mesh.push_node(_ndat);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
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
                        _mesh.push_tri4(_tdat);
            }
    
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
    
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._xgrid._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d._xpos.
                    push_tail(_hmsh.
                        _xgrid._data[_ipos]);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._ygrid._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d._ypos.
                    push_tail(_hmsh.
                        _ygrid._data[_ipos]);
            }
        
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._value._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_2d._hmat.
                    push_tail(_hmsh.
                        _value._data[_ipos]);
            }
            
            }
            else
            if (_hmsh._zgrid._size!= 0)
            {
    /*--------------------------------- euclidean-grid-2d */
            _hfun._kind 
                = jmsh_kind::euclidean_grid ;
            _hfun._ndim = +3;
    
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._xgrid._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d._xpos.
                    push_tail(_hmsh.
                        _xgrid._data[_ipos]);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._ygrid._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d._ypos.
                    push_tail(_hmsh.
                        _ygrid._data[_ipos]);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._zgrid._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d._zpos.
                    push_tail(_hmsh.
                        _zgrid._data[_ipos]);
            }
        
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._value._size ; 
                    ++_ipos )
            {
                _hfun._euclidean_grid_3d._hmat.
                    push_tail(_hmsh.
                        _value._data[_ipos]);
            }
            
            }
        }
        else
        if (_hmsh._flags == 
                JIGSAW_ELLIPSOID_GRID )
        {
            if (_hmsh._zgrid._size== 0)
            {
    /*--------------------------------- euclidean-grid-2d */
            _hfun._kind 
                = jmsh_kind::ellipsoid_grid ;
            _hfun._ndim = +3;
    
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._xgrid._size ; 
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d._xpos.
                    push_tail(_hmsh.
                        _xgrid._data[_ipos]);
            }
            
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._ygrid._size ; 
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d._ypos.
                    push_tail(_hmsh.
                        _ygrid._data[_ipos]);
            }
        
            for (auto _ipos = (iptr_type)+0 ;
                _ipos != _hmsh._value._size ; 
                    ++_ipos )
            {
                _hfun._ellipsoid_grid_3d._hmat.
                    push_tail(_hmsh.
                        _value._data[_ipos]);
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
            _jcfg, _jlog, _hfun, _hmsh) ;
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
            #define __max2 (iptr_type) \
        _hfun._euclidean_mesh_2d._mesh._set1.count()
        
            real_type _hmin = 
            std::numeric_limits<real_type>::max() ;
            
            iptr_type _imin = 
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax = 
            std::numeric_limits<iptr_type>::min() ;

            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh._set1.head() ;
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh._set1.tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;
            
                _hmin = std::min(
                    _hmin, _iter->hval ()) ;
            }

            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh._set3.head() ;
                      _iter != _hfun.
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

            if (_hmin <= (real_type) +0.)
            {
                _jlog.push (
    "  **input error: HFUN. values must be strictly +ve.\n") ;
        
                _errv = __invalid_argument ;
            }

            if (_imin < +0 || _imax >= __max2)
            {
                _jlog.push (
    "  **input error: HFUN. tria. indexing is incorrect.\n") ;
        
                _errv = __invalid_argument ;
            }
        }
        else
        if (_hfun._ndim == +2 && 
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-2d */
            real_type _hmin = 
            std::numeric_limits<real_type>::infinity();
  
            for (auto _iter  = _hfun.
            _euclidean_grid_2d._hmat.head();
                      _iter != _hfun.
            _euclidean_grid_2d._hmat.tend();
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
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
 
            if (_hnum != _xnum * _ynum)
            {
                _jlog.push (
    "  **input error: HFUN. matrix incorrect dimensions.\n") ;
        
                _errv = __invalid_argument ;
            }
 
            if (_hmin <= (real_type) +0.)
            {
                _jlog.push (
    "  **input error: HFUN. values must be strictly +ve.\n") ;
        
                _errv = __invalid_argument ;
            }
 
            if (!_mono)
            {
                _jlog.push (
    "  **input error: grid must be monotonic increasing.\n") ;
        
                _errv = __invalid_argument ;
            }        
        }
        else
        if (_hfun._ndim == +3 && 
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            #define __max3 (iptr_type) \
        _hfun._euclidean_mesh_3d._mesh._set1.count()
        
            real_type _hmin = 
            std::numeric_limits<real_type>::max() ;
            
            iptr_type _imin = 
            std::numeric_limits<iptr_type>::max() ;
            iptr_type _imax = 
            std::numeric_limits<iptr_type>::min() ;

            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh._set1.head() ;
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh._set1.tend() ;
                    ++_iter  )
            {
                if (_iter->mark() < 0) continue ;
            
                _hmin = std::min(
                    _hmin, _iter->hval ()) ;
            }

            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh._set4.head() ;
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh._set4.tend() ;
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

            if (_hmin <= (real_type) +0.)
            {
                _jlog.push (
    "  **input error: HFUN. values must be strictly +ve.\n") ;
        
                _errv = __invalid_argument ;
            }

            if (_imin < +0 || _imax >= __max3)
            {
                _jlog.push (
    "  **input error: HFUN. tria. indexing is incorrect.\n") ;
        
                _errv = __invalid_argument ;
            }
        }
        else
        if (_hfun._ndim == +3 && 
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-3d */
            real_type _hmin = 
            std::numeric_limits<real_type>::infinity();
  
            for (auto _iter  = _hfun.
            _euclidean_grid_3d._hmat.head();
                      _iter != _hfun.
            _euclidean_grid_3d._hmat.tend();
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
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
 
            if (_hnum!=_xnum*_ynum*_znum)
            {
                _jlog.push (
    "  **input error: HFUN. matrix incorrect dimensions.\n") ;
        
                _errv = __invalid_argument ;
            }
 
            if (_hmin <= (real_type) +0.)
            {
                _jlog.push (
    "  **input error: HFUN. values must be strictly +ve.\n") ;
        
                _errv = __invalid_argument ;
            }
 
            if (!_mono)
            {
                _jlog.push (
    "  **input error: grid must be monotonic increasing.\n") ;
        
                _errv = __invalid_argument ;
            }        
        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_grid)
        {
    /*--------------------------------- ellipsoid-grid-3d */
            real_type _hmin = 
            std::numeric_limits<real_type>::infinity();
  
            for (auto _iter  = _hfun.
            _ellipsoid_grid_3d._hmat.head();
                      _iter != _hfun.
            _ellipsoid_grid_3d._hmat.tend();
                    ++_iter  )
            {
                _hmin = std::min(_hmin, *_iter) ;
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
            }
 
            auto _xnum = _hfun.
                _ellipsoid_grid_3d._xpos.count();
                
            auto _ynum = _hfun.
                _ellipsoid_grid_3d._ypos.count();
                
            auto _hnum = _hfun.
                _ellipsoid_grid_3d._hmat.count();
 
            if (_hnum != _xnum * _ynum)
            {
                _jlog.push (
    "  **input error: HFUN. matrix incorrect dimensions.\n") ;
        
                _errv = __invalid_argument ;
            }
 
            if (_hmin <= (real_type) +0.)
            {
                _jlog.push (
    "  **input error: HFUN. values must be strictly +ve.\n") ;
        
                _errv = __invalid_argument ;
            }
 
            if (!_mono)
            {
                _jlog.push (
    "  **input error: grid must be monotonic increasing.\n") ;
        
                _errv = __invalid_argument ;
            }        
        }
        
        #undef  __max2
        #undef  __max3

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
            << std::setprecision(+2)        \
                  << __var << "\n" ;        \
            _jlog.push(_sstr.str());

        iptr_type _errv  = __no_error ;
        
        __unreferenced(_jcfg) ;

        real_type _hmin = 
            +std::numeric_limits
                <real_type>::infinity() ;
                
        real_type _hmax = 
            -std::numeric_limits
                <real_type>::infinity() ;
        
        if (_hfun._ndim == +0)
        {
    /*--------------------------------- constant-value-kd */
            _jlog.push(
                "  CONSTANT-VALUE\n\n") ;
            
            __dumpREAL(".VAL(H).", 
                _hfun._constant_value_kd._hval)
            
        }
        else
        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +2) ;
            
            _jlog.push("\n") ;
            
            iptr_type _num1 = +0 ;
            iptr_type _num3 = +0 ;
            
            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh._set1.head() ;
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh._set1.tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;
                
                _hmin = std::min (
                    _hmin, _iter->hval()) ;
                _hmax = std::max (
                    _hmax, _iter->hval()) ;
            }
            
            __dumpREAL(".MIN(H).", _hmin)
            __dumpREAL(".MAX(H).", _hmax)
            
            _jlog.push("  \n") ;
            
            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh._set1.head() ;
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh._set1.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }
            
            __dumpINTS("|COORD.|", _num1)
            
            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh._set3.head() ;
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh._set3.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num3 += +1 ;
            }
            
            __dumpINTS("|TRIA-3|", _num3)
            
        }
        else
        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-2d */
            _jlog.push(
                "  EUCLIDEAN-GRID\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +2) ;
            
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
            
            __dumpREAL(".MIN(H).", _hmin)
            __dumpREAL(".MAX(H).", _hmax)
            
            _jlog.push("  \n") ;
            
            auto _xnum = _hfun.
                _euclidean_grid_2d._xpos.count();
                
            auto _ynum = _hfun.
                _euclidean_grid_2d._ypos.count();
            
            __dumpINTS("|XGRID.|", _xnum)
            __dumpINTS("|YGRID.|", _ynum)
            
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            _jlog.push(
                "  EUCLIDEAN-MESH\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +3) ;
            
            _jlog.push("\n") ;
            
            iptr_type _num1 = +0 ;
            iptr_type _num4 = +0 ;
            
            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh._set1.head() ;
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh._set1.tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;
                
                _hmin = std::min (
                    _hmin, _iter->hval()) ;
                _hmax = std::max (
                    _hmax, _iter->hval()) ;
            }
            
            __dumpREAL(".MIN(H).", _hmin)
            __dumpREAL(".MAX(H).", _hmax)
            
            _jlog.push("  \n") ;
            
            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh._set1.head() ;
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh._set1.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num1 += +1 ;
            }
            
            __dumpINTS("|COORD.|", _num1)
            
            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh._set4.head() ;
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh._set4.tend() ;
                    ++_iter )
            {
            if (_iter->mark()>=+0) _num4 += +1 ;
            }
            
            __dumpINTS("|TRIA-4|", _num4)
            
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-3d */
            _jlog.push(
                "  EUCLIDEAN-GRID\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +3) ;
            
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
            
            __dumpREAL(".MIN(H).", _hmin)
            __dumpREAL(".MAX(H).", _hmax)
            
            _jlog.push("  \n") ;
            
            auto _xnum = _hfun.
                _euclidean_grid_3d._xpos.count();
                
            auto _ynum = _hfun.
                _euclidean_grid_3d._ypos.count();
                
            auto _znum = _hfun.
                _euclidean_grid_3d._zpos.count();
            
            __dumpINTS("|XGRID.|", _xnum)
            __dumpINTS("|YGRID.|", _ynum)
            __dumpINTS("|ZGRID.|", _znum)
            
        }
        else
        if (_hfun._kind ==
             jmsh_kind::ellipsoid_grid)
        {
    /*--------------------------------- ellipsoid-grid-3d */
            _jlog.push(
                "  ELLIPSOID-GRID\n\n") ;
            
            __dumpINTS("|NDIMS.|",  +3) ;
            
            _jlog.push("\n") ;
            
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
            
            __dumpREAL(".MIN(H).", _hmin)
            __dumpREAL(".MAX(H).", _hmax)
            
            _jlog.push("  \n") ;
            
            auto _xnum = _hfun.
                _ellipsoid_grid_3d._xpos.count();
                
            auto _ynum = _hfun.
                _ellipsoid_grid_3d._ypos.count();
            
            __dumpINTS("|XGRID.|", _xnum)
            __dumpINTS("|YGRID.|", _ynum)
            
        }

        _jlog.push("\n") ;

        #undef  __dumpINTS
        #undef  __dumpREAL


        return (  _errv) ;
    }
    

#   endif   //__HFN_LOAD__



