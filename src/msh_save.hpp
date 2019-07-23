
    /*
    --------------------------------------------------------
     * MSH-SAVE: parse MESH data into *.MSH / MSH_t.
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
     * Last updated: 20 July, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_SAVE__
#   define __MSH_SAVE__

    /*
    --------------------------------------------------------
     * SAVE-RDEL: save *.MSH output file.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_rdel (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        rdel_data &_rdel
        )
    {
        iptr_type _errv  = __no_error  ;

        __unreferenced(_jlog) ;

        try
        {
            containers::array<iptr_type> _nmap;
                
            std::ofstream _file ;

            std::string _path ;
            std::string _name ;
            std::string _fext ; 
            file_part(
                _jcfg._mesh_file, 
                    _path, _name, _fext);

            _file.open(
                _jcfg._mesh_file, 
                    std::ofstream::out );
                     
            if (_file.is_open())
            {
            if (_rdel._ndim == +2 &&
                _rdel._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 2-dim. mesh */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=2;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=2 \n" ;

                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_2d._tria._nset.count() , 
                        containers::tight_alloc , -1);

                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->
                          _data._node[0]] = 1 ;
                    _nmap[_item->
                          _data._node[1]] = 1 ;
                    }
                }

                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->
                          _data._node[0]] = 1 ;
                    _nmap[_item->
                          _data._node[1]] = 1 ;
                    _nmap[_item->
                          _data._node[2]] = 1 ;
                    }
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_2d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _last << "\n" ;
                
                iptr_type _npos  = +0 ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_jcfg._iter_opts.dual() )
                {
                if (_rdel._euclidean_rdel_2d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POWER data */   
                _file << "POWER=" 
                      << _last << ";1" << "\n" ;
                
                iptr_type _npos  = +0 ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << (real_type)0. << "\n" ;
                    }
                }
                }
                }
                
                if (_rdel._euclidean_rdel_2d.
                        _eset.count() > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                _file << "EDGE2=" << 
                    _rdel._euclidean_rdel_2d.
                        _eset.count() << "\n" ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _file << 
                    _nmap[_item->_data._node[0]]
                          << ";" << 
                    _nmap[_item->_data._node[1]]
                          << ";" << 
                    _item->_data._part << "\n" ;             
                    }
                }
                }
     
                if (_rdel._euclidean_rdel_2d.
                        _tset.count() > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << 
                    _rdel._euclidean_rdel_2d.
                        _tset.count() << "\n" ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _file << 
                    _nmap[_item->_data._node[0]]
                          << ";" << 
                    _nmap[_item->_data._node[1]]
                          << ";" << 
                    _nmap[_item->_data._node[2]]
                          << ";" << 
                    _item->_data._part << "\n" ;
                    }
                }
                }
     
            }
            else
            if (_rdel._ndim == +3 &&
                _rdel._kind == 
                jmsh_kind::euclidean_mesh) 
            {
            /*-------------------------- save 3-dim. mesh */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=2;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=3 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);
                    
            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_3d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->
                          _data._node[0]] = 1 ;
                    _nmap[_item->
                          _data._node[1]] = 1 ;                   
                    }
                }

                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._fset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._fset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->
                          _data._node[0]] = 1 ;
                    _nmap[_item->
                          _data._node[1]] = 1 ;
                    _nmap[_item->
                          _data._node[2]] = 1 ;                   
                    }
                }

                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->
                          _data._node[0]] = 1 ;
                    _nmap[_item->
                          _data._node[1]] = 1 ;
                    _nmap[_item->
                          _data._node[2]] = 1 ; 
                    _nmap[_item->
                          _data._node[3]] = 1 ;                   
                    }
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_3d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _last << "\n" ;
                
                iptr_type _npos  = +0 ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << _iter->pval(2) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_jcfg._iter_opts.dual() )
                {
                if (_rdel._euclidean_rdel_3d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POWER data */                
                _file << "POWER=" 
                      << _last << ";1" << "\n" ;
                
                iptr_type _npos  = +0 ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << (real_type)0. << "\n" ;
                    }
                }
                }
                }
                
                if (_rdel._euclidean_rdel_3d.
                        _eset.count() > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                _file << "EDGE2=" << 
                    _rdel._euclidean_rdel_3d.
                        _eset.count() << "\n" ;
                        
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _file << 
                    _nmap[_item->_data._node[0]]
                          << ";" << 
                    _nmap[_item->_data._node[1]]
                          << ";" << 
                    _item->_data._part << "\n" ;
                    }
                }
                }
                
                if (_rdel._euclidean_rdel_3d.
                        _fset.count() > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << 
                    _rdel._euclidean_rdel_3d.
                        _fset.count() << "\n" ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._fset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._fset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _file << 
                    _nmap[_item->_data._node[0]]
                          << ";" << 
                    _nmap[_item->_data._node[1]]
                          << ";" << 
                    _nmap[_item->_data._node[2]]
                          << ";" << 
                    _item->_data._part << "\n" ;
                    }
                }
                }
                
                if (_rdel._euclidean_rdel_3d.
                        _tset.count() > +0)
                {
            /*-------------------------- write TRIA4 data */ 
                _file << "TRIA4=" << 
                    _rdel._euclidean_rdel_3d.
                        _tset.count() << "\n" ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _file << 
                    _nmap[_item->_data._node[0]]
                          << ";" << 
                    _nmap[_item->_data._node[1]]
                          << ";" << 
                    _nmap[_item->_data._node[2]]
                          << ";" << 
                    _nmap[_item->_data._node[3]]
                          << ";" << 
                    _item->_data._part << "\n" ;
                    }
                }
                }
                
            }

            }
            else
            {
                _errv = __file_not_located ;
            }
            
            _file.close();

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }
    
    /*
    --------------------------------------------------------
     * SAVE-RDEL: save MSH_t output data.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_rdel (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        rdel_data &_rdel ,
        jigsaw_msh_t &_mmsh
        )
    {
        iptr_type _errv  = __no_error  ;

        try
        {
            containers::array<iptr_type> _nmap;
        
            __unreferenced (_jlog) ;
            __unreferenced (_jcfg) ;
    
            if (_rdel._ndim == +2 &&    // save 2-dim. mesh
                _rdel._kind ==
                jmsh_kind::euclidean_mesh)
            {
                _mmsh._flags = JIGSAW_EUCLIDEAN_MESH ;
            
            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_2d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->_data._node[0]]=1 ;
                    _nmap[_item->_data._node[1]]=1 ;                   
                    }
                }

                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->_data._node[0]]=1 ;
                    _nmap[_item->_data._node[1]]=1 ;
                    _nmap[_item->_data._node[2]]=1 ;                   
                    }
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_2d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                jigsaw_alloc_vert2 (
                     &_mmsh._vert2, _last) ;
                     
                jigsaw_alloc_reals (
                     &_mmsh._power, _last) ;
                
                iptr_type _npos  = +0 ;
                iptr_type _nout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0 )
                    {
                    _mmsh._vert2._data[_nout].
                        _ppos[0] = _iter->pval(0) ;
                    _mmsh._vert2._data[_nout].
                        _ppos[1] = _iter->pval(1) ;
                    
                    _mmsh._vert2.
                        _data[_nout]._itag = 0 ;
                        
                    _mmsh._power.
                    _data[_nout] = (real_type)+0. ;
                        
                    _nout = _nout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_2d._tria.
                    clear(containers::tight_alloc);

                }
                
                if (_rdel._euclidean_rdel_2d.
                        _eset.count() > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                jigsaw_alloc_edge2 (
                     &_mmsh._edge2,_rdel.
                _euclidean_rdel_2d._eset.count()) ;
                
                iptr_type _eout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _mmsh._edge2._data[_eout].
                        _node[0] = 
                        _nmap[_item->_data._node[0]] ;
                    _mmsh._edge2._data[_eout].
                        _node[1] = 
                        _nmap[_item->_data._node[1]] ;
                    
                    _mmsh._edge2._data[_eout].
                        _itag = _item->_data._part ;
                        
                    _eout = _eout + 1 ;               
                    }
                }

                _rdel._euclidean_rdel_2d._eset.
                    clear(containers::tight_alloc) ;

                }
     
                if (_rdel._euclidean_rdel_2d.
                        _tset.count() > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria3 (
                     &_mmsh._tria3,_rdel.
                _euclidean_rdel_2d._tset.count()) ;
                
                iptr_type _tout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _mmsh._tria3._data[_tout].
                        _node[0] = 
                        _nmap[_item->_data._node[0]] ;
                    _mmsh._tria3._data[_tout].
                        _node[1] = 
                        _nmap[_item->_data._node[1]] ;
                    _mmsh._tria3._data[_tout].
                        _node[2] = 
                        _nmap[_item->_data._node[2]] ;
                    
                    _mmsh._tria3._data[_tout].
                        _itag = _item->_data._part ;
                        
                    _tout = _tout + 1 ;               
                    }
                }

                _rdel._euclidean_rdel_2d._tset.
                    clear(containers::tight_alloc) ;

                }
     
            }
            else
            if (_rdel._ndim == +3 &&    // save 3-dim. mesh
                _rdel._kind == 
                jmsh_kind::euclidean_mesh) 
            {      
                _mmsh._flags = JIGSAW_EUCLIDEAN_MESH ;
            
            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_3d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->_data._node[0]]=1 ;
                    _nmap[_item->_data._node[1]]=1 ;                   
                    }
                }

                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._fset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._fset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->_data._node[0]]=1 ;
                    _nmap[_item->_data._node[1]]=1 ;
                    _nmap[_item->_data._node[2]]=1 ;                   
                    }
                }

                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _nmap[_item->_data._node[0]]=1 ;
                    _nmap[_item->_data._node[1]]=1 ;
                    _nmap[_item->_data._node[2]]=1 ; 
                    _nmap[_item->_data._node[3]]=1 ;                   
                    }
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_3d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                jigsaw_alloc_vert3 (
                     &_mmsh._vert3 , _last);
                     
                jigsaw_alloc_reals (
                     &_mmsh._power , _last);
                
                iptr_type _npos  = +0 ;
                iptr_type _nout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0 )
                    {
                    _mmsh._vert3._data[_nout].
                        _ppos[0] = _iter->pval(0) ;
                    _mmsh._vert3._data[_nout].
                        _ppos[1] = _iter->pval(1) ;
                    _mmsh._vert3._data[_nout].
                        _ppos[2] = _iter->pval(2) ;
                    
                    _mmsh._vert3.
                        _data[_nout]._itag = 0 ;
                        
                    _mmsh._power.
                    _data[_nout] = (real_type)+0. ;
                        
                    _nout = _nout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_3d._tria.
                    clear(containers::tight_alloc);

                }

                if (_rdel._euclidean_rdel_3d.
                        _eset.count() > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                jigsaw_alloc_edge2 (
                     &_mmsh._edge2,_rdel.
                _euclidean_rdel_3d._eset.count()) ;
   
                iptr_type _eout  = +0 ;                     
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._eset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._eset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _mmsh._edge2._data[_eout].
                        _node[0] = 
                        _nmap[_item->_data._node[0]] ;
                    _mmsh._edge2._data[_eout].
                        _node[1] = 
                        _nmap[_item->_data._node[1]] ;
                    
                    _mmsh._edge2._data[_eout].
                        _itag = _item->_data._part ;
                        
                    _eout = _eout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_3d._eset.
                    clear(containers::tight_alloc) ;

                }
                
                if (_rdel._euclidean_rdel_3d.
                        _fset.count() > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria3 (
                     &_mmsh._tria3,_rdel.
                _euclidean_rdel_3d._fset.count()) ;
 
                iptr_type _fout  = +0 ;               
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._fset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._fset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _mmsh._tria3._data[_fout].
                        _node[0] = 
                        _nmap[_item->_data._node[0]] ;
                    _mmsh._tria3._data[_fout].
                        _node[1] = 
                        _nmap[_item->_data._node[1]] ;
                    _mmsh._tria3._data[_fout].
                        _node[2] = 
                        _nmap[_item->_data._node[2]] ;
                    
                    _mmsh._tria3._data[_fout].
                        _itag = _item->_data._part ;
                        
                    _fout = _fout + 1 ;               
                    }
                }

                _rdel._euclidean_rdel_3d._fset.
                    clear(containers::tight_alloc) ;

                }
                
                if (_rdel._euclidean_rdel_3d.
                        _tset.count() > +0)
                {
            /*-------------------------- write TRIA4 data */ 
                jigsaw_alloc_tria4 (
                     &_mmsh._tria4,_rdel.
                _euclidean_rdel_3d._tset.count()) ;
                
                iptr_type _tout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tset._lptr.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tset._lptr.tend();
                        ++_iter  )
                {
                    if ( *_iter == nullptr) continue ;
                    
                    for (auto _item  = *_iter ;
                        _item != nullptr;
                        _item  = _item->_next )
                    {
                    _mmsh._tria4._data[_tout].
                        _node[0] = 
                        _nmap[_item->_data._node[0]] ;
                    _mmsh._tria4._data[_tout].
                        _node[1] = 
                        _nmap[_item->_data._node[1]] ;
                    _mmsh._tria4._data[_tout].
                        _node[2] = 
                        _nmap[_item->_data._node[2]] ;
                    _mmsh._tria4._data[_tout].
                        _node[3] = 
                        _nmap[_item->_data._node[3]] ;
                    
                    _mmsh._tria4._data[_tout].
                        _itag = _item->_data._part ;
                        
                    _tout = _tout + 1 ;               
                    }
                }

                _rdel._euclidean_rdel_3d._tset.
                    clear(containers::tight_alloc) ;

                }
                
            }
            
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }
    
    /*
    --------------------------------------------------------
     * SAVE-TRIA: save *.MSH output file.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_tria (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        rdel_data &_rdel
        )
    {
        iptr_type _errv  = __no_error  ;

        __unreferenced(_jlog) ;

        try
        {
            containers::array<iptr_type> _nmap;
                
            std::ofstream _file ;

            std::string _path ;
            std::string _name ;
            std::string _fext ; 
            file_part(
                _jcfg._tria_file, 
                    _path, _name, _fext);

            _file.open(
                _jcfg._tria_file, 
                    std::ofstream::out );
                     
            if (_file.is_open())
            {
            if (_rdel._ndim == +2 &&
                _rdel._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 2-dim. mesh */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=2;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=2 \n" ;

                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_2d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _ntri = +0;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;
                    
                    _ntri += +1 ;
                    
                    _nmap[_iter->node(0)]=1;
                    _nmap[_iter->node(1)]=1;
                    _nmap[_iter->node(2)]=1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_2d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _last << "\n" ;
                
                iptr_type _npos  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_rdel._euclidean_rdel_2d.
                        _tria._tset.count() > 0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << _ntri << "\n" ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;
                
                    _file << 
                    _nmap[_iter->node(0)] << ";" << 
                    _nmap[_iter->node(1)] << ";" <<
                    _nmap[_iter->node(2)] << ";" <<
                       +0 << "\n" ;
                }
                }
     
            }
            else
            if (_rdel._ndim == +3 &&
                _rdel._kind == 
                jmsh_kind::euclidean_mesh) 
            {
            /*-------------------------- save 3-dim. mesh */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=2;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=3 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);
                    
            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_3d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _ntri  = +0;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;
                
                    _ntri += +1 ;
                    
                    _nmap[_iter->node(0)]=1;
                    _nmap[_iter->node(1)]=1;
                    _nmap[_iter->node(2)]=1;
                    _nmap[_iter->node(3)]=1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_3d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _last << "\n" ;
                
                iptr_type _npos  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << _iter->pval(2) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_rdel._euclidean_rdel_3d.
                        _tria._tset.count() > 0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA4=" << _ntri << "\n" ;
                
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;
                
                    _file << 
                    _nmap[_iter->node(0)] << ";" << 
                    _nmap[_iter->node(1)] << ";" <<
                    _nmap[_iter->node(2)] << ";" <<
                    _nmap[_iter->node(3)] << ";" <<
                       +0 << "\n" ;
                }
                }
                
            }

            }
            else
            {
                _errv = __file_not_located ;
            }
            
            _file.close();

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    /*
    --------------------------------------------------------
     * SAVE-TRIA: save MSH_t output data.
    --------------------------------------------------------
     */
     
    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_tria (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        rdel_data &_rdel ,
        jigsaw_msh_t &_mmsh
        )
    {
        iptr_type _errv  = __no_error  ;

        try
        {
            containers::array<iptr_type> _nmap;
        
            __unreferenced (_jlog) ;
            __unreferenced (_jcfg) ;
    
            if (_rdel._ndim == +2 &&    // save 2-dim. mesh
                _rdel._kind ==
                jmsh_kind::euclidean_mesh)
            {
                _mmsh._flags = JIGSAW_EUCLIDEAN_MESH ;
            
            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_2d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _ntri  = +0;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;
                
                    _ntri +=  +1 ;
                    
                    _nmap[_iter->node(0)]=1;
                    _nmap[_iter->node(1)]=1;
                    _nmap[_iter->node(2)]=1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_2d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                jigsaw_alloc_vert2 (
                     &_mmsh._vert2, _last) ;
                     
                jigsaw_alloc_reals (
                     &_mmsh._power, _last) ;
                
                iptr_type _npos  = +0 ;
                iptr_type _nout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0 )
                    {
                    _mmsh._vert2._data[_nout].
                        _ppos[0] = _iter->pval(0) ;
                    _mmsh._vert2._data[_nout].
                        _ppos[1] = _iter->pval(1) ;
                    
                    _mmsh._vert2.
                        _data[_nout]._itag = 0 ;
                        
                    _mmsh._power.
                    _data[_nout] = (real_type)+0. ;
                        
                    _nout = _nout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_2d.
                    _tria._nset.
                    clear(containers::tight_alloc);

                }
     
                if (_rdel._euclidean_rdel_2d.
                        _tria._tset.count() > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria3 (
                     &_mmsh._tria3, _ntri) ;
                
                iptr_type _tout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() >= 0)
                    {
                    _mmsh._tria3._data[_tout].
                        _node[0] = 
                        _nmap[_iter->node(0)];
                    _mmsh._tria3._data[_tout].
                        _node[1] = 
                        _nmap[_iter->node(1)];
                    _mmsh._tria3._data[_tout].
                        _node[2] = 
                        _nmap[_iter->node(2)];

                    _mmsh._tria3.
                        _data[_tout]._itag = 0 ;
                     
                    _tout = _tout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_2d.
                    _tria._tset.
                    clear(containers::tight_alloc) ;

                }
     
            }
            else
            if (_rdel._ndim == +3 &&    // save 3-dim. mesh
                _rdel._kind == 
                jmsh_kind::euclidean_mesh) 
            {      
                _mmsh._flags = JIGSAW_EUCLIDEAN_MESH ;
            
            /*------------ index mapping for active nodes */
                _nmap.set_count(_rdel.
                    _euclidean_rdel_3d._tria._nset.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _ntri  = +0;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;
                
                    _ntri +=  +1 ;
                    
                    _nmap[_iter->node(0)]=1;
                    _nmap[_iter->node(1)]=1;
                    _nmap[_iter->node(2)]=1;
                    _nmap[_iter->node(3)]=1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _last ++ ;
                    }
                }

                if (_rdel._euclidean_rdel_3d.
                        _tria._nset.count() > 0)
                {
            /*-------------------------- write POINT data */                
                jigsaw_alloc_vert3 (
                     &_mmsh._vert3 , _last);
                     
                jigsaw_alloc_reals (
                     &_mmsh._power , _last);
                
                iptr_type _npos  = +0 ;
                iptr_type _nout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._nset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._nset.tend();
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0 )
                    {
                    _mmsh._vert3._data[_nout].
                        _ppos[0] = _iter->pval(0) ;
                    _mmsh._vert3._data[_nout].
                        _ppos[1] = _iter->pval(1) ;
                    _mmsh._vert3._data[_nout].
                        _ppos[2] = _iter->pval(2) ;
                    
                    _mmsh._vert3.
                        _data[_nout]._itag = 0 ;
                        
                    _mmsh._power.
                    _data[_nout] = (real_type)+0. ;
                        
                    _nout = _nout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_3d.
                    _tria._nset.
                    clear(containers::tight_alloc);

                }

                if (_rdel._euclidean_rdel_3d.
                        _tria._tset.count() > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria4 (
                     &_mmsh._tria4, _ntri) ;
                
                iptr_type _tout  = +0 ;
                for (auto _iter  = _rdel.
                _euclidean_rdel_3d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_3d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() >= 0)
                    {
                    _mmsh._tria4._data[_tout].
                        _node[0] = 
                        _nmap[_iter->node(0)];
                    _mmsh._tria4._data[_tout].
                        _node[1] = 
                        _nmap[_iter->node(1)];
                    _mmsh._tria4._data[_tout].
                        _node[2] = 
                        _nmap[_iter->node(2)];
                    _mmsh._tria4._data[_tout].
                        _node[3] = 
                        _nmap[_iter->node(3)];

                    _mmsh._tria4.
                        _data[_tout]._itag = 0 ;
                     
                    _tout = _tout + 1 ;
                    }
                }

                _rdel._euclidean_rdel_3d.
                    _tria._tset.
                    clear(containers::tight_alloc) ;

                }

            }
            
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    /*
    --------------------------------------------------------
     * SAVE-MESH: save *.MSH output file.
    --------------------------------------------------------
     */
    
    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_mesh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_mesh
        )
    {
        iptr_type _errv  = __no_error  ;

        __unreferenced(_jlog) ;

        try
        {
            containers::array<iptr_type> _nmap;
                
            std::ofstream _file ;

            std::string _path ;
            std::string _name ;
            std::string _fext ; 
            file_part(
                _jcfg._mesh_file, 
                    _path, _name, _fext);

            _file.open(
                _jcfg._mesh_file, 
                    std::ofstream::out );
                     
            if (_file.is_open())
            {
            if (_mesh._ndim == +2 &&
                _mesh._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 2-dim. mesh */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=2 \n" ;

                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*------------ index mapping for active nodes */
                _nmap.set_count(_mesh.
                    _euclidean_mesh_2d._mesh._set1.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _nnum = +0 ;
                iptr_type _enum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;                   
                    _enum +=  +1 ;
                    }
                }

                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;
                    _nmap[_iter->node(2)] = +1 ;                   
                    _tnum +=  +1 ;
                    }
                }

                for (auto _iter  = _nmap.head();
                          _iter != _nmap.tend();
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _nnum ++ ;
                    }
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _nnum << "\n" ;
                
                iptr_type _npos  = +0;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set1.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_nnum > +0)
                {
            /*-------------------------- write POWER data */                
                if (_jcfg._iter_opts.dual() )
                {
                _file << "POWER=" 
                      << _nnum << ";1" << "\n" ;
                
                iptr_type _npos  = +0;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set1.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(2) << "\n" ;
                    }
                }
                }
                }
                
                if (_enum > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                _file << "EDGE2=" << _enum << "\n" ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << 
                    _nmap[_iter->node(0)] << ";" 
                          << 
                    _nmap[_iter->node(1)] << ";"
                          << 
                    _iter->itag() << "\n"  ;
                    }
                }
                }
     
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << _tnum << "\n" ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << 
                    _nmap[_iter->node(0)] << ";" 
                          << 
                    _nmap[_iter->node(1)] << ";"
                          << 
                    _nmap[_iter->node(2)] << ";"
                          << 
                    _iter->itag() << "\n"  ;
                    }
                }
                }
     
            }
            else
            if (_mesh._ndim == +3 &&
                _mesh._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 3-dim. mesh */
                _file << "# " << _name << ".msh" 
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=3 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*------------ index mapping for active nodes */
                _nmap.set_count(_mesh.
                    _euclidean_mesh_3d._mesh._set1.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _nnum = +0 ;
                iptr_type _enum = +0 ;
                iptr_type _fnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;                   
                    _enum +=  +1 ;
                    }
                }

                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;
                    _nmap[_iter->node(2)] = +1 ;                   
                    _fnum +=  +1 ;
                    }
                }
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set4.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set4.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;
                    _nmap[_iter->node(2)] = +1 ;
                    _nmap[_iter->node(3)] = +1 ;                   
                    _tnum +=  +1 ;
                    }
                }

                for (auto _iter  = _nmap.head();
                          _iter != _nmap.tend();
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _nnum ++ ;
                    }
                }                  
             
                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _nnum << "\n" ;
                
                iptr_type _npos  = +0;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set1.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << _iter->pval(2) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_nnum > +0)
                {
            /*-------------------------- write POWER data */                
                if (_jcfg._iter_opts.dual() )
                {
                _file << "POWER=" 
                      << _nnum << ";1" << "\n" ;
                
                iptr_type _npos  = +0;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set1.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0)
                    {
                    _file << _iter->pval(3) << "\n" ;
                    }
                }
                }
                }
                
                if (_enum > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                _file << "EDGE2=" << _enum << "\n" ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << 
                    _nmap[_iter->node(0)] << ";" 
                          << 
                    _nmap[_iter->node(1)] << ";"
                          << 
                    _iter->itag() << "\n"  ;
                    }
                }
                }
                
                if (_fnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << _fnum << "\n" ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << 
                    _nmap[_iter->node(0)] << ";" 
                          << 
                    _nmap[_iter->node(1)] << ";"
                          << 
                    _nmap[_iter->node(2)] << ";"
                          << 
                    _iter->itag() << "\n"  ;
                    }
                }
                }
                
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA4=" << _tnum << "\n" ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set4.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set4.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << 
                    _nmap[_iter->node(0)] << ";" 
                          << 
                    _nmap[_iter->node(1)] << ";"
                          << 
                    _nmap[_iter->node(2)] << ";"
                          << 
                    _nmap[_iter->node(3)] << ";"
                          << 
                    _iter->itag() << "\n"  ;
                    }
                }
                }
            
            }

            }
            else
            {
                _errv = __file_not_located  ;
            }
            
            _file.close();

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    /*
    --------------------------------------------------------
     * SAVE-MESH: save MSH_t output data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_mesh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        mesh_data &_mesh ,
        jigsaw_msh_t &_mmsh
        )
    {
        iptr_type _errv  = __no_error  ;

        try
        {
            containers::array<iptr_type> _nmap ;
        
            __unreferenced (_jcfg) ;
            __unreferenced (_jlog) ;
     
            if (_mesh._ndim == +2 &&    // save 2-dim. mesh
                _mesh._kind ==
                jmsh_kind::euclidean_mesh)
            {
                _mmsh._flags = JIGSAW_EUCLIDEAN_MESH ;
    
            /*------------ index mapping for active nodes */
                _nmap.set_count(_mesh.
                    _euclidean_mesh_2d._mesh._set1.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _nnum = +0 ;
                iptr_type _enum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;                   
                    _enum +=  +1 ;
                    }
                }

                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;
                    _nmap[_iter->node(2)] = +1 ;                   
                    _tnum +=  +1 ;
                    }
                }

                for (auto _iter  = _nmap.head();
                          _iter != _nmap.tend();
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _nnum ++ ;
                    }
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                jigsaw_alloc_vert2 (
                     &_mmsh._vert2, _nnum) ;

                jigsaw_alloc_reals (
                     &_mmsh._power, _nnum) ;

                iptr_type _npos  = +0 ;
                iptr_type _nout  = +0 ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set1.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0 )
                    {
                    _mmsh._vert2._data[_nout].
                        _ppos[0] = _iter->pval(0) ;
                    _mmsh._vert2._data[_nout].
                        _ppos[1] = _iter->pval(1) ;
                    
                    _mmsh._vert2.
                        _data[_nout]._itag = 0 ;
                        
                    _mmsh._power.
                    _data[_nout] = _iter->pval(2) ;
                        
                    _nout = _nout + 1 ;
                    }
                }

                _mesh._euclidean_mesh_2d._mesh._set1.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_2d._mesh._adj1.
                    clear(containers::tight_alloc);

                }
                
                if (_enum > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                jigsaw_alloc_edge2 (
                     &_mmsh._edge2, _enum) ;
                
                iptr_type _eout  = +0 ;
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _mmsh._edge2.
                        _data[_eout]._node[0] 
                            = _nmap[_iter->node(0)] ;
                    _mmsh._edge2.
                        _data[_eout]._node[1] 
                            = _nmap[_iter->node(1)] ;
                            
                    _mmsh._edge2.
                        _data[_eout].
                            _itag = _iter->itag() ; 
                        
                    _eout = _eout + 1 ;                    
                    }
                }

                _mesh._euclidean_mesh_2d._mesh._set2.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_2d._mesh._adj2.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_2d._mesh._map2.
                    clear(containers::tight_alloc);

                }
     
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria3 (
                     &_mmsh._tria3, _tnum) ;
                
                iptr_type _tout  = +0 ;
                for (auto _iter  = _mesh.
                _euclidean_mesh_2d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_2d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _mmsh._tria3.
                        _data[_tout]._node[0] 
                            = _nmap[_iter->node(0)] ;
                    _mmsh._tria3.
                        _data[_tout]._node[1] 
                            = _nmap[_iter->node(1)] ;
                    _mmsh._tria3.
                        _data[_tout]._node[2] 
                            = _nmap[_iter->node(2)] ;
                            
                    _mmsh._tria3.
                        _data[_tout].
                            _itag = _iter->itag() ; 
                        
                    _tout = _tout + 1 ;
                    }
                }

                _mesh._euclidean_mesh_2d._mesh._set3.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_2d._mesh._adj3.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_2d._mesh._map3.
                    clear(containers::tight_alloc);

                }
     
            }
            else
            if (_mesh._ndim == +3 &&    // save 3-dim. mesh
                _mesh._kind ==
                jmsh_kind::euclidean_mesh)
            {
                _mmsh._flags = JIGSAW_EUCLIDEAN_MESH ;
            
            /*------------ index mapping for active nodes */
                _nmap.set_count(_mesh.
                    _euclidean_mesh_3d._mesh._set1.count() , 
                        containers::tight_alloc, -1) ;

                iptr_type _nnum = +0 ;
                iptr_type _enum = +0 ;
                iptr_type _fnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;                   
                    _enum +=  +1 ;
                    }
                }

                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;
                    _nmap[_iter->node(2)] = +1 ;                   
                    _fnum +=  +1 ;
                    }
                }
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set4.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set4.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                    _nmap[_iter->node(0)] = +1 ;
                    _nmap[_iter->node(1)] = +1 ;
                    _nmap[_iter->node(2)] = +1 ;
                    _nmap[_iter->node(3)] = +1 ;                   
                    _tnum +=  +1 ;
                    }
                }

                for (auto _iter  = _nmap.head();
                          _iter != _nmap.tend();
                        ++_iter  )
                {
                    if ( *_iter >= +0)
                    {
                         *_iter = _nnum ++ ;
                    }
                }
                
                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                jigsaw_alloc_vert3 (
                     &_mmsh._vert3, _nnum) ;
                     
                jigsaw_alloc_reals (
                     &_mmsh._power, _nnum) ;

                iptr_type _npos  = +0 ;
                iptr_type _nout  = +0 ;
                
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set1.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0 &&
                        _nmap[_npos ] >= 0 )
                    {
                    _mmsh._vert3._data[_nout].
                        _ppos[0] = _iter->pval(0) ;
                    _mmsh._vert3._data[_nout].
                        _ppos[1] = _iter->pval(1) ;
                    _mmsh._vert3._data[_nout].
                        _ppos[2] = _iter->pval(2) ;
                    
                    _mmsh._vert3.
                        _data[_nout]._itag = 0 ;
                        
                    _mmsh._power.
                    _data[_nout] = _iter->pval(3) ;
                        
                    _nout = _nout + 1 ;
                    }
                }

                _mesh._euclidean_mesh_3d._mesh._set1.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._adj1.
                    clear(containers::tight_alloc);

                }
                
                if (_enum > +0)
                {
            /*-------------------------- write EDGE2 data */ 
                jigsaw_alloc_edge2 (
                     &_mmsh._edge2, _enum) ;
                
                iptr_type _eout  = +0 ;
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set2.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set2.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _mmsh._edge2.
                        _data[_eout]._node[0] 
                            = _nmap[_iter->node(0)] ;
                    _mmsh._edge2.
                        _data[_eout]._node[1] 
                            = _nmap[_iter->node(1)] ;
                            
                    _mmsh._edge2.
                        _data[_eout].
                            _itag = _iter->itag() ;  
                        
                    _eout = _eout + 1 ;                    
                    }
                }

                _mesh._euclidean_mesh_3d._mesh._set2.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._adj2.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._map2.
                    clear(containers::tight_alloc);

                }
                
                if (_fnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria3 (
                     &_mmsh._tria3, _fnum) ;
                
                iptr_type _fout  = +0 ;
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set3.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _mmsh._tria3.
                        _data[_fout]._node[0] 
                            = _nmap[_iter->node(0)] ;
                    _mmsh._tria3.
                        _data[_fout]._node[1] 
                            = _nmap[_iter->node(1)] ;
                    _mmsh._tria3.
                        _data[_fout]._node[2] 
                            = _nmap[_iter->node(2)] ;
                            
                    _mmsh._tria3.
                        _data[_fout].
                            _itag = _iter->itag() ; 
                        
                    _fout = _fout + 1 ;
                    }
                }

                _mesh._euclidean_mesh_3d._mesh._set3.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._adj3.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._map3.
                    clear(containers::tight_alloc);

                }
                
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                jigsaw_alloc_tria4 (
                     &_mmsh._tria4, _tnum) ;
                
                iptr_type _tout  = +0 ;
                for (auto _iter  = _mesh.
                _euclidean_mesh_3d._mesh._set4.head() ;
                          _iter != _mesh.
                _euclidean_mesh_3d._mesh._set4.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _mmsh._tria4.
                        _data[_tout]._node[0] 
                            = _nmap[_iter->node(0)] ;
                    _mmsh._tria4.
                        _data[_tout]._node[1] 
                            = _nmap[_iter->node(1)] ;
                    _mmsh._tria4.
                        _data[_tout]._node[2] 
                            = _nmap[_iter->node(2)] ;
                    _mmsh._tria4.
                        _data[_tout]._node[3] 
                            = _nmap[_iter->node(3)] ;
                            
                    _mmsh._tria4.
                        _data[_tout].
                            _itag = _iter->itag() ; 
                        
                    _tout = _tout + 1 ;
                    }
                }

                _mesh._euclidean_mesh_3d._mesh._set4.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._adj4.
                    clear(containers::tight_alloc);

                _mesh._euclidean_mesh_3d._mesh._map4.
                    clear(containers::tight_alloc);

                }
 
            }
               
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    /*
    --------------------------------------------------------
     * SAVE-HFUN: save *.MSH output file.
    --------------------------------------------------------
     */
    
    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_ffun
        )
    {
        iptr_type _errv  = __no_error  ;

        __unreferenced(_jlog) ;

        try
        {
            std::ofstream _file ;

            std::string _path ;
            std::string _name ;
            std::string _fext ; 
            file_part(
                _jcfg._mesh_file, 
                    _path, _name, _fext);

            _file.open(
                _jcfg._mesh_file, 
                    std::ofstream::out );
                     
            if (_file.is_open())
            {
            if (_ffun._ndim == +2 &&
                _ffun._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 2-dim. mesh */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=2 \n" ;

                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*-------------------------- count mesh items */
                iptr_type _nnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _ffun.
                _euclidean_mesh_2d._mesh._set1.head() ;
                          _iter != _ffun.
                _euclidean_mesh_2d._mesh._set1.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                        _nnum +=  +1 ;
                    }
                }

                for (auto _iter  = _ffun.
                _euclidean_mesh_2d._mesh._set3.head() ;
                          _iter != _ffun.
                _euclidean_mesh_2d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                        _tnum +=  +1 ;
                    }
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _nnum << "\n" ;
                
                iptr_type _npos  = +0;
                
                for (auto _iter  = _ffun.
                _euclidean_mesh_2d._mesh._set1.head() ;
                          _iter != _ffun.
                _euclidean_mesh_2d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_nnum > +0)
                {
            /*-------------------------- write VALUE data */                
                _file << "VALUE=" 
                      << _nnum << ";1" << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_mesh_2d._hval.head() ;
                          _iter != _ffun.
                    _euclidean_mesh_2d._hval.tend() ;
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }
                
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << _tnum << "\n" ;
                
                for (auto _iter  = _ffun.
                _euclidean_mesh_2d._mesh._set3.head() ;
                          _iter != _ffun.
                _euclidean_mesh_2d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << _iter->node(0) << ";" 
                          << _iter->node(1) << ";"
                          << _iter->node(2) << ";"
                          <<    +0 << "\n" ;
                    }
                }
                }
     
            }
            else
            if (_ffun._ndim == +3 &&
                _ffun._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 3-dim. mesh */
                _file << "# " << _name << ".msh" 
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;EUCLIDEAN-MESH \n" ;
                _file << "NDIMS=3 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*-------------------------- count mesh items */
                iptr_type _nnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _ffun.
                _euclidean_mesh_3d._mesh._set1.head() ;
                          _iter != _ffun.
                _euclidean_mesh_3d._mesh._set1.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                        _nnum +=  +1 ;
                    }
                }

                for (auto _iter  = _ffun.
                _euclidean_mesh_3d._mesh._set4.head() ;
                          _iter != _ffun.
                _euclidean_mesh_3d._mesh._set4.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                        _tnum +=  +1 ;
                    }
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _nnum << "\n" ;
                
                iptr_type _npos  = +0;
                
                for (auto _iter  = _ffun.
                _euclidean_mesh_3d._mesh._set1.head() ;
                          _iter != _ffun.
                _euclidean_mesh_3d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= 0)
                    {
                    _file << _iter->pval(0) << ";"
                          << _iter->pval(1) << ";"
                          << _iter->pval(2) << ";"
                          << +0 << "\n" ;
                    }
                }
                }
                
                if (_nnum > +0)
                {
            /*-------------------------- write VALUE data */                
                _file << "VALUE=" 
                      << _nnum << ";1" << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_mesh_3d._hval.head() ;
                          _iter != _ffun.
                    _euclidean_mesh_3d._hval.tend() ;
                        ++_iter  )
                {   
                    _file << *_iter << "\n" ;
                }
                }
                
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA4=" << _tnum << "\n" ;
                
                for (auto _iter  = _ffun.
                _euclidean_mesh_3d._mesh._set4.head() ;
                          _iter != _ffun.
                _euclidean_mesh_3d._mesh._set4.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << _iter->node(0) << ";" 
                          << _iter->node(1) << ";"
                          << _iter->node(2) << ";"
                          << _iter->node(3) << ";"
                          <<    +0 << "\n" ;
                    }
                }
                }
            
            }
            else
            if (_ffun._ndim == +2 &&
                _ffun._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
            /*-------------------------- save 3-dim. mesh */
                _file << "# " << _name << ".msh" 
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;ELLIPSOID-MESH \n" ;
                _file << "NDIMS=2 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*-------------------------- write class data */

                _file << "RADII=" << 
                _ffun._ellipsoid_mesh_3d._radA 
                      << ";" <<
                _ffun._ellipsoid_mesh_3d._radB 
                      << ";" <<
                _ffun._ellipsoid_mesh_3d._radC
                      << "\n" ;

            /*-------------------------- count mesh items */
                iptr_type _nnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _ffun.
                _ellipsoid_mesh_3d._mesh._set1.head() ;
                          _iter != _ffun.
                _ellipsoid_mesh_3d._mesh._set1.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                        _nnum +=  +1 ;
                    }
                }

                for (auto _iter  = _ffun.
                _ellipsoid_mesh_3d._mesh._set3.head() ;
                          _iter != _ffun.
                _ellipsoid_mesh_3d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                    {
                        _tnum +=  +1 ;
                    }
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */                
                _file << "POINT=" << _nnum << "\n" ;
                
                iptr_type _npos = +0 ;
                iptr_type _itag = +0 ;

                real_type _ppos[3] ;
                real_type _apos[2] ;
                for (auto _iter  = _ffun.
                _ellipsoid_mesh_3d._mesh._set1.head() ;
                          _iter != _ffun.
                _ellipsoid_mesh_3d._mesh._set1.tend() ;
                        ++_iter, ++_npos)
                {
                    if (_iter->mark() >= +0)
                    {
                    _ppos[0] =_iter->pval(0);
                    _ppos[1] =_iter->pval(1);
                    _ppos[2] =_iter->pval(2);

                    _ffun._ellipsoid_mesh_3d.
                        toS2(_ppos, _apos);

                    _file << _apos[0] << ";"
                          << _apos[1] << ";"
                          << _itag    <<"\n";
                    }
                }
                }
                
                if (_nnum > +0)
                {
            /*-------------------------- write VALUE data */                
                _file << "VALUE=" 
                      << _nnum << ";1" << "\n" ;
                
                for (auto _iter  = _ffun.
                    _ellipsoid_mesh_3d._hval.head() ;
                          _iter != _ffun.
                    _ellipsoid_mesh_3d._hval.tend() ;
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }
                
                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */ 
                _file << "TRIA3=" << _tnum << "\n" ;
                
                for (auto _iter  = _ffun.
                _ellipsoid_mesh_3d._mesh._set3.head() ;
                          _iter != _ffun.
                _ellipsoid_mesh_3d._mesh._set3.tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= 0 &&
                        _iter->self() >= 1 )
                    {
                    _file << _iter->node(0) << ";" 
                          << _iter->node(1) << ";"
                          << _iter->node(2) << ";"
                          <<    +0 << "\n" ;
                    }
                }
                }

            }
            else
            if (_ffun._ndim == +2 &&
                _ffun._kind ==
                jmsh_kind::euclidean_grid)
            {
            /*-------------------------- save 2-dim. grid */
                _file << "# " << _name << ".msh" 
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;EUCLIDEAN-GRID \n" ;
                _file << "NDIMS=2 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*-------------------------- write class data */

                if (_ffun._euclidean_grid_2d.
                        _xpos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=1;" << 
                    _ffun._euclidean_grid_2d.
                        _xpos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_2d._xpos.head();
                          _iter != _ffun.
                    _euclidean_grid_2d._xpos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._euclidean_grid_2d.
                        _ypos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=2;" << 
                    _ffun._euclidean_grid_2d.
                        _ypos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_2d._ypos.head();
                          _iter != _ffun.
                    _euclidean_grid_2d._ypos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._euclidean_grid_2d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                _file << "VALUE=" << 
                    _ffun._euclidean_grid_2d.
                        _hmat.count() << ";1" "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_2d._hmat.head();
                          _iter != _ffun.
                    _euclidean_grid_2d._hmat.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

            }
            else
            if (_ffun._ndim == +3 &&
                _ffun._kind ==
                jmsh_kind::euclidean_grid)
            {
            /*-------------------------- save 3-dim. grid */
                _file << "# " << _name << ".msh" 
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;EUCLIDEAN-GRID \n" ;
                _file << "NDIMS=3 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*-------------------------- write class data */

                if (_ffun._euclidean_grid_3d.
                        _xpos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=1;" << 
                    _ffun._euclidean_grid_3d.
                        _xpos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_3d._xpos.head();
                          _iter != _ffun.
                    _euclidean_grid_3d._xpos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._euclidean_grid_3d.
                        _ypos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=2;" << 
                    _ffun._euclidean_grid_3d.
                        _ypos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_3d._ypos.head();
                          _iter != _ffun.
                    _euclidean_grid_3d._ypos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._euclidean_grid_3d.
                        _zpos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=3;" << 
                    _ffun._euclidean_grid_3d.
                        _zpos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_3d._zpos.head();
                          _iter != _ffun.
                    _euclidean_grid_3d._zpos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._euclidean_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                _file << "VALUE=" << 
                    _ffun._euclidean_grid_3d.
                        _hmat.count() << ";1" "\n" ;
                
                for (auto _iter  = _ffun.
                    _euclidean_grid_3d._hmat.head();
                          _iter != _ffun.
                    _euclidean_grid_3d._hmat.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

            }
            else
            if (_ffun._ndim == +2 &&
                _ffun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
            /*-------------------------- save 3-dim. mesh */
                _file << "# " << _name << ".msh" 
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;ELLIPSOID-GRID \n" ;
                _file << "NDIMS=2 \n" ;
                
                _file << std::scientific ;
                _file << 
                    std::setprecision(16);

            /*-------------------------- write class data */

                _file << "RADII=" << 
                _ffun._ellipsoid_grid_3d._radA 
                      << ";" <<
                _ffun._ellipsoid_grid_3d._radB 
                      << ";" <<
                _ffun._ellipsoid_grid_3d._radC
                      << "\n" ;

                if (_ffun._ellipsoid_grid_3d.
                        _xpos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=1;" << 
                    _ffun._ellipsoid_grid_3d.
                        _xpos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _ellipsoid_grid_3d._xpos.head();
                          _iter != _ffun.
                    _ellipsoid_grid_3d._xpos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._ellipsoid_grid_3d.
                        _ypos.count() > +0)
                {
            /*-------------------------- write COORD data */                
                _file << "COORD=2;" << 
                    _ffun._ellipsoid_grid_3d.
                        _ypos.count() << "\n" ;
                
                for (auto _iter  = _ffun.
                    _ellipsoid_grid_3d._ypos.head();
                          _iter != _ffun.
                    _ellipsoid_grid_3d._ypos.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

                if (_ffun._ellipsoid_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                _file << "VALUE=" << 
                    _ffun._ellipsoid_grid_3d.
                        _hmat.count() << ";1" "\n" ;
                
                for (auto _iter  = _ffun.
                    _ellipsoid_grid_3d._hmat.head();
                          _iter != _ffun.
                    _ellipsoid_grid_3d._hmat.tend();
                        ++_iter  )
                {
                    _file << *_iter << "\n" ;
                }
                }

            }

            }
            else
            {
                _errv = __file_not_located  ;
            }
            
            _file.close();

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    /*
    --------------------------------------------------------
     * SAVE-MESH: save MSH_t output data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        hfun_data &_ffun ,
        jigsaw_msh_t &_fmsh
        )
    {
        iptr_type _errv  = __no_error  ;

    //  NB. this routine handles output for MARCHE --- only 
    //  *overwriting* the "value" arrays in FMSH, since all
    //  else is passed read-only...

        try
        {
            __unreferenced (_jcfg) ;
            __unreferenced (_jlog) ;
     
            if (_ffun._ndim == +2 &&    // save 2-dim. mesh
                _ffun._kind ==
                jmsh_kind::euclidean_mesh)
            {
                _fmsh._flags = JIGSAW_EUCLIDEAN_MESH ;

                if (_ffun._euclidean_mesh_2d.
                        _hval.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                __assert( _fmsh._value._size == 
                    _ffun._euclidean_mesh_2d.
                               _hval.count() &&
                "savemsh: inconsistent VALUE array") ;                
    
                iptr_type _ipos  = +0 ;
                for (auto _iter  = _ffun.
                    _euclidean_mesh_2d._hval.head();
                          _iter != _ffun.
                    _euclidean_mesh_2d._hval.tend();
                        ++_iter, ++_ipos)
                {
                    _fmsh._value.
                        _data[_ipos] = *_iter ;
                }
                }

            }
            else
            if (_ffun._ndim == +3 &&    // save 3-dim. mesh
                _ffun._kind ==
                jmsh_kind::euclidean_mesh)
            {
                _fmsh._flags = JIGSAW_EUCLIDEAN_MESH ;

                if (_ffun._euclidean_mesh_3d.
                        _hval.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                __assert( _fmsh._value._size == 
                    _ffun._euclidean_mesh_3d.
                               _hval.count() &&
                "savemsh: inconsistent VALUE array") ;                
    
                iptr_type _ipos  = +0 ;
                for (auto _iter  = _ffun.
                    _euclidean_mesh_3d._hval.head();
                          _iter != _ffun.
                    _euclidean_mesh_3d._hval.tend();
                        ++_iter, ++_ipos)
                {
                    _fmsh._value.
                        _data[_ipos] = *_iter ;
                }
                }

            }
            else
            if (_ffun._ndim == +2 &&    // save 2-sph. mesh
                _ffun._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
                _fmsh._flags = JIGSAW_ELLIPSOID_MESH ;

                if (_ffun._ellipsoid_mesh_3d.
                        _hval.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                __assert( _fmsh._value._size == 
                    _ffun._ellipsoid_mesh_3d.
                               _hval.count() &&
                "savemsh: inconsistent VALUE array") ;                
    
                iptr_type _ipos  = +0 ;
                for (auto _iter  = _ffun.
                    _ellipsoid_mesh_3d._hval.head();
                          _iter != _ffun.
                    _ellipsoid_mesh_3d._hval.tend();
                        ++_iter, ++_ipos)
                {
                    _fmsh._value.
                        _data[_ipos] = *_iter ;
                }
                }

            }
            else    
            if (_ffun._ndim == +2 &&    // save 2-dim. grid
                _ffun._kind ==
                jmsh_kind::euclidean_grid)
            {
                _fmsh._flags = JIGSAW_EUCLIDEAN_GRID ;

                if (_ffun._euclidean_grid_2d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                __assert( _fmsh._value._size == 
                    _ffun._euclidean_grid_2d.
                               _hmat.count() &&
                "savemsh: inconsistent VALUE array") ;                
    
                iptr_type _ipos  = +0 ;
                for (auto _iter  = _ffun.
                    _euclidean_grid_2d._hmat.head();
                          _iter != _ffun.
                    _euclidean_grid_2d._hmat.tend();
                        ++_iter, ++_ipos)
                {
                    _fmsh._value.
                        _data[_ipos] = *_iter ;
                }
                }

            }
            else
            if (_ffun._ndim == +3 &&    // save 3-dim. grid
                _ffun._kind ==
                jmsh_kind::euclidean_grid)
            {
                _fmsh._flags = JIGSAW_EUCLIDEAN_GRID ;

                if (_ffun._euclidean_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                __assert( _fmsh._value._size == 
                    _ffun._euclidean_grid_3d.
                               _hmat.count() &&
                "savemsh: inconsistent VALUE array") ;                
    
                iptr_type _ipos  = +0 ;
                for (auto _iter  = _ffun.
                    _euclidean_grid_3d._hmat.head();
                          _iter != _ffun.
                    _euclidean_grid_3d._hmat.tend();
                        ++_iter, ++_ipos)
                {
                    _fmsh._value.
                        _data[_ipos] = *_iter ;
                }
                }

            }
            else
            if (_ffun._ndim == +2 &&    // save 2-sph. grid
                _ffun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
                _fmsh._flags = JIGSAW_ELLIPSOID_GRID ;

                if (_ffun._ellipsoid_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                __assert( _fmsh._value._size == 
                    _ffun._ellipsoid_grid_3d.
                               _hmat.count() &&
                "savemsh: inconsistent VALUE array") ;                
    
                iptr_type _ipos  = +0 ;
                for (auto _iter  = _ffun.
                    _ellipsoid_grid_3d._hmat.head();
                          _iter != _ffun.
                    _ellipsoid_grid_3d._hmat.tend();
                        ++_iter, ++_ipos)
                {
                    _fmsh._value.
                        _data[_ipos] = *_iter ;
                }
                }

            }

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }


#   endif   //__MSH_SAVE__



