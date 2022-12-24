
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
     * Last updated: 30 May, 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from msh_save.hpp


    /*
    --------------------------------------------------------
     * SAVE-TRIA: save *.MSH output file.
    --------------------------------------------------------
     */

    #define save_vert2(_mesh)       \
        _file << "POINT=" << _last << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0,       \
              _npos = +0;                   \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh._nset.head();      \
                  _iter != _mesh._nset.tend();      \
                ++_iter, ++_npos)   \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _nmap[_npos ] >= 0)         \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next],      \
                PRINTCHUNK,         \
                "%.17g;%.17g;+0\n", \
                _iter->pval(0) + _xoff[0] ,         \
                _iter->pval(1) + _xoff[1]),         \
                VERT2CHUNK)         \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_vert3(_mesh)       \
        _file << "POINT=" << _last << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0,       \
              _npos = +0;                   \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh._nset.head();      \
                  _iter != _mesh._nset.tend();      \
                ++_iter, ++_npos)           \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _nmap[_npos ] >= 0)         \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next],      \
                PRINTCHUNK,         \
                "%.17g;%.17g;%.17g;+0\n",   \
                _iter->pval(0) + _xoff[0] , \
                _iter->pval(1) + _xoff[1] , \
                _iter->pval(2) + _xoff[2]), \
                VERT3CHUNK)         \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_vpwrk(_koff, _mesh)        \
        _file << "POWER=" << _last << ";1" << "\n" ;    \
                                    \
         auto _next = +0, _roll = +0,       \
              _npos = +0;                   \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh._nset.head();      \
                  _iter != _mesh._nset.tend();      \
                ++_iter, ++_npos)           \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _nmap[_npos ] >= 0)         \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next] ,     \
                PRINTCHUNK, "%.17g\n" ,     \
            _iter->pval(_koff) ), VALUECHUNK)       \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_tria3(_mesh)       \
        _file << "TRIA3=" << _ntri << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0;       \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh._tset.head();      \
                  _iter != _mesh._tset.tend();      \
                ++_iter  )          \
        {                           \
            if (_iter->mark() < +0) continue ;      \
                                    \
            PRINTCHARS(snprintf(&_fbuf[_next],      \
                PRINTCHUNK, "%u;%u;%u;+0\n",        \
            _nmap[_iter->node(0)] , \
            _nmap[_iter->node(1)] , \
            _nmap[_iter->node(2)]),TRIA3CHUNK)      \
        }                           \
            PRINTFINAL;             \


    #define save_tria4(_mesh)       \
        _file << "TRIA4=" << _ntri << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0;       \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh._tset.head();      \
                  _iter != _mesh._tset.tend();      \
                ++_iter  )          \
        {                           \
            if (_iter->mark() < +0) continue ;      \
                                    \
            PRINTCHARS(snprintf(&_fbuf[_next],      \
                PRINTCHUNK,         \
                "%u;%u;%u;%u;+0\n", \
            _nmap[_iter->node(0)] , \
            _nmap[_iter->node(1)] , \
            _nmap[_iter->node(2)] , \
            _nmap[_iter->node(3)]),TRIA4CHUNK)      \
        }                           \
            PRINTFINAL;             \


    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_tria (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float     *_xoff ,
        mesh_data &_rdel ,
        mshfile_tag const&
        )
    {
        iptr_type _errv  = __no_error  ;

        try
        {
            containers::array<iptr_type> _nmap;

            __unreferenced(_jlog) ;

            std::string _path, _name, _fext;
            file_part(
                _jcfg._tria_file, _path, _name, _fext) ;

            std::ofstream  _file;
            _file.open(
                _jcfg._tria_file,
            std::ofstream::out | std::ofstream::trunc) ;

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

                    _ntri += +1  ;

                    _nmap[_iter->node(0)] = 1;
                    _nmap[_iter->node(1)] = 1;
                    _nmap[_iter->node(2)] = 1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0) *_iter = _last ++  ;
                }

                if (_last > +0)
                {
            /*-------------------------- write POINT data */
                    save_vert2(
                        _rdel._euclidean_rdel_2d._tria)
                }
                if (_last > +0)
                {
            /*-------------------------- write POWER data */
                    save_vpwrk(
                    +2, _rdel._euclidean_rdel_2d._tria)
                }
                if (_ntri > +0)
                {
            /*-------------------------- write TRIA3 data */
                    save_tria3(
                        _rdel._euclidean_rdel_2d._tria)
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

                    _ntri += +1  ;

                    _nmap[_iter->node(0)] = 1;
                    _nmap[_iter->node(1)] = 1;
                    _nmap[_iter->node(2)] = 1;
                    _nmap[_iter->node(3)] = 1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0) *_iter = _last ++  ;
                }

                if (_last > +0)
                {
            /*-------------------------- write POINT data */
                    save_vert3(
                        _rdel._euclidean_rdel_3d._tria)
                }
                if (_last > +0)
                {
            /*-------------------------- write POWER data */
                    save_vpwrk(
                    +3, _rdel._euclidean_rdel_3d._tria)
                }
                if (_ntri > +0)
                {
            /*-------------------------- write TRIA4 data */
                    save_tria4(
                        _rdel._euclidean_rdel_3d._tria)
                }

            }

            }
            else
            {
                _errv = __file_not_created ;
            }

            _file.close();

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    #undef save_vert2
    #undef save_vert3
    #undef save_vpwrk
    #undef save_tria3
    #undef save_tria4

    /*
    --------------------------------------------------------
     * SAVE-TRIA: save *.NC output file.
    --------------------------------------------------------
     */

    #ifdef __use_netcdf

    #define save_vert2(_mesh)   \
        ncutil::def_dim(        \
            _jcfg._tria_file, "vert2_row", _last);  \
        ncutil::def_dim(        \
            _jcfg._tria_file, "vert2_col", +2) ;    \
                                \
        ncutil::def_var(        \
            _jcfg._tria_file, "vert2_pos", NC_DOUBLE,   \
                "2-vertex coordinates" ,            \
                {"vert2_row", "vert2_col"} ) ;      \
        ncutil::def_var(        \
            _jcfg._tria_file, "vert2_tag", NC_INT,  \
                "2-vertex id-tags" ,                \
                {"vert2_row", "singleton"} ) ;      \
        ncutil::def_var(        \
            _jcfg._tria_file, "power_val", NC_DOUBLE,   \
                "vertex power values" ,             \
                {"vert2_row", "singleton"} ) ;      \
                                \
        _pos.clear(); _pos.set_alloc(_last * 2);    \
        _tag.clear(); _tag.set_alloc(_last * 1);    \
        _val.clear(); _val.set_alloc(_last * 1);    \
                                \
        auto _npos = +0;        \
        for (auto _iter  = _mesh._nset.head() ;     \
                  _iter != _mesh._nset.tend() ;     \
                ++_iter, ++_npos)       \
        {                       \
            if (_iter->mark() >= 0 &&   \
                _nmap[_npos ] >= 0 )    \
            {                   \
                _pos.push_tail( \
                    _iter->pval(0) + _xoff[0]) ;    \
                _pos.push_tail( \
                    _iter->pval(1) + _xoff[1]) ;    \
                                \
                _tag.push_tail (0) ;    \
                                \
                _val.push_tail(_iter->pval(2)) ;    \
            }                   \
        }                       \
                                \
        ncutil::put_var(        \
        _jcfg._tria_file, "vert2_pos", &_pos[0]) ;  \
        ncutil::put_var(        \
        _jcfg._tria_file, "vert2_tag", &_tag[0]) ;  \
        ncutil::put_var(        \
        _jcfg._tria_file, "power_val", &_val[0]) ;  \


    #define save_vert3(_mesh)   \
        ncutil::def_dim(        \
            _jcfg._tria_file, "vert3_row", _last);  \
        ncutil::def_dim(        \
            _jcfg._tria_file, "vert3_col", +3) ;    \
                                \
        ncutil::def_var(        \
            _jcfg._tria_file, "vert3_pos", NC_DOUBLE,   \
                "3-vertex coordinates" ,            \
                {"vert3_row", "vert3_col"} ) ;      \
        ncutil::def_var(        \
            _jcfg._tria_file, "vert3_tag", NC_INT,  \
                "3-vertex id-tags" ,                \
                {"vert3_row", "singleton"} ) ;      \
        ncutil::def_var(        \
            _jcfg._tria_file, "power_val", NC_DOUBLE,   \
                "vertex power values" ,             \
                {"vert3_row", "singleton"} ) ;      \
                                \
        _pos.clear(); _pos.set_alloc(_last * 3);    \
        _tag.clear(); _tag.set_alloc(_last * 1);    \
        _val.clear(); _val.set_alloc(_last * 1);    \
                                \
        auto _npos = +0;        \
        for (auto _iter  = _mesh._nset.head() ;     \
                  _iter != _mesh._nset.tend() ;     \
                ++_iter, ++_npos)       \
        {                       \
            if (_iter->mark() >= 0 &&   \
                _nmap[_npos ] >= 0 )    \
            {                   \
                _pos.push_tail( \
                    _iter->pval(0) + _xoff[0]) ;    \
                _pos.push_tail( \
                    _iter->pval(1) + _xoff[1]) ;    \
                _pos.push_tail( \
                    _iter->pval(2) + _xoff[2]) ;    \
                                \
                _tag.push_tail (0) ;    \
                                \
                _val.push_tail(_iter->pval(3)) ;    \
            }                   \
        }                       \
                                \
        ncutil::put_var(        \
        _jcfg._tria_file, "vert3_pos", &_pos[0]) ;  \
        ncutil::put_var(        \
        _jcfg._tria_file, "vert3_tag", &_tag[0]) ;  \
        ncutil::put_var(        \
        _jcfg._tria_file, "power_val", &_val[0]) ;  \


    #define save_tria3(_mesh)   \
        ncutil::def_dim(        \
            _jcfg._tria_file, "tria3_row", _ntri);  \
        ncutil::def_dim(        \
            _jcfg._tria_file, "tria3_col", +3) ;    \
                                \
        ncutil::def_var(        \
            _jcfg._tria_file, "tria3_idx", NC_INT,  \
                "3-node tria indexes" ,             \
                {"tria3_row", "tria3_col"} ) ;      \
        ncutil::def_var(        \
            _jcfg._tria_file, "tria3_tag", NC_INT,  \
                "3-node tria id-tags" ,             \
                {"tria3_row", "singleton"} ) ;      \
                                \
        _idx.clear(); _idx.set_alloc(_ntri * 3);    \
        _tag.clear(); _tag.set_alloc(_ntri * 1);    \
                                \
        for (auto _iter  = _mesh._tset.head() ;     \
                  _iter != _mesh._tset.tend() ;     \
                ++_iter  )      \
        {                       \
            if (_iter->mark() >= 0)     \
            {                   \
                _idx.push_tail( \
                    _nmap[_iter->node(0)]) ;        \
                _idx.push_tail( \
                    _nmap[_iter->node(1)]) ;        \
                _idx.push_tail( \
                    _nmap[_iter->node(2)]) ;        \
                                \
                _tag.push_tail(0) ;     \
            }                   \
        }                       \
                                \
        ncutil::put_var(        \
        _jcfg._tria_file, "tria3_idx", &_idx[0]) ;  \
        ncutil::put_var(        \
        _jcfg._tria_file, "tria3_tag", &_tag[0]) ;  \


    #define save_tria4(_mesh)   \
        ncutil::def_dim(        \
            _jcfg._tria_file, "tria4_row", _ntri);  \
        ncutil::def_dim(        \
            _jcfg._tria_file, "tria4_col", +4) ;    \
                                \
        ncutil::def_var(        \
            _jcfg._tria_file, "tria4_idx", NC_INT,  \
                "4-node tria indexes" ,             \
                {"tria4_row", "tria4_col"} ) ;      \
        ncutil::def_var(        \
            _jcfg._tria_file, "tria4_tag", NC_INT,  \
                "4-node tria id-tags" ,             \
                {"tria4_row", "singleton"} ) ;      \
                                \
        _idx.clear(); _idx.set_alloc(_ntri * 4);    \
        _tag.clear(); _tag.set_alloc(_ntri * 1);    \
                                \
        for (auto _iter  = _mesh._tset.head() ;     \
                  _iter != _mesh._tset.tend() ;     \
                ++_iter  )      \
        {                       \
            if (_iter->mark() >= 0)     \
            {                   \
                _idx.push_tail( \
                    _nmap[_iter->node(0)]) ;        \
                _idx.push_tail( \
                    _nmap[_iter->node(1)]) ;        \
                _idx.push_tail( \
                    _nmap[_iter->node(2)]) ;        \
                _idx.push_tail( \
                    _nmap[_iter->node(3)]) ;        \
                                \
                _tag.push_tail(0) ;     \
            }                   \
        }                       \
                                \
        ncutil::put_var(        \
        _jcfg._tria_file, "tria4_idx", &_idx[0]) ;  \
        ncutil::put_var(        \
        _jcfg._tria_file, "tria4_tag", &_tag[0]) ;  \


    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_tria (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float     *_xoff ,
        mesh_data &_rdel ,
        _netcdf_tag const&
        )
    {
        iptr_type _errv  = __no_error  ;

        try
        {
            containers::array<iptr_type> _nmap;

            __unreferenced(_jlog) ;

            std::string _path, _name, _fext;
            file_part(
                _jcfg._tria_file, _path, _name, _fext) ;

            int _retv, _ncid;
            if ((_retv = nc_create(
                    _jcfg._tria_file.c_str(),
                        NC_CLOBBER|NC_NETCDF4, &_ncid)))
            {
                return __file_not_created;
            }

            if (_rdel._ndim == +2 &&
                _rdel._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 2-dim. mesh */
                ncutil::put_str(_jcfg._tria_file,
                    "title",
                    _name + "; created by " + __JGSWVSTR);

                ncutil::put_str(_jcfg._tria_file,
                    "MSHID", "EUCLIDEAN-MESH" ) ;

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

                    _ntri += +1  ;

                    _nmap[_iter->node(0)] = 1;
                    _nmap[_iter->node(1)] = 1;
                    _nmap[_iter->node(2)] = 1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0) *_iter = _last ++  ;
                }

                ncutil::def_dim(
                    _jcfg._tria_file, "singleton", +1);

                containers::array<real_type> _pos;
                containers::array<iptr_type> _idx;
                containers::array<iptr_type> _tag;
                containers::array<real_type> _val;

                if (_last > +0)
                {
            /*-------------------------- write POINT data */
                    save_vert2(
                        _rdel._euclidean_rdel_2d._tria)
                }
                if (_ntri > +0)
                {
            /*-------------------------- write TRIA3 data */
                    save_tria3(
                        _rdel._euclidean_rdel_2d._tria)
                }

            }
            else
            if (_rdel._ndim == +3 &&
                _rdel._kind ==
                jmsh_kind::euclidean_mesh)
            {
            /*-------------------------- save 3-dim. mesh */
                ncutil::put_str(_jcfg._tria_file,
                    "title",
                    _name + "; created by " + __JGSWVSTR);

                ncutil::put_str(_jcfg._tria_file,
                    "MSHID", "EUCLIDEAN-MESH" ) ;

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

                    _ntri += +1  ;

                    _nmap[_iter->node(0)] = 1;
                    _nmap[_iter->node(1)] = 1;
                    _nmap[_iter->node(2)] = 1;
                    _nmap[_iter->node(3)] = 1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0) *_iter = _last ++  ;
                }

                ncutil::def_dim(
                    _jcfg._tria_file, "singleton", +1);

                containers::array<real_type> _pos;
                containers::array<iptr_type> _idx;
                containers::array<iptr_type> _tag;
                containers::array<real_type> _val;

                if (_last > +0)
                {
            /*-------------------------- write POINT data */
                    save_vert3(
                        _rdel._euclidean_rdel_3d._tria)
                }
                if (_ntri > +0)
                {
            /*-------------------------- write TRIA3 data */
                    save_tria4(
                        _rdel._euclidean_rdel_3d._tria)
                }

            }
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    #undef save_vert2
    #undef save_vert3
    #undef save_tria3
    #undef save_tria4

    #endif  //__use_netcdf

    /*
    --------------------------------------------------------
     * SAVE-TRIA: save MSH_t output data.
    --------------------------------------------------------
     */

    #define save_vert2(_mesh)       \
        jigsaw_alloc_vert2(&_mmsh._vert2, _last) ;  \
        jigsaw_alloc_reals(&_mmsh._power, _last) ;  \
                                    \
        iptr_type _npos  = +0 ;     \
        iptr_type _nout  = +0 ;     \
                                    \
        for (auto _iter  = _mesh._nset.head();      \
                  _iter != _mesh._nset.tend();      \
                ++_iter, ++_npos)   \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _nmap[_npos ] >= 0 )        \
            {                       \
            _mmsh._vert2._data[_nout].      \
            _ppos[0] = _iter->pval(0) + _xoff[0] ;  \
            _mmsh._vert2._data[_nout].      \
            _ppos[1] = _iter->pval(1) + _xoff[1] ;  \
                                    \
            _mmsh._vert2. _data[_nout]._itag = 0 ;  \
                                    \
            _mmsh._power.           \
            _data[_nout] = _iter->pval(2) ; \
                                    \
            _nout = _nout + 1 ;     \
            }                       \
        }                           \
                                    \
        _mesh._nset.clear(containers::tight_alloc) ;\


    #define save_vert3(_mesh)       \
        jigsaw_alloc_vert3(&_mmsh._vert3, _last) ;  \
        jigsaw_alloc_reals(&_mmsh._power, _last) ;  \
                                    \
        iptr_type _npos  = +0 ;     \
        iptr_type _nout  = +0 ;     \
                                    \
        for (auto _iter  = _mesh._nset.head();      \
                  _iter != _mesh._nset.tend();      \
                ++_iter, ++_npos)   \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _nmap[_npos ] >= 0 )        \
            {                       \
            _mmsh._vert3._data[_nout].      \
            _ppos[0] = _iter->pval(0) + _xoff[0] ;  \
            _mmsh._vert3._data[_nout].      \
            _ppos[1] = _iter->pval(1) + _xoff[1] ;  \
            _mmsh._vert3._data[_nout].      \
            _ppos[2] = _iter->pval(2) + _xoff[2] ;  \
                                    \
            _mmsh._vert3. _data[_nout]._itag = 0 ;  \
                                    \
            _mmsh._power.           \
            _data[_nout] = _iter->pval(3) ; \
                                    \
            _nout = _nout + 1 ;     \
            }                       \
        }                           \
                                    \
        _mesh._nset.clear(containers::tight_alloc) ;\


    #define save_tria3(_mesh)       \
        jigsaw_alloc_tria3(&_mmsh._tria3, _ntri) ;  \
                                    \
        iptr_type _tout  = +0 ;     \
                                    \
        for (auto _iter  = _mesh._tset.head();      \
                  _iter != _mesh._tset.tend();      \
                ++_iter  )          \
        {                           \
            if (_iter->mark() >= 0) \
            {                       \
            _mmsh._tria3._data[_tout].      \
                _node[0] = _nmap[_iter->node(0)];   \
            _mmsh._tria3._data[_tout].      \
                _node[1] = _nmap[_iter->node(1)];   \
            _mmsh._tria3._data[_tout].      \
                _node[2] = _nmap[_iter->node(2)];   \
                                    \
            _mmsh._tria3._data[_tout]._itag = 0 ;   \
                                    \
            _tout = _tout + 1 ;     \
            }                       \
        }                           \
                                    \
        _mesh._tset.clear(containers::tight_alloc) ;\


    #define save_tria4(_mesh)       \
        jigsaw_alloc_tria4(&_mmsh._tria4, _ntri) ;  \
                                    \
        iptr_type _tout  = +0 ;     \
                                    \
        for (auto _iter  = _mesh._tset.head();      \
                  _iter != _mesh._tset.tend();      \
                ++_iter  )          \
        {                           \
            if (_iter->mark() >= 0) \
            {                       \
            _mmsh._tria4._data[_tout].      \
                _node[0] = _nmap[_iter->node(0)];   \
            _mmsh._tria4._data[_tout].      \
                _node[1] = _nmap[_iter->node(1)];   \
            _mmsh._tria4._data[_tout].      \
                _node[2] = _nmap[_iter->node(2)];   \
            _mmsh._tria4._data[_tout].      \
                _node[3] = _nmap[_iter->node(3)];   \
                                    \
            _mmsh._tria4._data[_tout]._itag = 0 ;   \
                                    \
            _tout = _tout + 1 ;     \
            }                       \
        }                           \
                                    \
        _mesh._tset.clear(containers::tight_alloc) ;\


    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_tria (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float     *_xoff ,
        mesh_data &_rdel ,
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

                iptr_type _ntri = +0;
                for (auto _iter  = _rdel.
                _euclidean_rdel_2d._tria._tset.head();
                          _iter != _rdel.
                _euclidean_rdel_2d._tria._tset.tend();
                        ++_iter  )
                {
                    if (_iter->mark() < +0) continue ;

                    _ntri += +1  ;

                    _nmap[_iter->node(0)] = 1;
                    _nmap[_iter->node(1)] = 1;
                    _nmap[_iter->node(2)] = 1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0) *_iter = _last ++  ;
                }

                if (_last > +0)
                {
            /*-------------------------- write POINT data */
                    save_vert2(
                        _rdel._euclidean_rdel_2d._tria)
                }
                if (_ntri > +0)
                {
            /*-------------------------- write TRIA3 data */
                    save_tria3(
                        _rdel._euclidean_rdel_2d._tria)
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

                    _ntri += +1  ;

                    _nmap[_iter->node(0)] = 1;
                    _nmap[_iter->node(1)] = 1;
                    _nmap[_iter->node(2)] = 1;
                    _nmap[_iter->node(3)] = 1;
                }

                iptr_type _last  = +0;
                for (auto _iter  = _nmap.head() ;
                          _iter != _nmap.tend() ;
                        ++_iter  )
                {
                    if ( *_iter >= +0) *_iter = _last ++  ;
                }

                if (_last > +0)
                {
            /*-------------------------- write POINT data */
                    save_vert3(
                        _rdel._euclidean_rdel_3d._tria)
                }
                if (_ntri > +0)
                {
            /*-------------------------- write TRIA3 data */
                    save_tria4(
                        _rdel._euclidean_rdel_3d._tria)
                }

            }

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    #undef save_vert2
    #undef save_vert3
    #undef save_tria3
    #undef save_tria4

    /*
    --------------------------------------------------------
     * SAVE-TRIA: save to *.MSH or *.NC
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_tria (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float     *_xoff ,
        mesh_data &_rdel
        )
    {
        iptr_type _errv  = __no_error ;

        try
        {
    /*------------------------- disambiguate via file-ext */
            std::string _path, _name, _fext;
            file_part(
                _jcfg._tria_file, _path, _name, _fext) ;

            if(_fext.find("msh") != std::string::npos)
            {
                return save_tria(_jcfg, _jlog,
                    _xoff, _rdel, mshfile_tag()) ;
            }
            else
            if(_fext.find( "nc") != std::string::npos)
            {
        #       ifdef __use_netcdf
                return save_tria(_jcfg, _jlog,
                    _xoff, _rdel, _netcdf_tag()) ;
        #       else
                return __netcdf_not_available;
        #       endif
            }
        }
        catch (...)
        {
            _errv = __unknown_error ;   // can't get here??
        }

        return ( _errv ) ;
    }



