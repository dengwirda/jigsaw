
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
     * SAVE-HFUN: save *.MSH output file.
    --------------------------------------------------------
     */

    template <
    typename      file_type ,
    typename      list_type
             >
    __normal_call void_type save_vals (
        file_type &_file,
            typename
        list_type::data_type const _xoff ,
        list_type &_list
        )
    {
        auto constexpr VALUEBLOCK = VALUECHUNK / 4 ;

    /*--------------------- save msh_t REALS data to file */
        auto _next = +0, _roll = +0 ;
        char _fbuf[PRINTCHUNK] ;

        auto _this =  _list.head () ;

        for (auto _iter  = _this,
                  _last  =(_list.tend() - 4);
                  _iter  < _last;
                  _iter += +4 , _this += +4)
        {
    /*--------------------- (apparently) faster to unroll */
            PRINTCHARS(snprintf(&_fbuf[_next] ,
            PRINTCHUNK,
            "%.9g\n%.9g\n%.9g\n%.9g\n" ,
            *(_iter + 0) + _xoff,
            *(_iter + 1) + _xoff,
            *(_iter + 2) + _xoff,
            *(_iter + 3) + _xoff) , VALUEBLOCK)
        }

        for (auto _iter  = _this,
                  _last  =(_list.tend() - 0);
                  _iter  < _last;
                  _iter += +1 )
        {
    /*--------------------- (apparently) faster to unroll */
            PRINTCHARS(snprintf(&_fbuf[_next] ,
            PRINTCHUNK, "%.9g\n",
            *(_iter + 0) + _xoff) , VALUECHUNK)
        }

        {
    /*--------------------- ensure last buffer is written */
            PRINTFINAL ;
        }
    }

    #define save_vert2(_mesh)       \
        _file << "POINT=" << _nnum << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0,       \
              _npos = +0;           \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh.node().head() ;    \
                  _iter != _mesh.node().tend() ;    \
                ++_iter, ++_npos)   \
        {                           \
            if (_iter->mark() >= 0 )        \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next]  ,    \
            PRINTCHUNK, "%.17g;%.17g;+0\n",         \
                _iter->pval(0) + _xoff[0] , \
                _iter->pval(1) + _xoff[1]), \
            VERT2CHUNK)             \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_vert3(_mesh)       \
        _file << "POINT=" << _nnum << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0,       \
              _npos = +0;           \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh.node().head() ;    \
                  _iter != _mesh.node().tend() ;    \
                ++_iter, ++_npos)   \
        {                           \
            if (_iter->mark() >= 0 )        \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next]  ,    \
            PRINTCHUNK,"%.17g;%.17g;%.17g;+0\n",    \
                _iter->pval(0) + _xoff[0] , \
                _iter->pval(1) + _xoff[1] , \
                _iter->pval(2) + _xoff[2]), \
            VERT3CHUNK)             \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_vertS(_mesh)       \
        _file << "POINT=" << _nnum << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0,       \
              _npos = +0;           \
         char _fbuf[PRINTCHUNK] ;   \
                                    \
        for (auto _iter  = _mesh.node().head() ;    \
                  _iter != _mesh.node().tend() ;    \
                ++_iter, ++_npos)   \
        {                           \
            if (_iter->mark() >= +0)        \
            {                       \
            real_type _ppos[3], _apos[2] ;  \
            _ppos[0] =_iter->pval(0);       \
            _ppos[1] =_iter->pval(1);       \
            _ppos[2] =_iter->pval(2);       \
                                    \
            _ffun._ellipsoid_mesh_3d.       \
                toS2( _ppos, _apos );       \
                                    \
            PRINTCHARS(snprintf(&_fbuf[_next] ,     \
            PRINTCHUNK ,"%.17g;%.17g;+0\n" ,        \
            _apos[ 0], _apos[ 1]) , VERT2CHUNK)     \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_tria3(_mesh)       \
        _file << "TRIA3=" << _tnum << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0;       \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh.tri3().head() ;    \
                  _iter != _mesh.tri3().tend() ;    \
                ++_iter  )          \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _iter->self() >= 1 )        \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next] ,     \
            PRINTCHUNK, "%u;%u;%u;+0\n" ,   \
                _iter->node(0) ,    \
                _iter->node(1) ,    \
                _iter->node(2) ), TRIA3CHUNK)       \
            }                       \
        }                           \
            PRINTFINAL;             \


    #define save_tria4(_mesh)       \
        _file << "TRIA4=" << _tnum << "\n" ;        \
                                    \
         auto _next = +0, _roll = +0;       \
         char _fbuf[PRINTCHUNK] ;           \
                                    \
        for (auto _iter  = _mesh.tri4().head() ;    \
                  _iter != _mesh.tri4().tend() ;    \
                ++_iter  )          \
        {                           \
            if (_iter->mark() >= 0 &&       \
                _iter->self() >= 1 )        \
            {                       \
            PRINTCHARS(snprintf(&_fbuf[_next] ,     \
            PRINTCHUNK, "%u;%u;%u;%u;+0\n" ,        \
                _iter->node(0) ,    \
                _iter->node(1) ,    \
                _iter->node(2) ,    \
                _iter->node(3) ), TRIA4CHUNK)       \
            }                       \
        }                           \
            PRINTFINAL;             \


    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float     *_xoff ,
        hfun_data &_ffun ,
        mshfile_tag const&
        )
    {
        iptr_type _errv  = __no_error  ;

    //  NB. this routine handles output for MARCHE --- this
    //  will *overwrite* the hfun-file input...

        __unreferenced(_jlog) ;

        try
        {
            std::string _path, _name, _fext;
            file_part(
                _jcfg._hfun_file, _path, _name, _fext) ;

            std::ofstream  _file;
            _file.open(
                _jcfg._hfun_file,
            std::ofstream::out | std::ofstream::trunc) ;

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

            /*-------------------------- count mesh items */
                iptr_type _nnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _ffun.
                _euclidean_mesh_2d._mesh.node().head() ;
                          _iter != _ffun.
                _euclidean_mesh_2d._mesh.node().tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                        _nnum +=  +1 ;
                }

                for (auto _iter  = _ffun.
                _euclidean_mesh_2d._mesh.tri3().head() ;
                          _iter != _ffun.
                _euclidean_mesh_2d._mesh.tri3().tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                        _tnum +=  +1 ;
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */
                save_vert2(
                    _ffun._euclidean_mesh_2d._mesh);
                }

                if (_nnum > +0)
                {
            /*-------------------------- write VALUE data */
                _file << "VALUE="
                      << _nnum << ";1" << "\n" ;

                save_vals(_file, +0.f,
                    _ffun._euclidean_mesh_2d._hval);
                }

                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */
                save_tria3(
                    _ffun._euclidean_mesh_2d._mesh);
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

            /*-------------------------- count mesh items */
                iptr_type _nnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _ffun.
                _euclidean_mesh_3d._mesh.node().head() ;
                          _iter != _ffun.
                _euclidean_mesh_3d._mesh.node().tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                        _nnum +=  +1 ;
                }

                for (auto _iter  = _ffun.
                _euclidean_mesh_3d._mesh.tri4().head() ;
                          _iter != _ffun.
                _euclidean_mesh_3d._mesh.tri4().tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                        _tnum +=  +1 ;
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */
                save_vert2(
                    _ffun._euclidean_mesh_3d._mesh);
                }

                if (_nnum > +0)
                {
            /*-------------------------- write VALUE data */
                _file << "VALUE="
                      << _nnum << ";1" << "\n" ;

                save_vals(_file, +0.f,
                    _ffun._euclidean_mesh_3d._hval);
                }

                if (_tnum > +0)
                {
            /*-------------------------- write TRIA4 data */
                save_tria4(
                    _ffun._euclidean_mesh_3d._mesh);
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

            /*-------------------------- write class data */
                {
                 char _fbuf [PRINTCHUNK] ;
                 auto _char = snprintf(&_fbuf[ 0 ] ,
                    PRINTCHUNK,
                    "%.17g;%.17g;%.17g\n",
                _ffun._ellipsoid_mesh_3d._radA ,
                _ffun._ellipsoid_mesh_3d._radB ,
                _ffun._ellipsoid_mesh_3d._radC) ;

                _file << "RADII="
                      << std::string(_fbuf, _char) ;
                }

            /*-------------------------- count mesh items */
                iptr_type _nnum = +0 ;
                iptr_type _tnum = +0 ;

                for (auto _iter  = _ffun.
                _ellipsoid_mesh_3d._mesh.node().head() ;
                          _iter != _ffun.
                _ellipsoid_mesh_3d._mesh.node().tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                        _nnum +=  +1 ;
                }

                for (auto _iter  = _ffun.
                _ellipsoid_mesh_3d._mesh.tri3().head() ;
                          _iter != _ffun.
                _ellipsoid_mesh_3d._mesh.tri3().tend() ;
                        ++_iter  )
                {
                    if (_iter->mark() >= +0 &&
                        _iter->self() >= +1 )
                        _tnum +=  +1 ;
                }

                if (_nnum > +0)
                {
            /*-------------------------- write POINT data */
                save_vertS(
                    _ffun._ellipsoid_mesh_3d._mesh);
                }

                if (_nnum > +0)
                {
            /*-------------------------- write VALUE data */
                _file << "VALUE="
                      << _nnum << ";1" << "\n" ;

                save_vals(_file, +0.f,
                    _ffun._ellipsoid_mesh_3d._hval);
                }

                if (_tnum > +0)
                {
            /*-------------------------- write TRIA3 data */
                save_tria3(
                    _ffun._ellipsoid_mesh_3d._mesh);
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

            /*-------------------------- write class data */

                if (_ffun._euclidean_grid_2d.
                        _xpos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=1;" <<
                    _ffun._euclidean_grid_2d.
                        _xpos.count() << "\n" ;

                save_vals(_file, _xoff[0],
                    _ffun._euclidean_grid_2d._xpos);
                }

                if (_ffun._euclidean_grid_2d.
                        _ypos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=2;" <<
                    _ffun._euclidean_grid_2d.
                        _ypos.count() << "\n" ;

                save_vals(_file, _xoff[1],
                    _ffun._euclidean_grid_2d._ypos);
                }

                if (_ffun._euclidean_grid_2d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */
                _file << "VALUE=" <<
                    _ffun._euclidean_grid_2d.
                        _hmat.count() << ";1" "\n" ;

                save_vals(_file, +0.f ,
                    _ffun._euclidean_grid_2d._hmat);
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

            /*-------------------------- write class data */

                if (_ffun._euclidean_grid_3d.
                        _xpos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=1;" <<
                    _ffun._euclidean_grid_3d.
                        _xpos.count() << "\n" ;

                save_vals(_file, _xoff[0],
                    _ffun._euclidean_grid_3d._xpos);
                }

                if (_ffun._euclidean_grid_3d.
                        _ypos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=2;" <<
                    _ffun._euclidean_grid_3d.
                        _ypos.count() << "\n" ;

                save_vals(_file, _xoff[1],
                    _ffun._euclidean_grid_3d._ypos);
                }

                if (_ffun._euclidean_grid_3d.
                        _zpos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=3;" <<
                    _ffun._euclidean_grid_3d.
                        _zpos.count() << "\n" ;

                save_vals(_file, _xoff[2],
                    _ffun._euclidean_grid_3d._zpos);
                }

                if (_ffun._euclidean_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */
                _file << "VALUE=" <<
                    _ffun._euclidean_grid_3d.
                        _hmat.count() << ";1" "\n" ;

                save_vals(_file, +0.f ,
                    _ffun._euclidean_grid_3d._hmat);
                }

            }
            else
            if (_ffun._ndim == +2 &&
                _ffun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
            /*-------------------------- save 2-dim. grid */
                _file << "# " << _name << ".msh"
                      << "; created by " ;
                _file << __JGSWVSTR "\n" ;
                _file << "MSHID=3;ELLIPSOID-GRID \n" ;
                _file << "NDIMS=2 \n" ;

            /*-------------------------- write class data */
                {
                 char _fbuf [PRINTCHUNK] ;
                 auto _char = snprintf(&_fbuf[ 0 ] ,
                    PRINTCHUNK,
                    "%.17g;%.17g;%.17g\n",
                _ffun._ellipsoid_grid_3d._radA ,
                _ffun._ellipsoid_grid_3d._radB ,
                _ffun._ellipsoid_grid_3d._radC) ;

                _file << "RADII="
                      << std::string(_fbuf, _char) ;
                }

                if (_ffun._ellipsoid_grid_3d.
                        _xpos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=1;" <<
                    _ffun._ellipsoid_grid_3d.
                        _xpos.count() << "\n" ;

                save_vals(_file, +0.0,
                    _ffun._ellipsoid_grid_3d._xpos);
                }

                if (_ffun._ellipsoid_grid_3d.
                        _ypos.count() > +0)
                {
            /*-------------------------- write COORD data */
                _file << "COORD=2;" <<
                    _ffun._ellipsoid_grid_3d.
                        _ypos.count() << "\n" ;

                save_vals(_file, +0.0,
                    _ffun._ellipsoid_grid_3d._ypos);
                }

                if (_ffun._ellipsoid_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */
                _file << "VALUE=" <<
                    _ffun._ellipsoid_grid_3d.
                        _hmat.count() << ";1" "\n" ;

                save_vals(_file, +0.f ,
                    _ffun._ellipsoid_grid_3d._hmat);
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
    #undef save_vertS
    #undef save_tria3
    #undef save_tria4

    /*
    --------------------------------------------------------
     * SAVE-MESH: save *.NC output file.
    --------------------------------------------------------
     */

    #ifdef __use_netcdf

    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float*   /*_xoff*/ ,
        hfun_data &_ffun ,
        _netcdf_tag const&
        )
    {
        iptr_type _errv  = __no_error  ;

    //  NB. this routine handles output for MARCHE --- only
    //  *overwriting* the "value" arrays in *.NC - all else
    //  is passed read-only...

        try
        {
            std::string _path, _name, _fext;
            file_part(
                _jcfg._hfun_file, _path, _name, _fext) ;

            int _retv, _ncid;
            if ((_retv = nc_create(
                    _jcfg._hfun_file.c_str(),
                        NC_CLOBBER|NC_NETCDF4, &_ncid)))
            {
                return __file_not_created;
            }

            __unreferenced (_jcfg) ;
            __unreferenced (_jlog) ;

            if (_ffun._ndim == +2 &&    // save 2-dim. mesh
                _ffun._kind ==
                jmsh_kind::euclidean_mesh)
            {
                ncutil::put_str(_jcfg._hfun_file,
                    "TITLE" ,
                _name + "; created by " + __JGSWVSTR);

                if (_ffun._euclidean_mesh_2d.
                        _hval.count() > +0)
                {
            /*-------------------------- write VALUE data */
                ncutil::put_var(
                    _jcfg._hfun_file, "value_mat",
                &_ffun._euclidean_mesh_2d._hval [0]) ;
                }

            }
            else
            if (_ffun._ndim == +3 &&    // save 3-dim. mesh
                _ffun._kind ==
                jmsh_kind::euclidean_mesh)
            {
                ncutil::put_str(_jcfg._hfun_file,
                    "TITLE" ,
                _name + "; created by " + __JGSWVSTR);

                if (_ffun._euclidean_mesh_3d.
                        _hval.count() > +0)
                {
            /*-------------------------- write VALUE data */
                ncutil::put_var(
                    _jcfg._hfun_file, "value_mat",
                &_ffun._euclidean_mesh_3d._hval [0]) ;
                }

            }
            else
            if (_ffun._ndim == +2 &&    // save 2-sph. mesh
                _ffun._kind ==
                jmsh_kind::ellipsoid_mesh)
            {
                ncutil::put_str(_jcfg._hfun_file,
                    "TITLE" ,
                _name + "; created by " + __JGSWVSTR);

                if (_ffun._ellipsoid_mesh_3d.
                        _hval.count() > +0)
                {
            /*-------------------------- write VALUE data */
                ncutil::put_var(
                    _jcfg._hfun_file, "value_mat",
                &_ffun._ellipsoid_mesh_3d._hval [0]) ;
                }

            }
            else
            if (_ffun._ndim == +2 &&    // save 2-dim. grid
                _ffun._kind ==
                jmsh_kind::euclidean_grid)
            {
                ncutil::put_str(_jcfg._hfun_file,
                    "TITLE" ,
                _name + "; created by " + __JGSWVSTR);

                if (_ffun._euclidean_grid_2d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */
                ncutil::put_var(
                    _jcfg._hfun_file, "value_mat",
                &_ffun._euclidean_grid_2d._hmat [0]) ;
                }

            }
            else
            if (_ffun._ndim == +3 &&    // save 3-dim. grid
                _ffun._kind ==
                jmsh_kind::euclidean_grid)
            {
                ncutil::put_str(_jcfg._hfun_file,
                    "TITLE" ,
                _name + "; created by " + __JGSWVSTR);

                if (_ffun._euclidean_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */
                ncutil::put_var(
                    _jcfg._hfun_file, "value_mat",
                &_ffun._euclidean_grid_3d._hmat [0]) ;
                }

            }
            else
            if (_ffun._ndim == +2 &&    // save 2-sph. grid
                _ffun._kind ==
                jmsh_kind::ellipsoid_grid)
            {
                ncutil::put_str(_jcfg._hfun_file,
                    "TITLE" ,
                _name + "; created by " + __JGSWVSTR);

                if (_ffun._ellipsoid_grid_3d.
                        _hmat.count() > +0)
                {
            /*-------------------------- write VALUE data */                
                ncutil::put_var(
                    _jcfg._hfun_file, "value_mat",
                &_ffun._ellipsoid_grid_3d._hmat [0]) ;
                }

            }

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return ( _errv ) ;
    }

    #endif  //__use_netcdf

    /*
    --------------------------------------------------------
     * SAVE-MESH: save MSH_t output data.
    --------------------------------------------------------
     */

    template <
    typename      _dst_list ,
    typename      _src_list
             >
    __inline_call void_type copy_vals (
        _dst_list &_ldst,
        _src_list &_lsrc
        )
    {
    /*--------------------- copy msh_t REALS data to data */
        std::memcpy(&_ldst._data[0] ,
           &_lsrc[0], _ldst._size * sizeof(fp32_type)) ;
    }

    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float*   /*_xoff*/ ,
        hfun_data &_ffun ,
        jigsaw_msh_t &_fmsh
        )
    {
        iptr_type _errv  = __no_error ;

    //  NB. this routine handles output for MARCHE --- only
    //  *overwriting* the "value" arrays in FMSH - all else
    //  is passed read-only...

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

                copy_vals(_fmsh._value,
                    _ffun._euclidean_mesh_2d._hval ) ;
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

                copy_vals(_fmsh._value,
                    _ffun._euclidean_mesh_3d._hval ) ;
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

                copy_vals(_fmsh._value,
                    _ffun._ellipsoid_mesh_3d._hval ) ;
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

                copy_vals(_fmsh._value,
                    _ffun._euclidean_grid_2d._hmat ) ;
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

                copy_vals(_fmsh._value,
                    _ffun._euclidean_grid_3d._hmat ) ;
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

                copy_vals(_fmsh._value,
                    _ffun._ellipsoid_grid_3d._hmat ) ;
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
     * SAVE-MESH: save to *.MSH or *.NC
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type save_hfun (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        float     *_xoff ,
        hfun_data &_ffun
        )
    {
        iptr_type _errv  = __no_error ;

        try
        {
    /*------------------------- disambiguate via file-ext */
            std::string _path, _name, _fext;
            file_part(
                _jcfg._hfun_file, _path, _name, _fext) ;

            if(_fext.find("msh") != std::string::npos)
            {
                return save_hfun(_jcfg, _jlog,
                    _xoff, _ffun, mshfile_tag()) ;
            }
            else
            if(_fext.find( "nc") != std::string::npos)
            {
        #       ifdef __use_netcdf
                return save_hfun(_jcfg, _jlog,
                    _xoff, _ffun, _netcdf_tag()) ;
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



