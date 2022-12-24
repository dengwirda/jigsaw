
    /*
    --------------------------------------------------------
     * MSH-READ: parse *.MSH file into mesh data.
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
     * Last updated: 08 Feb., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_READ__
#   define __MSH_READ__

    class jmsh_reader_base
    {
/*------------------------ inherit from base to read file */
    public  :
    __normal_call void_type initialise() {}
    __normal_call void_type postscript() {}
    __normal_call void_type push_mshid (
        std::int32_t  /*_ftag*/,
        jmsh_kind::enum_data  /*_kind*/
        ) { }
    __normal_call void_type push_ndims (
        std:: size_t  /*_ndim*/
        ) { }
    __normal_call void_type push_radii (
        double      * /*_erad*/
        ) { }
    __normal_call void_type open_point (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_point (
        std:: size_t  /*_irow*/,
        double     *  /*_pval*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_seeds (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_seeds (
        std:: size_t  /*_irow*/,
        double     *  /*_pval*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_power (
        std:: size_t  /*_nrow*/,
        std:: size_t  /*_npwr*/
        ) { }
    __normal_call void_type push_power (
        std:: size_t  /*_irow*/,
        double     *  /*_xpwr*/
        ) { }
    __normal_call void_type open_edge2 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_edge2 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_tria3 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_tria3 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_quad4 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_quad4 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_tria4 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_tria4 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_hexa8 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_hexa8 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_wedg6 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_wedg6 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_pyra5 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_pyra5 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_value (
        std:: size_t  /*_nrow*/,
        std:: size_t  /*_nval*/
        ) { }
    __normal_call void_type push_value (
        std:: size_t  /*_irow*/,
        float       * /*_vdat*/
        ) { }
    __normal_call void_type open_slope (
        std:: size_t  /*_nrow*/,
        std:: size_t  /*_nval*/
        ) { }
    __normal_call void_type push_slope (
        std:: size_t  /*_irow*/,
        float       * /*_vdat*/
        ) { }
    __normal_call void_type open_bound (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_bound (
        std:: size_t  /*_irow*/,
        std::int32_t  /*_itag*/,
        std::int32_t  /*_inum*/,
        std::int32_t  /*_kind*/
        ) { }
    __normal_call void_type open_coord (
        std:: size_t  /*_idim*/,
        std:: size_t  /*_ndim*/
        ) { }
    __normal_call void_type push_coord (
        std:: size_t  /*_idim*/,
        std:: size_t  /*_irow*/,
        double        /*_ppos*/
        ) { }
    } ;

    /*
    --------------------------------------------------------
     * JMSH-READER: read *.MSH mesh files
    --------------------------------------------------------
     */

    class jmsh_reader
    {
    public  :

    typedef containers::array <
            std::string > string_tokens;

    string_tokens           _errs ;

    std::int32_t            _ftag ;
    std:: size_t            _ndim ;

    jmsh_kind::enum_data    _kind ;

    public  :

    /*
    --------------------------------------------------------
     * READ-MSHID: read MSHID data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_mshid (
        std::ifstream &_ffid,
        string_tokens &_stok,
        dest_type     &_dest
        )
    {
        __unreferenced (_ffid) ;

        if (_stok.count() == +2)
        {
        this->_ftag = std::stol(_stok[1]) ;
        this->_kind =
                jmsh_kind::euclidean_mesh ;

        _dest.push_mshid(this->_ftag ,
                         this->_kind
                         ) ;
        }
        else
        if (_stok.count() == +3)
        {

        this->_ftag = std::stol(_stok[1]) ;

        if (_stok[2] == "EUCLIDEAN-MESH")
        {
            this->_kind =
                jmsh_kind::euclidean_mesh ;
        }
        else
        if (_stok[2] == "EUCLIDEAN-GRID")
        {
            this->_kind =
                jmsh_kind::euclidean_grid ;
        }
        else
        if (_stok[2] == "EUCLIDEAN-DUAL")
        {
            this->_kind =
                jmsh_kind::euclidean_dual ;
        }
        else
        if (_stok[2] == "ELLIPSOID-MESH")
        {
            this->_kind =
                jmsh_kind::ellipsoid_mesh ;
        }
        else
        if (_stok[2] == "ELLIPSOID-GRID")
        {
            this->_kind =
                jmsh_kind::ellipsoid_grid ;
        }
        else
        if (_stok[2] == "ELLIPSOID-DUAL")
        {
            this->_kind =
                jmsh_kind::ellipsoid_dual ;
        }
        else
        {
            this->_errs.
            push_tail("Invalid MSHID!") ;
        }

        _dest.push_mshid(this->_ftag ,
                         this->_kind
                         ) ;

        }
        else
        {
            this->_errs.
            push_tail("Invalid MSHID!") ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-NDIMS: read NDIMS data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_ndims (
        std::ifstream &_ffid,
        string_tokens &_stok,
        dest_type     &_dest
        )
    {
        __unreferenced (_ffid) ;

        if (_stok.count() == +2)
        {
            this->
           _ndim = std::stoi(_stok[ 1]);

           _dest.
            push_ndims (this->_ndim) ;
        }
        else
        {
            this->_errs.
            push_tail("Invalid NDIMS!");
        }
    }

    /*
    --------------------------------------------------------
     * READ-RADII: read RADII data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_radii (
        std::ifstream &_ffid,
        string_tokens &_stok,
        dest_type     &_dest
        )
    {
        __unreferenced (_ffid) ;

        double _erad[ +3] ;
        if (_stok.count() == +4)
        {
            _erad[0] = std::stod(_stok[1]) ;
            _erad[1] = std::stod(_stok[2]) ;
            _erad[2] = std::stod(_stok[3]) ;

            _dest.push_radii(_erad);
        }
        else
        {
            this->_errs.
            push_tail("Invalid RADII!");
        }
    }

    /*
    --------------------------------------------------------
     * READ-POINT: read POINT data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_point (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid POINT!");
        }

        _dest.open_point(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                std:: size_t static
                    constexpr _VMAX = +256 ;

                double _pval [_VMAX] ;
                char  *_next =(char*)_line.data() ;

                for (size_t _ipos  = +0;
                    _ipos < this->_ndim; ++_ipos)
                {
                    _pval[_ipos] =
                     std::strtod(_next, &_next) ;

                    _next += +1;
                }

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_point (
                   _irow, _pval, _itag) ;

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-SEEDS: read SEEDS data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_seeds (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid SEEDS!");
        }

        _dest.open_seeds(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                std:: size_t static
                    constexpr _VMAX = +256 ;

                double _pval [_VMAX] ;
                char  *_next =(char*)_line.data() ;

                for (size_t _ipos  = +0;
                    _ipos < this->_ndim; ++_ipos)
                {
                    _pval[_ipos] =
                     std::strtod(_next, &_next) ;

                    _next += +1;
                }

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_seeds (
                   _irow, _pval, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-COORD: read COORD data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_coord (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _idim = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +3)
        {
            _idim = std::stol(_stok[1]);
            _nrow = std::stol(_stok[2]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid COORD!");
        }

        _dest.open_coord(_idim, _nrow) ;

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                _dest.push_coord ( _idim,
                    _irow, std::stod(_line)) ;
            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-POWER: read POWER data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_power (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _npwr = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +3)
        {
            _nrow = std::stol(_stok[1]);
            _npwr = std::stol(_stok[2]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid POWER!");
        }

        _dest.open_power(_nrow , _npwr);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                std:: size_t static
                    constexpr _VMAX = +256 ;

                real_t _vpwr [_VMAX] ;
                char  *_next =(char*)_line.data() ;

                for (size_t _ipos  = +0;
                    _ipos < _npwr; ++_ipos)
                {
                    _vpwr[_ipos] =
                     std::strtod(_next, &_next) ;

                    _next += +1 ;
                }

                _dest.push_power(_irow,  _vpwr) ;

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-VALUE: read VALUE data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_value (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _nval = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +3)
        {
            _nrow = std::stol(_stok[1]);
            _nval = std::stol(_stok[2]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid VALUE!");
        }

        _dest.open_value (_nrow, _nval);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                std:: size_t static
                    constexpr _VMAX = +256 ;

                fp32_t _vval [_VMAX] ;
                char  *_next =(char*)_line.data() ;

                for (size_t _ipos  = +0;
                    _ipos < _nval; ++_ipos)
                {
                    _vval[_ipos] =
                     std::strtof(_next, &_next) ;

                    _next += +1 ;
                }

                _dest.push_value(_irow,  _vval) ;

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-SLOPE: read SLOPE data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_slope (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _nval = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +3)
        {
            _nrow = std::stol(_stok[1]);
            _nval = std::stol(_stok[2]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid SLOPE!");
        }

        _dest.open_slope (_nrow, _nval);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                std:: size_t static
                    constexpr _VMAX = +256 ;

                fp32_t _vval [_VMAX] ;
                char  *_next =(char*)_line.data() ;

                for (size_t _ipos  = +0;
                    _ipos < _nval; ++_ipos)
                {
                    _vval[_ipos] =
                     std::strtof(_next, &_next) ;

                    _next += +1 ;
                }

                _dest.push_slope(_irow,  _vval) ;

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-EDGE2: read EDGE2 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_edge2 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid EDGE2!");
        }

        _dest.open_edge2(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[2];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_edge2 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-TRIA3: read TRIA3 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_tria3 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid TRIA3!");
        }

        _dest.open_tria3(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[3];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[2] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_tria3 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-QUAD4: read QUAD4 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_quad4 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid QUAD4!");
        }

        _dest.open_quad4(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[4];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[2] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[3] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_quad4 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-TRIA4: read TRIA4 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_tria4 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid TRIA4!");
        }

        _dest.open_tria4(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[4];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[2] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[3] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_tria4 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-HEXA8: read HEXA8 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_hexa8 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid HEXA8!");
        }

        _dest.open_hexa8(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[8];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[2] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[3] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[4] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[5] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[6] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[7] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_hexa8 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-WEDG6: read WEDG6 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_wedg6 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid WEDG6!");
        }

        _dest.open_wedg6(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[6];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[2] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[3] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[4] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[5] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_wedg6 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-PYRA5: read PYRA5 data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_pyra5 (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid PYRA5!");
        }

        _dest.open_pyra5(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _node[5];
                _node[0] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[1] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[2] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[3] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;
                _node[4] =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;

                _dest.push_pyra5 (
                   _irow, _node, _itag);

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-BOUND: read BOUND data section
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_bound (
        std::ifstream&_ffid ,
        string_tokens&_stok ,
        dest_type    &_dest
        )
    {
    /*----------------------------------------- read head */
        std:: size_t _nrow = +0;
        std:: size_t _irow = +0;
        if (_stok.count() == +2)
        {
            _nrow = std::stol(_stok[1]);
        }
        else
        {
            this->_errs.
            push_tail("Invalid BOUND!");
        }

        _dest.open_bound(_nrow);

    /*----------------------------------------- read data */
        std::string _line ;
        while (std::getline(_ffid, _line))
        {
            try
            {
                char *_next =(char *)_line.data() ;

                std::int32_t _itag =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _inum =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                std::int32_t _KIND =
                std::strtol( _next, &_next, 10) ;
                _next += +1 ;

                _dest.push_bound (
                _irow,  _itag, _inum, _KIND) ;

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }

            _irow += +1 ;

            if (--_nrow == +0) break ;
        }
    }

    /*
    --------------------------------------------------------
     * READ-FILE: read *.MSH file into MESH
    --------------------------------------------------------
     */

    template <
        typename  dest_type
             >
    __normal_call void_type read_file (
        std::ifstream&_ffid ,
        dest_type   &&_dest
        )
    {
        this->_ftag = -1 ;
        this->_ndim = +0 ;

       _dest.initialise();

        std::string _line;
        while (std::getline(_ffid, _line))
        {
            _line = trim(_line) ;

            if (_line.size() <= 0) continue ;
            if (_line[ +0] == '#') continue ;

            try
            {
            std::transform(_line.begin() ,
                           _line.  end() ,
                           _line.begin() ,
            [](unsigned char c){ return
              (unsigned char)::toupper(c); } ) ;

            containers::
                array<std::string> _stok ;

            find_toks(_line, "=;", _stok);

            for (auto _iter  = _stok.head() ;
                      _iter != _stok.tend() ;
                    ++_iter  )
        /*---------------------------- trim on each token */
            *_iter = trim( *_iter ) ;

            if (_stok[0] == "MSHID")
                {
            read_mshid(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "NDIMS")
                {
            read_ndims(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "RADII")
                {
            read_radii(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "POINT")
                {
            read_point(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "SEEDS")
                {
            read_seeds(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "COORD")
                {
            read_coord(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "POWER")
                {
            read_power(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "VALUE")
                {
            read_value(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "SLOPE")
                {
            read_slope(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "EDGE2")
                {
            read_edge2(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "TRIA3")
                {
            read_tria3(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "QUAD4")
                {
            read_quad4(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "TRIA4")
                {
            read_tria4(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "HEXA8")
                {
            read_hexa8(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "WEDG6")
                {
            read_wedg6(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "PYRA5")
                {
            read_pyra5(_ffid, _stok, _dest) ;
                }
            else
            if (_stok[0] == "BOUND")
                {
            read_bound (_ffid, _stok, _dest) ;
                }

            }
            catch (...)
            {
                this->_errs.push_tail(_line) ;
            }
        }
        {
           _dest.postscript();
        }
    }

    } ;

#   endif   //__MSH_READ__



