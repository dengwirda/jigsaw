    /*
    --------------------------------------------------------
     * MSH-READ-NC: parse NetCDF into mesh data.
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
     * Last updated: 12 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

     // from msh_read.hpp


    /*
    --------------------------------------------------------
     * MSH-READER-NC: read *.nc mesh files
    --------------------------------------------------------
     */

    class msh_reader_nc_
    {
    public  :

    typedef containers::array <
            std::string > string_tokens;

    string_tokens           _errs ;

    std::int32_t            _ftag = +3 ;
    std:: size_t            _ndim = +0 ;

    jmsh_kind::enum_data    _kind =
               jmsh_kind::euclidean_mesh ;

    private :

#   ifdef __use_netcdf

    std::string             _file ;

    __normal_call size_t checked_size (
        size_t const*_lens,
        int _rank
        ) const
    {
        size_t constexpr _imax =
            std::numeric_limits<size_t>::max();

        size_t _size = 1;
        for (int _idim = 0; _idim < _rank; ++_idim)
        {
            if (_lens[_idim] != 0 &&
                    _size > _imax / _lens[_idim])
                throw std::invalid_argument(
                "NetCDF variable shape overflow") ;

            _size *= _lens [_idim] ;
        }
        return _size ;
    }

    __normal_call bool has_var (
        std::string const&_name
        ) const
    {
        return ncutil::has_var(this->_file, _name) ;
    }

    __normal_call void get_var_len (
        std::string const&_name,
        int &_rank,
        size_t *_lens
        ) const
    {
        ncutil::get_var_len(this->_file,
            _name, _rank , _lens, NC_MAX_VAR_DIMS) ;
    }

    template <
        typename  data_type
             >
    __normal_call void read_var (
        std::string const&_name,
        size_t _size,
        containers::array<data_type> &_data
        ) const
    {
        _data.set_count(_size, containers::tight_alloc) ;
        if (_size == 0) return ;

        ncutil::get_var( this->_file, _name, &_data[0]) ;
    }

    __normal_call void read_header (
        )
    {
        if (!ncutil::has_att(this->_file, "MSHID"))
            throw std::invalid_argument(
                "Missing NetCDF MSHID attribute");

        std::string _kind_name ;
        ncutil::get_str(this->_file, "MSHID", _kind_name) ;
        while (!_kind_name.empty() &&
                   _kind_name.back() == '\0')
            _kind_name.pop_back() ;

        _kind_name = trim(_kind_name) ;

        std::transform(_kind_name.begin(),
                       _kind_name.  end(),
                       _kind_name.begin(),
        [](unsigned char c){ return
          (unsigned char)::toupper(c); } ) ;

        if (_kind_name == "EUCLIDEAN-MESH")
            this->_kind = jmsh_kind::euclidean_mesh ;
        else
        if (_kind_name == "EUCLIDEAN-GRID")
            this->_kind = jmsh_kind::euclidean_grid ;
        else
        if (_kind_name == "EUCLIDEAN-DUAL")
            this->_kind = jmsh_kind::euclidean_dual ;
        else
        if (_kind_name == "ELLIPSOID-MESH")
            this->_kind = jmsh_kind::ellipsoid_mesh ;
        else
        if (_kind_name == "ELLIPSOID-GRID")
            this->_kind = jmsh_kind::ellipsoid_grid ;
        else
        if (_kind_name == "ELLIPSOID-DUAL")
            this->_kind = jmsh_kind::ellipsoid_dual ;
        else
            throw std::invalid_argument(
                "Invalid NetCDF MSHID attribute");

        if (ncutil::has_att(this->_file, "NDIMS"))
        {
            int _ndim_val = +0 ;
            ncutil::get_att(
                this->_file, "NDIMS", &_ndim_val) ;
            if (_ndim_val < +1 || _ndim_val > +3)
                throw std::invalid_argument(
                    "Invalid NetCDF NDIMS attribute") ;

            this->_ndim = (size_t)_ndim_val ;
        }
        else
        {
            if (has_var("xgrid_val")) this->_ndim = 1 ;
            if (has_var("ygrid_val")) this->_ndim = 2 ;
            if (has_var("zgrid_val")) this->_ndim = 3 ;

            if (has_var("vert2_pos")) this->_ndim = 2 ;
            if (has_var("seed2_pos")) this->_ndim = 2 ;
            if (has_var("vert3_pos")) this->_ndim = 3 ;
            if (has_var("seed3_pos")) this->_ndim = 3 ;

            if (this->_ndim <= 0)
                throw std::invalid_argument(
                    "Missing NetCDF NDIMS attribute") ;
        }
    }

    template <
        typename  dest_type
             >
    __normal_call void read_radii (
        dest_type &_dest
        )
    {
        if(!ncutil::has_att(this->_file, "RADII"))
            return ;

        size_t _nrad = +0 ;
        ncutil::get_att_len(this->_file, "RADII",
                           _nrad) ;
        if (_nrad != +3)
            throw std::invalid_argument(
                "Invalid NetCDF RADII attribute");

        double _radii[3] = {1., 1., 1.} ;
        ncutil::get_att(this->_file, "RADII", _radii) ;
        _dest.push_radii (_radii) ;
    }

    template <
        typename  dest_type,
        typename  open_type,
        typename  push_type
             >
    __normal_call void read_point (
        std::string const&_pos_name,
        std::string const&_tag_name,
        dest_type &_dest,
        open_type  _open,
        push_type  _push
        )
    {
        bool _has_pos = has_var(_pos_name) ;
        bool _has_tag = has_var(_tag_name) ;
        if (!_has_pos && !_has_tag)
            return ;
        if (!_has_pos || !_has_tag)
            throw std::invalid_argument(
            "Missing NetCDF " + _tag_name + " variable") ;

        int _rank = +0 ;
        size_t _lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_pos_name, _rank, _lens) ;
        if (_rank != +2 || _lens[1] != this->_ndim)
            throw std::invalid_argument(
            "Invalid NetCDF " + _pos_name + " dimensions") ;

        size_t _nrow = _lens[0] ;
        int _tag_rank = +0 ;
        size_t _tag_lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_tag_name, _tag_rank, _tag_lens) ;
        if ((_tag_rank != +1 && _tag_rank != +2) ||
            _tag_lens[0] != _nrow ||
            (_tag_rank == +2 && _tag_lens[1] != +1))
            throw std::invalid_argument(
            "Invalid NetCDF " + _tag_name + " dimensions") ;

        size_t _size = checked_size(_lens, _rank) ;

        containers::array<double> _pos ;
        containers::array<std::int32_t> _tag ;
        read_var(_pos_name, _size, _pos) ;
        read_var(_tag_name, _nrow, _tag) ;

        _open(_dest, _nrow) ;
        for (size_t _irow = +0; _irow < _nrow; ++_irow)
            _push(_dest, _irow,
                &_pos[_irow * this->_ndim], _tag[_irow]
                ) ;
    }

    template <
        typename  dest_type,
        typename  open_type,
        typename  push_type
             >
    __normal_call void read_coord (
        size_t _idim,
        dest_type &_dest,
        open_type  _open,
        push_type  _push
        )
    {
        static char const*_names[3] = {
            "xgrid_val", "ygrid_val", "zgrid_val"
        } ;

        std::string _name = _names[_idim - 1] ;
        if (!has_var(_name))
            return ;
        if (_idim > this->_ndim)
            throw std::invalid_argument(
                "NetCDF coord dimension exceeds NDIMS");

        int _rank = +0 ;
        size_t _lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_name, _rank, _lens) ;
        if (_rank != +1 && _rank != +2)
            throw std::invalid_argument(
            "Invalid NetCDF " + _name + " dimensions");
        if (_rank == +2 && _lens[1] != +1)
            throw std::invalid_argument(
            "Invalid NetCDF " + _name + " dimensions");

        size_t _nrow = _lens[0] ;
        containers::array<double> _data ;
        read_var(_name, _nrow, _data) ;

        _open(_dest, _nrow, _idim) ;
        for (size_t _irow = +0; _irow < _nrow; ++_irow)
            _push(_dest,
                _idim, _irow, _data[_irow]) ;
    }

    template <
        typename  data_type,
        typename  dest_type,
        typename  open_type,
        typename  push_type
             >
    __normal_call void read_value (
        std::string const&_name,
        bool _flatten,
        dest_type &_dest,
        open_type  _open,
        push_type  _push
        )
    {
        if (!has_var(_name)) return ;

        int _rank = +0 ;
        size_t _lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_name, _rank, _lens) ;
        if (_rank < +1 || (_rank > +2 && !_flatten) ||
            (_rank == +2 && _lens[1] == +0))
            throw std::invalid_argument(
            "Invalid NetCDF " + _name + " dimensions") ;

        size_t _size = checked_size(_lens, _rank) ;
        size_t _nrow = _lens[0] ;
        size_t _nval = _rank == +1 ? +1 : _lens[1] ;
        if (_flatten &&
            (this->_kind == jmsh_kind::euclidean_grid ||
             this->_kind == jmsh_kind::ellipsoid_grid ||
             _rank > +2))
        {
            _nrow = _size ;
            _nval = +1 ;
        }
        containers::array<data_type> _data ;
        read_var(_name, _size, _data) ;

        _open(_dest, _nrow, _nval) ;
        for (size_t _irow = +0; _irow < _nrow; ++_irow)
            _push(_dest, _irow,
                 &_data[ _irow * _nval]) ;
    }

    template <
        typename  dest_type,
        typename  open_type,
        typename  push_type
             >
    __normal_call void read_elements (
        std::string const&_name,
        size_t _nnod,
        dest_type &_dest,
        open_type  _open,
        push_type  _push
        )
    {
        std::string _idx_name = _name + "_idx" ;
        std::string _tag_name = _name + "_tag" ;
        bool _has_idx = has_var(_idx_name) ;
        bool _has_tag = has_var(_tag_name) ;
        if (!_has_idx && !_has_tag)
            return ;
        if (!_has_idx || !_has_tag)
            throw std::invalid_argument(
            "Missing NetCDF " + _tag_name + " variable") ;

        int _rank = +0 ;
        size_t _lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_idx_name, _rank, _lens) ;
        if (_rank != +2 || _lens[1] != _nnod)
            throw std::invalid_argument(
            "Invalid NetCDF " + _idx_name + " dimensions") ;

        size_t _nrow = _lens[0] ;
        int _tag_rank = +0 ;
        size_t _tag_lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_tag_name, _tag_rank, _tag_lens) ;
        if ((_tag_rank != +1 && _tag_rank != +2) ||
            _tag_lens[0] != _nrow ||
            (_tag_rank == +2 && _tag_lens[1] != +1))
            throw std::invalid_argument(
            "Invalid NetCDF " + _tag_name + " dimensions") ;

        size_t _size = checked_size(_lens, _rank) ;

        containers::array<std::int32_t> _idx ;
        containers::array<std::int32_t> _tag ;
        read_var(_idx_name, _size, _idx) ;
        read_var(_tag_name, _nrow, _tag) ;

        _open(_dest, _nrow) ;
        std::int32_t _node[8] = {} ;
        for (size_t _irow = +0; _irow < _nrow; ++_irow)
        {
            for (size_t _inod = +0; _inod < _nnod ;
                      ++_inod )
                _node[_inod] = _idx[_irow * _nnod + _inod] ;

            _push(_dest, _irow, _node, _tag[_irow]) ;
        }
    }

    template <
        typename  dest_type
             >
    __normal_call void read_bound (
        std::string const&_name,
        dest_type &_dest
        )
    {
        if (!has_var(_name)) return ;

        int _rank = +0 ;
        size_t _lens[NC_MAX_VAR_DIMS] = {} ;
        get_var_len(_name, _rank, _lens) ;
        if (_rank != +2 || _lens[1] != +3)
            throw std::invalid_argument(
                "Invalid NetCDF bound dimensions") ;

        size_t _nrow = _lens[0] ;
        containers::array<std::int32_t> _data ;
        read_var(_name,
            checked_size(_lens, _rank), _data) ;
        _dest.open_bound(_nrow) ;
        for (size_t _irow = +0; _irow < _nrow; ++_irow)
            _dest.push_bound(_irow,
                _data[_irow * 3 + 0],
                _data[_irow * 3 + 1],
                _data[_irow * 3 + 2]) ;
    }

    template <
        typename  dest_type
             >
    __normal_call void read_all (
        dest_type &_dest
        )
    {
        if (this->_ndim >= +2)
        {
        read_point (
            this->_ndim == +2 ? "vert2_pos" :
                                "vert3_pos" ,
            this->_ndim == +2 ? "vert2_tag" :
                                "vert3_tag" , _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_point(_n); },
            [](dest_type &_d, size_t _i, double *_p,
               std::int32_t _t)
                { _d.push_point(_i, _p, _t); }) ;

        read_point (
            this->_ndim == +2 ? "seed2_pos" :
                                "seed3_pos" ,
            this->_ndim == +2 ? "seed2_tag" :
                                "seed3_tag" , _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_seeds(_n); },
            [](dest_type &_d, size_t _i, double *_p,
               std::int32_t _t)
                { _d.push_seeds(_i, _p, _t); }) ;
        }

        read_radii(_dest) ;

        read_coord(+1, _dest,
            [](dest_type &_d, size_t _n, size_t _i)
                { _d.open_coord(_n, _i); },
            [](dest_type &_d, size_t _i, size_t _r, 
                              double _v)
                { _d.push_coord(_i, _r, _v); }) ;
        read_coord(+2, _dest,
            [](dest_type &_d, size_t _n, size_t _i)
                { _d.open_coord(_n, _i); },
            [](dest_type &_d, size_t _i, size_t _r,
                              double _v)
                { _d.push_coord(_i, _r, _v); }) ;
        read_coord(+3, _dest,
            [](dest_type &_d, size_t _n, size_t _i)
                { _d.open_coord(_n, _i); },
            [](dest_type &_d, size_t _i, size_t _r, 
                              double _v)
                { _d.push_coord(_i, _r, _v); }) ;

        read_value<double>("power_val", false, _dest,
            [](dest_type &_d, size_t _n, size_t _v)
                { _d.open_power(_n, _v); },
            [](dest_type &_d, size_t _i, double *_v)
                { _d.push_power(_i, _v); }) ;
        read_value<fp32_t>("value_mat", true, _dest,
            [](dest_type &_d, size_t _n, size_t _v)
                { _d.open_value(_n, _v); },
            [](dest_type &_d, size_t _i, fp32_t*_v)
                { _d.push_value(_i, _v); }) ;
        read_value<fp32_t>("slope_mat", true, _dest,
            [](dest_type &_d, size_t _n, size_t _v)
                { _d.open_slope(_n, _v); },
            [](dest_type &_d, size_t _i, fp32_t*_v)
                { _d.push_slope(_i, _v); }) ;

        read_elements("edge2", +2, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_edge2(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_edge2(_i, _n, _t); }) ;
        read_elements("tria3", +3, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_tria3(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_tria3(_i, _n, _t); }) ;
        read_elements("quad4", +4, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_quad4(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_quad4(_i, _n, _t); }) ;
        read_elements("tria4", +4, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_tria4(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_tria4(_i, _n, _t); }) ;
        read_elements("hexa8", +8, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_hexa8(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_hexa8(_i, _n, _t); }) ;
        read_elements("wedg6", +6, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_wedg6(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_wedg6(_i, _n, _t); }) ;
        read_elements("pyra5", +5, _dest,
            [](dest_type &_d, size_t _n)
                { _d.open_pyra5(_n); },
            [](dest_type &_d, size_t _i, std::int32_t *_n,
               std::int32_t _t)
                { _d.push_pyra5(_i, _n, _t); }) ;

        read_bound("bound_idx", _dest) ;
    }

#   endif

    public  :

    template <
        typename  dest_type
             >
    __normal_call void_type read_file (
        std::string const&_path,
        dest_type &_dest
        )
    {
        this->_errs.clear() ;
        this->_ftag = +3 ;
        this->_ndim = +0 ;
#   ifdef __use_netcdf
        this->_file = _path ;
        try
        {
            _dest.initialise() ;
             read_header() ;
            _dest.push_mshid(this->_ftag, this->_kind) ;
            _dest.push_ndims(this->_ndim) ;
             read_all(_dest) ;
            _dest.postscript() ;
        }
        catch (std::exception const&_error)
        {
            this->_errs.push_tail(_error.what()) ;
        }
#   else
        __unreferenced(_path) ;
        __unreferenced(_dest) ;
        this->_errs.push_tail(
            "NetCDF support is not available") ;
#   endif
    }

    } ;



