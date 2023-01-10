
    /*
    --------------------------------------------------------
     * NCUTIL: simple interface to the NetCDF library.
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
     * Last updated: 20 May, 2022
     *
     * Copyright 2022--
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   include <string>
#   include <stdexcept>
#   include <initializer_list>

#   include "netcdf.h"

#   pragma once

#   ifndef __NCUTIL__
#   define __NCUTIL__

    namespace ncutil {

    /*
    --------------------------------------------------------
     * GET/PUT-DIM: handlers for NetCDF dimensions
    --------------------------------------------------------
     */

    __inline_call void def_dim (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of dimension
        size_t _dlen  // length of dimension
        )
    {
        int _retv, _ncid, _dtag;

        if ((_retv = nc_open(_file.c_str(), NC_WRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for write: " + std::to_string(_retv));

        if ((_retv = nc_def_dim(
                _ncid, _name.c_str(), _dlen, &_dtag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting dimension " +
                    _name + " :" + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    __inline_call void get_dim (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of dimension
        size_t &_dlen  // length of dimension
        )
    {
        int _retv, _ncid, _dtag;

        if ((_retv = nc_open(_file.c_str(), NC_NOWRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for read.: " + std::to_string(_retv));

        if ((_retv = nc_inq_dimid(_ncid, _name.c_str(), &_dtag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting dimension " +
                    _name + " : " + std::to_string(_retv));
        }

        if ((_retv = nc_inq_dimlen(_ncid, _dtag, &_dlen)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting dimension " +
                    _name + " : " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    /*
    --------------------------------------------------------
     * GET/PUT-ATT: handlers for NetCDF attributes
    --------------------------------------------------------
     */

    template <
    typename      data_type
             >
    __inline_call void get_att (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of attribute
      __write_ptr(data_type) _vals,  // ptr to attr. value array
        int _vtag = NC_GLOBAL  // NetCDF varid; default global
        )
    {
        int _retv, _ncid;

        if ((_retv = nc_open(_file.c_str(), NC_NOWRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for read.: " + std::to_string(_retv));

        if ((_retv = nc_get_att(_ncid,
                _vtag, _name.c_str (), (void *) _vals)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    __inline_call void get_str (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of attribute
        std::string &_vals,  // attr. string
        int _vtag = NC_GLOBAL  // NetCDF varid; default global
        )
    {
        int _retv, _ncid;
        size_t _slen;

        if ((_retv = nc_open(_file.c_str(), NC_NOWRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for read.: " + std::to_string(_retv));

        if ((_retv = nc_inq_attlen(
                _ncid, _vtag, _name.c_str(), &_slen)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }

        _vals.resize(_slen + 1) ;  // note: init to null
        if ((_retv = nc_get_att_text(
                _ncid, _vtag, _name.c_str(), &_vals[0])))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    template <
    typename      data_type
             >
    __inline_call void put_att (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of attribute
        nc_type  _type,  // NetCDF data-type
        data_type const&_vals,  // single attr. value
        int _vtag = NC_GLOBAL  // NetCDF varid; default global
        )
    {
        int _retv, _ncid;

        if ((_retv = nc_open(_file.c_str(), NC_WRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for write: " + std::to_string(_retv));

        if ((_retv = nc_put_att(
                _ncid, _vtag,
                _name.c_str(), _type, +1, (void *) &_vals)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting attribute " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    __inline_call void put_str (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of attribute
        std::string const&_vals,  // attr. string
        int _vtag = NC_GLOBAL  // NetCDF varid; default global
        )
    {
        int _retv, _ncid;

        if ((_retv = nc_open(_file.c_str(), NC_WRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for read.: " + std::to_string(_retv));

        if ((_retv = nc_put_att_text(_ncid, _vtag,
                _name.c_str (), _vals.size(), _vals.c_str())))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    template <
    typename      data_type
             >
    __inline_call void put_att (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of attribute
        nc_type  _type,  // NetCDF data-type
        size_t   _nval,  // length of attr. value array
      __const_ptr(data_type) _vals,  // ptr to att. value array
        int _vtag = NC_GLOBAL  // NetCDF varid; default global
        )
    {
        int _retv, _ncid;

        if ((_retv = nc_open(_file.c_str(), NC_WRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for write: " + std::to_string(_retv));

        if ((_retv = nc_put_att(
                _ncid, _vtag ,
                _name.c_str(), _type, _nval, (void *) _vals)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting attribute " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    /*
    --------------------------------------------------------
     * GET/PUT-VAR: handlers for NetCDF variables
    --------------------------------------------------------
     */

    __inline_call void def_var (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of variable
        nc_type  _type,  // NetCDF data-type
        std::string const&_long,  // NetCDF long_name
    std::initializer_list<std::string> _dims  // dim. name list
        )
    {
        int _retv, _ncid, _vtag, _dtag[256];
        int _idim = +0;

        if (_dims.size() > +256)
            throw std::invalid_argument(
                "Error putting variable " +
                    _name + " : too many dimensions!");

        if ((_retv = nc_open(_file.c_str(), NC_WRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for write: " + std::to_string(_retv));

        for (auto _elem : _dims)
        {
        if ((_retv = nc_inq_dimid(
                _ncid, _elem.c_str(), &_dtag[ _idim++ ])))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting variable " +
                    _name + " : " + std::to_string(_retv));
        }
        }

        if ((_retv = nc_def_var(
                _ncid, _name.c_str(),
                    _type, _dims.size (), _dtag, &_vtag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting variable " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_put_att_text(_ncid, _vtag,
                "long_name", _long.size(), _long.c_str())))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting variable " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    template <
    typename      data_type
             >
    __inline_call void get_var (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of variable
      __write_ptr(data_type) _vals  // ptr to var. value array
        )
    {
        int _retv, _ncid, _vtag;

        if ((_retv = nc_open(_file.c_str(), NC_NOWRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for read.: " + std::to_string(_retv));

        if ((_retv = nc_inq_varid(_ncid, _name.c_str(), &_vtag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting variable " +
                    _name + " : " + std::to_string(_retv));
        }

        if ((_retv = nc_get_var(_ncid, _vtag, (void *) _vals)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting variable " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    template <
    typename      data_type
             >
    __inline_call void put_var (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of variable
      __const_ptr(data_type) _vals  // ptr to var. value array
        )
    {
        int _retv, _ncid, _vtag;

        if ((_retv = nc_open(_file.c_str(), NC_WRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for write: " + std::to_string(_retv));

        if ((_retv = nc_inq_varid(_ncid, _name.c_str(), &_vtag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting variable " +
                    _name + " : " + std::to_string(_retv));
        }

        if ((_retv = nc_put_var(_ncid, _vtag, (void *) _vals)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error putting variable " +
                    _name + ": " + std::to_string(_retv));
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    /*
    --------------------------------------------------------
     * GET/PUT-TAG: handlers for NetCDF ID tags.
    --------------------------------------------------------
     */

    int static const NC_VARIABLE  = 100;
    int static const NC_ATTRIBUTE = 101;
    int static const NC_DIMENSION = 102;

    __inline_call void get_tag (
        std::string const&_file,  // name of file
        std::string const&_name,  // name of dim./att./var.
        int _kind, int &_itag,  // kind and idtag of name
        int _vtag = NC_GLOBAL   // NetCDF varid; default global
        )
    {
        int _retv, _ncid;

        if ((_retv = nc_open(_file.c_str(), NC_NOWRITE, &_ncid)))
            throw std::invalid_argument(
                "Can't open " + _file +
                    " for read.: " + std::to_string(_retv));

        if (_kind == NC_DIMENSION)
        {
        if ((_retv = nc_inq_dimid(
                _ncid, _name.c_str(), &_itag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }
        }
        else
        if (_kind == NC_VARIABLE)
        {
        if ((_retv = nc_inq_varid(
                _ncid, _name.c_str(), &_itag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }
        }
        else
        if (_kind == NC_ATTRIBUTE)
        {
        if ((_retv = nc_inq_attid(
                _ncid, _vtag, _name.c_str(), &_itag)))
        {
            nc_close(_ncid) ;
            throw std::invalid_argument(
                "Error getting attribute " +
                    _name + ": " + std::to_string(_retv));
        }
        }

        if ((_retv = nc_close(_ncid)))
            throw std::invalid_argument(
                "Error handling " +
                    _file + " close: " + std::to_string(_retv));
    }

    }

#   endif   //__NCUTIL__



