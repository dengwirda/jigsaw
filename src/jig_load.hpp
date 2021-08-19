
    /*
    --------------------------------------------------------
     * JIG-LOAD: parse *.JIG file into JCFG data.
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
     * Last updated: 07 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   include "jig_read.hpp"

#   ifndef __JIG_LOAD__
#   define __JIG_LOAD__

    /*
    --------------------------------------------------------
     * READ-JCFG: read *.JCFG input file.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type read_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog
        )
    {
        class jcfg_loader:
            public jcfg_reader_base
        {
        public  :
        jcfg_data             *_jjig;

        public  :
    /*------------------------------------- create loader */
        __normal_call jcfg_loader (
            jcfg_data *_jsrc = nullptr
            ) : _jjig( _jsrc ) {}

    /*------------------------------------- MISC keywords */
        __normal_call void_type push_verbosity (
            std::int32_t  _verb
            )
        {
            this->_jjig->_verbosity = _verb ;
            this->_jjig->
           _mesh_opts.verb() = _verb;
            this->_jjig->
           _iter_opts.verb() = _verb;
        }

    /*------------------------------------- GEOM keywords */
        __normal_call void_type push_geom_file (
            std::string   _file
            )
        {
            this->
           _jjig->_geom_file = _file;
        }
        __normal_call void_type push_geom_seed (
            std::int32_t  _seed
            )
        {
            this->_jjig->
           _mesh_opts.seed() = _seed;
        }
        __normal_call void_type push_geom_feat (
            bool          _feat
            )
        {
            this->_jjig->
           _mesh_opts.feat() = _feat;
        }
        __normal_call void_type push_geom_phi1 (
            double        _phi1
            )
        {
            this->_jjig->
           _mesh_opts.phi1() = _phi1;
        }
        __normal_call void_type push_geom_phi2 (
            double        _phi2
            )
        {
            this->_jjig->
           _mesh_opts.phi2() = _phi2;
        }
        __normal_call void_type push_geom_eta1 (
            double        _eta1
            )
        {
            this->_jjig->
           _mesh_opts.eta1() = _eta1;
        }
        __normal_call void_type push_geom_eta2 (
            double        _eta2
            )
        {
            this->_jjig->
           _mesh_opts.eta2() = _eta2;
        }

    /*------------------------------------- INIT keywords */
        __normal_call void_type push_init_near (
            double        _near
            )
        {
            this->_jjig->
           _mesh_opts.near() = _near;
        }

    /*------------------------------------- HFUN keywords */
        __normal_call void_type push_hfun_file (
            std::string   _file
            )
        {
            this->
           _jjig->_hfun_file = _file;
        }
        __normal_call void_type push_hfun_scal (
            std::int32_t  _scal
            )
        {
            this->_jjig->_hfun_scal =
        (jcfg_data::hfun_scal::enum_data)_scal ;
        }

        __normal_call void_type push_hfun_hmax (
            double        _hmax
            )
        {
            this->
           _jjig->_hfun_hmax = _hmax;
        }
        __normal_call void_type push_hfun_hmin (
            double        _hmin
            )
        {
            this->
           _jjig->_hfun_hmin = _hmin;
        }

    /*------------------------------------- INIT keywords */
        __normal_call void_type push_init_file (
            std::string   _file
            )
        {
            this->
           _jjig->_init_file = _file;
        }

    /*------------------------------------- BNDS keywords */
        __normal_call void_type push_bnds_kern (
            std::int32_t  _kern
            )
        {
            this->_jjig->_bnds_pred =
        (jcfg_data::bnds_pred::enum_data)_kern ;
        }

    /*------------------------------------- MESH keywords */
        __normal_call void_type push_mesh_kern (
            std::int32_t  _kern
            )
        {
            this->_jjig->_mesh_pred =
        (jcfg_data::mesh_pred::enum_data)_kern ;
        }

        __normal_call void_type push_mesh_file (
            std::string   _file
            )
        {
            this->
           _jjig->_mesh_file = _file;
        }
        __normal_call void_type push_tria_file (
            std::string   _file
            )
        {
            this->
           _jjig->_tria_file = _file;
        }
        __normal_call void_type push_bnds_file (
            std::string   _file
            )
        {
            this->
           _jjig->_bnds_file = _file;
        }

        __normal_call void_type push_mesh_dims (
            std::int32_t  _dims
            )
        {
            this->_jjig->
           _mesh_opts.dims() = _dims;
        }
        __normal_call void_type push_mesh_iter (
            std::int32_t  _iter
            )
        {
            this->_jjig->
           _mesh_opts.iter() = _iter;
        }
        __normal_call void_type push_mesh_rule (
            std::int32_t  _rule
            )
        {
            this->_jjig->
           _mesh_opts.rule() = _rule;
        }
        __normal_call void_type push_mesh_siz1 (
            double        _siz1
            )
        {
            this->_jjig->
           _mesh_opts.siz1() = _siz1;
        }
        __normal_call void_type push_mesh_siz2 (
            double        _siz2
            )
        {
            this->_jjig->
           _mesh_opts.siz2() = _siz2;
        }
        __normal_call void_type push_mesh_siz3 (
            double        _siz3
            )
        {
            this->_jjig->
           _mesh_opts.siz3() = _siz3;
        }
        __normal_call void_type push_mesh_top1 (
            bool          _top1
            )
        {
            this->_jjig->
           _mesh_opts.top1() = _top1;
        }
        __normal_call void_type push_mesh_top2 (
            bool          _top2
            )
        {
            this->_jjig->
           _mesh_opts.top2() = _top2;
        }
        __normal_call void_type push_mesh_rad2 (
            double        _rad2
            )
        {
            this->_jjig->
           _mesh_opts.rad2() = _rad2;
        }
        __normal_call void_type push_mesh_rad3 (
            double        _rad3
            )
        {
            this->_jjig->
           _mesh_opts.rad3() = _rad3;
        }
        __normal_call void_type push_mesh_off2 (
            double        _off2
            )
        {
            this->_jjig->
           _mesh_opts.off2() = _off2;
        }
        __normal_call void_type push_mesh_off3 (
            double        _off3
            )
        {
            this->_jjig->
           _mesh_opts.off3() = _off3;
        }
        __normal_call void_type push_mesh_snk2 (
            double        _snk2
            )
        {
            this->_jjig->
           _mesh_opts.snk2() = _snk2;
        }
        __normal_call void_type push_mesh_snk3 (
            double        _snk3
            )
        {
            this->_jjig->
           _mesh_opts.snk3() = _snk3;
        }
        __normal_call void_type push_mesh_eps1 (
            double        _eps1
            )
        {
            this->_jjig->
           _mesh_opts.eps1() = _eps1;
        }
        __normal_call void_type push_mesh_eps2 (
            double        _eps2
            )
        {
            this->_jjig->
           _mesh_opts.eps2() = _eps2;
        }
        __normal_call void_type push_mesh_vol3 (
            double        _vol3
            )
        {
            this->_jjig->
           _mesh_opts.vol3() = _vol3;
        }

    /*------------------------------------- OPTM keywords */
        __normal_call void_type push_optm_kern (
            std::int32_t  _kern
            )
        {
            this->_jjig->_iter_pred =
        (jcfg_data::iter_pred::enum_data)_kern ;
        }

        __normal_call void_type push_optm_iter (
            std::int32_t  _iter
            )
        {
            this->_jjig->
           _iter_opts.iter() = _iter;
        }
        __normal_call void_type push_optm_qtol (
            double        _qtol
            )
        {
            this->_jjig->
           _iter_opts.qtol() = _qtol;
        }
        __normal_call void_type push_optm_qlim (
            double        _qlim
            )
        {
            this->_jjig->
           _iter_opts.qlim() = _qlim;
        }
        __normal_call void_type push_optm_tria (
            bool          _flag
            )
        {
            this->_jjig->
           _iter_opts.tria() = _flag;
        }
        __normal_call void_type push_optm_dual (
            bool          _flag
            )
        {
            this->_jjig->
           _iter_opts.dual() = _flag;
        }
        __normal_call void_type push_optm_div_ (
            bool          _flag
            )
        {
            this->_jjig->
           _iter_opts.div_() = _flag;
        }
        __normal_call void_type push_optm_zip_ (
            bool          _flag
            )
        {
            this->_jjig->
           _iter_opts.zip_() = _flag;
        }

        } ;

    /*---------------------------------- parse JCFG. file */
        iptr_type _errv  = __no_error ;

        try
        {
            jcfg_reader   _read;
            std::ifstream _file;
            _file. open(
            _jcfg._jcfg_file, std::ifstream::in) ;

            if (_file.is_open())
            {
                _read.read_file (
                _file, jcfg_loader(&_jcfg));
            }
            else
            {
                _jlog.push(
            "**parse error: file not found!\n" ) ;

                _errv = __file_not_located ;
            }
            _file.close ();

            for (auto _iter  =
                      _read._errs.head();
                      _iter !=
                      _read._errs.tend();
                    ++_iter  )
            {
                _jlog.push(
            "**parse error: " + * _iter + "\n" ) ;

                _errv = __invalid_useropts ;
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
     * COPY-JCFG: read *.JCFG input data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type copy_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        jigsaw_jig_t const&_jjig
        )
    {
        iptr_type _errv  = __no_error ;

        try
        {
            __unreferenced(_jlog) ;

    /*------------------------------------- MISC keywords */
            _jcfg._verbosity =
                         _jjig._verbosity ;
            _jcfg._mesh_opts.
                verb() = _jjig._verbosity ;
            _jcfg._iter_opts.
                verb() = _jjig._verbosity ;

    /*------------------------------------- BNDS keywords */
            if (_jjig._bnds_kern ==
                    JIGSAW_BNDS_TRIACELL)
            _jcfg._bnds_pred =
                jcfg_data::bnds_pred::bnd_tria ;
            else
            if (_jjig._bnds_kern ==
                    JIGSAW_BNDS_DUALCELL)
            _jcfg._bnds_pred =
                jcfg_data::bnds_pred::bnd_dual ;
            else
            _errv = __invalid_useropts ;

    /*------------------------------------- GEOM keywords */
            _jcfg._mesh_opts.
                seed() = _jjig._geom_seed ;
            _jcfg._mesh_opts.
                feat() = _jjig._geom_feat ;
            _jcfg._mesh_opts.
                eta1() = _jjig._geom_eta1 ;
            _jcfg._mesh_opts.
                eta2() = _jjig._geom_eta2 ;

    /*------------------------------------- INIT keywords */
            _jcfg._mesh_opts.
                near() = _jjig._init_near ;

    /*------------------------------------- HFUN keywords */
            if (_jjig._hfun_scal ==
                    JIGSAW_HFUN_RELATIVE)
            _jcfg._hfun_scal =
                jcfg_data::hfun_scal::relative ;
            else
            if (_jjig._hfun_scal ==
                    JIGSAW_HFUN_ABSOLUTE)
            _jcfg._hfun_scal =
                jcfg_data::hfun_scal::absolute ;
            else
            _errv = __invalid_useropts ;

            _jcfg.
            _hfun_hmax = _jjig._hfun_hmax ;
            _jcfg.
            _hfun_hmin = _jjig._hfun_hmin ;

    /*------------------------------------- RDEL keywords */
            if (_jjig._mesh_kern ==
                    JIGSAW_KERN_DELFRONT)
            _jcfg._mesh_pred =
                jcfg_data::mesh_pred::delfront ;
            else
            if (_jjig._mesh_kern ==
                    JIGSAW_KERN_DELAUNAY)
            _jcfg._mesh_pred =
                jcfg_data::mesh_pred::delaunay ;
            else
            if (_jjig._mesh_kern ==
                    JIGSAW_KERN_BISECTOR)
            _jcfg._mesh_pred =
                jcfg_data::mesh_pred::bisector ;
            else
            _errv = __invalid_useropts ;

            _jcfg._mesh_opts.
                dims() = _jjig._mesh_dims ;
            _jcfg._mesh_opts.
                iter() = _jjig._mesh_iter ;

            _jcfg._mesh_opts.
                siz1() = _jjig._mesh_siz1 ;
            _jcfg._mesh_opts.
                siz2() = _jjig._mesh_siz2 ;
            _jcfg._mesh_opts.
                siz3() = _jjig._mesh_siz3 ;

            _jcfg._mesh_opts.
                top1() = _jjig._mesh_top1 ;
            _jcfg._mesh_opts.
                top2() = _jjig._mesh_top2 ;

            _jcfg._mesh_opts.
                rad2() = _jjig._mesh_rad2 ;
            _jcfg._mesh_opts.
                rad3() = _jjig._mesh_rad3 ;
            _jcfg._mesh_opts.
                off2() = _jjig._mesh_off2 ;
            _jcfg._mesh_opts.
                off3() = _jjig._mesh_off3 ;
            _jcfg._mesh_opts.
                snk2() = _jjig._mesh_snk2 ;
            _jcfg._mesh_opts.
                snk3() = _jjig._mesh_snk3 ;

            _jcfg._mesh_opts.
                eps1() = _jjig._mesh_eps1 ;
            _jcfg._mesh_opts.
                eps2() = _jjig._mesh_eps2 ;

            _jcfg._mesh_opts.
                vol3() = _jjig._mesh_vol3 ;

    /*------------------------------------- ITER keywords */
            if (_jjig._optm_kern ==
                    JIGSAW_KERN_ODT_DQDX)
            _jcfg._iter_pred =
                jcfg_data::iter_pred::odt_dqdx ;
            else
            if (_jjig._optm_kern ==
                    JIGSAW_KERN_CVT_DQDX)
            _jcfg._iter_pred =
                jcfg_data::iter_pred::cvt_dqdx ;
            else
            if (_jjig._optm_kern ==
                    JIGSAW_KERN_H95_DQDX)
            _jcfg._iter_pred =
                jcfg_data::iter_pred::h95_dqdx ;
            else
            _errv = __invalid_useropts ;

            _jcfg._iter_opts.
                iter() = _jjig._optm_iter ;

            _jcfg._iter_opts.
                qtol() = _jjig._optm_qtol ;
            _jcfg._iter_opts.
                qlim() = _jjig._optm_qlim ;

            _jcfg._iter_opts.
                tria() = _jjig._optm_tria ;
            _jcfg._iter_opts.
                dual() = _jjig._optm_dual ;
            _jcfg._iter_opts.
                div_() = _jjig._optm_div_ ;
            _jcfg._iter_opts.
                zip_() = _jjig._optm_zip_ ;

        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * TEST-JCFG: test the JCFG data for validity.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type test_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog
        )
    {
        std::stringstream  _sstr;

    /*---------------------------------- read "ints" data */
        #define __testINTS(__tag, __var, __vlo, __vhi) \
            if ( (__var < __vlo) ||             \
                 (__var > __vhi) )              \
            {                                   \
                _errv = __invalid_useropts ;    \
                _sstr.str("");                  \
                _sstr.clear();                  \
                _sstr <<                        \
                "**input error: " __tag " = "   \
                      << std::setw(+9)          \
                      << std::setfill(' ')      \
                      << __var << "\n";         \
                _jlog.push(_sstr.str()) ;       \
            }                                   \

    /*---------------------------------- test "real" data */
        #define __testREAL(__tag, __var, __vlo, __vhi) \
            if ( (__var < __vlo) ||             \
                 (__var > __vhi) )              \
            {                                   \
                _errv = __invalid_useropts ;    \
                _sstr.str("");                  \
                _sstr.clear();                  \
                _sstr <<                        \
                "**input error: " __tag " = "   \
                      << std::scientific        \
                      << std::setprecision(2)   \
                      << __var << "\n";         \
                _jlog.push(_sstr.str()) ;       \
            }                                   \

    /*---------------------------------- warn "real" data */
        #define __warnREAL(__tag, __var, __vlo, __vhi) \
            if ( (__var < __vlo) ||             \
                 (__var > __vhi) )              \
            {                                   \
                _sstr.str("");                  \
                _sstr.clear();                  \
                _sstr <<                        \
              "**input warning: " __tag " = "   \
                      << std::scientific        \
                      << std::setprecision(2)   \
                      << __var << "\n";         \
                _jlog.push(_sstr.str()) ;       \
            }                                   \


        iptr_type _errv = __no_error ;

    /*---------------------------- test GEOM keywords */
        __testINTS("GEOM-SEED",
            _jcfg._mesh_opts.seed(),
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())

        __testREAL("GEOM-PHI1",
            _jcfg._mesh_opts.phi1(),
            (real_type)  0.,
            (real_type)180.)
        __testREAL("GEOM-PHI2",
            _jcfg._mesh_opts.phi2(),
            (real_type)  0.,
            (real_type)180.)

        __testREAL("GEOM-ETA1",
            _jcfg._mesh_opts.eta1(),
            (real_type)  0.,
            (real_type)180.)
        __testREAL("GEOM-ETA2",
            _jcfg._mesh_opts.eta2(),
            (real_type)  0.,
            (real_type)180.)

    /*---------------------------- test INIT keywords */
        __testREAL("INIT-NEAR",
            _jcfg._mesh_opts.near(),
            (real_type)  0.,
            (real_type)  1.)

    /*---------------------------- test HFUN keywords */
        __testREAL("HFUN-HMAX",
            _jcfg ._hfun_hmax ,
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())
        __testREAL("HFUN-HMIN",
            _jcfg ._hfun_hmin ,
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())

    /*---------------------------- test MESH keywords */
        __testINTS("MESH-ITER",
            _jcfg._mesh_opts.iter(),
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())

        __testINTS("MESH-DIMS",
            _jcfg._mesh_opts.dims(),
            (iptr_type) + 0,
            (iptr_type) + 3)

        __testREAL("MESH-SIZ1",
            _jcfg._mesh_opts.siz1(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-SIZ2",
            _jcfg._mesh_opts.siz2(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-SIZ3",
            _jcfg._mesh_opts.siz3(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-EPS1",
            _jcfg._mesh_opts.eps1(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-EPS2",
            _jcfg._mesh_opts.eps2(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-RAD2",
            _jcfg._mesh_opts.rad2(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-RAD3",
            _jcfg._mesh_opts.rad3(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())

        __warnREAL("MESH-RAD2",
            _jcfg._mesh_opts.rad2(),
            (real_type)  1.,
        std::numeric_limits<real_type>::infinity())
        __warnREAL("MESH-RAD3",
            _jcfg._mesh_opts.rad3(),
            (real_type)  2.,
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-OFF2",
            _jcfg._mesh_opts.off2(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())
        __testREAL("MESH-OFF3",
            _jcfg._mesh_opts.off3(),
            (real_type)  0.,
        std::numeric_limits<real_type>::infinity())

        __testREAL("MESH-SNK2",
            _jcfg._mesh_opts.snk2(),
            (real_type)  0.,
            (real_type)  1.)
        __testREAL("MESH-SNK3",
            _jcfg._mesh_opts.snk3(),
            (real_type)  0.,
            (real_type)  1.)

        __testREAL("MESH-VOL3",
            _jcfg._mesh_opts.vol3(),
            (real_type)  0.,
            (real_type)  1.)
        __warnREAL("MESH-VOL3",
            _jcfg._mesh_opts.vol3(),
            (real_type)  0.,
            (real_type)  .3)

    /*---------------------------- test OPTM keywords */
        __testINTS("OPTM-ITER",
            _jcfg._iter_opts.iter(),
            (iptr_type) + 0,
        std::numeric_limits<iptr_type>::     max())

        __testREAL("OPTM-QTOL",
            _jcfg._iter_opts.qtol(),
            (real_type)  0.,
            (real_type)  1.)

        __testREAL("OPTM-QLIM",
            _jcfg._iter_opts.qlim(),
            (real_type)  0.,
            (real_type)  1.)

        #undef  __testINTS
        #undef  __testREAL
        #undef  __warnREAL

        return (  _errv ) ;
    }

    /*
    --------------------------------------------------------
     * ECHO-JCFG: print a summary of the JCFG data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type echo_jcfg (
        jcfg_data &_jcfg ,
        jlog_data &_jlog
        )
    {
        iptr_type _errv  = __no_error ;

    /*---------------------------------- push "bool" data */
        #define __pushBVAL(__var)           \
        std::string(__var ? "TRUE" : "FALSE") + "\n"

    /*---------------------------------- push "ints" data */
        #define __pushIVAL(__var)           \
           (__var !=                        \
        std::numeric_limits<iptr_type>::max() ? \
            std::to_string(__var) : "MAXINT"\
            ) + "\n"

    /*---------------------------------- push "real" data */
        #define __pushRVAL(__var)           \
           (__var !=                        \
        std::numeric_limits<real_type>::max() ? \
        to_string_prec(__var, +2) : "MAXFLT"\
            ) + "\n"

        _jlog.push("  GEOM-FILE = " +
                    _jcfg._geom_file + "\n" ) ;
        _jlog.push("  MESH-FILE = " +
                    _jcfg._mesh_file + "\n" ) ;
        _jlog.push("  HFUN-FILE = " +
                    _jcfg._hfun_file + "\n" ) ;
        _jlog.push("  INIT-FILE = " +
                    _jcfg._init_file + "\n" ) ;
        _jlog.push("  TRIA-FILE = " +
                    _jcfg._tria_file + "\n" ) ;
        _jlog.push("  BNDS-FILE = " +
                    _jcfg._bnds_file + "\n" ) ;

        _jlog.push("\n") ;

        if (_jcfg._verbosity > +0)
        {
    /*---------------------------- push GEOM keywords */
        _jlog.push("  GEOM-SEED = " +
            __pushIVAL(_jcfg._mesh_opts.seed()));

        _jlog.push("  GEOM-PHI1 = " +
            __pushRVAL(_jcfg._mesh_opts.phi1()));
        _jlog.push("  GEOM-PHI2 = " +
            __pushRVAL(_jcfg._mesh_opts.phi2()));

        _jlog.push("  GEOM-ETA1 = " +
            __pushRVAL(_jcfg._mesh_opts.eta1()));
        _jlog.push("  GEOM-ETA2 = " +
            __pushRVAL(_jcfg._mesh_opts.eta2()));

        _jlog.push("  GEOM-FEAT = " +
            __pushBVAL(_jcfg._mesh_opts.feat()));

        _jlog.push("\n") ;

    /*---------------------------- push INIT keywords */
        _jlog.push("  INIT-NEAR = " +
            __pushRVAL(_jcfg._mesh_opts.near()));

        _jlog.push("\n") ;

    /*---------------------------- push HFUN keywords */
        if(_jcfg._hfun_scal ==
         jcfg_data::hfun_scal::absolute)
        _jlog.push("  HFUN-SCAL = ABSOLUTE \n") ;
        else
        if(_jcfg._hfun_scal ==
         jcfg_data::hfun_scal::relative)
        _jlog.push("  HFUN-SCAL = RELATIVE \n") ;

        _jlog.push("  HFUN-HMAX = " +
            __pushRVAL( _jcfg._hfun_hmax )) ;
        _jlog.push("  HFUN-HMIN = " +
            __pushRVAL( _jcfg._hfun_hmin )) ;

        _jlog.push("\n") ;

    /*---------------------------- push BNDS keywords */
        if(_jcfg._bnds_pred ==
         jcfg_data::bnds_pred::bnd_tria)
        _jlog.push("  BNDS-KERN = BND-TRIA \n") ;
        else
        if(_jcfg._bnds_pred ==
         jcfg_data::bnds_pred::bnd_dual)
        _jlog.push("  BNDS-KERN = BND-DUAL \n") ;

    /*---------------------------- push MESH keywords */
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::delaunay)
        _jlog.push("  MESH-KERN = DELAUNAY \n") ;
        else
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::delfront)
        _jlog.push("  MESH-KERN = DELFRONT \n") ;
        else
        if(_jcfg._mesh_pred ==
         jcfg_data::mesh_pred::bisector)
        _jlog.push("  MESH-KERN = BISECTOR \n") ;

        _jlog.push("  MESH-TOP1 = " +
            __pushBVAL(_jcfg._mesh_opts.top1()));
        _jlog.push("  MESH-TOP2 = " +
            __pushBVAL(_jcfg._mesh_opts.top2()));

        _jlog.push("  MESH-ITER = " +
            __pushIVAL(_jcfg._mesh_opts.iter()));
        _jlog.push("  MESH-DIMS = " +
            __pushIVAL(_jcfg._mesh_opts.dims()));

        _jlog.push("  MESH-SIZ1 = " +
            __pushRVAL(_jcfg._mesh_opts.siz1()));
        _jlog.push("  MESH-SIZ2 = " +
            __pushRVAL(_jcfg._mesh_opts.siz2()));
        _jlog.push("  MESH-SIZ2 = " +
            __pushRVAL(_jcfg._mesh_opts.siz3()));

        _jlog.push("  MESH-EPS1 = " +
            __pushRVAL(_jcfg._mesh_opts.eps1()));
        _jlog.push("  MESH-EPS2 = " +
            __pushRVAL(_jcfg._mesh_opts.eps2()));

        _jlog.push("  MESH-RAD2 = " +
            __pushRVAL(_jcfg._mesh_opts.rad2()));
        _jlog.push("  MESH-RAD3 = " +
            __pushRVAL(_jcfg._mesh_opts.rad3()));

        _jlog.push("  MESH-OFF2 = " +
            __pushRVAL(_jcfg._mesh_opts.off2()));
        _jlog.push("  MESH-OFF3 = " +
            __pushRVAL(_jcfg._mesh_opts.off3()));

        _jlog.push("  MESH-SNK2 = " +
            __pushRVAL(_jcfg._mesh_opts.snk2()));
        _jlog.push("  MESH-SNK3 = " +
            __pushRVAL(_jcfg._mesh_opts.snk3()));

        _jlog.push("  MESH-VOL3 = " +
            __pushRVAL(_jcfg._mesh_opts.vol3()));

        _jlog.push("\n") ;

    /*---------------------------- push OPTM keywords */
        if(_jcfg._iter_pred ==
         jcfg_data::iter_pred::odt_dqdx)
        _jlog.push("  OPTM-KERN = ODT+DQDX \n") ;
        else
        if(_jcfg._iter_pred ==
         jcfg_data::iter_pred::cvt_dqdx)
        _jlog.push("  OPTM-KERN = CVT+DQDX \n") ;
        else
        if(_jcfg._iter_pred ==
         jcfg_data::iter_pred::h95_dqdx)
        _jlog.push("  OPTM-KERN = H95+DQDX \n") ;

        _jlog.push("  OPTM-ITER = " +
            __pushIVAL(_jcfg._iter_opts.iter()));

        _jlog.push("  OPTM-QTOL = " +
            __pushRVAL(_jcfg._iter_opts.qtol()));
        _jlog.push("  OPTM-QLIM = " +
            __pushRVAL(_jcfg._iter_opts.qlim()));

        _jlog.push("  OPTM-ZIP_ = " +
            __pushBVAL(_jcfg._iter_opts.zip_()));
        _jlog.push("  OPTM-DIV_ = " +
            __pushBVAL(_jcfg._iter_opts.div_()));
        _jlog.push("  OPTM-TRIA = " +
            __pushBVAL(_jcfg._iter_opts.tria()));
        _jlog.push("  OPTM-DUAL = " +
            __pushBVAL(_jcfg._iter_opts.dual()));

        _jlog.push("\n") ;

        }

        #undef  __pushRVAL
        #undef  __pushIVAL
        #undef  __pushBVAL

        return (  _errv) ;
    }


#   endif   //__JIG_LOAD__



