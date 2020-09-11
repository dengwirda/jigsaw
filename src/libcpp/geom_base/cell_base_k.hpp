
    /*
    --------------------------------------------------------
     * CELL-BASE-K: operations on linear cells in R^k.
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
     * Last updated: 26 July, 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __CELL_BASE_K__
#   define __CELL_BASE_K__

    namespace geometry {

    /*
    --------------------------------------------------------
     * tria. area/volume/normals/etc.
    --------------------------------------------------------
     */

    template <
    typename      data_type
             >
    __inline_call data_type tria_area_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ev12[2];
        vector_2d(_p1, _p2, _ev12);

        data_type _ev13[2];
        vector_2d(_p1, _p3, _ev13);

        data_type _aval =
            _ev12[0] * _ev13[1] -
            _ev12[1] * _ev13[0] ;

        _aval *= (data_type)+.5 ;

        return ( _aval )  ;
    }

    template <
    typename      data_type
             >
    __inline_call data_type tria_area_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ev12[3];
        vector_3d(_p1, _p2, _ev12);

        data_type _ev13[3];
        vector_3d(_p1, _p3, _ev13);

        data_type  _avec[3] = {
        _ev12[1] * _ev13[2] -
        _ev12[2] * _ev13[1] ,
        _ev12[2] * _ev13[0] -
        _ev12[0] * _ev13[2] ,
        _ev12[0] * _ev13[1] -
        _ev12[1] * _ev13[0] } ;

        data_type _aval =
            geometry::length_3d(_avec) ;

        return (data_type)+.5 * _aval  ;
    }

    template <
    typename      data_type
             >
    __inline_call void_type quad_axes_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4,
    __write_ptr  (data_type) _x1,
    __write_ptr  (data_type) _x2
        ) // see VERDICT ref. manual
    {
        _x1[0] = (_p2[0] - _p1[0])
               + (_p4[0] - _p3[0]);
        _x1[1] = (_p2[1] - _p1[1])
               + (_p4[1] - _p3[1]);

        _x2[0] = (_p3[0] - _p2[0])
               + (_p4[0] - _p1[0]);
        _x2[1] = (_p3[1] - _p2[1])
               + (_p4[1] - _p1[1]);
    }

    template <
    typename      data_type
             >
    __inline_call void_type quad_axes_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4,
    __write_ptr  (data_type) _x1,
    __write_ptr  (data_type) _x2
        ) // see VERDICT ref. manual
    {
        _x1[0] = (_p2[0] - _p1[0])
               + (_p4[0] - _p3[0]);
        _x1[1] = (_p2[1] - _p1[1])
               + (_p4[1] - _p3[1]);
        _x1[2] = (_p2[2] - _p1[2])
               + (_p4[2] - _p3[2]);

        _x2[0] = (_p3[0] - _p2[0])
               + (_p4[0] - _p1[0]);
        _x2[1] = (_p3[1] - _p2[1])
               + (_p4[1] - _p1[1]);
        _x2[2] = (_p3[2] - _p2[2])
               + (_p4[2] - _p1[2]);
    }

    template <
    typename      data_type
             >
    __inline_call data_type quad_area_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        ) // see VERDICT ref. manual
    {
        data_type _a1 =
        tria_area_2d(_p1, _p2, _p3);
        data_type _a2 =
        tria_area_2d(_p1, _p3, _p4);
        data_type _a3 =
        tria_area_2d(_p1, _p2, _p4);
        data_type _a4 =
        tria_area_2d(_p2, _p3, _p4);

        return ((_a1+_a2) + (_a3+_a4)) / +2.;

        /*
    //  data_type _x1[2], _x2[2], _nc;
    //  quad_axes_2d(_p1, _p2, _p3, _p4, _x1, _x2);

        data_type _v1[2], _a1;
        data_type _v2[2], _a2;
        data_type _v3[2], _a3;
        data_type _v4[2], _a4;
        vector_2d(_p1, _p2, _v1);
        vector_2d(_p2, _p3, _v2);
        vector_2d(_p3, _p4, _v3);
        vector_2d(_p4, _p1, _v4);

        cross_2d (_v4, _v1, _a1);
        cross_2d (_v1, _v2, _a2);
        cross_2d (_v2, _v3, _a3);
        cross_2d (_v3, _v4, _a4);

        return (_a1 + _a2 + _a3 + _a4);
        */
    }

    template <
    typename      data_type
             >
    __inline_call data_type quad_area_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        ) // see VERDICT ref. manual
    {
        data_type _a1 =
        tria_area_3d(_p1, _p2, _p3);
        data_type _a2 =
        tria_area_3d(_p1, _p3, _p4);
        data_type _a3 =
        tria_area_3d(_p1, _p2, _p4);
        data_type _a4 =
        tria_area_3d(_p2, _p3, _p4);

        return ((_a1+_a2) + (_a3+_a4)) / +2.;

        /*
        data_type _x1[3], _x2[3], _nc[3];
        quad_axes_3d(_p1, _p2, _p3, _p4, _x1, _x2);

        data_type _v1[3], _n1[3];
        data_type _v2[3], _n2[3];
        data_type _v3[3], _n3[3];
        data_type _v4[3], _n4[3];
        vector_3d(_p1, _p2, _v1);
        vector_3d(_p2, _p3, _v2);
        vector_3d(_p3, _p4, _v3);
        vector_3d(_p4, _p1, _v4);

        cross_3d (_v4, _v1, _n1);
        cross_3d (_v1, _v2, _n2);
        cross_3d (_v2, _v3, _n3);
        cross_3d (_v3, _v4, _n4);

        cross_3d (_x1, _x2, _nc);

        data_type _lc = length_3d(_nc);

        data_type _a1 =
            dot_3d(_nc, _n1) / _lc;
        data_type _a2 =
            dot_3d(_nc, _n2) / _lc;
        data_type _a3 =
            dot_3d(_nc, _n3) / _lc;
        data_type _a4 =
            dot_3d(_nc, _n4) / _lc;

        return (_a1 + _a2 + _a3 + _a4);
        */
    }

    template <
    typename      data_type
             >
    __inline_call void_type line_norm_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __write_ptr  (data_type) _nv
         )
    {
        _nv[0] = _p1[1] - _p2[1];
        _nv[1] = _p2[0] - _p1[0];
    }

    template <
    typename      data_type
             >
    __inline_call void_type tria_norm_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __write_ptr  (data_type) _nv
         )
    {
        data_type _ev12[3], _ev13[3] ;
        vector_3d(_p1, _p2, _ev12);
        vector_3d(_p1, _p3, _ev13);

        _nv[0] = _ev12[1] * _ev13[2] -
                 _ev12[2] * _ev13[1] ;
        _nv[1] = _ev12[2] * _ev13[0] -
                 _ev12[0] * _ev13[2] ;
        _nv[2] = _ev12[0] * _ev13[1] -
                 _ev12[1] * _ev13[0] ;
    }

    template <
    typename      data_type
             >
    __inline_call void_type quad_norm_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4,
    __write_ptr  (data_type) _nv
         )
    {
        data_type _xx11[3], _xx22[3] ;
        quad_axes_3d(
        _p1, _p2, _p3, _p4, _xx11, _xx22);

        _nv[0] = _xx11[1] * _xx22[2] -
                 _xx11[2] * _xx22[1] ;
        _nv[1] = _xx11[2] * _xx22[0] -
                 _xx11[0] * _xx22[2] ;
        _nv[2] = _xx11[0] * _xx22[1] -
                 _xx11[1] * _xx22[0] ;
    }

    template <
    typename      data_type
             >
    __inline_call data_type tetra_vol_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4    // +ve if CCW from _p4
        )
    {
        data_type _v14[3], _v24[3],
                  _v34[3];
        vector_3d(_p1, _p4, _v14);
        vector_3d(_p2, _p4, _v24);
        vector_3d(_p3, _p4, _v34);

        data_type _vdet =
      + _v14[0] * (_v24[1] * _v34[2] -
                   _v24[2] * _v34[1] )
      - _v14[1] * (_v24[0] * _v34[2] -
                   _v24[2] * _v34[0] )
      + _v14[2] * (_v24[0] * _v34[1] -
                   _v24[1] * _v34[0]);

        return _vdet / (data_type)6. ;
    }

    /*
    --------------------------------------------------------
     * cell "quality" scores.
    --------------------------------------------------------
     */

    template <
    typename      data_type
             >
    __inline_call
        data_type tria_quality_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {   // "volume-length" metric, see Shewchuk

        data_type static // +4. * std::sqrt(3.)
            constexpr _mul =
       (data_type)+6.928203230275509 ;

        data_type _len1 =
            lensqr_2d(_p1, _p2) ;
        data_type _len2 =
            lensqr_2d(_p2, _p3) ;
        data_type _len3 =
            lensqr_2d(_p3, _p1) ;

        data_type _lbar =
            _len1+_len2+_len3 ;

        data_type _area =
        tria_area_2d(_p1, _p2, _p3) ;

        return _mul * _area / _lbar ;
    }

    template <
    typename      data_type
             >
    __inline_call
        data_type tria_quality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {   // "volume-length" metric, see Shewchuk

        data_type static // +4. * std::sqrt(3.)
            constexpr _mul =
       (data_type)+6.928203230275509 ;

        data_type _len1 =
            lensqr_3d(_p1, _p2) ;
        data_type _len2 =
            lensqr_3d(_p2, _p3) ;
        data_type _len3 =
            lensqr_3d(_p3, _p1) ;

        data_type _lbar =
            _len1+_len2+_len3 ;

        data_type _area =
        tria_area_3d(_p1, _p2, _p3) ;

        return _mul * _area / _lbar ;
    }

    template <
    typename      data_type
             >
    __inline_call
        data_type tria_quality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        )
    {   // "volume-length" metric, see Shewchuk

        data_type static // +6. * std::sqrt(2.)
            constexpr _scal =
       (data_type)+8.485281374238571 ;

        data_type _tvol = tetra_vol_3d (
            _p1, _p2, _p3, _p4);

        data_type _lrms =
            lensqr_3d(_p1, _p2) +
            lensqr_3d(_p2, _p3) +
            lensqr_3d(_p3, _p1) +
            lensqr_3d(_p1, _p4) +
            lensqr_3d(_p2, _p4) +
            lensqr_3d(_p3, _p4) ;

        _lrms =  _lrms/ (data_type)6.0 ;
        _lrms =
        std::pow(_lrms, (data_type)1.5);

        return _scal * _tvol / _lrms ;
    }

    template <
    typename      data_type
             >
    __inline_call
        data_type quad_quality_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        )
    {   // mean sub-tria "volume-length" metrics

        data_type static // +2. / std::sqrt(+3.)
            constexpr _scal =
       (data_type)+1.154700538379252 ;

        data_type _val1 =
        tria_quality_2d(_p1, _p2, _p3) ;

        data_type _val2 =
        tria_quality_2d(_p1, _p3, _p4) ;

        data_type _val3 =
        tria_quality_2d(_p1, _p2, _p4) ;

        data_type _val4 =
        tria_quality_2d(_p2, _p3, _p4) ;

        data_type _vmin =
        std::min(std::min(_val1, _val2),
                 std::min(_val3, _val4)) ;

        data_type _prod =
            _val1*_val2*_val3*_val4 ;

        data_type _mean = _scal *
        std::pow(std::abs(_prod), 1./4.) ;

        return
        (_vmin > 0.) ? +_mean : -_mean ;
    }

    template <
    typename      data_type
             >
    __inline_call
        data_type quad_quality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        )
    {   // mean sub-tria "volume-length" metrics

        data_type static // +2. / std::sqrt(+3.)
            constexpr _scal =
       (data_type)+1.154700538379252 ;

        data_type _val1 =
        tria_quality_3d(_p1, _p2, _p3) ;

        data_type _val2 =
        tria_quality_3d(_p1, _p3, _p4) ;

        data_type _val3 =
        tria_quality_3d(_p1, _p2, _p4) ;

        data_type _val4 =
        tria_quality_3d(_p2, _p3, _p4) ;

        data_type _vmin =
        std::min(std::min(_val1, _val2),
                 std::min(_val3, _val4)) ;

        data_type _prod =
            _val1*_val2*_val3*_val4 ;

        data_type _mean = _scal *
        std::pow(std::abs(_prod), 1./4.) ;

        return
        (_vmin > 0.) ? +_mean : -_mean ;
    }

    /*
    --------------------------------------------------------
     * dual "quality" scores.
    --------------------------------------------------------
     */

    template <
    typename      data_type
             >
    __normal_call
        data_type tria_duality_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ob[3];
        perp_ball_2d(_ob, _p1, _p2, _p3);

        data_type _o1[3], _o2[3], _o3[3];
        perp_ball_2d(_o1, _p1, _p2) ;
        perp_ball_2d(_o2, _p2, _p3) ;
        perp_ball_2d(_o3, _p3, _p1) ;

        data_type _mb[3];
        mass_ball_2d(_mb, _p1, _p2, _p3);

        data_type _m1[3], _m2[3], _m3[3];
        mass_ball_2d(_m1, _p1, _p2) ;
        mass_ball_2d(_m2, _p2, _p3) ;
        mass_ball_2d(_m3, _p3, _p1) ;

        data_type _lb =
        geometry::lensqr_2d(_ob, _mb) ;

        data_type _l1 =
        geometry::lensqr_2d(_o1, _m1) ;
        data_type _l2 =
        geometry::lensqr_2d(_o2, _m2) ;
        data_type _l3 =
        geometry::lensqr_2d(_o3, _m3) ;

        data_type _r1 = _m1[2] ;
        data_type _r2 = _m2[2] ;
        data_type _r3 = _m3[2] ;

        data_type _rb =             // chara.-length
       (_r1+_r2+_r3) / (data_type)+3. ;

        data_type _qb = _lb / _rb ;

        data_type _q1 = _l1 / _r1 ;
        data_type _q2 = _l2 / _r2 ;
        data_type _q3 = _l3 / _r3 ;

        data_type _qe =             // straight mean
       (_q1+_q2+_q3) / (data_type)+3. ;

        data_type _qq =
      ((data_type)+1.-.33) * _qb +
      ((data_type)+0.+.33) * _qe ;

        return (data_type)1.-_qq ;
    }

    template <
    typename      data_type
             >
    __normal_call
        data_type tria_duality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ob[4];
        perp_ball_3d(_ob, _p1, _p2, _p3);

        data_type _o1[4], _o2[4], _o3[4];
        perp_ball_3d(_o1, _p1, _p2) ;
        perp_ball_3d(_o2, _p2, _p3) ;
        perp_ball_3d(_o3, _p3, _p1) ;

        data_type _mb[4];
        mass_ball_3d(_mb, _p1, _p2, _p3);

        data_type _m1[4], _m2[4], _m3[4];
        mass_ball_3d(_m1, _p1, _p2) ;
        mass_ball_3d(_m2, _p2, _p3) ;
        mass_ball_3d(_m3, _p3, _p1) ;

        data_type _lb =
        geometry::lensqr_3d(_ob, _mb) ;

        data_type _l1 =
        geometry::lensqr_3d(_o1, _m1) ;
        data_type _l2 =
        geometry::lensqr_3d(_o2, _m2) ;
        data_type _l3 =
        geometry::lensqr_3d(_o3, _m3) ;

        data_type _r1 = _m1[3] ;
        data_type _r2 = _m2[3] ;
        data_type _r3 = _m3[3] ;

        data_type _rb =             // chara.-length
       (_r1+_r2+_r3) / (data_type)+3. ;

        data_type _qb = _lb / _rb ;

        data_type _q1 = _l1 / _r1 ;
        data_type _q2 = _l2 / _r2 ;
        data_type _q3 = _l3 / _r3 ;

        data_type _qe =             // straight mean
       (_q1+_q2+_q3) / (data_type)+3. ;

        data_type _qq =
      ((data_type)+1.-.33) * _qb +
      ((data_type)+0.+.33) * _qe ;

        return (data_type)1.-_qq ;
    }

    /*
    template <
    typename      data_type
             >
    __normal_call
        data_type tria_duality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        )
    {
        //!! to-do...

        return    _qq ;
    }
    */

    }

#   endif//__CELL_BASE_K__


