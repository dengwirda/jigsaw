
    /*
    --------------------------------------------------------
     * TRIA-ELEM-K: operations on simplexes in R^k. 
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
     * Last updated: 06 June, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __TRIA_ELEM_K__
#   define __TRIA_ELEM_K__

    namespace geometry {

    /*
    --------------------------------------------------------
     * tria. area/volume/normals/etc. 
    --------------------------------------------------------
     */

    template <
    typename      data_type
             > 
    __normal_call data_type tria_area_2d (
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
    __normal_call data_type tria_area_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ev12[3], _ev13[3] ;
        vector_3d(_p1, _p2, _ev12);
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
     * tria. "quality" scores. 
    --------------------------------------------------------
     */

    template <
    typename      data_type
             > 
    __normal_call 
        data_type tria_quality_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        // mean of condition no. + gradient-error metrics
        // see Shewchuk
    
        // 4. * std::sqrt(3.)
        data_type static 
            constexpr _mulA = 
       (data_type)+6.928203230275509 ;
        
        // 4. / std::sqrt(3.)
        data_type static 
            constexpr _mulB = 
       (data_type)+2.309401076758503 ;
       
        data_type _len1 = 
            lensqr_2d(_p1, _p2) ;
        data_type _len2 = 
            lensqr_2d(_p2, _p3) ;
        data_type _len3 =  
            lensqr_2d(_p3, _p1) ;

        data_type _barA = 
            _len1+_len2+_len3 ;
            
        data_type _barB = 
            _len1*_len2*_len3 ;

        _barB = std::pow(
            _barB, (data_type)+1./3.);

        data_type _area = 
        tria_area_2d(_p1, _p2, _p3);

        data_type _scrA = 
            _mulA * _area / _barA ;
        
        data_type _scrB = 
            _mulB * _area / _barB ;

        return (data_type).67* _scrA +
               (data_type).33* _scrB ;
                
        /*
        // 4. * std::sqrt(3.)
        data_type static 
            constexpr _scal = 
       (data_type)+6.928203230275509 ;
    
        data_type _elen = 
            lensqr_2d(_p1, _p2) + 
            lensqr_2d(_p2, _p3) + 
            lensqr_2d(_p3, _p1) ;

        data_type _area = 
        tria_area_2d(_p1, _p2, _p3);

        return _scal * _area / _elen ;
        */
    }

    template <
    typename      data_type
             > 
    __normal_call 
        data_type tria_quality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        // mean of condition no. + gradient-error metrics
        // see Shewchuk
    
        // 4. * std::sqrt(3.)
        data_type static 
            constexpr _mulA = 
       (data_type)+6.928203230275509 ;
        
        // 4. / std::sqrt(3.)
        data_type static 
            constexpr _mulB = 
       (data_type)+2.309401076758503 ;
       
        data_type _len1 = 
            lensqr_3d(_p1, _p2) ;
        data_type _len2 = 
            lensqr_3d(_p2, _p3) ;
        data_type _len3 =  
            lensqr_3d(_p3, _p1) ;
            
        data_type _barA = 
            _len1+_len2+_len3 ;
            
        data_type _barB = 
            _len1*_len2*_len3 ;

        _barB = std::pow(
            _barB, (data_type)+1./3.);

        data_type _area = 
        tria_area_3d(_p1, _p2, _p3);
        
        data_type _scrA = 
            _mulA * _area / _barA ;
        
        data_type _scrB = 
            _mulB * _area / _barB ;

        return (data_type).67* _scrA +
               (data_type).33* _scrB ;
               
        /*
        // 4. * std::sqrt(3.)
        data_type static 
            constexpr _scal = 
       (data_type)+6.928203230275509 ;
    
        data_type _elen = 
            lensqr_3d(_p1, _p2) + 
            lensqr_3d(_p2, _p3) + 
            lensqr_3d(_p3, _p1) ;

        data_type _area = 
        tria_area_3d(_p1, _p2, _p3);

        return _scal * _area / _elen ;
        */
    }

    template <
    typename      data_type
             > 
    __normal_call 
        data_type tria_quality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        )
    {
        // 6. * std::sqrt(2.)
        data_type static 
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

    /*
    --------------------------------------------------------
     * voro. "quality" scores. 
    --------------------------------------------------------
     */
    
    template <
    typename      data_type
             > 
    __normal_call 
        data_type dual_quality_2d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ob[ +3] ;
        geometry::perp_ball_2d(_ob, 
            _p1 , _p2, _p3);
            
        data_type _o1[ +3] ;
        geometry::perp_ball_2d(_o1, 
            _p1 , _p2) ;
        data_type _o2[ +3] ;
        geometry::perp_ball_2d(_o2, 
            _p2 , _p3) ;
        data_type _o3[ +3] ;
        geometry::perp_ball_2d(_o3, 
            _p3 , _p1) ;
            
        data_type _mb[ +3] ;
        geometry::mass_ball_2d(_mb, 
            _p1 , _p2, _p3);
            
        data_type _m1[ +3] ;
        geometry::mass_ball_2d(_m1, 
            _p1 , _p2) ;
        data_type _m2[ +3] ;
        geometry::mass_ball_2d(_m2, 
            _p2 , _p3) ;       
        data_type _m3[ +3] ;
        geometry::mass_ball_2d(_m3, 
            _p3 , _p1) ;

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
              
        data_type _qb = 
       (data_type)+1. - _lb / _rb ;
       
        data_type _q1 = 
       (data_type)+1. - _l1 / _r1 ;
        data_type _q2 = 
       (data_type)+1. - _l2 / _r2 ;
        data_type _q3 = 
       (data_type)+1. - _l3 / _r3 ;
       
        data_type _qe =             // straight mean
       (_q1+_q2+_q3) / (data_type)+3. ;
         
        /*
        data_type _qe =             // harmonic mean
            (data_type) +1. / _q1 + 
            (data_type) +1. / _q2 + 
            (data_type) +1. / _q3 ;
        _qe = 
            (data_type) +3. / _qe ;
         */
        
        data_type _qq = 
            (data_type)+.67 * _qb + 
                (data_type)+.33 * _qe ;

        return    _qq ;
    }
    
    template <
    typename      data_type
             > 
    __normal_call 
        data_type dual_quality_3d (
    __const_ptr  (data_type) _p1,
    __const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        data_type _ob[ +4] ;
        geometry::perp_ball_3d(_ob, 
            _p1 , _p2, _p3);
            
        data_type _o1[ +4] ;
        geometry::perp_ball_3d(_o1, 
            _p1 , _p2) ;          
        data_type _o2[ +4] ;
        geometry::perp_ball_3d(_o2, 
            _p2 , _p3) ;
        data_type _o3[ +4] ;
        geometry::perp_ball_3d(_o3, 
            _p3 , _p1) ;
            
        data_type _mb[ +4] ;
        geometry::mass_ball_3d(_mb, 
            _p1 , _p2, _p3);
            
        data_type _m1[ +4] ;
        geometry::mass_ball_3d(_m1, 
            _p1 , _p2) ;
        data_type _m2[ +4] ;
        geometry::mass_ball_3d(_m2, 
            _p2 , _p3) ;
        data_type _m3[ +4] ;
        geometry::mass_ball_3d(_m3, 
            _p3 , _p1) ;

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
              
        data_type _qb = 
       (data_type)+1. - _lb / _rb ;
       
        data_type _q1 = 
       (data_type)+1. - _l1 / _r1 ;
        data_type _q2 = 
       (data_type)+1. - _l2 / _r2 ;
        data_type _q3 = 
       (data_type)+1. - _l3 / _r3 ;
       
        data_type _qe =             // straight mean
       (_q1+_q2+_q3) / (data_type)+3. ;
         
        /*
        data_type _qe =             // harmonic mean
            (data_type) +1. / _q1 + 
            (data_type) +1. / _q2 + 
            (data_type) +1. / _q3 ;
        _qe = 
            (data_type) +3. / _qe ;
         */
         
        data_type _qq = 
            (data_type)+.67 * _qb + 
                (data_type)+.33 * _qe ;

        return    _qq ;
    }
    
    /*
    template <
    typename      data_type
             > 
    __normal_call 
        data_type dual_quality_3d (
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

#   endif//__TRIA_ELEM_K__


