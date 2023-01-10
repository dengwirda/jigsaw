
    /*
    --------------------------------------------------------
     * PREDICATE-k: robust geometric predicates in E^k.
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
     * Last updated: 30 Apr., 2020
     *
     * Copyright 2020--
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from predicate_k.hpp...


    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the unweighted "bisect" predicate in E^2.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect2d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect2d predicate, "exact" version */
        mp::expansion< 2 > _ab_xx_, _ab_yy_,
                           _ac_xx_, _ac_yy_,
                           _bc_xx_, _bc_yy_;
        mp::expansion< 4 > _tt_xx_, _tt_yy_;
        mp::expansion< 32> _absum_;

        _OK = true;

    /*----------------------------------- compute: d(p,q) */
        _ab_xx_.from_sub(_pa[0], _pb[0]);
        _ab_yy_.from_sub(_pa[1], _pb[1]);

        _ac_xx_.from_sub(_pa[0], _pc[0]);
        _ac_yy_.from_sub(_pa[1], _pc[1]);

        _bc_xx_.from_sub(_pb[0], _pc[0]);
        _bc_yy_.from_sub(_pb[1], _pc[1]);

        mp::expansion_add(_ac_xx_, _bc_xx_,
                          _tt_xx_) ;
        mp::expansion_add(_ac_yy_, _bc_yy_,
                          _tt_yy_) ;

        mp::expansion_dot(_ab_xx_, _tt_xx_,
                          _ab_yy_, _tt_yy_,
                          _absum_) ;

    /*----------------------------------- return signed d */
        return mp::expansion_est(_absum_) ;
    }

    __normal_call REAL_TYPE bisect2d_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect2d predicate, "bound" version */
        ia_flt    _acx, _acy,
                  _bcx, _bcy,
                  _abx, _aby, _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _abx.from_sub(_pa[0], _pb[0]) ;   // coord. diff.
        _aby.from_sub(_pa[1], _pb[1]) ;

        _acx.from_sub(_pa[0], _pc[0]) ;
        _acy.from_sub(_pa[1], _pc[1]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;

        _sgn = (_abx * (_acx + _bcx))
             + (_aby * (_acy + _bcy)) ;

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid () ) ;
    }

    __normal_call REAL_TYPE bisect2d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect2d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  5. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy;
        REAL_TYPE _bcx, _bcy;
        REAL_TYPE _acsqr, _bcsqr ;

        REAL_TYPE _sgn, _FT ;

        REAL_TYPE _ACSQR, _BCSQR ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;

        _acsqr = _acx * _acx +
                 _acy * _acy ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy ;

        _ACSQR = std::abs(_acsqr);
        _BCSQR = std::abs(_bcsqr);

        _FT  = _ACSQR + _BCSQR ;          // roundoff tol
        _FT *= _ER  ;

        _sgn = _acsqr - _bcsqr ;          // d_ab - d_bc

        _OK  = _sgn > +_FT || _sgn < -_FT ;

        return _sgn ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the weighted "bisect" predicate in E^2.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect2w_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect2w predicate, "exact" version */
        mp::expansion< 2 > _ab_xx_, _ab_yy_,
                           _ab_ww_,
                           _ac_xx_, _ac_yy_,
                           _bc_xx_, _bc_yy_;
        mp::expansion< 4 > _tt_xx_, _tt_yy_;
        mp::expansion< 32> _ttsum_;
        mp::expansion< 34> _absum_;

        _OK = true;

    /*----------------------------------- compute: d(p,q) */
        _ab_xx_.from_sub(_pa[0], _pb[0]);
        _ab_yy_.from_sub(_pa[1], _pb[1]);
        _ab_ww_.from_sub(_pa[2], _pb[2]);

        _ac_xx_.from_sub(_pa[0], _pc[0]);
        _ac_yy_.from_sub(_pa[1], _pc[1]);

        _bc_xx_.from_sub(_pb[0], _pc[0]);
        _bc_yy_.from_sub(_pb[1], _pc[1]);

        mp::expansion_add(_ac_xx_, _bc_xx_,
                          _tt_xx_) ;
        mp::expansion_add(_ac_yy_, _bc_yy_,
                          _tt_yy_) ;

        mp::expansion_dot(_ab_xx_, _tt_xx_,
                          _ab_yy_, _tt_yy_,
                          _ttsum_) ;

        mp::expansion_sub(_ttsum_, _ab_ww_,
                          _absum_) ;

    /*----------------------------------- return signed d */
        return mp::expansion_est(_absum_) ;
    }

    __normal_call REAL_TYPE bisect2w_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect2w predicate, "bound" version */
        ia_flt    _acx, _acy, _abw,
                  _bcx, _bcy,
                  _abx, _aby, _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _abx.from_sub(_pa[0], _pb[0]) ;   // coord. diff.
        _aby.from_sub(_pa[1], _pb[1]) ;
        _abw.from_sub(_pa[2], _pb[2]) ;

        _acx.from_sub(_pa[0], _pc[0]) ;
        _acy.from_sub(_pa[1], _pc[1]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;

        _sgn = (_abx * (_acx + _bcx))
             + (_aby * (_acy + _bcy)) ;

        _sgn-=  _abw ;

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid () ) ;
    }

    __normal_call REAL_TYPE bisect2w_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect2w predicate, "float" version */
        REAL_TYPE static const _ER =
        +  6. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy;
        REAL_TYPE _bcx, _bcy;
        REAL_TYPE _acsqr, _bcsqr ;
        REAL_TYPE _a_sum, _b_sum ;

        REAL_TYPE _A_SUM, _B_SUM ;

        REAL_TYPE _sgn, _FT;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;

        _acsqr = _acx * _acx +
                 _acy * _acy ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy ;

        _a_sum = _acsqr - _pa[2] ;
        _b_sum = _bcsqr - _pb[2] ;

        _A_SUM = std::abs(_acsqr)
               + std::abs(_pa[2]);
        _B_SUM = std::abs(_bcsqr)
               + std::abs(_pb[2]);

        _FT  = _A_SUM + _B_SUM ;          // roundoff tol
        _FT *= _ER ;

        _sgn = _a_sum - _b_sum ;          // d_ab - d_bc

        _OK  = _sgn > +_FT || _sgn < -_FT ;

        return _sgn ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the unweighted "bisect" predicate in E^3.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect3d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect3d predicate, "exact" version */
        mp::expansion< 2 > _ab_xx_, _ab_yy_,
                           _ab_zz_,
                           _ac_xx_, _ac_yy_,
                           _ac_zz_,
                           _bc_xx_, _bc_yy_,
                           _bc_zz_;
        mp::expansion< 4 > _tt_xx_, _tt_yy_,
                           _tt_zz_;
        mp::expansion< 48> _absum_;

        _OK = true;

    /*----------------------------------- compute: d(p,q) */
        _ab_xx_.from_sub(_pa[0], _pb[0]);
        _ab_yy_.from_sub(_pa[1], _pb[1]);
        _ab_zz_.from_sub(_pa[2], _pb[2]);

        _ac_xx_.from_sub(_pa[0], _pc[0]);
        _ac_yy_.from_sub(_pa[1], _pc[1]);
        _ac_zz_.from_sub(_pa[2], _pc[2]);

        _bc_xx_.from_sub(_pb[0], _pc[0]);
        _bc_yy_.from_sub(_pb[1], _pc[1]);
        _bc_zz_.from_sub(_pb[2], _pc[2]);

        mp::expansion_add(_ac_xx_, _bc_xx_,
                          _tt_xx_) ;
        mp::expansion_add(_ac_yy_, _bc_yy_,
                          _tt_yy_) ;
        mp::expansion_add(_ac_zz_, _bc_zz_,
                          _tt_zz_) ;

        mp::expansion_dot(_ab_xx_, _tt_xx_,
                          _ab_yy_, _tt_yy_,
                          _ab_zz_, _tt_zz_,
                          _absum_) ;

    /*----------------------------------- return signed d */
        return mp::expansion_est(_absum_) ;
    }

    __normal_call REAL_TYPE bisect3d_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect3d predicate, "bound" version */
        ia_flt    _acx, _acy, _acz ,
                  _bcx, _bcy, _bcz ,
                  _abx, _aby, _abz ;
        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _abx.from_sub(_pa[0], _pb[0]) ;   // coord. diff.
        _aby.from_sub(_pa[1], _pb[1]) ;
        _abz.from_sub(_pa[2], _pb[2]) ;

        _acx.from_sub(_pa[0], _pc[0]) ;
        _acy.from_sub(_pa[1], _pc[1]) ;
        _acz.from_sub(_pa[2], _pc[2]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;
        _bcz.from_sub(_pb[2], _pc[2]) ;

        _sgn = (_abx * (_acx + _bcx))
             + (_aby * (_acy + _bcy))
             + (_abz * (_acz + _bcz)) ;

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid () ) ;
    }

    __normal_call REAL_TYPE bisect3d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect3d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  6. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy, _acz ;
        REAL_TYPE _bcx, _bcy, _bcz ;
        REAL_TYPE _acsqr, _bcsqr ;

        REAL_TYPE _sgn, _FT ;

        REAL_TYPE _ACSQR, _BCSQR ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;
        _acz = _pa [2] - _pc [2] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;
        _bcz = _pb [2] - _pc [2] ;

        _acsqr = _acx * _acx +
                 _acy * _acy +
                 _acz * _acz ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy +
                 _bcz * _bcz ;

        _ACSQR = std::abs(_acsqr);
        _BCSQR = std::abs(_bcsqr);

        _FT  = _ACSQR + _BCSQR ;          // roundoff tol
        _FT *= _ER  ;

        _sgn = _acsqr - _bcsqr ;          // d_ab - d_bc

        _OK  = _sgn > +_FT || _sgn < -_FT ;

        return _sgn ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the weighted "bisect" predicate in E^3.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect3w_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect3w predicate, "exact" version */
        mp::expansion< 2 > _ab_xx_, _ab_yy_,
                           _ab_zz_, _ab_ww_,
                           _ac_xx_, _ac_yy_,
                           _ac_zz_,
                           _bc_xx_, _bc_yy_,
                           _bc_zz_;
        mp::expansion< 4 > _tt_xx_, _tt_yy_,
                           _tt_zz_;
        mp::expansion< 48> _ttsum_;
        mp::expansion< 50> _absum_;

        _OK = true;

    /*----------------------------------- compute: d(p,q) */
        _ab_xx_.from_sub(_pa[0], _pb[0]);
        _ab_yy_.from_sub(_pa[1], _pb[1]);
        _ab_zz_.from_sub(_pa[2], _pb[2]);
        _ab_ww_.from_sub(_pa[3], _pb[3]);

        _ac_xx_.from_sub(_pa[0], _pc[0]);
        _ac_yy_.from_sub(_pa[1], _pc[1]);
        _ac_zz_.from_sub(_pa[2], _pc[2]);

        _bc_xx_.from_sub(_pb[0], _pc[0]);
        _bc_yy_.from_sub(_pb[1], _pc[1]);
        _bc_zz_.from_sub(_pb[2], _pc[2]);

        mp::expansion_add(_ac_xx_, _bc_xx_,
                          _tt_xx_) ;
        mp::expansion_add(_ac_yy_, _bc_yy_,
                          _tt_yy_) ;
        mp::expansion_add(_ac_zz_, _bc_zz_,
                          _tt_zz_) ;

        mp::expansion_dot(_ab_xx_, _tt_xx_,
                          _ab_yy_, _tt_yy_,
                          _ab_zz_, _tt_zz_,
                          _ttsum_) ;

        mp::expansion_sub(_ttsum_, _ab_ww_,
                          _absum_) ;

    /*----------------------------------- return signed d */
        return mp::expansion_est(_absum_) ;
    }

    __normal_call REAL_TYPE bisect3w_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect3w predicate, "bound" version */
        ia_flt    _acx, _acy, _acz ,
                  _abw,
                  _bcx, _bcy, _bcz ,
                  _abx, _aby, _abz ;
        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _abx.from_sub(_pa[0], _pb[0]) ;   // coord. diff.
        _aby.from_sub(_pa[1], _pb[1]) ;
        _abz.from_sub(_pa[2], _pb[2]) ;
        _abw.from_sub(_pa[3], _pb[3]) ;

        _acx.from_sub(_pa[0], _pc[0]) ;
        _acy.from_sub(_pa[1], _pc[1]) ;
        _acz.from_sub(_pa[2], _pc[2]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;
        _bcz.from_sub(_pb[2], _pc[2]) ;

        _sgn = (_abx * (_acx + _bcx))
             + (_aby * (_acy + _bcy))
             + (_abz * (_acz + _bcz)) ;

        _sgn-=  _abw ;

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid () ) ;
    }

    __normal_call REAL_TYPE bisect3w_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect3w predicate, "float" version */
        REAL_TYPE static const _ER =
        +  7. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy, _acz ;
        REAL_TYPE _bcx, _bcy, _bcz ;
        REAL_TYPE _acsqr, _bcsqr ;
        REAL_TYPE _a_sum, _b_sum ;

        REAL_TYPE _A_SUM, _B_SUM ;

        REAL_TYPE _sgn, _FT;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;
        _acz = _pa [2] - _pc [2] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;
        _bcz = _pb [2] - _pc [2] ;

        _acsqr = _acx * _acx +
                 _acy * _acy +
                 _acz * _acz ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy +
                 _bcz * _bcz ;

        _a_sum = _acsqr - _pa[3] ;
        _b_sum = _bcsqr - _pb[3] ;

        _A_SUM = std::abs(_acsqr)
               + std::abs(_pa[3]);
        _B_SUM = std::abs(_bcsqr)
               + std::abs(_pb[3]);

        _FT  = _A_SUM + _B_SUM ;          // roundoff tol
        _FT *= _ER ;

        _sgn = _a_sum - _b_sum ;          // d_ab - d_bc

        _OK  = _sgn > +_FT || _sgn < -_FT ;

        return _sgn ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the unweighted "bisect" predicate in E^4.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect4d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect4d predicate, "exact" version */
        mp::expansion< 2 > _ab_xx_, _ab_yy_,
                           _ab_zz_, _ab_qq_,
                           _ac_xx_, _ac_yy_,
                           _ac_zz_, _ac_qq_,
                           _bc_xx_, _bc_yy_,
                           _bc_zz_, _bc_qq_;
        mp::expansion< 4 > _tt_xx_, _tt_yy_,
                           _tt_zz_, _tt_qq_;
        mp::expansion< 64> _absum_;

        _OK = true;

    /*----------------------------------- compute: d(p,q) */
        _ab_xx_.from_sub(_pa[0], _pb[0]);
        _ab_yy_.from_sub(_pa[1], _pb[1]);
        _ab_zz_.from_sub(_pa[2], _pb[2]);
        _ab_qq_.from_sub(_pa[3], _pb[3]);

        _ac_xx_.from_sub(_pa[0], _pc[0]);
        _ac_yy_.from_sub(_pa[1], _pc[1]);
        _ac_zz_.from_sub(_pa[2], _pc[2]);
        _ac_qq_.from_sub(_pa[3], _pc[3]);

        _bc_xx_.from_sub(_pb[0], _pc[0]);
        _bc_yy_.from_sub(_pb[1], _pc[1]);
        _bc_zz_.from_sub(_pb[2], _pc[2]);
        _bc_qq_.from_sub(_pb[3], _pc[3]);

        mp::expansion_add(_ac_xx_, _bc_xx_,
                          _tt_xx_) ;
        mp::expansion_add(_ac_yy_, _bc_yy_,
                          _tt_yy_) ;
        mp::expansion_add(_ac_zz_, _bc_zz_,
                          _tt_zz_) ;
        mp::expansion_add(_ac_qq_, _bc_qq_,
                          _tt_qq_) ;

        mp::expansion_dot(_ab_xx_, _tt_xx_,
                          _ab_yy_, _tt_yy_,
                          _ab_zz_, _tt_zz_,
                          _ab_qq_, _tt_qq_,
                          _absum_) ;

    /*----------------------------------- return signed d */
        return mp::expansion_est(_absum_) ;
    }

    __normal_call REAL_TYPE bisect4d_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect4d predicate, "bound" version */
        ia_flt    _acx, _acy, _acz, _acq ,
                  _bcx, _bcy, _bcz, _bcq ,
                  _abx, _aby, _abz, _abq ;
        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _abx.from_sub(_pa[0], _pb[0]) ;   // coord. diff.
        _aby.from_sub(_pa[1], _pb[1]) ;
        _abz.from_sub(_pa[2], _pb[2]) ;
        _abq.from_sub(_pa[3], _pb[3]) ;

        _acx.from_sub(_pa[0], _pc[0]) ;
        _acy.from_sub(_pa[1], _pc[1]) ;
        _acz.from_sub(_pa[2], _pc[2]) ;
        _acq.from_sub(_pa[3], _pc[3]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;
        _bcz.from_sub(_pb[2], _pc[2]) ;
        _bcq.from_sub(_pb[3], _pc[3]) ;

        _sgn = (_abx * (_acx + _bcx))
             + (_aby * (_acy + _bcy))
             + (_abz * (_acz + _bcz))
             + (_abq * (_acq + _bcq)) ;

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid () ) ;
    }

    __normal_call REAL_TYPE bisect4d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect4d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  7. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy, _acz, _acq;
        REAL_TYPE _bcx, _bcy, _bcz, _bcq;
        REAL_TYPE _acsqr, _bcsqr ;

        REAL_TYPE _sgn, _FT ;

        REAL_TYPE _ACSQR, _BCSQR ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;
        _acz = _pa [2] - _pc [2] ;
        _acq = _pa [3] - _pc [3] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;
        _bcz = _pb [2] - _pc [2] ;
        _bcq = _pb [3] - _pc [3] ;

        _acsqr = _acx * _acx +
                 _acy * _acy +
                 _acz * _acz +
                 _acq * _acq ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy +
                 _bcz * _bcz +
                 _bcq * _bcq ;

        _ACSQR = std::abs(_acsqr);
        _BCSQR = std::abs(_bcsqr);

        _FT  = _ACSQR + _BCSQR ;          // roundoff tol
        _FT *= _ER  ;

        _sgn = _acsqr - _bcsqr ;          // d_ab - d_bc

        _OK  = _sgn > +_FT || _sgn < -_FT ;

        return _sgn ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the weighted "bisect" predicate in E^4.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect4w_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect4w predicate, "exact" version */
        mp::expansion< 2 > _ab_xx_, _ab_yy_,
                           _ab_zz_, _ab_qq_,
                           _ab_ww_,
                           _ac_xx_, _ac_yy_,
                           _ac_zz_, _ac_qq_,
                           _bc_xx_, _bc_yy_,
                           _bc_zz_, _bc_qq_;
        mp::expansion< 4 > _tt_xx_, _tt_yy_,
                           _tt_zz_, _tt_qq_;
        mp::expansion< 64> _ttsum_;
        mp::expansion< 66> _absum_;

        _OK = true;

    /*----------------------------------- compute: d(p,q) */
        _ab_xx_.from_sub(_pa[0], _pb[0]);
        _ab_yy_.from_sub(_pa[1], _pb[1]);
        _ab_zz_.from_sub(_pa[2], _pb[2]);
        _ab_qq_.from_sub(_pa[3], _pb[3]);
        _ab_ww_.from_sub(_pa[4], _pb[4]);

        _ac_xx_.from_sub(_pa[0], _pc[0]);
        _ac_yy_.from_sub(_pa[1], _pc[1]);
        _ac_zz_.from_sub(_pa[2], _pc[2]);
        _ac_qq_.from_sub(_pa[3], _pc[3]);

        _bc_xx_.from_sub(_pb[0], _pc[0]);
        _bc_yy_.from_sub(_pb[1], _pc[1]);
        _bc_zz_.from_sub(_pb[2], _pc[2]);
        _bc_qq_.from_sub(_pb[3], _pc[3]);

        mp::expansion_add(_ac_xx_, _bc_xx_,
                          _tt_xx_) ;
        mp::expansion_add(_ac_yy_, _bc_yy_,
                          _tt_yy_) ;
        mp::expansion_add(_ac_zz_, _bc_zz_,
                          _tt_zz_) ;
        mp::expansion_add(_ac_qq_, _bc_qq_,
                          _tt_qq_) ;

        mp::expansion_dot(_ab_xx_, _tt_xx_,
                          _ab_yy_, _tt_yy_,
                          _ab_zz_, _tt_zz_,
                          _ab_qq_, _tt_qq_,
                          _ttsum_) ;

        mp::expansion_sub(_ttsum_, _ab_ww_,
                          _absum_) ;

    /*----------------------------------- return signed d */
        return mp::expansion_est(_absum_) ;
    }

    __normal_call REAL_TYPE bisect4w_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect4w predicate, "bound" version */
        ia_flt    _acx, _acy, _acz, _acq,
                  _abw,
                  _bcx, _bcy, _bcz, _bcq,
                  _abx, _aby, _abz, _abq;
        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _abx.from_sub(_pa[0], _pb[0]) ;   // coord. diff.
        _aby.from_sub(_pa[1], _pb[1]) ;
        _abz.from_sub(_pa[2], _pb[2]) ;
        _abq.from_sub(_pa[3], _pb[3]) ;
        _abw.from_sub(_pa[4], _pb[4]) ;

        _acx.from_sub(_pa[0], _pc[0]) ;
        _acy.from_sub(_pa[1], _pc[1]) ;
        _acz.from_sub(_pa[2], _pc[2]) ;
        _acq.from_sub(_pa[3], _pc[3]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;
        _bcz.from_sub(_pb[2], _pc[2]) ;
        _bcq.from_sub(_pb[3], _pc[3]) ;

        _sgn = (_abx * (_acx + _bcx))
             + (_aby * (_acy + _bcy))
             + (_abz * (_acz + _bcz))
             + (_abq * (_acq + _bcq)) ;

        _sgn-=  _abw ;

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid () ) ;
    }

    __normal_call REAL_TYPE bisect4w_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- bisect4w predicate, "float" version */
        REAL_TYPE static const _ER =
        +  8. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy, _acz, _acq;
        REAL_TYPE _bcx, _bcy, _bcz, _bcq;
        REAL_TYPE _acsqr, _bcsqr ;
        REAL_TYPE _a_sum, _b_sum ;

        REAL_TYPE _A_SUM, _B_SUM ;

        REAL_TYPE _sgn, _FT;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;
        _acz = _pa [2] - _pc [2] ;
        _acq = _pa [3] - _pc [3] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;
        _bcz = _pb [2] - _pc [2] ;
        _bcq = _pb [3] - _pc [3] ;

        _acsqr = _acx * _acx +
                 _acy * _acy +
                 _acz * _acz +
                 _acq * _acq ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy +
                 _bcz * _bcz +
                 _bcq * _bcq ;

        _a_sum = _acsqr - _pa[4] ;
        _b_sum = _bcsqr - _pb[4] ;

        _A_SUM = std::abs(_acsqr)
               + std::abs(_pa[4]);
        _B_SUM = std::abs(_bcsqr)
               + std::abs(_pb[4]);

        _FT  = _A_SUM + _B_SUM ;          // roundoff tol
        _FT *= _ER ;

        _sgn = _a_sum - _b_sum ;          // d_ab - d_bc

        _OK  = _sgn > +_FT || _sgn < -_FT ;

        return _sgn ;
    }



