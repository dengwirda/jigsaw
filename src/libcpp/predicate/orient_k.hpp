
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
     * Compute an exact determinant using multi-precision
     * expansions, a'la shewchuk
     *
     *   | ax  ay  +1. |
     *   | bx  by  +1. |
     *   | cx  cy  +1. |
     *
     * This is the planar "orientation" predicate in E^2.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE orient2d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- orient2d predicate, "exact" version */
        mp::expansion< 4 > _d2_ab_, _d2_ac_,
                           _d2_bc_;
        mp::expansion< 12> _d3full;

        _OK = true;

    /*-------------------------------------- 2 x 2 minors */
        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pb[ 0], _pb[ 1],
                        _d2_ab_ ) ;

        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pc[ 0], _pc[ 1],
                        _d2_ac_ ) ;

        compute_det_2x2(_pb[ 0], _pb[ 1],
                        _pc[ 0], _pc[ 1],
                        _d2_bc_ ) ;

    /*-------------------------------------- 3 x 3 result */
        unitary_det_3x3(_d2_bc_, _d2_ac_,
                        _d2_ab_,
                        _d3full, +3) ;

    /*-------------------------------------- leading det. */
        return mp::expansion_est(_d3full) ;
    }

    __normal_call REAL_TYPE orient2d_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- orient2d predicate, "bound" version */
        ia_flt    _acx, _acy ;
        ia_flt    _bcx, _bcy ;
        ia_flt    _acxbcy, _acybcx ;

        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _acx.from_sub(_pa[0], _pc[0]) ;   // coord. diff.
        _acy.from_sub(_pa[1], _pc[1]) ;

        _bcx.from_sub(_pb[0], _pc[0]) ;
        _bcy.from_sub(_pb[1], _pc[1]) ;

        _acxbcy = _acx * _bcy ;
        _acybcx = _acy * _bcx ;

        _sgn = _acxbcy - _acybcx ;        // 2 x 2 result

        _OK  =
           _sgn.lo() >= (REAL_TYPE)0.
        || _sgn.up() <= (REAL_TYPE)0. ;

        return ( _sgn.mid() ) ;
    }

    __normal_call REAL_TYPE orient2d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        bool_type &_OK
        )
    {
    /*--------------- orient2d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  4. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy ;
        REAL_TYPE _bcx, _bcy ;
        REAL_TYPE _acxbcy, _acybcx ;

        REAL_TYPE _ACXBCY, _ACYBCX ;

        REAL_TYPE _sgn, _FT;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;

        _acxbcy = _acx * _bcy ;
        _acybcx = _acy * _bcx ;

        _ACXBCY = std::abs(_acxbcy);
        _ACYBCX = std::abs(_acybcx);

        _FT  = _ACXBCY + _ACYBCX ;        // roundoff tol
        _FT *= _ER ;

        _sgn = _acxbcy - _acybcx ;        // 2 x 2 result

        _OK  =
          _sgn > +_FT || _sgn < -_FT ;

        return ( _sgn ) ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact determinant using multi-precision
     * expansions, a'la shewchuk
     *
     *   | ax  ay  az  +1. |
     *   | bx  by  bz  +1. |
     *   | cx  cy  cz  +1. |
     *   | dx  dy  dz  +1. |
     *
     * This is the planar "orientation" predicate in E^3.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE orient3d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
        bool_type &_OK
        )
    {
    /*--------------- orient3d predicate, "exact" version */
        mp::expansion< 4 > _d2_ab_, _d2_ac_,
                           _d2_ad_,
                           _d2_bc_, _d2_bd_,
                           _d2_cd_;
        mp::expansion< 12> _d3_abc, _d3_abd,
                           _d3_acd, _d3_bcd;
        mp::expansion< 96> _d4full;

        _OK = true;

        mp::expansion< 1 > _pa_zz_(_pa[ 2]);
        mp::expansion< 1 > _pb_zz_(_pb[ 2]);
        mp::expansion< 1 > _pc_zz_(_pc[ 2]);
        mp::expansion< 1 > _pd_zz_(_pd[ 2]);

    /*-------------------------------------- 2 x 2 minors */
        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pb[ 0], _pb[ 1],
                        _d2_ab_ ) ;

        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pc[ 0], _pc[ 1],
                        _d2_ac_ ) ;

        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pd[ 0], _pd[ 1],
                        _d2_ad_ ) ;

        compute_det_2x2(_pb[ 0], _pb[ 1],
                        _pc[ 0], _pc[ 1],
                        _d2_bc_ ) ;

        compute_det_2x2(_pb[ 0], _pb[ 1],
                        _pd[ 0], _pd[ 1],
                        _d2_bd_ ) ;

        compute_det_2x2(_pc[ 0], _pc[ 1],
                        _pd[ 0], _pd[ 1],
                        _d2_cd_ ) ;

    /*-------------------------------------- 3 x 3 minors */
        unitary_det_3x3(_d2_cd_, _d2_bd_,
                        _d2_bc_,
                        _d3_bcd, +3) ;

        unitary_det_3x3(_d2_cd_, _d2_ad_,
                        _d2_ac_,
                        _d3_acd, +3) ;

        unitary_det_3x3(_d2_bd_, _d2_ad_,
                        _d2_ab_,
                        _d3_abd, +3) ;

        unitary_det_3x3(_d2_bc_, _d2_ac_,
                        _d2_ab_,
                        _d3_abc, +3) ;

    /*-------------------------------------- 4 x 4 result */
        compute_det_4x4(_d3_bcd, _pa_zz_,
                        _d3_acd, _pb_zz_,
                        _d3_abd, _pc_zz_,
                        _d3_abc, _pd_zz_,
                        _d4full, +3) ;

    /*-------------------------------------- leading det. */
        return mp::expansion_est(_d4full) ;
    }

    __normal_call REAL_TYPE orient3d_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
        bool_type &_OK
        )
    {
    /*--------------- orient3d predicate, "bound" version */
        ia_flt    _adx, _ady, _adz ,
                  _bdx, _bdy, _bdz ,
                  _cdx, _cdy, _cdz ;
        ia_flt    _bdxcdy, _cdxbdy ,
                  _cdxady, _adxcdy ,
                  _adxbdy, _bdxady ;

        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _adx.from_sub(_pa[0], _pd[0]) ;   // coord. diff.
        _ady.from_sub(_pa[1], _pd[1]) ;
        _adz.from_sub(_pa[2], _pd[2]) ;

        _bdx.from_sub(_pb[0], _pd[0]) ;
        _bdy.from_sub(_pb[1], _pd[1]) ;
        _bdz.from_sub(_pb[2], _pd[2]) ;

        _cdx.from_sub(_pc[0], _pd[0]) ;
        _cdy.from_sub(_pc[1], _pd[1]) ;
        _cdz.from_sub(_pc[2], _pd[2]) ;

        _bdxcdy = _bdx * _cdy ;           // 2 x 2 minors
        _cdxbdy = _cdx * _bdy ;
        _cdxady = _cdx * _ady ;
        _adxcdy = _adx * _cdy ;
        _adxbdy = _adx * _bdy ;
        _bdxady = _bdx * _ady ;

        _sgn =                            // 3 x 3 result
          _adz * (_bdxcdy - _cdxbdy)
        + _bdz * (_cdxady - _adxcdy)
        + _cdz * (_adxbdy - _bdxady);

        _OK  =
          _sgn.lo() >= (REAL_TYPE)0.
        ||_sgn.up() <= (REAL_TYPE)0.;

        return ( _sgn.mid() ) ;
    }

    __normal_call REAL_TYPE orient3d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
        bool_type &_OK
        )
    {
    /*--------------- orient3d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  8. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _adx, _ady, _adz ,
                  _bdx, _bdy, _bdz ,
                  _cdx, _cdy, _cdz ;
        REAL_TYPE _bdxcdy, _cdxbdy ,
                  _cdxady, _adxcdy ,
                  _adxbdy, _bdxady ;

        REAL_TYPE _ADZ, _BDZ, _CDZ ;
        REAL_TYPE _BDXCDY, _CDXBDY ,
                  _CDXADY, _ADXCDY ,
                  _ADXBDY, _BDXADY ;

        REAL_TYPE _sgn, _FT;

        _adx = _pa [0] - _pd [0] ;        // coord. diff.
        _ady = _pa [1] - _pd [1] ;
        _adz = _pa [2] - _pd [2] ;

        _ADZ = std::abs (_adz) ;

        _bdx = _pb [0] - _pd [0] ;
        _bdy = _pb [1] - _pd [1] ;
        _bdz = _pb [2] - _pd [2] ;

        _BDZ = std::abs (_bdz) ;

        _cdx = _pc [0] - _pd [0] ;
        _cdy = _pc [1] - _pd [1] ;
        _cdz = _pc [2] - _pd [2] ;

        _CDZ = std::abs (_cdz) ;

        _bdxcdy = _bdx * _cdy ;           // 2 x 2 minors
        _cdxbdy = _cdx * _bdy ;
        _cdxady = _cdx * _ady ;
        _adxcdy = _adx * _cdy ;
        _adxbdy = _adx * _bdy ;
        _bdxady = _bdx * _ady ;

        _BDXCDY = std::abs (_bdxcdy) ;
        _CDXBDY = std::abs (_cdxbdy) ;
        _CDXADY = std::abs (_cdxady) ;
        _ADXCDY = std::abs (_adxcdy) ;
        _ADXBDY = std::abs (_adxbdy) ;
        _BDXADY = std::abs (_bdxady) ;

        _FT  =                            // roundoff tol
          _ADZ * (_BDXCDY + _CDXBDY)
        + _BDZ * (_CDXADY + _ADXCDY)
        + _CDZ * (_ADXBDY + _BDXADY) ;

        _FT *= _ER ;

        _sgn =                            // 3 x 3 result
          _adz * (_bdxcdy - _cdxbdy)
        + _bdz * (_cdxady - _adxcdy)
        + _cdz * (_adxbdy - _bdxady) ;

        _OK  =
          _sgn > +_FT || _sgn < -_FT ;

        return ( _sgn ) ;
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact determinant using multi-precision
     * expansions, a'la shewchuk
     *
     *   | ax  ay  az  aq  +1. |
     *   | bx  by  bz  bq  +1. |
     *   | cx  cy  cz  cq  +1. |
     *   | dx  dy  dz  dq  +1. |
     *   | ex  ey  ez  eq  +1. |
     *
     * This is the planar "orientation" predicate in E^4.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE orient4d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
      __const_ptr(REAL_TYPE) _pe ,
        bool_type &_OK
        )
    {
    /*--------------- orient4d predicate, "exact" version */
        mp::expansion< 4 > _d2_ab_, _d2_ac_,
                           _d2_ad_, _d2_ae_,
                           _d2_bc_, _d2_bd_,
                           _d2_be_,
                           _d2_cd_, _d2_ce_,
                           _d2_de_;
        mp::expansion< 24> _d3_abc, _d3_abd,
                           _d3_abe,
                           _d3_acd, _d3_ace,
                           _d3_ade,
                           _d3_bcd, _d3_bce,
                           _d3_bde, _d3_cde;
        mp::expansion< 96> _d4abcd, _d4abce,
                           _d4abde, _d4acde,
                           _d4bcde;
        mp::expansion<960> _d5full;

        _OK = true;

        mp::expansion< 1 > _pa_zz_(_pa[ 2]);
        mp::expansion< 1 > _pb_zz_(_pb[ 2]);
        mp::expansion< 1 > _pc_zz_(_pc[ 2]);
        mp::expansion< 1 > _pd_zz_(_pd[ 2]);
        mp::expansion< 1 > _pe_zz_(_pe[ 2]);

        mp::expansion< 1 > _pa_qq_(_pa[ 3]);
        mp::expansion< 1 > _pb_qq_(_pb[ 3]);
        mp::expansion< 1 > _pc_qq_(_pc[ 3]);
        mp::expansion< 1 > _pd_qq_(_pd[ 3]);
        mp::expansion< 1 > _pe_qq_(_pe[ 3]);

    /*-------------------------------------- 2 x 2 minors */
        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pb[ 0], _pb[ 1],
                        _d2_ab_ ) ;

        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pc[ 0], _pc[ 1],
                        _d2_ac_ ) ;

        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pd[ 0], _pd[ 1],
                        _d2_ad_ ) ;

        compute_det_2x2(_pa[ 0], _pa[ 1],
                        _pe[ 0], _pe[ 1],
                        _d2_ae_ ) ;

        compute_det_2x2(_pb[ 0], _pb[ 1],
                        _pc[ 0], _pc[ 1],
                        _d2_bc_ ) ;

        compute_det_2x2(_pb[ 0], _pb[ 1],
                        _pd[ 0], _pd[ 1],
                        _d2_bd_ ) ;

        compute_det_2x2(_pb[ 0], _pb[ 1],
                        _pe[ 0], _pe[ 1],
                        _d2_be_ ) ;

        compute_det_2x2(_pc[ 0], _pc[ 1],
                        _pd[ 0], _pd[ 1],
                        _d2_cd_ ) ;

        compute_det_2x2(_pc[ 0], _pc[ 1],
                        _pe[ 0], _pe[ 1],
                        _d2_ce_ ) ;

        compute_det_2x2(_pd[ 0], _pd[ 1],
                        _pe[ 0], _pe[ 1],
                        _d2_de_ ) ;

    /*-------------------------------------- 3 x 3 minors */
        compute_det_3x3(_d2_bc_, _pa_zz_,
                        _d2_ac_, _pb_zz_,
                        _d2_ab_, _pc_zz_,
                        _d3_abc, +3) ;

        compute_det_3x3(_d2_bd_, _pa_zz_,
                        _d2_ad_, _pb_zz_,
                        _d2_ab_, _pd_zz_,
                        _d3_abd, +3) ;

        compute_det_3x3(_d2_be_, _pa_zz_,
                        _d2_ae_, _pb_zz_,
                        _d2_ab_, _pe_zz_,
                        _d3_abe, +3) ;

        compute_det_3x3(_d2_cd_, _pa_zz_,
                        _d2_ad_, _pc_zz_,
                        _d2_ac_, _pd_zz_,
                        _d3_acd, +3) ;

        compute_det_3x3(_d2_ce_, _pa_zz_,
                        _d2_ae_, _pc_zz_,
                        _d2_ac_, _pe_zz_,
                        _d3_ace, +3) ;

        compute_det_3x3(_d2_de_, _pa_zz_,
                        _d2_ae_, _pd_zz_,
                        _d2_ad_, _pe_zz_,
                        _d3_ade, +3) ;

        compute_det_3x3(_d2_cd_, _pb_zz_,
                        _d2_bd_, _pc_zz_,
                        _d2_bc_, _pd_zz_,
                        _d3_bcd, +3) ;

        compute_det_3x3(_d2_ce_, _pb_zz_,
                        _d2_be_, _pc_zz_,
                        _d2_bc_, _pe_zz_,
                        _d3_bce, +3) ;

        compute_det_3x3(_d2_de_, _pb_zz_,
                        _d2_be_, _pd_zz_,
                        _d2_bd_, _pe_zz_,
                        _d3_bde, +3) ;

        compute_det_3x3(_d2_de_, _pc_zz_,
                        _d2_ce_, _pd_zz_,
                        _d2_cd_, _pe_zz_,
                        _d3_cde, +3) ;

    /*-------------------------------------- 4 x 4 minors */
        unitary_det_4x4(_d3_cde, _d3_bde,
                        _d3_bce, _d3_bcd,
                        _d4bcde, +4) ;

        unitary_det_4x4(_d3_cde, _d3_ade,
                        _d3_ace, _d3_acd,
                        _d4acde, +4) ;

        unitary_det_4x4(_d3_bde, _d3_ade,
                        _d3_abe, _d3_abd,
                        _d4abde, +4) ;

        unitary_det_4x4(_d3_bce, _d3_ace,
                        _d3_abe, _d3_abc,
                        _d4abce, +4) ;

        unitary_det_4x4(_d3_bcd, _d3_acd,
                        _d3_abd, _d3_abc,
                        _d4abcd, +4) ;

    /*-------------------------------------- 5 x 5 result */
        compute_det_5x5(_d4bcde, _pa_qq_,
                        _d4acde, _pb_qq_,
                        _d4abde, _pc_qq_,
                        _d4abce, _pd_qq_,
                        _d4abcd, _pe_qq_,
                        _d5full, +4) ;

    /*-------------------------------------- leading det. */
        return mp::expansion_est(_d5full) ;
    }

    __normal_call REAL_TYPE orient4d_i (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
      __const_ptr(REAL_TYPE) _pe ,
        bool_type &_OK
        )
    {
    /*--------------- orient4d predicate, "bound" version */
        ia_flt    _aex, _aey, _aez ,
                  _aeq,
                  _bex, _bey, _bez ,
                  _beq,
                  _cex, _cey, _cez ,
                  _ceq,
                  _dex, _dey, _dez ,
                  _deq;
        ia_flt    _aexbey, _bexaey ,
                  _aexcey, _cexaey ,
                  _bexcey, _cexbey ,
                  _cexdey, _dexcey ,
                  _dexaey, _aexdey ,
                  _bexdey, _dexbey ;
        ia_flt    _ab_, _bc_, _cd_, _da_,
                  _ac_, _bd_;
        ia_flt    _abc, _bcd, _cda, _dab;
        ia_flt    _sgn;

        ia_rnd    _rnd;                   // up rounding!

        _aex.from_sub(_pa[0], _pe[0]) ;   // coord. diff.
        _aey.from_sub(_pa[1], _pe[1]) ;
        _aez.from_sub(_pa[2], _pe[2]) ;
        _aeq.from_sub(_pa[3], _pe[3]) ;

        _bex.from_sub(_pb[0], _pe[0]) ;
        _bey.from_sub(_pb[1], _pe[1]) ;
        _bez.from_sub(_pb[2], _pe[2]) ;
        _beq.from_sub(_pb[3], _pe[3]) ;

        _cex.from_sub(_pc[0], _pe[0]) ;
        _cey.from_sub(_pc[1], _pe[1]) ;
        _cez.from_sub(_pc[2], _pe[2]) ;
        _ceq.from_sub(_pc[3], _pe[3]) ;

        _dex.from_sub(_pd[0], _pe[0]) ;
        _dey.from_sub(_pd[1], _pe[1]) ;
        _dez.from_sub(_pd[2], _pe[2]) ;
        _deq.from_sub(_pd[3], _pe[3]) ;

        _aexbey = _aex * _bey ;           // 2 x 2 minors
        _bexaey = _bex * _aey ;
        _ab_ = _aexbey - _bexaey ;

        _bexcey = _bex * _cey;
        _cexbey = _cex * _bey;
        _bc_ = _bexcey - _cexbey ;

        _cexdey = _cex * _dey;
        _dexcey = _dex * _cey;
        _cd_ = _cexdey - _dexcey ;

        _dexaey = _dex * _aey;
        _aexdey = _aex * _dey;
        _da_ = _dexaey - _aexdey ;

        _aexcey = _aex * _cey;
        _cexaey = _cex * _aey;
        _ac_ = _aexcey - _cexaey ;

        _bexdey = _bex * _dey;
        _dexbey = _dex * _bey;
        _bd_ = _bexdey - _dexbey ;

        _abc =                            // 3 x 3 minors
          _aez * _bc_ - _bez * _ac_
        + _cez * _ab_ ;

        _bcd =
          _bez * _cd_ - _cez * _bd_
        + _dez * _bc_ ;

        _cda =
          _cez * _da_ + _dez * _ac_
        + _aez * _cd_ ;

        _dab =
          _dez * _ab_ + _aez * _bd_
        + _bez * _da_ ;

        _sgn =                            // 4 x 4 result
          _deq * _abc - _ceq * _dab
        + _beq * _cda - _aeq * _bcd ;

        _OK =
          _sgn.lo() >= (REAL_TYPE)0.
        ||_sgn.up() <= (REAL_TYPE)0.;

        return ( _sgn.mid() ) ;
    }

    __normal_call REAL_TYPE orient4d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
      __const_ptr(REAL_TYPE) _pe ,
        bool_type &_OK
        )
    {
    /*--------------- orient4d predicate, "float" version */
        REAL_TYPE static const _ER =
        + 13. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _aex, _aey, _aez ,
                  _aeq,
                  _bex, _bey, _bez ,
                  _beq,
                  _cex, _cey, _cez ,
                  _ceq,
                  _dex, _dey, _dez ,
                  _deq;
        REAL_TYPE _aexbey, _bexaey ,
                  _aexcey, _cexaey ,
                  _bexcey, _cexbey ,
                  _cexdey, _dexcey ,
                  _dexaey, _aexdey ,
                  _bexdey, _dexbey ;
        REAL_TYPE _ab_, _bc_, _cd_, _da_,
                  _ac_, _bd_;
        REAL_TYPE _abc, _bcd, _cda, _dab;

        REAL_TYPE _AEZ, _BEZ, _CEZ, _DEZ;
        REAL_TYPE _AEQ, _BEQ, _CEQ, _DEQ;
        REAL_TYPE _AEXBEY, _BEXAEY ,
                  _CEXAEY, _AEXCEY ,
                  _BEXCEY, _CEXBEY ,
                  _CEXDEY, _DEXCEY ,
                  _DEXAEY, _AEXDEY ,
                  _BEXDEY, _DEXBEY ;
        REAL_TYPE _AB_, _BC_, _CD_, _DA_,
                  _AC_, _BD_;
        REAL_TYPE _ABC, _BCD, _CDA, _DAB;

        REAL_TYPE _sgn, _FT ;

        _aex = _pa [0] - _pe [0] ;        // coord. diff.
        _aey = _pa [1] - _pe [1] ;
        _aez = _pa [2] - _pe [2] ;
        _aeq = _pa [3] - _pe [3] ;

        _AEZ = std::abs (_aez) ;
        _AEQ = std::abs (_aeq) ;

        _bex = _pb [0] - _pe [0] ;
        _bey = _pb [1] - _pe [1] ;
        _bez = _pb [2] - _pe [2] ;
        _beq = _pb [3] - _pe [3] ;

        _BEZ = std::abs (_bez) ;
        _BEQ = std::abs (_beq) ;

        _cex = _pc [0] - _pe [0] ;
        _cey = _pc [1] - _pe [1] ;
        _cez = _pc [2] - _pe [2] ;
        _ceq = _pc [3] - _pe [3] ;

        _CEZ = std::abs (_cez) ;
        _CEQ = std::abs (_ceq) ;

        _dex = _pd [0] - _pe [0] ;
        _dey = _pd [1] - _pe [1] ;
        _dez = _pd [2] - _pe [2] ;
        _deq = _pd [3] - _pe [3] ;

        _DEZ = std::abs (_dez) ;
        _DEQ = std::abs (_deq) ;

        _aexbey = _aex * _bey;            // 2 x 2 minors
        _bexaey = _bex * _aey;
        _ab_ = _aexbey - _bexaey ;

        _AEXBEY = std::abs (_aexbey) ;
        _BEXAEY = std::abs (_bexaey) ;
        _AB_ = _AEXBEY + _BEXAEY ;

        _bexcey = _bex * _cey;
        _cexbey = _cex * _bey;
        _bc_ = _bexcey - _cexbey ;

        _BEXCEY = std::abs (_bexcey) ;
        _CEXBEY = std::abs (_cexbey) ;
        _BC_ = _BEXCEY + _CEXBEY ;

        _cexdey = _cex * _dey;
        _dexcey = _dex * _cey;
        _cd_ = _cexdey - _dexcey ;

        _CEXDEY = std::abs (_cexdey) ;
        _DEXCEY = std::abs (_dexcey) ;
        _CD_ = _CEXDEY + _DEXCEY ;

        _dexaey = _dex * _aey;
        _aexdey = _aex * _dey;
        _da_ = _dexaey - _aexdey ;

        _DEXAEY = std::abs (_dexaey) ;
        _AEXDEY = std::abs (_aexdey) ;
        _DA_ = _DEXAEY + _AEXDEY ;

        _aexcey = _aex * _cey;
        _cexaey = _cex * _aey;
        _ac_ = _aexcey - _cexaey ;

        _AEXCEY = std::abs (_aexcey) ;
        _CEXAEY = std::abs (_cexaey) ;
        _AC_ = _AEXCEY + _CEXAEY ;

        _bexdey = _bex * _dey;
        _dexbey = _dex * _bey;
        _bd_ = _bexdey - _dexbey ;

        _BEXDEY = std::abs (_bexdey) ;
        _DEXBEY = std::abs (_dexbey) ;
        _BD_ = _BEXDEY + _DEXBEY ;

        _abc =                            // 3 x 3 minors
          _aez * _bc_ - _bez * _ac_
        + _cez * _ab_ ;
        _ABC =
          _AEZ * _BC_ + _BEZ * _AC_
        + _CEZ * _AB_ ;

        _bcd =
          _bez * _cd_ - _cez * _bd_
        + _dez * _bc_ ;
        _BCD =
          _BEZ * _CD_ + _CEZ * _BD_
        + _DEZ * _BC_ ;

        _cda =
          _cez * _da_ + _dez * _ac_
        + _aez * _cd_ ;
        _CDA =
          _CEZ * _DA_ + _DEZ * _AC_
        + _AEZ * _CD_ ;

        _dab =
          _dez * _ab_ + _aez * _bd_
        + _bez * _da_ ;
        _DAB =
          _DEZ * _AB_ + _AEZ * _BD_
        + _BEZ * _DA_ ;

        _FT  =                            // roundoff tol
          _DEQ * _ABC + _CEQ * _DAB
        + _BEQ * _CDA + _AEQ * _BCD ;

        _FT *= _ER ;

        _sgn =                            // 4 x 4 result
          _deq * _abc - _ceq * _dab
        + _beq * _cda - _aeq * _bcd ;

        _OK  =
          _sgn > _FT || _sgn < -_FT ;

        return ( _sgn ) ;
    }


