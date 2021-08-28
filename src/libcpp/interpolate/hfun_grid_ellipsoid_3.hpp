
    /*
    --------------------------------------------------------
     * HFUN-GRID-ELLIPSOID-3D: ellipsoidal H(x) in R^3.
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
     * Last updated: 20 Apr., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_GRID_ELLIPSOID_3__
#   define __HFUN_GRID_ELLIPSOID_3__

    namespace mesh {

    template <
    typename R ,
    typename V ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_grid_ellipsoid_3d
        : public hfun_base_kd <I, R>
    {
    public :

    /*-------------- lat.-lon. spheroidal size-fun in R^3 */

    typedef R                       real_type ;
    typedef V                       vals_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef hfun_grid_ellipsoid_3d  <
            real_type ,
            vals_type ,
            iptr_type >             hfun_type ;

    typedef typename  hfun_base_kd  <
            iptr_type ,
            real_type >::hint_type  hint_type ;

    typedef containers::array   <
            real_type ,
            allocator >             real_list ;

    typedef containers::array   <
            vals_type ,
            allocator >             vals_list ;

    public  :

    /*--------------- (x/a)**2 + (y/b)**2 + (z/c)**2 = 1. */

    real_type                      _radA = 1. ;
    real_type                      _radB = 1. ;
    real_type                      _radC = 1. ;

    containers::array <
        real_type, allocator>      _xpos ;
    containers::array <
        real_type, allocator>      _ypos ;

    containers::array <
        real_type, allocator>      _sinx ;  // to precompute
    containers::array <
        real_type, allocator>      _siny ;
    containers::array <
        real_type, allocator>      _cosx ;
    containers::array <
        real_type, allocator>      _cosy ;

    containers::array <
        vals_type, allocator>      _hmat ;

    vals_type                      _htop ;  // mean at poles
    vals_type                      _hbot ;

    containers::array <
        vals_type, allocator>      _dhdx ;

    bool_type                      _xvar ;
    bool_type                      _yvar ;

    bool_type                      _wrap ;

    public  :

    __inline_call void_type indx_from_subs (
        iptr_type _ipos,
        iptr_type _jpos,
        iptr_type&_indx
        ) const
    {
    /*------------ helper: convert into "un-rolled" index */
        iptr_type _ynum =
       (iptr_type)this->_ypos.count() ;

        _indx = _jpos * _ynum + _ipos ;
    }

    __inline_call void_type subs_from_indx (
        iptr_type  _indx,
        iptr_type &_ipos,
        iptr_type &_jpos
        ) const
    {
    /*------------ helper: convert from "un-rolled" index */
        iptr_type _ynum =
       (iptr_type)this->_ypos.count() ;

        _ipos = _indx % _ynum ;
        _jpos =(_indx - _ipos )/_ynum ;
    }

    __inline_call void_type toR3 (
      __const_ptr(real_type) _apos ,
      __write_ptr(real_type) _ppos
        ) const
    {
    /*------------ helper: convert from S^2 to R^3 coord. */
        _ppos[0] = this->_radA *
            std::cos( _apos[0] ) *
            std::cos( _apos[1] ) ;

        _ppos[1] = this->_radB *
            std::sin( _apos[0] ) *
            std::cos( _apos[1] ) ;

        _ppos[2] = this->_radC *
            std::sin( _apos[1] ) ;
    }

    __inline_call void_type toS2 (
      __const_ptr(real_type) _ppos ,
      __write_ptr(real_type) _apos
        ) const
    {
    /*------------ helper: convert from R^3 to S^2 coord. */
        real_type _xmul =
            _ppos[0] * this->_radB ;
        real_type _ymul =
            _ppos[1] * this->_radA ;
        real_type _zrat =
            _ppos[2] / this->_radC ;

        _zrat = std::min(+1.,_zrat);
        _zrat = std::max(-1.,_zrat);

        _apos[0]= std::atan2(_ymul,
                             _xmul);
        _apos[1]= std::asin (_zrat);
    }

    /*
    --------------------------------------------------------
     * INIT: init. size-fun. class.
    --------------------------------------------------------
     */

    __inline_call void_type init (
        )
    {
        real_type static const _FTOL =
            std::pow (
            std::numeric_limits<real_type>
            ::epsilon(), (real_type) +.75) ;

        this->_xvar = false ;
        this->_yvar = false ;
        this->_wrap = false ;

        if (this->_xpos.empty()) return ;
        if (this->_ypos.empty()) return ;

    /*-------------------------- test for uniform spacing */
        real_type _xbar, _xmin, _xmax ;
        _xbar = *this->_xpos.tail() -
                *this->_xpos.head() ;

        _xbar /=(this->_xpos.count () - 1) ;

        _xmin = _xbar - _FTOL * _xbar ;
        _xmax = _xbar + _FTOL * _xbar ;

        real_type _ybar, _ymin, _ymax ;
        _ybar = *this->_ypos.tail() -
                *this->_ypos.head() ;

        _ybar /=(this->_ypos.count () - 1) ;

        _ymin = _ybar - _FTOL * _ybar ;
        _ymax = _ybar + _FTOL * _ybar ;

        for (auto _iter  = _xpos.head() ;
                  _iter != _xpos.tail() ;
                ++_iter  )
        {
            real_type _xdel =
                *(_iter+1)-*(_iter+0) ;

            if (_xdel < _xmin ||
                _xdel > _xmax )
            {
                _xvar =  true ; break ;
            }
        }

        for (auto _iter  = _ypos.head() ;
                  _iter != _ypos.tail() ;
                ++_iter  )
        {
            real_type _ydel =
                *(_iter+1)-*(_iter+0) ;

            if (_ydel < _ymin ||
                _ydel > _ymax )
            {
                _yvar =  true ; break ;
            }
        }

    /*-------------------------- test for "x" periodicity */
        real_type _xdel =
            std::cos(*_xpos.tail()) -
            std::cos(*_xpos.head()) ;

        this->_wrap =
            std::abs(_xdel) < _FTOL ;

    /*-------------------------- sum for extrap. at poles */
        double _xtop = 0., _xbot = 0. ;

        iptr_type _indx, _jpos = +0 ;
        for (auto _iter  = _xpos.head() ;
                  _iter != _xpos.tail() ;
                ++_iter, ++_jpos)
        {
            iptr_type _n = -1 +
           (iptr_type)this->_ypos.count() ;

            indx_from_subs(
                +0, _jpos, _indx) ;

            _xtop += this->_hmat[_indx] ;

            indx_from_subs(
                _n, _jpos, _indx) ;

            _xbot += this->_hmat[_indx] ;
        }

        this->_htop = (vals_type)(
            _xtop / this->_xpos.count() ) ;
        this->_hbot = (vals_type)(
            _xbot / this->_xpos.count() ) ;
    }

    /*
    --------------------------------------------------------
     * CLIP-HFUN: impose |dh/dx| limits.
    --------------------------------------------------------
     */

    __normal_call void_type clip (
        )
    {
        class  less_than
        {
    /*-------------------- "LESS-THAN" operator for queue */
        public  :
            typename
            vals_list::_write_it _hptr;

        public  :
        __inline_call less_than  (
            typename
            vals_list::_write_it _hsrc
            ) : _hptr(_hsrc) {}

        __inline_call
            bool_type operator() (
            iptr_type _ipos,
            iptr_type _jpos
            )
        {   return *(this->_hptr+_ipos) <
                   *(this->_hptr+_jpos) ;
        }
        } ;

        typedef typename
            allocator:: size_type   uint_type ;

        uint_type static constexpr
            _null =
        std::numeric_limits<uint_type>::max() ;

        containers::prioritymap <
            iptr_type ,
            less_than ,
            allocator >
        _sort((less_than(this->_hmat.head())));

        containers:: array      <
            typename
            allocator:: size_type,
            allocator >     _keys;

    /*-------------------- init. values for periodic bc's */
        iptr_type IBEG = +0;
        iptr_type IEND =
       (iptr_type)this->_ypos.count() - 1 ;

        iptr_type JBEG = +0;
        iptr_type JEND =
       (iptr_type)this->_xpos.count() - 1 ;

        if (this->_wrap)
        {
        iptr_type _inum  = +0;
        for (auto _iter  =
                   this->_ypos.head() ;
                  _iter !=
                   this->_ypos.tend() ;
                ++_iter , ++_inum)
        {
            iptr_type  _left, _pair ;
            indx_from_subs(
                _inum,  JBEG, _left);
            indx_from_subs(
                _inum,  JEND, _pair);

            vals_type _hmin = std::min(
                this->_hmat[_left],
                this->_hmat[_pair]) ;

            this->_hmat[_left] = _hmin;
            this->_hmat[_pair] = _hmin;
        }
        }

    /*-------------------- push nodes onto priority queue */
        _keys.set_count (
            _hmat.count(),
        containers::tight_alloc, _null) ;

        {
        iptr_type _inum  = +0;
        for (auto _iter  =
                   this->_hmat.head() ;
                  _iter !=
                   this->_hmat.tend() ;
                ++_iter , ++_inum)
        {
            _keys[_inum] = _sort.push(_inum);
        }
        }

    /*-------------------- prebuild lon-lat to R^3 coeff. */
        this->_sinx.set_count(_xpos.count (),
            containers::tight_alloc, 0.) ;
        this->_cosx.set_count(_xpos.count (),
            containers::tight_alloc, 0.) ;

        {
        iptr_type _inum  = +0;
        for (auto _iter  =
                   this->_xpos.head() ;
                  _iter !=
                   this->_xpos.tend() ;
                ++_iter , ++_inum)
        {
            this->_sinx[_inum] =
                std::sin(this->_xpos[_inum]);
            this->_cosx[_inum] =
                std::cos(this->_xpos[_inum]);
        }
        }

        this->_siny.set_count(_ypos.count (),
            containers::tight_alloc, 0.) ;
        this->_cosy.set_count(_ypos.count (),
            containers::tight_alloc, 0.) ;

        {
        iptr_type _inum  = +0;
        for (auto _iter  =
                   this->_ypos.head() ;
                  _iter !=
                   this->_ypos.tend() ;
                ++_iter , ++_inum)
        {
            this->_siny[_inum] =
                std::sin(this->_ypos[_inum]);
            this->_cosy[_inum] =
                std::cos(this->_ypos[_inum]);
        }
        }

    /*-------------------- compute h(x) via fast-marching */
        for ( ; !_sort.empty() ; )
        {
            iptr_type _base ;
            _sort._pop_root(_base) ;

            _keys[_base] = _null ;

            iptr_type  _ipos, _jpos ;
            subs_from_indx(
                _base, _ipos, _jpos);

            vals_type _hnow = _hmat[_base] ;

            for (auto _IPOS = _ipos - 1 ;
                      _IPOS < _ipos + 1 ;
                    ++_IPOS )
            for (auto _JPOS = _jpos - 1 ;
                      _JPOS < _jpos + 1 ;
                    ++_JPOS )
            {
                if (_IPOS >= IBEG && _IPOS < IEND)
                if (_JPOS >= JBEG && _JPOS < JEND)
                {
    /*-------------------- un-pack implicit cell indexing */
                 auto _ipii = _IPOS + 0 ;
                 auto _ipjj = _JPOS + 0 ;

                iptr_type  _inod;
                indx_from_subs(
                    _ipii, _ipjj, _inod);

                 auto _jpii = _IPOS + 1 ;
                 auto _jpjj = _JPOS + 0 ;

                iptr_type  _jnod;
                indx_from_subs(
                    _jpii, _jpjj, _jnod);

                 auto _kpii = _IPOS + 1 ;
                 auto _kpjj = _JPOS + 1 ;

                iptr_type  _knod;
                indx_from_subs(
                    _kpii, _kpjj, _knod);

                 auto _lpii = _IPOS + 0 ;
                 auto _lpjj = _JPOS + 1 ;

                iptr_type  _lnod;
                indx_from_subs(
                    _lpii, _lpjj, _lnod);

    /*-------------------- skip any cells with null nodes */
                if (_keys[_inod] == _null &&
                    _inod != _base) continue ;
                if (_keys[_jnod] == _null &&
                    _jnod != _base) continue ;
                if (_keys[_knod] == _null &&
                    _knod != _base) continue ;
                if (_keys[_lnod] == _null &&
                    _lnod != _base) continue ;

    /*-------------------- skip cells due to sorted order */
                vals_type _hmax;
                _hmax = this->_hmat[_inod] ;
                _hmax = std::max(
                _hmax , this->_hmat[_jnod]);
                _hmax = std::max(
                _hmax , this->_hmat[_knod]);
                _hmax = std::max(
                _hmax , this->_hmat[_lnod]);

                if (_hmax <= _hnow) continue ;

    /*-------------------- set-up cell vertex coordinates */
                real_type  _IXYZ[3] ;
                _IXYZ[0] =  this->_radA *
                    this->_cosx[_ipjj] *
                    this->_cosy[_ipii] ;
                _IXYZ[1] =  this->_radB *
                    this->_sinx[_ipjj] *
                    this->_cosy[_ipii] ;
                _IXYZ[2] =  this->_radC *
                    this->_siny[_ipii] ;

                real_type  _JXYZ[3] ;
                _JXYZ[0] =  this->_radA *
                    this->_cosx[_jpjj] *
                    this->_cosy[_jpii] ;
                _JXYZ[1] =  this->_radB *
                    this->_sinx[_jpjj] *
                    this->_cosy[_jpii] ;
                _JXYZ[2] =  this->_radC *
                    this->_siny[_jpii] ;

                real_type  _KXYZ[3] ;
                _KXYZ[0] =  this->_radA *
                    this->_cosx[_kpjj] *
                    this->_cosy[_kpii] ;
                _KXYZ[1] =  this->_radB *
                    this->_sinx[_kpjj] *
                    this->_cosy[_kpii] ;
                _KXYZ[2] =  this->_radC *
                    this->_siny[_kpii] ;

                real_type  _LXYZ[3] ;
                _LXYZ[0] =  this->_radA *
                    this->_cosx[_lpjj] *
                    this->_cosy[_lpii] ;
                _LXYZ[1] =  this->_radB *
                    this->_sinx[_lpjj] *
                    this->_cosy[_lpii] ;
                _LXYZ[2] =  this->_radC *
                    this->_siny[_lpii] ;

    /*-------------------- solve for local |dh/dx| limits */
                vals_type _iold =
                     this->_hmat[_inod] ;
                vals_type _jold =
                     this->_hmat[_jnod] ;
                vals_type _kold =
                     this->_hmat[_knod] ;
                vals_type _lold =
                     this->_hmat[_lnod] ;

                if (this->_dhdx.count() >1)
                {
    /*-------------------- update adj. set, g = g(x) case */
                if (eikonal_quad_3d (
                   _IXYZ , _JXYZ ,
                   _KXYZ , _LXYZ ,
                    this->_hmat[_inod],
                    this->_hmat[_jnod],
                    this->_hmat[_knod],
                    this->_hmat[_lnod],
                    this->_dhdx[_inod],
                    this->_dhdx[_jnod],
                    this->_dhdx[_knod],
                    this->_dhdx[_lnod])  )
                {

                if (_keys[_inod] != _null)
                if (_hmat[_inod] != _iold)
                    _sort.update(
                    _keys[_inod] ,  _inod) ;

                if (_keys[_jnod] != _null)
                if (_hmat[_jnod] != _jold)
                    _sort.update(
                    _keys[_jnod] ,  _jnod) ;

                if (_keys[_knod] != _null)
                if (_hmat[_knod] != _kold)
                    _sort.update(
                    _keys[_knod] ,  _knod) ;

                if (_keys[_lnod] != _null)
                if (_hmat[_lnod] != _lold)
                    _sort.update(
                    _keys[_lnod] ,  _lnod) ;

                if (this->_wrap)
                {
    /*-------------------- update adj. set, periodic bc's */
                if (_ipjj==JBEG)
                {
                iptr_type _pair;
                indx_from_subs(
                    _ipii, JEND  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_inod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_inod] != _iold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }

                if (_jpjj==JEND)
                {
                iptr_type _pair;
                indx_from_subs(
                    _jpii, JBEG  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_jnod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_jnod] != _jold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }

                if (_kpjj==JEND)
                {
                iptr_type _pair;
                indx_from_subs(
                    _kpii, JBEG  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_knod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_knod] != _kold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }

                if (_lpjj==JBEG)
                {
                iptr_type _pair;
                indx_from_subs(
                    _lpii, JEND  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_lnod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_lnod] != _lold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }
                }

                }
                }
                else
                if (this->_dhdx.count()==1)
                {
    /*-------------------- update adj. set, const. g case */
                if (eikonal_quad_3d (
                   _IXYZ , _JXYZ ,
                   _KXYZ , _LXYZ ,
                    this->_hmat[_inod],
                    this->_hmat[_jnod],
                    this->_hmat[_knod],
                    this->_hmat[_lnod],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ])  )
                {

                if (_keys[_inod] != _null)
                if (_hmat[_inod] != _iold)
                    _sort.update(
                    _keys[_inod] ,  _inod) ;

                if (_keys[_jnod] != _null)
                if (_hmat[_jnod] != _jold)
                    _sort.update(
                    _keys[_jnod] ,  _jnod) ;

                if (_keys[_knod] != _null)
                if (_hmat[_knod] != _kold)
                    _sort.update(
                    _keys[_knod] ,  _knod) ;

                if (_keys[_lnod] != _null)
                if (_hmat[_lnod] != _lold)
                    _sort.update(
                    _keys[_lnod] ,  _lnod) ;

                if (this->_wrap)
                {
    /*-------------------- update adj. set, periodic bc's */
                if (_ipjj==JBEG)
                {
                iptr_type _pair;
                indx_from_subs(
                    _ipii, JEND  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_inod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_inod] != _iold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }

                if (_jpjj==JEND)
                {
                iptr_type _pair;
                indx_from_subs(
                    _jpii, JBEG  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_jnod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_jnod] != _jold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }

                if (_kpjj==JEND)
                {
                iptr_type _pair;
                indx_from_subs(
                    _kpii, JBEG  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_knod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_knod] != _kold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }

                if (_lpjj==JBEG)
                {
                iptr_type _pair;
                indx_from_subs(
                    _lpii, JEND  ,  _pair) ;

                this->_hmat [_pair] =
                    this->_hmat [_lnod] ;

                if (_keys[_pair] != _null)
                if (_hmat[_lnod] != _lold)
                    _sort.update(
                    _keys[_pair] ,  _pair) ;
                }
                }

                }
                }

                }
            }
        }

    }

    /*
    --------------------------------------------------------
     * EVAL: eval. size-fun. value.
    --------------------------------------------------------
     */

    __normal_call real_type eval (
        real_type *_ppos ,
        hint_type &_hint
        )
    {
        real_type static const FT =
            std::pow (
        +std::numeric_limits <real_type>
            ::epsilon(),(real_type)+.8);

        real_type _hval =
            +std::numeric_limits
                <real_type>::infinity();

        __unreferenced(_hint) ;

        if (this->_xpos.count() == +0)
            return _hval ;

        if (this->_ypos.count() == +0)
            return _hval ;

    /*---------------------------- compute xyz to lat-lon */
        real_type _apos[2] ;
        if (this->_radA == +1. &&
            this->_radB == +1. &&
            this->_radC == +1.)
        {
    /*---------------------------- "historical" workflow! */
        real_type _radius = std::sqrt (
            _ppos[0]*_ppos[0] +
            _ppos[1]*_ppos[1] +
            _ppos[2]*_ppos[2] ) ;

        if (_radius < FT)
        _radius  = _radius + FT ;

        _apos[1] =
        std::asin (_ppos[2]/ _radius ) ;

        _apos[0] =
        std::atan2(_ppos[1], _ppos[0]) ;
        }
        else
        {
    /*---------------------------- for a "full" ellipsoid */
            toS2(_ppos, _apos);
        }

        real_type _alon = _apos[ 0] ;
        real_type _alat = _apos[ 1] ;

        iptr_type _pole = (iptr_type)0 ;
        real_type _bias = (real_type)0.;

        real_type static const PI =
       (real_type)std::atan(+1.0) * +4.;

        real_type static const PI_h =
       (real_type)+.5 * PI ;
        real_type static const PI_1 =
       (real_type)+1. * PI ;
        real_type static const PI_2 =
       (real_type)+2. * PI ;

        if (_alon<-PI_1) _alon += PI_2 ;
        if (_alon>=PI_1) _alon -= PI_2 ;

        if (_alon < *this->_xpos.head())
            _alon = *this->_xpos.head();
        if (_alon > *this->_xpos.tail())
            _alon = *this->_xpos.tail();

        if (_alat<-PI_h) _alat  =-PI_h ;
        if (_alat>=PI_h) _alat  = PI_h ;

        if (_alat < *this->_ypos.head())
        {
            _bias =             // deal w interp. at pole
           (_alat - *this->_ypos.head()) /
           (-PI_h - *this->_ypos.head()) ;
            _pole = -1 ;
            _alat = *this->_ypos.head();

            _bias = std::sin(_bias*PI_h) ;
        }

        if (_alat > *this->_ypos.tail())
        {
            _bias =             // deal w interp. at pole
           (_alat - *this->_ypos.tail()) /
           ( PI_h - *this->_ypos.tail()) ;
            _pole = +1 ;
            _alat = *this->_ypos.tail();

            _bias = std::sin(_bias*PI_h) ;
        }

    /*---------------------------- find enclosing x-range */
        iptr_type _ipos = (iptr_type) -1 ;
        iptr_type _jpos = (iptr_type) -1 ;

        if (this->_xvar)
        {
            auto _joff =
            algorithms::upper_bound (
                this->_xpos.head(),
                this->_xpos.tend(),
            _alon,std::less<real_type>() ) ;

            _jpos = (iptr_type) (
            _joff - this->_xpos.head()- 1) ;
        }
        else
        {
            real_type _xmin, _xmax, _xdel;
            _xmin = *this->_xpos.head();
            _xmax = *this->_xpos.tail();

            _xdel = (_xmax - _xmin) /
                (this->_xpos.count () - 1) ;

            _jpos = (iptr_type)
                ((_alon-_xmin)/_xdel) ;
        }

    /*---------------------------- find enclosing y-range */
        if (this->_yvar)
        {
            auto _ioff =
            algorithms::upper_bound (
                this->_ypos.head(),
                this->_ypos.tend(),
            _alat,std::less<real_type>() ) ;

            _ipos = (iptr_type) (
            _ioff - this->_ypos.head()- 1) ;
        }
        else
        {
            real_type _ymin, _ymax, _ydel;
            _ymin = *this->_ypos.head();
            _ymax = *this->_ypos.tail();

            _ydel = (_ymax - _ymin) /
                (this->_ypos.count () - 1) ;

            _ipos = (iptr_type)
                ((_alat-_ymin)/_ydel) ;
        }

        if (_ipos ==
       (iptr_type)this->_ypos.count() - 1)
            _ipos = _ipos - 1 ;

        if (_jpos ==
       (iptr_type)this->_xpos.count() - 1)
            _jpos = _jpos - 1 ;

    /*---------------------------- a linear interpolation */
        real_type _xp11 =
            this->_xpos[_jpos + 0] ;
        real_type _xp22 =
            this->_xpos[_jpos + 1] ;

        real_type _yp11 =
            this->_ypos[_ipos + 0] ;
        real_type _yp22 =
            this->_ypos[_ipos + 1] ;

        real_type _xval = _alon ;
        real_type _yval = _alat ;

        real_type _aa22 =
           (_yval-_yp11) * (_xval-_xp11) ;
        real_type _aa21 =
           (_yval-_yp11) * (_xp22-_xval) ;
        real_type _aa12 =
           (_yp22-_yval) * (_xval-_xp11) ;
        real_type _aa11 =
           (_yp22-_yval) * (_xp22-_xval) ;

        iptr_type _kk11 ;
        indx_from_subs(
            _ipos + 0, _jpos + 0, _kk11) ;
        iptr_type _kk12 ;
        indx_from_subs(
            _ipos + 0, _jpos + 1, _kk12) ;
        iptr_type _kk21 ;
        indx_from_subs(
            _ipos + 1, _jpos + 0, _kk21) ;
        iptr_type _kk22 ;
        indx_from_subs(
            _ipos + 1, _jpos + 1, _kk22) ;

        real_type _hbar =
          ( _aa11*this->_hmat[_kk11]
          + _aa12*this->_hmat[_kk12]
          + _aa21*this->_hmat[_kk21]
          + _aa22*this->_hmat[_kk22] )
        / ( _aa11+_aa12+_aa21+_aa22) ;

        if (_pole >= +1)
        {
            return              // deal w interp. at pole
        ((real_type)1. - _bias) * _hbar +
        ((real_type)0. + _bias) * _htop ;
        }
        else
        if (_pole <= -1)
        {
            return              // deal w interp. at pole
        ((real_type)1. - _bias) * _hbar +
        ((real_type)0. + _bias) * _hbot ;
        }
        else { return _hbar ; } // return interior interp
    }

    } ;


    }

#   endif   //__HFUN_GRID_ELLIPSOID_3__



