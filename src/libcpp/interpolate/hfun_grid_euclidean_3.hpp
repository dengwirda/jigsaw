
    /*
    --------------------------------------------------------
     * HFUN-GRID-EUCLIDEAN-kD: structured H(X) in R^k.
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
     * Last updated: 15 Jun., 2022
     *
     * Copyright 2013-2022
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_GRID_EUCLIDEAN_3__
#   define __HFUN_GRID_EUCLIDEAN_3__

    namespace mesh {

    template <
    typename R ,
    typename V ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_grid_euclidean_3d
        : public hfun_base_kd <I, R>
    {
    public  :

    /*---------------------- "grid"-based size-fun in R^3 */

    typedef R                       real_type ;
    typedef V                       vals_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef hfun_grid_euclidean_3d  <
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


    containers::array <
        real_type, allocator>      _xpos;
    containers::array <
        real_type, allocator>      _ypos;
    containers::array <
        real_type, allocator>      _zpos;

    containers::array <
        vals_type, allocator>      _hmat;

    containers::array <
        vals_type, allocator>      _dhdx;

    bool_type                      _xvar;
    bool_type                      _yvar;
    bool_type                      _zvar;

    public  :

    __inline_call void_type indx_from_subs (
        iptr_type _ipos,
        iptr_type _jpos,
        iptr_type _kpos,
        iptr_type&_indx
        )
    {
    /*------------ helper: convert into "un-rolled" index */
        iptr_type _ynum =
       (iptr_type)this->_ypos.count() ;

        iptr_type _xnum =
       (iptr_type)this->_xpos.count() ;

        _indx = _kpos * _xnum * _ynum
              + _jpos * _ynum + _ipos ;
    }

    /*
    --------------------------------------------------------
     * INIT: init. size-fun. class.
    --------------------------------------------------------
     */

    __inline_call void_type init (
        float _xoff = + 0.f ,
        float _yoff = + 0.f ,
        float _zoff = + 0.f
        )
    {
        real_type static const _FTOL =
            std::pow (
            std::numeric_limits<real_type>
            ::epsilon(), (real_type) +.75) ;

        this->_xvar = false ;
        this->_yvar = false ;
        this->_zvar = false ;

        if (this->_xpos.empty()) return ;
        if (this->_ypos.empty()) return ;
        if (this->_zpos.empty()) return ;

        for (auto
            _iter  = this->_xpos.head() ;
            _iter != this->_xpos.tend() ;
          ++_iter  )
        {
           *_iter -=  _xoff ;
        }

        real_type _xbar, _xmin, _xmax ;
        _xbar = *this->_xpos.tail() -
                *this->_xpos.head() ;

        _xbar /=(this->_xpos.count () - 1) ;

        _xmin = _xbar - _FTOL * _xbar ;
        _xmax = _xbar + _FTOL * _xbar ;

        for (auto
            _iter  = this->_ypos.head() ;
            _iter != this->_ypos.tend() ;
          ++_iter  )
        {
           *_iter -=  _yoff ;
        }

        real_type _ybar, _ymin, _ymax ;
        _ybar = *this->_ypos.tail() -
                *this->_ypos.head() ;

        _ybar /=(this->_ypos.count () - 1) ;

        _ymin = _ybar - _FTOL * _ybar ;
        _ymax = _ybar + _FTOL * _ybar ;

        for (auto
            _iter  = this->_zpos.head() ;
            _iter != this->_zpos.tend() ;
          ++_iter  )
        {
           *_iter -=  _zoff ;
        }

        real_type _zbar, _zmin, _zmax ;
        _zbar = *this->_zpos.tail() -
                *this->_zpos.head() ;

        _zbar /=(this->_zpos.count () - 1) ;

        _zmin = _zbar - _FTOL * _zbar ;
        _zmax = _zbar + _FTOL * _zbar ;

        for (auto
            _iter  = this->_xpos.head() ;
            _iter != this->_xpos.tail() ;
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

        for (auto
            _iter  = this->_ypos.head() ;
            _iter != this->_ypos.tail() ;
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

        for (auto
            _iter  = this->_zpos.head() ;
            _iter != this->_zpos.tail() ;
          ++_iter  )
        {
            real_type _zdel =
                *(_iter+1)-*(_iter+0) ;

            if (_zdel < _zmin ||
                _zdel > _zmax )
            {
                _zvar =  true ; break ;
            }
        }

    }

    /*
    --------------------------------------------------------
     * EVAL: eval. size-fun. value.
    --------------------------------------------------------
     */

    __inline_call real_type eval (
        real_type *_ppos,
        hint_type &_hint
        )
    {
        real_type _hOUT =
    +std::numeric_limits<real_type>::infinity();

        __unreferenced (_hint) ;

        if (this->_xpos.count() == +0)
            return _hOUT ;

        real_type _XPOS = _ppos[0] ;

        if (_XPOS < *this->_xpos.head() )
            _XPOS = *this->_xpos.head() ;
        if (_XPOS > *this->_xpos.tail() )
            _XPOS = *this->_xpos.tail() ;

        if (this->_ypos.count() == +0)
            return _hOUT ;

        real_type _YPOS = _ppos[1] ;

        if (_YPOS < *this->_ypos.head() )
            _YPOS = *this->_ypos.head() ;
        if (_YPOS > *this->_ypos.tail() )
            _YPOS = *this->_ypos.tail() ;

        if (this->_zpos.count() == +0)
            return _hOUT ;

        real_type _ZPOS = _ppos[2] ;

        if (_ZPOS < *this->_zpos.head() )
            _ZPOS = *this->_zpos.head() ;
        if (_ZPOS > *this->_zpos.tail() )
            _ZPOS = *this->_zpos.tail() ;

    /*---------------------------- find enclosing x-range */
        auto _ipos = this-> null_hint() ;
        auto _jpos = this-> null_hint() ;
        auto _kpos = this-> null_hint() ;

        if (this->_xvar == true)
        {
            auto _joff =
            algorithms::upper_bound (
                this->_xpos.head(),
                this->_xpos.tend(),
            _XPOS,std::less<real_type>());

            _jpos = (iptr_type) (
            _joff-this->_xpos.head() - 1);
        }
        else
        {
            real_type _xmin, _xmax, _xdel;
            _xmin = *this->_xpos.head();
            _xmax = *this->_xpos.tail();

            _xdel = (_xmax - _xmin) /
                (this->_xpos.count() - 1);

            _jpos = (iptr_type)
              ( (_XPOS - _xmin) / _xdel );
        }

    /*---------------------------- find enclosing y-range */
        if (this->_yvar == true)
        {
            auto _ioff =
            algorithms::upper_bound (
                this->_ypos.head(),
                this->_ypos.tend(),
            _YPOS,std::less<real_type>());

            _ipos = (iptr_type) (
            _ioff-this->_ypos.head() - 1);
        }
        else
        {
            real_type _ymin, _ymax, _ydel;
            _ymin = *this->_ypos.head();
            _ymax = *this->_ypos.tail();

            _ydel = (_ymax - _ymin) /
                (this->_ypos.count() - 1);

            _ipos = (iptr_type)
              ( (_YPOS - _ymin) / _ydel );
        }

    /*---------------------------- find enclosing z-range */
        if (this->_zvar == true )
        {
            auto _koff =
            algorithms::upper_bound (
                this->_zpos.head(),
                this->_zpos.tend(),
            _ZPOS,std::less<real_type>());

            _kpos = (iptr_type) (
            _koff-this->_zpos.head() - 1);
        }
        else
        {
            real_type _zmin, _zmax, _zdel;
            _zmin = *this->_zpos.head();
            _zmax = *this->_zpos.tail();

            _zdel = (_zmax - _zmin) /
                (this->_zpos.count() - 1);

            _kpos = (iptr_type)
              ( (_ZPOS - _zmin) / _zdel );
        }

        if (_ipos ==
       (iptr_type)this->_ypos.count() - 1)
            _ipos = _ipos - 1 ;

        if (_jpos ==
       (iptr_type)this->_xpos.count() - 1)
            _jpos = _jpos - 1 ;

        if (_kpos ==
       (iptr_type)this->_zpos.count() - 1)
            _kpos = _kpos - 1 ;

    /*---------------------------- a linear interpolation */
        real_type _xx11 =
            this->_xpos[_jpos + 0] ;
        real_type _xx22 =
            this->_xpos[_jpos + 1] ;

        real_type _yy11 =
            this->_ypos[_ipos + 0] ;
        real_type _yy22 =
            this->_ypos[_ipos + 1] ;

        real_type _zz11 =
            this->_zpos[_kpos + 0] ;
        real_type _zz22 =
            this->_zpos[_kpos + 1] ;

        real_type _v222 =(_ZPOS-_zz11) *
           (_YPOS-_yy11)*(_XPOS-_xx11) ;
        real_type _v212 =(_ZPOS-_zz11) *
           (_YPOS-_yy11)*(_xx22-_XPOS) ;
        real_type _v122 =(_ZPOS-_zz11) *
           (_yy22-_YPOS)*(_XPOS-_xx11) ;
        real_type _v112 =(_ZPOS-_zz11) *
           (_yy22-_YPOS)*(_xx22-_XPOS) ;

        real_type _v221 =(_zz22-_ZPOS) *
           (_YPOS-_yy11)*(_XPOS-_xx11) ;
        real_type _v211 =(_zz22-_ZPOS) *
           (_YPOS-_yy11)*(_xx22-_XPOS) ;
        real_type _v121 =(_zz22-_ZPOS) *
           (_yy22-_YPOS)*(_XPOS-_xx11) ;
        real_type _v111 =(_zz22-_ZPOS) *
           (_yy22-_YPOS)*(_xx22-_XPOS) ;

        iptr_type _k111, _k112, _k121,
                  _k211, _k122, _k212,
                  _k221, _k222 ;
        indx_from_subs(
        _ipos+0, _jpos+0, _kpos+0, _k111);
        indx_from_subs(
        _ipos+1, _jpos+0, _kpos+0, _k211);
        indx_from_subs(
        _ipos+0, _jpos+1, _kpos+0, _k121);
        indx_from_subs(
        _ipos+1, _jpos+1, _kpos+0, _k221);

        indx_from_subs(
        _ipos+0, _jpos+0, _kpos+1, _k112);
        indx_from_subs(
        _ipos+1, _jpos+0, _kpos+1, _k212);
        indx_from_subs(
        _ipos+0, _jpos+1, _kpos+1, _k122);
        indx_from_subs(
        _ipos+1, _jpos+1, _kpos+1, _k222);

        real_type _hBAR =
          ( _v111*this->_hmat[_k111]
          + _v112*this->_hmat[_k112]
          + _v121*this->_hmat[_k121]
          + _v211*this->_hmat[_k211]
          + _v122*this->_hmat[_k122]
          + _v212*this->_hmat[_k212]
          + _v221*this->_hmat[_k221]
          + _v222*this->_hmat[_k222] )
        / ( _v111+_v112+_v121+_v211+
            _v122+_v212+_v221+_v222) ;

        return (  _hBAR ) ;
    }


    } ;


    }

#   endif   //__HFUN_GRID_EUCLIDEAN_3__



