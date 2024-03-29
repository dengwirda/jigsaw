
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

#   ifndef __HFUN_GRID_EUCLIDEAN_2__
#   define __HFUN_GRID_EUCLIDEAN_2__

    namespace mesh {

    template <
    typename R ,
    typename V ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_grid_euclidean_2d
        : public hfun_base_kd <I, R>
    {
    public  :

    /*---------------------- "grid"-based size-fun in R^2 */

    typedef R                       real_type ;
    typedef V                       vals_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef hfun_grid_euclidean_2d  <
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
        vals_type, allocator>      _hmat;

    containers::array <
        vals_type, allocator>      _dhdx;

    bool_type                      _xvar;
    bool_type                      _yvar;

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
        typename vals_list::_write_it _hptr ;

        public  :
        __inline_call less_than  (
        typename vals_list::_write_it _hsrc
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

        containers::priorityidx <
            iptr_type ,
            iptr_type ,
            less_than ,
            allocator >
        _sort((less_than(this->_hmat.head())));

    /*-------------------- push nodes onto priority queue */
        iptr_type _inum  = +0;
        for (auto _iter  =
                   this->_hmat.head() ;
                  _iter !=
                   this->_hmat.tend() ;
                ++_iter , ++_inum)
        {
            _sort.push(_inum, _inum ) ;
        }

    /*-------------------- compute h(x) via fast-marching */
        vals_type static _FTOL =
            (vals_type)std::pow(
        std::numeric_limits<vals_type>::epsilon(), .75) ;

        iptr_type IBEG = +0;
        iptr_type IEND =
       (iptr_type)this->_ypos.count() - 1 ;

        iptr_type JBEG = +0;
        iptr_type JEND =
       (iptr_type)this->_xpos.count() - 1 ;

    #   define ISALIVE(__idx)           \
           (_sort.keys(__idx) != _sort.null())

    #   define UPDATED(__new, __old)    \
        std::abs(__new - __old) > _FTOL*std::abs(__new)

        for ( ; !_sort.empty() ; )
        {
            iptr_type _base, _bidx;
            _sort._pop_root( _bidx, _base) ;

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

    /*-------------------- skip cells due to sorted order */
                if (_inod != _base &&
                   !ISALIVE(_inod)) continue ;
                if (_jnod != _base &&
                   !ISALIVE(_jnod)) continue ;
                if (_knod != _base &&
                   !ISALIVE(_knod)) continue ;
                if (_lnod != _base &&
                   !ISALIVE(_lnod)) continue ;

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
                real_type _IXYZ[2];
                _IXYZ[0] = this->_xpos[_ipjj];
                _IXYZ[1] = this->_ypos[_ipii];

                real_type _JXYZ[2];
                _JXYZ[0] = this->_xpos[_jpjj];
                _JXYZ[1] = this->_ypos[_jpii];

                real_type _KXYZ[2];
                _KXYZ[0] = this->_xpos[_kpjj];
                _KXYZ[1] = this->_ypos[_kpii];

                real_type _LXYZ[2];
                _LXYZ[0] = this->_xpos[_lpjj];
                _LXYZ[1] = this->_ypos[_lpii];

    /*-------------------- solve for local |dh/dx| limits */
                vals_type _iold =
                     this->_hmat[_inod] ;
                vals_type _jold =
                     this->_hmat[_jnod] ;
                vals_type _kold =
                     this->_hmat[_knod] ;
                vals_type _lold =
                     this->_hmat[_lnod] ;

                vals_type _inew =
                     this->_hmat[_inod] ;
                vals_type _jnew =
                     this->_hmat[_jnod] ;
                vals_type _knew =
                     this->_hmat[_knod] ;
                vals_type _lnew =
                     this->_hmat[_lnod] ;

                if (this->_dhdx.count() >1)
                {
    /*-------------------- update adj. set, g = g(x) case */
                if (eikonal_grid_2d (
                   _IXYZ , _JXYZ ,
                   _KXYZ , _LXYZ , _hnow ,
                   _inew , _jnew ,
                   _knew , _lnew ,
                    this->_dhdx[_inod],
                    this->_dhdx[_jnod],
                    this->_dhdx[_knod],
                    this->_dhdx[_lnod])  )
                {

            //  push updates one-at-a-time to ensure heap
            //  maintains its sorted order

                if (_sort.
                     keys(_inod) != _sort.null())
                if ( UPDATED(_inew, _iold) )
                {
                    _hmat[_inod]  = _inew;
                    _sort.reduce(_inod , _inod) ;
                }

                if (_sort.
                     keys(_jnod) != _sort.null())
                if ( UPDATED(_jnew, _jold) )
                {
                    _hmat[_jnod]  = _jnew;
                    _sort.reduce(_jnod , _jnod) ;
                }

                if (_sort.
                     keys(_knod) != _sort.null())
                if ( UPDATED(_knew, _kold) )
                {
                    _hmat[_knod]  = _knew;
                    _sort.reduce(_knod , _knod) ;
                }

                if (_sort.
                     keys(_lnod) != _sort.null())
                if ( UPDATED(_lnew, _lold) )
                {
                    _hmat[_lnod]  = _lnew;
                    _sort.reduce(_lnod , _lnod) ;
                }

                }
                }
                else
                if (this->_dhdx.count()==1)
                {
    /*-------------------- update adj. set, const. g case */
                if (eikonal_grid_2d (
                   _IXYZ , _JXYZ ,
                   _KXYZ , _LXYZ , _hnow ,
                   _inew , _jnew ,
                   _knew , _lnew ,
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ])  )
                {

            //  push updates one-at-a-time to ensure heap
            //  maintains its sorted order

                if (_sort.
                     keys(_inod) != _sort.null())
                if ( UPDATED(_inew, _iold) )
                {
                    _hmat[_inod]  = _inew;
                    _sort.reduce(_inod , _inod) ;
                }

                if (_sort.
                     keys(_jnod) != _sort.null())
                if ( UPDATED(_jnew, _jold) )
                {
                    _hmat[_jnod]  = _jnew;
                    _sort.reduce(_jnod , _jnod) ;
                }

                if (_sort.
                     keys(_knod) != _sort.null())
                if ( UPDATED(_knew, _kold) )
                {
                    _hmat[_knod]  = _knew;
                    _sort.reduce(_knod , _knod) ;
                }

                if (_sort.
                     keys(_lnod) != _sort.null())
                if ( UPDATED(_lnew, _lold) )
                {
                    _hmat[_lnod]  = _lnew;
                    _sort.reduce(_lnod , _lnod) ;
                }

                }
                }

                }
            }
        }

    #   undef ISALIVE
    #   undef UPDATED
    }

    /*
    --------------------------------------------------------
     * INIT: init. size-fun. class.
    --------------------------------------------------------
     */

    __inline_call void_type init (
        float _xoff = + 0.f ,
        float _yoff = + 0.f
        )
    {
        real_type static const _FTOL =
            std::pow (
            std::numeric_limits<real_type>
            ::epsilon(), (real_type) +.75) ;

        this->_xvar = false ;
        this->_yvar = false ;

        if (this->_xpos.empty()) return ;
        if (this->_ypos.empty()) return ;

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

    /*---------------------------- find enclosing x-range */
        auto _ipos = this-> null_hint() ;
        auto _jpos = this-> null_hint() ;

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

        if (_ipos ==
       (iptr_type)this->_ypos.count() - 1)
            _ipos = _ipos - 1 ;

        if (_jpos ==
       (iptr_type)this->_xpos.count() - 1)
            _jpos = _jpos - 1 ;

    /*---------------------------- a linear interpolation */
        real_type _xx11 =
            this->_xpos[_jpos + 0] ;
        real_type _xx22 =
            this->_xpos[_jpos + 1] ;

        real_type _yy11 =
            this->_ypos[_ipos + 0] ;
        real_type _yy22 =
            this->_ypos[_ipos + 1] ;

        real_type _aa22 =
           (_YPOS-_yy11) * (_XPOS-_xx11) ;
        real_type _aa21 =
           (_YPOS-_yy11) * (_xx22-_XPOS) ;
        real_type _aa12 =
           (_yy22-_YPOS) * (_XPOS-_xx11) ;
        real_type _aa11 =
           (_yy22-_YPOS) * (_xx22-_XPOS) ;

        iptr_type _kk11, _kk12 ,
                  _kk21, _kk22 ;
        indx_from_subs(
            _ipos + 0, _jpos + 0, _kk11) ;
        indx_from_subs(
            _ipos + 0, _jpos + 1, _kk12) ;
        indx_from_subs(
            _ipos + 1, _jpos + 0, _kk21) ;
        indx_from_subs(
            _ipos + 1, _jpos + 1, _kk22) ;

        real_type _hBAR =
          ( _aa11*this->_hmat[_kk11]
          + _aa12*this->_hmat[_kk12]
          + _aa21*this->_hmat[_kk21]
          + _aa22*this->_hmat[_kk22] )
        / ( _aa11+_aa12+_aa21+_aa22) ;

        return (  _hBAR ) ;
    }

    } ;


    }

#   endif   //__HFUN_GRID_EUCLIDEAN_2__



