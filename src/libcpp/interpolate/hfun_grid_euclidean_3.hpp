
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
     * Last updated: 20 Oct., 2024
     *
     * Copyright 2013-2024
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
        ) const
    {
    /*------------ helper: convert into "un-rolled" index */
        iptr_type _ynum =
       (iptr_type)this->_ypos.count() ;

        iptr_type _xnum =
       (iptr_type)this->_xpos.count() ;

        _indx = _kpos * _xnum * _ynum
              + _jpos * _ynum + _ipos ;
    }
    
    __inline_call void_type subs_from_indx (
        iptr_type _indx,
        iptr_type&_ipos,
        iptr_type&_jpos,
        iptr_type&_kpos
        ) const
    {
    /*------------ helper: convert from "un-rolled" index */
        iptr_type _ynum =
       (iptr_type)this->_ypos.count() ;

        iptr_type _xnum =
       (iptr_type)this->_xpos.count() ;
       
        _kpos = _indx /(_xnum * _ynum);
        _indx-= _kpos * _xnum * _ynum ;
        _jpos = _indx / _ynum ;
        _ipos = _indx - _jpos * _ynum ;
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
        for (auto _iter  = this->_hmat.head() ;
                  _iter != this->_hmat.tend() ;
                ++_iter, ++_inum)
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
       
        iptr_type KBEG = +0;
        iptr_type KEND =
       (iptr_type)this->_zpos.count() - 1 ;

    #   define ISALIVE(__idx)           \
           (_sort.keys(__idx) != _sort.null())

    #   define UPDATED(__new, __old)    \
        std::abs(__new - __old) > _FTOL*std::abs(__new)

    #   define REQUEUE(__new, __old, __idx)     \
            if(UPDATED(__new, __old)) {         \
            if(ISALIVE(__idx))          \
            {                           \
                _hmat[__idx]  = __new;  \
                _sort.reduce(__idx , __idx) ;   \
            }                           \
            else                        \
            {                           \
                _hmat[__idx]  = __new;  \
                _sort.push  (__idx , __idx) ;   \
            } }

        for ( ; !_sort.empty() ; )
        {
            iptr_type  _base, _bidx;
            _sort._pop_root ( _bidx, _base) ;

            iptr_type  _ipos, _jpos, _kpos;
            subs_from_indx(
                _base, _ipos, _jpos, _kpos) ;

            vals_type _hnow = _hmat [_base] ;
            
            for (auto _IPOS = _ipos - 1 ;
                      _IPOS < _ipos + 1 ;
                    ++_IPOS )
            for (auto _JPOS = _jpos - 1 ;
                      _JPOS < _jpos + 1 ;
                    ++_JPOS )
            for (auto _KPOS = _kpos - 1 ;
                      _KPOS < _kpos + 1 ;
                    ++_KPOS )
            {
                if (_IPOS >= IBEG && _IPOS < IEND)
                if (_JPOS >= JBEG && _JPOS < JEND)
                if (_KPOS >= KBEG && _KPOS < KEND)
                {
    /*-------------------- un-pack implicit cell indexing */
                 auto _ipii = _IPOS + 0 ;
                 auto _ipjj = _JPOS + 0 ;
                 auto _ipkk = _KPOS + 0 ;

                iptr_type  _inod;
                indx_from_subs(
                    _ipii, _ipjj, _ipkk, _inod);
                    
                 auto _jpii = _IPOS + 1 ;
                 auto _jpjj = _JPOS + 0 ;
                 auto _jpkk = _KPOS + 0 ;

                iptr_type  _jnod;
                indx_from_subs(
                    _jpii, _jpjj, _jpkk, _jnod);
                    
                 auto _kpii = _IPOS + 1 ;
                 auto _kpjj = _JPOS + 1 ;
                 auto _kpkk = _KPOS + 0 ;

                iptr_type  _knod;
                indx_from_subs(
                    _kpii, _kpjj, _kpkk, _knod);
                    
                 auto _lpii = _IPOS + 0 ;
                 auto _lpjj = _JPOS + 1 ;
                 auto _lpkk = _KPOS + 0 ;

                iptr_type  _lnod;
                indx_from_subs(
                    _lpii, _lpjj, _lpkk, _lnod);
                     
                 auto _mpii = _IPOS + 0 ;
                 auto _mpjj = _JPOS + 0 ;
                 auto _mpkk = _KPOS + 1 ;

                iptr_type  _mnod;
                indx_from_subs(
                    _mpii, _mpjj, _mpkk, _mnod);
                    
                 auto _npii = _IPOS + 1 ;
                 auto _npjj = _JPOS + 0 ;
                 auto _npkk = _KPOS + 1 ;

                iptr_type  _nnod;
                indx_from_subs(
                    _npii, _npjj, _npkk, _nnod);
                    
                 auto _opii = _IPOS + 1 ;
                 auto _opjj = _JPOS + 1 ;
                 auto _opkk = _KPOS + 1 ;

                iptr_type  _onod;
                indx_from_subs(
                    _opii, _opjj, _opkk, _onod);
                    
                 auto _ppii = _IPOS + 0 ;
                 auto _ppjj = _JPOS + 1 ;
                 auto _ppkk = _KPOS + 1 ;

                iptr_type  _pnod;
                indx_from_subs(
                    _ppii, _ppjj, _ppkk, _pnod);
                    
    /*-------------------- skip cells due to sorted order */
                iptr_type _near = +0;
                if (_inod != _base &&
                    ISALIVE(_inod)) _near++;
                if (_jnod != _base &&
                    ISALIVE(_jnod)) _near++;
                if (_knod != _base &&
                    ISALIVE(_knod)) _near++;
                if (_lnod != _base &&
                    ISALIVE(_lnod)) _near++;
                if (_mnod != _base &&
                    ISALIVE(_mnod)) _near++;
                if (_nnod != _base &&
                    ISALIVE(_nnod)) _near++;
                if (_onod != _base &&
                    ISALIVE(_onod)) _near++;
                if (_pnod != _base &&
                    ISALIVE(_pnod)) _near++;
                
                if (_near == 0)     continue ;

                vals_type _hmax = .0;
                _hmax = this->_hmat[_inod] ;
                _hmax = std::max(
                _hmax , this->_hmat[_jnod]);
                _hmax = std::max(
                _hmax , this->_hmat[_knod]);
                _hmax = std::max(
                _hmax , this->_hmat[_lnod]);
                _hmax = std::max(
                _hmax , this->_hmat[_mnod]);
                _hmax = std::max(
                _hmax , this->_hmat[_nnod]);
                _hmax = std::max(
                _hmax , this->_hmat[_onod]);
                _hmax = std::max(
                _hmax , this->_hmat[_pnod]);

                if (_hmax <= _hnow) continue ;
                
    /*-------------------- set-up cell vertex coordinates */
                real_type _IXYZ[3];
                _IXYZ[0] = this->_xpos[_ipjj];
                _IXYZ[1] = this->_ypos[_ipii];
                _IXYZ[2] = this->_zpos[_ipkk];

                real_type _JXYZ[3];
                _JXYZ[0] = this->_xpos[_jpjj];
                _JXYZ[1] = this->_ypos[_jpii];
                _JXYZ[2] = this->_zpos[_jpkk];

                real_type _KXYZ[3];
                _KXYZ[0] = this->_xpos[_kpjj];
                _KXYZ[1] = this->_ypos[_kpii];
                _KXYZ[2] = this->_zpos[_kpkk];

                real_type _LXYZ[3];
                _LXYZ[0] = this->_xpos[_lpjj];
                _LXYZ[1] = this->_ypos[_lpii];
                _LXYZ[2] = this->_zpos[_lpkk];
                
                real_type _MXYZ[3];
                _MXYZ[0] = this->_xpos[_mpjj];
                _MXYZ[1] = this->_ypos[_mpii];
                _MXYZ[2] = this->_zpos[_mpkk];

                real_type _NXYZ[3];
                _NXYZ[0] = this->_xpos[_npjj];
                _NXYZ[1] = this->_ypos[_npii];
                _NXYZ[2] = this->_zpos[_npkk];      
                
                real_type _OXYZ[3];
                _OXYZ[0] = this->_xpos[_opjj];
                _OXYZ[1] = this->_ypos[_opii];
                _OXYZ[2] = this->_zpos[_opkk];
                
                real_type _PXYZ[3];
                _PXYZ[0] = this->_xpos[_ppjj];
                _PXYZ[1] = this->_ypos[_ppii];
                _PXYZ[2] = this->_zpos[_ppkk];   
                
    /*-------------------- solve for local |dh/dx| limits */
                vals_type _iold =
                     this->_hmat[_inod] ;
                vals_type _jold =
                     this->_hmat[_jnod] ;
                vals_type _kold =
                     this->_hmat[_knod] ;
                vals_type _lold =
                     this->_hmat[_lnod] ;
                vals_type _mold =
                     this->_hmat[_mnod] ;
                vals_type _nold =
                     this->_hmat[_nnod] ;
                vals_type _oold =
                     this->_hmat[_onod] ;
                vals_type _pold =
                     this->_hmat[_pnod] ;

                vals_type _inew =
                     this->_hmat[_inod] ;
                vals_type _jnew =
                     this->_hmat[_jnod] ;
                vals_type _knew =
                     this->_hmat[_knod] ;
                vals_type _lnew =
                     this->_hmat[_lnod] ;    
                vals_type _mnew =
                     this->_hmat[_mnod] ;
                vals_type _nnew =
                     this->_hmat[_nnod] ;        
                vals_type _onew =
                     this->_hmat[_onod] ;
                vals_type _pnew =
                     this->_hmat[_pnod] ;
                     
                if (this->_dhdx.count() > +1)
                {
    /*-------------------- update adj. set, g = g(x) case */
                if (eikonal_hexa_3d (
                   _IXYZ , _JXYZ , _KXYZ , _LXYZ , 
                   _MXYZ , _NXYZ , _OXYZ , _PXYZ ,  
                   _hnow ,
                   _inew , _jnew , _knew , _lnew ,
                   _mnew , _nnew , _onew , _pnew ,
                    this->_dhdx[_inod],
                    this->_dhdx[_jnod],
                    this->_dhdx[_knod],
                    this->_dhdx[_lnod],
                    this->_dhdx[_mnod],
                    this->_dhdx[_nnod],
                    this->_dhdx[_onod],
                    this->_dhdx[_pnod])  )
                {

            //  push updates one-at-a-time to ensure heap
            //  maintains its sorted order
            
                REQUEUE (_inew, _iold, _inod)
                REQUEUE (_jnew, _jold, _jnod)
                REQUEUE (_knew, _kold, _knod)
                REQUEUE (_lnew, _lold, _lnod)
                REQUEUE (_mnew, _mold, _mnod)
                REQUEUE (_nnew, _nold, _nnod)
                REQUEUE (_onew, _oold, _onod)
                REQUEUE (_pnew, _pold, _pnod)
                
                }
                }
                else
                if (this->_dhdx.count() == 1)
                {
    /*-------------------- update adj. set, const. g case */
                if (eikonal_hexa_3d (
                   _IXYZ , _JXYZ , _KXYZ , _LXYZ , 
                   _MXYZ , _NXYZ , _OXYZ , _PXYZ ,  
                   _hnow ,
                   _inew , _jnew , _knew , _lnew ,
                   _mnew , _nnew , _onew , _pnew ,
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ])  )
                {
    
            //  push updates one-at-a-time to ensure heap
            //  maintains its sorted order
            
                REQUEUE (_inew, _iold, _inod)
                REQUEUE (_jnew, _jold, _jnod)
                REQUEUE (_knew, _kold, _knod)
                REQUEUE (_lnew, _lold, _lnod)
                REQUEUE (_mnew, _mold, _mnod)
                REQUEUE (_nnew, _nold, _nnod)
                REQUEUE (_onew, _oold, _onod)
                REQUEUE (_pnew, _pold, _pnod)
            
                }
                }
    
                }
            }
        }      
        
    #   undef ISALIVE
    #   undef UPDATED
    #   undef REQUEUE
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
        real_type *_ppos
        )
    {
        auto _hint = this->null_hint();
        return eval(_ppos, _hint);
    }

    __normal_call real_type eval (
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



