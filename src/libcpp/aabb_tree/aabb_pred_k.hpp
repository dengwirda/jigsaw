
/*
------------------------------------------------------------
 * AABB-PRED-K: search pred. for AABB-trees in R^k.
------------------------------------------------------------
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
------------------------------------------------------------
 *
 * Last updated: 08 December, 2019
 *
 * Copyright 2013-2019
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __AABB_PRED_K__
#   define __AABB_PRED_K__

    namespace geom_tree {


    template <
    typename R ,
    typename I
             >
    class aabb_pred_node_2
    {
/*---------------------- node-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 2;

    private :

    real_type                  _ppos [ 2] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_node_2 (
    __const_ptr  (real_type) _psrc
        )
    {
        this->_ppos[0] = _psrc[0];
        this->_ppos[1] = _psrc[1];
    }

/*---------------------------------- TRUE if intersection */
    __inline_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::node_rect_2d (
             this->_ppos ,
            _bmin, _bmax))
        {
    /*------------------------------ bbox can't intersect */
            return false ;
        }
        else
        {
            return  true ;
        }
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_node_3
    {
/*---------------------- node-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 3;

    private :

    real_type                  _ppos [ 3] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_node_3 (
    __const_ptr  (real_type) _psrc
        )
    {
        this->_ppos[0] = _psrc[0];
        this->_ppos[1] = _psrc[1];
        this->_ppos[2] = _psrc[2];
    }

/*---------------------------------- TRUE if intersection */
    __inline_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::node_rect_3d (
             this->_ppos ,
            _bmin, _bmax))
        {
    /*------------------------------ bbox can't intersect */
            return false ;
        }
        else
        {
            return  true ;
        }
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_rect_2
    {
/*---------------------- rect-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 2;

    private :

    real_type                  _rmin [ 2] ;
    real_type                  _rmax [ 2] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_rect_2 (
    __const_ptr  (real_type) _asrc ,
    __const_ptr  (real_type) _bsrc
        )
    {
        this->_rmin[0] = _asrc[0];
        this->_rmin[1] = _asrc[1];

        this->_rmax[0] = _bsrc[0];
        this->_rmax[1] = _bsrc[1];
    }

/*---------------------------------- TRUE if intersection */
    __inline_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::rect_rect_2d (
             this->_rmin ,
             this->_rmax ,
            _bmin, _bmax))
        {
    /*------------------------------ bbox can't intersect */
            return false ;
        }
        else
        {
            return  true ;
        }
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_rect_3
    {
/*---------------------- rect-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 3;

    private :

    real_type                  _rmin [ 3] ;
    real_type                  _rmax [ 3] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_rect_3 (
    __const_ptr  (real_type) _asrc ,
    __const_ptr  (real_type) _bsrc
        )
    {
        this->_rmin[0] = _asrc[0];
        this->_rmin[1] = _asrc[1];
        this->_rmin[2] = _asrc[2];

        this->_rmax[0] = _bsrc[0];
        this->_rmax[1] = _bsrc[1];
        this->_rmax[2] = _bsrc[2];
    }

/*---------------------------------- TRUE if intersection */
    __inline_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::rect_rect_3d (
             this->_rmin ,
             this->_rmax ,
            _bmin, _bmax))
        {
    /*------------------------------ bbox can't intersect */
            return false ;
        }
        else
        {
            return  true ;
        }
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_line_2
    {
/*---------------------- line-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 2;

    private :

    real_type                  _ipos [ 2] ;
    real_type                  _jpos [ 2] ;

    real_type                  _xmul [ 2] ;

    real_type                  _rmin [ 2] ;
    real_type                  _rmax [ 2] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_line_2 (
    __const_ptr  (real_type) _isrc ,
    __const_ptr  (real_type) _jsrc
        )
    {
        this->_ipos[0] = _isrc[0];
        this->_ipos[1] = _isrc[1];

        this->_jpos[0] = _jsrc[0];
        this->_jpos[1] = _jsrc[1];

        this->_rmin[0] = std::min(
              _isrc[0] , _jsrc[0]) ;
        this->_rmin[1] = std::min(
              _isrc[1] , _jsrc[1]) ;

        this->_rmax[0] = std::max(
              _isrc[0] , _jsrc[0]) ;
        this->_rmax[1] = std::max(
              _isrc[1] , _jsrc[1]) ;

        this->_xmul[0] = (real_type)1. /
            ( _jsrc[0] - _isrc[0]) ;
        this->_xmul[1] = (real_type)1. /
            ( _jsrc[1] - _isrc[1]) ;
    }

/*---------------------------------- TRUE if intersection */
    __inline_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::rect_rect_2d (
             this->_rmin ,
             this->_rmax ,
            _bmin, _bmax))
    /*------------------------------ bbox can't intersect */
            return false ;

    /*------------------------------ test if line overlap */
        real_type _aval, _bval ;
        _aval = (_bmin[0]-_ipos[0])
                *_xmul[0];
        _bval = (_bmax[0]-_ipos[0])
                *_xmul[0];

        real_type _tmin, _tmax ;
        _tmin =
        std::min( _aval, _bval )  ;
        _tmax =
        std::max( _aval, _bval )  ;

        if (_tmax<_tmin) return false ;

        _aval = (_bmin[1]-_ipos[1])
                *_xmul[1];
        _bval = (_bmax[1]-_ipos[1])
                *_xmul[1];

        _tmin = std::max(_tmin,
        std::min( _aval, _bval )) ;
        _tmax = std::min(_tmax,
        std::max( _aval, _bval )) ;

        if (_tmax<_tmin) return false ;

        if (_tmax<  +0.) return false ;

        return true ;
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_line_3
    {
/*---------------------- line-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 3;

    private :

    real_type                  _ipos [ 3] ;
    real_type                  _jpos [ 3] ;

    real_type                  _xmul [ 3] ;

    real_type                  _rmin [ 3] ;
    real_type                  _rmax [ 3] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_line_3 (
    __const_ptr  (real_type) _isrc ,
    __const_ptr  (real_type) _jsrc
        )
    {
        this->_ipos[0] = _isrc[0];
        this->_ipos[1] = _isrc[1];
        this->_ipos[2] = _isrc[2];

        this->_jpos[0] = _jsrc[0];
        this->_jpos[1] = _jsrc[1];
        this->_jpos[2] = _jsrc[2];

        this->_rmin[0] = std::min(
              _isrc[0] , _jsrc[0]) ;
        this->_rmin[1] = std::min(
              _isrc[1] , _jsrc[1]) ;
        this->_rmin[2] = std::min(
              _isrc[2] , _jsrc[2]) ;

        this->_rmax[0] = std::max(
              _isrc[0] , _jsrc[0]) ;
        this->_rmax[1] = std::max(
              _isrc[1] , _jsrc[1]) ;
        this->_rmax[2] = std::max(
              _isrc[2] , _jsrc[2]) ;

        this->_xmul[0] = (real_type)1. /
            ( _jsrc[0] - _isrc[0]) ;
        this->_xmul[1] = (real_type)1. /
            ( _jsrc[1] - _isrc[1]) ;
        this->_xmul[2] = (real_type)1. /
            ( _jsrc[2] - _isrc[2]) ;
    }

/*---------------------------------- TRUE if intersection */
    __inline_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::rect_rect_3d (
             this->_rmin ,
             this->_rmax ,
            _bmin, _bmax))
    /*------------------------------ bbox can't intersect */
            return false ;

    /*------------------------------ test if line overlap */
        real_type _aval, _bval ;
        _aval = (_bmin[0]-_ipos[0])
                *_xmul[0];
        _bval = (_bmax[0]-_ipos[0])
                *_xmul[0];

        real_type _tmin, _tmax ;
        _tmin =
        std::min( _aval, _bval )  ;
        _tmax =
        std::max( _aval, _bval )  ;

        if (_tmax<_tmin) return false ;

        _aval = (_bmin[1]-_ipos[1])
                *_xmul[1];
        _bval = (_bmax[1]-_ipos[1])
                *_xmul[1];

        _tmin = std::max(_tmin,
        std::min( _aval, _bval )) ;
        _tmax = std::min(_tmax,
        std::max( _aval, _bval )) ;

        if (_tmax<_tmin) return false ;

        _aval = (_bmin[2]-_ipos[2])
                *_xmul[2];
        _bval = (_bmax[2]-_ipos[2])
                *_xmul[2];

        _tmin = std::max(_tmin,
        std::min( _aval, _bval )) ;
        _tmax = std::min(_tmax,
        std::max( _aval, _bval )) ;

        if (_tmax<_tmin) return false ;

        if (_tmax<  +0.) return false ;

        return true ;
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_flat_3
    {
/*---------------------- flat-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 3;

    private :

    real_type                  _ppos [ 3] ;
    real_type                  _vnrm [ 3] ;
    real_type                  _vabs [ 3] ;
    real_type                  _dval ;

    real_type                  _rmin [ 3] ;
    real_type                  _rmax [ 3] ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_flat_3 (
    __const_ptr  (real_type) _psrc ,
    __const_ptr  (real_type) _vsrc ,
    __const_ptr  (real_type) _asrc ,
    __const_ptr  (real_type) _bsrc
        )
    {
        this->_dval =
        geometry::dot_3d(_vsrc, _psrc) ;

        this->_ppos[0] = _psrc[0];
        this->_ppos[1] = _psrc[1];
        this->_ppos[2] = _psrc[2];

        this->_vnrm[0] = _vsrc[0];
        this->_vnrm[1] = _vsrc[1];
        this->_vnrm[2] = _vsrc[2];

        this->_vabs[0] =
               std::abs (_vsrc[0]) ;
        this->_vabs[1] =
               std::abs (_vsrc[1]) ;
        this->_vabs[2] =
               std::abs (_vsrc[2]) ;

        this->_rmin[0] = _asrc[0];
        this->_rmin[1] = _asrc[1];
        this->_rmin[2] = _asrc[2];

        this->_rmax[0] = _bsrc[0];
        this->_rmax[1] = _bsrc[1];
        this->_rmax[2] = _bsrc[2];
    }

/*---------------------------------- TRUE if intersection */
    __normal_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
    /*------------------------------ test if bbox overlap */
        if(!geometry::rect_rect_3d (
             this->_rmin ,
             this->_rmax ,
            _bmin, _bmax))
    /*------------------------------ bbox can't intersect */
            return false ;

    /*------------------------------ check flat intersect */
        real_type _bmid[3] = {
       (real_type) +.5 * _bmin [0] +
       (real_type) +.5 * _bmax [0] ,
       (real_type) +.5 * _bmin [1] +
       (real_type) +.5 * _bmax [1] ,
       (real_type) +.5 * _bmin [2] +
       (real_type) +.5 * _bmax [2] ,
            } ;
        real_type _blen[3] = {
       (real_type) +1. * _bmax [0] -
       (real_type) +1. * _bmid [0] ,
       (real_type) +1. * _bmax [1] -
       (real_type) +1. * _bmid [1] ,
       (real_type) +1. * _bmax [2] -
       (real_type) +1. * _bmid [2] ,
            } ;

        real_type _rval =
        _blen[0] * this->_vabs [0] +
        _blen[1] * this->_vabs [1] +
        _blen[2] * this->_vabs [2] ;

        real_type _sval =
        geometry::dot_3d(_vnrm, _bmid) ;

        _sval -=   this->_dval;

        return std::abs(_sval)<=_rval;
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_ball_2
    {
/*---------------------- ball-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 2;

    private :

    real_type                  _ppos [ 2] ;
    real_type                  _rsqr ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_ball_2 (
    __const_ptr  (real_type) _psrc ,
                  real_type  _rsrc
        )
    {
        this->_ppos[0] = _psrc[0];
        this->_ppos[1] = _psrc[1];

        this->_rsqr    = _rsrc *
                         _rsrc ;
    }

/*---------------------------------- TRUE if intersection */
    __normal_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
        real_type _dsqr = (real_type)+0. ;
        real_type _diff = (real_type)+0. ;

    /*------------------------------ check ball intersect */
        if (this->_ppos[0] < _bmin[0])
        {
            _diff  = _bmin[0] -
               this->_ppos[0] ;

            _dsqr += _diff * _diff ;
        }
        else
        if (this->_ppos[0] > _bmax[0])
        {
            _diff  = _bmax[0] -
               this->_ppos[0] ;

            _dsqr += _diff * _diff ;
        }

        if (this->_ppos[1] < _bmin[1])
        {
            _diff  = _bmin[1] -
               this->_ppos[1] ;

            _dsqr += _diff * _diff ;
        }
        else
        if (this->_ppos[1] > _bmax[1])
        {
            _diff  = _bmax[1] -
               this->_ppos[1] ;

            _dsqr += _diff * _diff ;
        }

        return _dsqr<= this->_rsqr ;
    }

    } ;

    template <
    typename R ,
    typename I
             >
    class aabb_pred_ball_3
    {
/*---------------------- ball-tree intersection predicate */
    public  :

    typedef R                   real_type ;
    typedef I                   iptr_type ;

    iptr_type static constexpr _dims = 3;

    private :

    real_type                  _ppos [ 3] ;
    real_type                  _rsqr ;

    public  :
/*---------------------------------- construct from _src. */
    __normal_call aabb_pred_ball_3 (
    __const_ptr  (real_type) _psrc ,
                  real_type  _rsrc
        )
    {
        this->_ppos[0] = _psrc[0];
        this->_ppos[1] = _psrc[1];
        this->_ppos[2] = _psrc[2];

        this->_rsqr    = _rsrc *
                         _rsrc ;
    }

/*---------------------------------- TRUE if intersection */
    __normal_call bool_type operator() (
    __const_ptr  (real_type) _bmin ,
    __const_ptr  (real_type) _bmax
        ) const
    {
        real_type _dsqr = (real_type)+0. ;
        real_type _diff = (real_type)+0. ;

    /*------------------------------ check ball intersect */
        if (this->_ppos[0] < _bmin[0])
        {
            _diff  = _bmin[0] -
               this->_ppos[0] ;

            _dsqr += _diff * _diff ;
        }
        else
        if (this->_ppos[0] > _bmax[0])
        {
            _diff  = _bmax[0] -
               this->_ppos[0] ;

            _dsqr += _diff * _diff ;
        }

        if (this->_ppos[1] < _bmin[1])
        {
            _diff  = _bmin[1] -
               this->_ppos[1] ;

            _dsqr += _diff * _diff ;
        }
        else
        if (this->_ppos[1] > _bmax[1])
        {
            _diff  = _bmax[1] -
               this->_ppos[1] ;

            _dsqr += _diff * _diff ;
        }

        if (this->_ppos[2] < _bmin[2])
        {
            _diff  = _bmin[2] -
               this->_ppos[2] ;

            _dsqr += _diff * _diff ;
        }
        else
        if (this->_ppos[2] > _bmax[2])
        {
            _diff  = _bmax[2] -
               this->_ppos[2] ;

            _dsqr += _diff * _diff ;
        }

        return _dsqr<= this->_rsqr ;
    }

    } ;


    }

#   endif   //__AABB_PRED_K__



