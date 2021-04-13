
    /*
    --------------------------------------------------------
     * INTERSECT-K: various (robust) intersection tests.
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
     * Last updated: 12 August, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __INTERSECT_K__
#   define __INTERSECT_K__

    namespace geometry
    {

    typedef char_type hits_type ;

    hits_type null_hits = +0 ;
    hits_type node_hits = +1 ;
    hits_type edge_hits = +2 ;
    hits_type face_hits = +3 ;
    hits_type tria_hits = +4 ;

    /*
    --------------------------------------------------------
     * BALL-LINE-KD: ball-line intersections
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call size_t ball_line_2d (
    __const_ptr  (real_type) _pc, // ball BB(pc,rc)
        real_type            _rc,
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qa, // intersections
    __write_ptr  (real_type) _qb
        )
    {
        real_type _vv[2] ;
        _vv[0] = _pb[0] -_pa[0] ;
        _vv[1] = _pb[1] -_pa[1] ;

        real_type _ca[2] ;
        _ca[0] = _pa[0] -_pc[0] ;
        _ca[1] = _pa[1] -_pc[1] ;

        real_type _aa =
            geometry::dot_2d(_vv, _vv) ;
        real_type _bb =
            geometry::dot_2d(_vv, _ca) *
                       (real_type) +2. ;
        real_type _cc =
            geometry::dot_2d(_ca, _ca) ;
        _cc -= _rc * _rc ;

        size_t _nn = +0  ;
        real_type _tt[2] ;
        if ( math::
        polyroots(_aa, _bb, _cc, _tt))
        {
        if (_tt[0] >= (real_type)+0. &&
            _tt[0] <= (real_type)+1. )
        {
    /*----------------------- compute 1st-root intersect. */
            real_type _WB =
           (real_type)+0.+_tt[0] ;

            real_type _WA =
           (real_type)+1.-_tt[0] ;

            _nn += +1 ;

            dd_flt _PA[2] ;
            _PA[0]=_pa[0] ;
            _PA[1]=_pa[1] ;

            dd_flt _PB[2] ;
            _PB[0]=_pb[0] ;
            _PB[1]=_pb[1] ;

            dd_flt _QQ[2] ;
            _QQ[0]=_PA[0] * _WA  +
                   _PB[0] * _WB  ;
            _QQ[1]=_PA[1] * _WA  +
                   _PB[1] * _WB  ;

            if (_nn == 1)
            {
            _qa[0]=_QQ[0] ;
            _qa[1]=_QQ[1] ;
            }
            else
            {
            _qb[0]=_QQ[0] ;
            _qb[1]=_QQ[1] ;
            }
        }
        if (_tt[1] >= (real_type)+0. &&
            _tt[1] <= (real_type)+1. )
        {
    /*----------------------- compute 2nd-root intersect. */
            real_type _WB =
           (real_type)+0.+_tt[1] ;

            real_type _WA =
           (real_type)+1.-_tt[1] ;

            _nn += +1 ;

            dd_flt _PA[2] ;
            _PA[0]=_pa[0] ;
            _PA[1]=_pa[1] ;

            dd_flt _PB[2] ;
            _PB[0]=_pb[0] ;
            _PB[1]=_pb[1] ;

            dd_flt _QQ[2] ;
            _QQ[0]=_PA[0] * _WA  +
                   _PB[0] * _WB  ;
            _QQ[1]=_PA[1] * _WA  +
                   _PB[1] * _WB  ;

            if (_nn == 1)
            {
            _qa[0]=_QQ[0] ;
            _qa[1]=_QQ[1] ;
            }
            else
            {
            _qb[0]=_QQ[0] ;
            _qb[1]=_QQ[1] ;
            }
        }
        }

        return ( _nn ) ; // return num roots
    }

    template <
    typename      real_type
             >
    __normal_call size_t ball_line_3d (
    __const_ptr  (real_type) _pc, // ball BB(pc,rc)
        real_type            _rc,
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qa, // intersections
    __write_ptr  (real_type) _qb
        )
    {
        real_type _vv[3] ;
        _vv[0] = _pb[0] -_pa[0] ;
        _vv[1] = _pb[1] -_pa[1] ;
        _vv[2] = _pb[2] -_pa[2] ;

        real_type _ca[3] ;
        _ca[0] = _pa[0] -_pc[0] ;
        _ca[1] = _pa[1] -_pc[1] ;
        _ca[2] = _pa[2] -_pc[2] ;

        real_type _aa =
            geometry::dot_3d(_vv, _vv) ;
        real_type _bb =
            geometry::dot_3d(_vv, _ca) *
                       (real_type) +2. ;
        real_type _cc =
            geometry::dot_3d(_ca, _ca) ;
        _cc -= _rc * _rc ;

        size_t _nn = +0  ;
        real_type _tt[2] ;
        if ( math::
        polyroots(_aa, _bb, _cc, _tt))
        {
        if (_tt[0] >= (real_type)+0. &&
            _tt[0] <= (real_type)+1. )
        {
    /*----------------------- compute 1st-root intersect. */
            real_type _WB =
           (real_type)+0.+_tt[0] ;

            real_type _WA =
           (real_type)+1.-_tt[0] ;

            _nn += +1 ;

            dd_flt _PA[3] ;
            _PA[0]=_pa[0] ;
            _PA[1]=_pa[1] ;
            _PA[2]=_pa[2] ;

            dd_flt _PB[3] ;
            _PB[0]=_pb[0] ;
            _PB[1]=_pb[1] ;
            _PB[2]=_pb[2] ;

            dd_flt _QQ[3] ;
            _QQ[0]=_PA[0] * _WA  +
                   _PB[0] * _WB  ;
            _QQ[1]=_PA[1] * _WA  +
                   _PB[1] * _WB  ;
            _QQ[2]=_PA[2] * _WA  +
                   _PB[2] * _WB  ;

            if (_nn == 1)
            {
            _qa[0]=_QQ[0] ;
            _qa[1]=_QQ[1] ;
            _qa[2]=_QQ[2] ;
            }
            else
            {
            _qb[0]=_QQ[0] ;
            _qb[1]=_QQ[1] ;
            _qb[2]=_QQ[2] ;
            }
        }
        if (_tt[1] >= (real_type)+0. &&
            _tt[1] <= (real_type)+1. )
        {
    /*----------------------- compute 2nd-root intersect. */
            real_type _WB =
           (real_type)+0.+_tt[1] ;

            real_type _WA =
           (real_type)+1.-_tt[1] ;

            _nn += +1 ;

            dd_flt _PA[3] ;
            _PA[0]=_pa[0] ;
            _PA[1]=_pa[1] ;
            _PA[2]=_pa[2] ;

            dd_flt _PB[3] ;
            _PB[0]=_pb[0] ;
            _PB[1]=_pb[1] ;
            _PB[2]=_pb[2] ;

            dd_flt _QQ[3] ;
            _QQ[0]=_PA[0] * _WA  +
                   _PB[0] * _WB  ;
            _QQ[1]=_PA[1] * _WA  +
                   _PB[1] * _WB  ;
            _QQ[2]=_PA[2] * _WA  +
                   _PB[2] * _WB  ;

            if (_nn == 1)
            {
            _qa[0]=_QQ[0] ;
            _qa[1]=_QQ[1] ;
            _qa[2]=_QQ[2] ;
            }
            else
            {
            _qb[0]=_QQ[0] ;
            _qb[1]=_QQ[1] ;
            _qb[2]=_QQ[2] ;
            }
        }
        }

        return ( _nn ) ; // return num roots
    }

    /*
    --------------------------------------------------------
     * RECT-RECT-KD: rect-rect intersections
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __inline_call bool node_rect_2d (
    __const_ptr  (real_type) _pp, // node
    __const_ptr  (real_type) _b0, // min. rect. B
    __const_ptr  (real_type) _b1  // max.
        )
    {
        if (_pp[0] >= _b0[0] &&
            _pp[0] <= _b1[0] &&
            _pp[1] >= _b0[1] &&
            _pp[1] <= _b1[1] )
        {
            return (  true ) ;
        }
        else
        {
            return ( false ) ;
        }
    }

    template <
    typename      real_type
             >
    __inline_call bool node_rect_3d (
    __const_ptr  (real_type) _pp, // node
    __const_ptr  (real_type) _b0, // min. rect. B
    __const_ptr  (real_type) _b1  // max.
        )
    {
        if (_pp[0] >= _b0[0] &&
            _pp[0] <= _b1[0] &&
            _pp[1] >= _b0[1] &&
            _pp[1] <= _b1[1] &&
            _pp[2] >= _b0[2] &&
            _pp[2] <= _b1[2] )
        {
            return (  true ) ;
        }
        else
        {
            return ( false ) ;
        }
    }

    template <
    typename      real_type
             >
    __inline_call bool rect_rect_2d (
    __const_ptr  (real_type) _a0, // min. rect. A
    __const_ptr  (real_type) _a1, // max.
    __const_ptr  (real_type) _b0, // min. rect. B
    __const_ptr  (real_type) _b1  // max.
        )
    {
        if (_a0[0] <= _b1[0] &&
            _b0[0] <= _a1[0] &&
            _a0[1] <= _b1[1] &&
            _b0[1] <= _a1[1] )
        {
            return (  true ) ;
        }
        else
        {
            return ( false ) ;
        }
    }

    template <
    typename      real_type
             >
    __inline_call bool rect_rect_3d (
    __const_ptr  (real_type) _a0, // min. rect. A
    __const_ptr  (real_type) _a1, // max.
    __const_ptr  (real_type) _b0, // min. rect. B
    __const_ptr  (real_type) _b1  // max.
        )
    {
        if (_a0[0] <= _b1[0] &&
            _b0[0] <= _a1[0] &&
            _a0[1] <= _b1[1] &&
            _b0[1] <= _a1[1] &&
            _a0[2] <= _b1[2] &&
            _b0[2] <= _a1[2] )
        {
            return (  true ) ;
        }
        else
        {
            return ( false ) ;
        }
    }

    /*
    --------------------------------------------------------
     * COPY-NODE-KD: node-node copy helpers.
    --------------------------------------------------------
     */

    template <
    typename      data_one_,
    typename      data_two_
             >
    __inline_call void copy_node_2d (
    __write_ptr  (data_one_) _pa, // copy: A := B
    __const_ptr  (data_two_) _pb
        )
    {
        _pa[0] = (data_one_) _pb[0] ;
        _pa[1] = (data_one_) _pb[1] ;
    }

    template <
    typename      data_one_,
    typename      data_two_
             >
    __inline_call void copy_node_3d (
    __write_ptr  (data_one_) _pa, // copy: A := B
    __const_ptr  (data_two_) _pb
        )
    {
        _pa[0] = (data_one_) _pb[0] ;
        _pa[1] = (data_one_) _pb[1] ;
        _pa[2] = (data_one_) _pb[2] ;
    }

    template <
    typename      data_one_,
    typename      data_two_
             >
    __inline_call void copy_node_xy (
    __write_ptr  (data_one_) _pa, // copy: A := B
    __const_ptr  (data_two_) _pb
        )
    {
        _pa[0] = (data_one_) _pb[0] ;
        _pa[1] = (data_one_) _pb[1] ;
    }

    template <
    typename      data_one_,
    typename      data_two_
             >
    __inline_call void copy_node_yz (
    __write_ptr  (data_one_) _pa, // copy: A := B
    __const_ptr  (data_two_) _pb
        )
    {
        _pa[0] = (data_one_) _pb[1] ;
        _pa[1] = (data_one_) _pb[2] ;
    }

    template <
    typename      data_one_,
    typename      data_two_
             >
    __inline_call void copy_node_xz (
    __write_ptr  (data_one_) _pa, // copy: A := B
    __const_ptr  (data_two_) _pb
        )
    {
        _pa[0] = (data_one_) _pb[0] ;
        _pa[1] = (data_one_) _pb[2] ;
    }

    /*
    --------------------------------------------------------
     * LINE-FLAT-KD: line-flat intersections
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call bool line_flat_3d (
    __const_ptr  (real_type) _pp, // (xx-pp).nv=0
    __const_ptr  (real_type) _nv,
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qq, // intersection
        bool _bind =   true
        )
    {
        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] &&
            _pa[2] == _pb[2] )
        {
          //std::cout << "node_flat_3d" << std::endl;

            return ( false ) ;
        }
        else
        if (_pp[0] == _pa[0] &&
            _pp[1] == _pa[1] &&
            _pp[2] == _pa[2] )
        {

            _qq[0] =  _pa[0] ;
            _qq[1] =  _pa[1] ;
            _qq[2] =  _pa[2] ;

            return (  true ) ;
        }
        else
        if (_pp[0] == _pb[0] &&
            _pp[1] == _pb[1] &&
            _pp[2] == _pb[2] )
        {

            _qq[0] =  _pb[0] ;
            _qq[1] =  _pb[1] ;
            _qq[2] =  _pb[2] ;

            return (  true ) ;
        }
        else
        {

        real_type _ab[3];
        _ab[0] = _pb[0] - _pa[0] ;
        _ab[1] = _pb[1] - _pa[1] ;
        _ab[2] = _pb[2] - _pa[2] ;

        real_type _ap[3];
        _ap[0] = _pp[0] - _pa[0] ;
        _ap[1] = _pp[1] - _pa[1] ;
        _ap[2] = _pp[2] - _pa[2] ;

        real_type _ep = + 1.0 *
        std::numeric_limits<real_type>::epsilon() ;

        real_type _d1 =
        geometry::dot_3d(_ap, _nv) ;
        real_type _d2 =
        geometry::dot_3d(_ab, _nv) ;

        if (std::abs(_d2) <= _ep * std::abs(_d1))
            return ( false ) ;

        real_type _tt =  _d1 / _d2 ;

        if (_bind)
        {
        if (_tt  < (real_type)+0.)
            return ( false ) ;
        if (_tt  > (real_type)+1.)
            return ( false ) ;
        }

        if (_tt == (real_type)+0.)
        {
            _qq[0] = _pa[0] ;
            _qq[1] = _pa[1] ;
            _qq[2] = _pa[2] ;
        }
        else
        if (_tt == (real_type)+1.)
        {
            _qq[0] = _pb[0] ;
            _qq[1] = _pb[1] ;
            _qq[2] = _pb[2] ;
        }
        else
        {
            dd_flt   _AB[3] ;
            _AB[0] = _pb[0] ;
            _AB[1] = _pb[1] ;
            _AB[2] = _pb[2] ;
            _AB[0]-= _pa[0] ;
            _AB[1]-= _pa[1] ;
            _AB[2]-= _pa[2] ;

            dd_flt   _QQ[3] ;
            _QQ[0] = _pa[0] ;
            _QQ[1] = _pa[1] ;
            _QQ[2] = _pa[2] ;
            _QQ[0]+= _AB[0] * _tt ;
            _QQ[1]+= _AB[1] * _tt ;
            _QQ[2]+= _AB[2] * _tt ;

            _qq[0] = _QQ[0] ;
            _qq[1] = _QQ[1] ;
            _qq[2] = _QQ[2] ;
        }

        return ( true ) ;

        }
    }

    /*
    --------------------------------------------------------
     * TRIA-FLAT-KD: tria-flat intersections
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call size_t tria_flat_3d (
    __const_ptr  (real_type) _pp, // (xx-pp).nv=0
    __const_ptr  (real_type) _nv,
    __const_ptr  (real_type) _pa, // tria
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc,
    __write_ptr  (real_type) _qa, // intersection
    __write_ptr  (real_type) _qb
        )
    {
        size_t _ni = +0;

        if (line_flat_3d (
            _pp, _nv,
            _pa, _pb,
           (_ni == +0) ? _qa : _qb))
            _ni += +1  ;
        if (line_flat_3d (
            _pp, _nv,
            _pb, _pc,
           (_ni == +0) ? _qa : _qb))
            _ni += +1  ;
        if (line_flat_3d (
            _pp, _nv,
            _pc, _pa,
           (_ni == +0) ? _qa : _qb))
            _ni += +1  ;

        return ( _ni ) ;
    }

    /*
    --------------------------------------------------------
     * PROJ-LINE-KD: node-line projections
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call bool proj_line_2d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // node on line
    __const_ptr  (real_type) _va, // vec. of line
                  real_type &_tt
        )
    {
        real_type _ap[2] ;
        _ap[0] = _pp[0] - _pa[0] ;
        _ap[1] = _pp[1] - _pa[1] ;

        real_type _d1 =
        geometry::dot_2d(_ap, _va) ;
        real_type _d2 =
        geometry::dot_2d(_va, _va) ;

        _tt =  _d1 / _d2 ;

        return (  true ) ;
    }

    template <
    typename      real_type
             >
    __normal_call bool proj_line_2d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // line nodes
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qq, // projection
                  hits_type &_ht
        )
    {
        real_type _vv[2], _tt;
        _vv[0] = _pb[0] - _pa[0] ;
        _vv[1] = _pb[1] - _pa[1] ;

        if (proj_line_2d(_pp, _pa,
                         _vv, _tt) )
        {
        if (_tt <= (real_type)+0.)
        {
            _ht  =  node_hits;

            _qq[0] = _pa[0] ;
            _qq[1] = _pa[1] ;

            return (  true ) ;
        }
        else
        if (_tt >= (real_type)+1.)
        {
            _ht  =  node_hits;

            _qq[0] = _pb[0] ;
            _qq[1] = _pb[1] ;

            return (  true ) ;
        }
        else
        {
            _ht  =  edge_hits;

            _qq[0] =
            _pa[0] + _tt*_vv[ 0] ;
            _qq[1] =
            _pa[1] + _tt*_vv[ 1] ;

            return (  true ) ;
        }
        }
        else
        {
            return ( false ) ;
        }
    }

    template <
    typename      real_type
             >
    __normal_call bool proj_line_3d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // node on line
    __const_ptr  (real_type) _va, // vec. of line
                  real_type &_tt
        )
    {
        real_type _ap[3] ;
        _ap[0] = _pp[0] - _pa[0] ;
        _ap[1] = _pp[1] - _pa[1] ;
        _ap[2] = _pp[2] - _pa[2] ;

        real_type _d1 =
        geometry::dot_3d(_ap, _va) ;
        real_type _d2 =
        geometry::dot_3d(_va, _va) ;

        _tt =  _d1 / _d2 ;

        return (  true ) ;
    }

    template <
    typename      real_type
             >
    __normal_call bool proj_line_3d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // line nodes
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qq, // projection
                  hits_type &_ht
        )
    {
        real_type _vv[3], _tt;
        _vv[0] = _pb[0] - _pa[0] ;
        _vv[1] = _pb[1] - _pa[1] ;
        _vv[2] = _pb[2] - _pa[2] ;

        if (proj_line_3d(_pp, _pa,
                         _vv, _tt) )
        {
        if (_tt <= (real_type)+0.)
        {
            _ht  =  node_hits;

            _qq[0] = _pa[0] ;
            _qq[1] = _pa[1] ;
            _qq[2] = _pa[2] ;

            return (  true ) ;
        }
        else
        if (_tt >= (real_type)+1.)
        {
            _ht  =  node_hits;

            _qq[0] = _pb[0] ;
            _qq[1] = _pb[1] ;
            _qq[2] = _pb[2] ;

            return (  true ) ;
        }
        else
        {
            _ht  =  edge_hits;

            _qq[0] =
            _pa[0] + _tt*_vv[ 0] ;
            _qq[1] =
            _pa[1] + _tt*_vv[ 1] ;
            _qq[2] =
            _pa[2] + _tt*_vv[ 2] ;

            return (  true ) ;
        }
        }
        else
        {
            return ( false ) ;
        }
    }

    /*
    --------------------------------------------------------
     * PROJ-TRIA-KD: node-tria projections
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __normal_call bool proj_flat_3d (
    __const_ptr  (real_type) _pi,
    __const_ptr  (real_type) _pp, // node on flat
    __const_ptr  (real_type) _nv, // nrm. of flat
    __write_ptr  (real_type) _qq
        )
    {
        real_type _ip[3];
        _ip[0] = _pp[0] - _pi[0] ;
        _ip[1] = _pp[1] - _pi[1] ;
        _ip[2] = _pp[2] - _pi[2] ;

        real_type _tt ;
        real_type _d1 =
        geometry::dot_3d(_ip, _nv) ;
        real_type _d2 =
        geometry::dot_3d(_nv, _nv) ;

        _tt    = _d1/_d2 ;

        _qq[0] =
        _pi[0] + _tt * _nv[0] ;
        _qq[1] =
        _pi[1] + _tt * _nv[1] ;
        _qq[2] =
        _pi[2] + _tt * _nv[2] ;

        return (  true ) ;
    }

    template <
    typename      real_type
             >
    __normal_call bool proj_tria_2d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // tria nodes
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc,
    __write_ptr  (real_type) _qq, // projection
                  hits_type &_ht
        )
    {

    /*----------------------- projector for low-dim. face */
    #   define __projface2d(_pi, _pj)       \
        if (proj_line_2d (              \
            _pp , _pi, _pj, _qt, _ht) ) \
        {                               \
            real_type _dd =             \
        geometry::lensqr_2d(_pp, _qt) ; \
                                        \
            if (_dd < _dm)              \
            {                           \
                _qq[0] =  _qt[0] ;      \
                _qq[1] =  _qt[1] ;      \
                                        \
                _dm = _dd ;             \
            }                           \
        }

        {
    /*----------------------- test node-tria intersection */
            double _PA[2] ;
            copy_node_2d(_PA, _pa) ;

            double _PB[2] ;
            copy_node_2d(_PB, _pb) ;

            double _PC[2] ;
            copy_node_2d(_PC, _pc) ;

            double _PP[2] ;
            copy_node_2d(_PP, _pp) ;

            double _s1 =
                geompred::orient2d (
                (double*)_PA ,
                (double*)_PB ,
                (double*)_PP )   ;

            double _s2 =
                geompred::orient2d (
                (double*)_PB ,
                (double*)_PC ,
                (double*)_PP )   ;

            double _s3 =
                geompred::orient2d (
                (double*)_PC ,
                (double*)_PA ,
                (double*)_PP )   ;

            if (_s1 * _s2 < +0.0 ||
                _s2 * _s3 < +0.0 ||
                _s3 * _s1 < +0.0 )
            {

            real_type _dm = +std::
                numeric_limits
           <real_type>::infinity() ;

            real_type _qt[2] ;

            _ht =  null_hits ;

            __projface2d(_pa, _pb) ;
            __projface2d(_pb, _pc) ;
            __projface2d(_pc, _pa) ;

            return ( true )  ;
            }
            else
            {

            _ht =  tria_hits ;

            copy_node_2d(_qq, _pp) ;

            return ( true )  ;
            }
        }

    #   undef  __projface2d
    }

    template <
    typename      real_type
             >
    __normal_call void tria_norm_3d (
    __const_ptr  (real_type) _pa, // tria nodes
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc,
    __write_ptr  (real_type) _nv
        ) ;

    template <
    typename      real_type
             >
    __normal_call bool proj_tria_3d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // tria nodes
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc,
    __write_ptr  (real_type) _qq, // projection
                  hits_type &_ht
        )
    {

    /*----------------------- projector for low-dim. face */
    #   define __projface3d(_pi, _pj)       \
        if (proj_line_3d (              \
            _pp , _pi, _pj, _qt, _ht) ) \
        {                               \
            real_type _dd =             \
        geometry::lensqr_3d(_pp, _qt) ; \
                                        \
            if (_dd < _dm)              \
            {                           \
                _qq[0] =  _qt[0] ;      \
                _qq[1] =  _qt[1] ;      \
                _qq[2] =  _qt[2] ;      \
                                        \
                _dm = _dd ;             \
            }                           \
        }

    /*----------------------- test node-tria intersection */
        double _s1[3], _s2[3], _s3[3];

        real_type _nv[3];
        tria_norm_3d(_pa, _pb, _pc, _nv);

        real_type _pt[3];
        proj_flat_3d(_pp, _pa, _nv, _pt);

        {
    /*----------------------- get orientation in xy-plane */
            double _TA[2] ;
            copy_node_xy(_TA, _pa) ;

            double _TB[2] ;
            copy_node_xy(_TB, _pb) ;

            double _TC[2] ;
            copy_node_xy(_TC, _pc) ;

            double _TP[2] ;
            copy_node_xy(_TP, _pt) ;

            _s1[0] =
                geompred::orient2d (
                (double*)_TA ,
                (double*)_TB ,
                (double*)_TP )   ;

            _s2[0] =
                geompred::orient2d (
                (double*)_TB ,
                (double*)_TC ,
                (double*)_TP )   ;

            _s3[0] =
                geompred::orient2d (
                (double*)_TC ,
                (double*)_TA ,
                (double*)_TP )   ;
        }

        {
    /*----------------------- get orientation in xz-plane */
            double _TA[2] ;
            copy_node_xz(_TA, _pa) ;

            double _TB[2] ;
            copy_node_xz(_TB, _pb) ;

            double _TC[2] ;
            copy_node_xz(_TC, _pc) ;

            double _TP[2] ;
            copy_node_xz(_TP, _pt) ;

            _s1[1] =
                geompred::orient2d (
                (double*)_TA ,
                (double*)_TB ,
                (double*)_TP )   ;

            _s2[1] =
                geompred::orient2d (
                (double*)_TB ,
                (double*)_TC ,
                (double*)_TP )   ;

            _s3[1] =
                geompred::orient2d (
                (double*)_TC ,
                (double*)_TA ,
                (double*)_TP )   ;
        }

        {
    /*----------------------- get orientation in yz-plane */
            double _TA[2] ;
            copy_node_yz(_TA, _pa) ;

            double _TB[2] ;
            copy_node_yz(_TB, _pb) ;

            double _TC[2] ;
            copy_node_yz(_TC, _pc) ;

            double _TP[2] ;
            copy_node_yz(_TP, _pt) ;

            _s1[2] =
                geompred::orient2d (
                (double*)_TA ,
                (double*)_TB ,
                (double*)_TP )   ;

            _s2[2] =
                geompred::orient2d (
                (double*)_TB ,
                (double*)_TC ,
                (double*)_TP )   ;

            _s3[2] =
                geompred::orient2d (
                (double*)_TC ,
                (double*)_TA ,
                (double*)_TP )   ;
        }

        {
    /*----------------------- test node-tria intersection */
            if (
            _s1[0] * _s2[0] < +0.0 ||
            _s2[0] * _s3[0] < +0.0 ||
            _s3[0] * _s1[0] < +0.0 ||
            _s1[1] * _s2[1] < +0.0 ||
            _s2[1] * _s3[1] < +0.0 ||
            _s3[1] * _s1[1] < +0.0 ||
            _s1[2] * _s2[2] < +0.0 ||
            _s2[2] * _s3[2] < +0.0 ||
            _s3[2] * _s1[2] < +0.0 )
            {

            real_type _dm = +std::
                numeric_limits
           <real_type>::infinity() ;

            real_type _qt[3] ;

            _ht =  null_hits ;

            __projface3d(_pa, _pb) ;
            __projface3d(_pb, _pc) ;
            __projface3d(_pc, _pa) ;

            return ( true )  ;
            }
            else
            {

            _ht =  face_hits ;

            copy_node_3d(_qq, _pt) ;

            return ( true )  ;
            }
        }

    #   undef  __projface3d
    }

    template <
    typename      real_type
             >
    __normal_call bool proj_tria_3d (
    __const_ptr  (real_type) _pp,
    __const_ptr  (real_type) _pa, // tria nodes
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc,
    __const_ptr  (real_type) _pd,
    __write_ptr  (real_type) _qq, // projection
                  hits_type &_ht
        )
    {

    /*----------------------- projector for low-dim. face */
    #   define __projface3d(_pi, _pj, _pk)  \
        if (proj_tria_3d (              \
            _pp , _pi, _pj, _pk, _qt, _ht)) \
        {                               \
            real_type _dd =             \
        geometry::lensqr_3d(_pp, _qt) ; \
                                        \
            if (_dd < _dm)              \
            {                           \
                _qq[0] =  _qt[0] ;      \
                _qq[1] =  _qt[1] ;      \
                _qq[2] =  _qt[2] ;      \
                                        \
                _dm = _dd ;             \
            }                           \
        }

        {
    /*----------------------- test node-tria intersection */
            double _PA[3] ;
            copy_node_3d(_PA, _pa) ;

            double _PB[3] ;
            copy_node_3d(_PB, _pb) ;

            double _PC[3] ;
            copy_node_3d(_PC, _pc) ;

            double _PD[3] ;
            copy_node_3d(_PD, _pd) ;

            double _PP[3] ;
            copy_node_3d(_PP, _pp) ;

            double _s1 =
                geompred::orient3d (
                (double*)_PB ,
                (double*)_PD ,
                (double*)_PC ,
                (double*)_PP )   ;

            double _s2 =
                geompred::orient3d (
                (double*)_PC ,
                (double*)_PD ,
                (double*)_PA ,
                (double*)_PP )   ;

            double _s3 =
                geompred::orient3d (
                (double*)_PA ,
                (double*)_PD ,
                (double*)_PB ,
                (double*)_PP )   ;

            double _s4 =
                geompred::orient3d (
                (double*)_PA ,
                (double*)_PB ,
                (double*)_PC ,
                (double*)_PP )   ;

            if (_s1 * _s2 < +0.0 ||
                _s2 * _s3 < +0.0 ||
                _s3 * _s1 < +0.0 ||
                _s1 * _s4 < +0.0 ||
                _s2 * _s4 < +0.0 ||
                _s3 * _s4 < +0.0 )
            {

            real_type _dm = +std::
                numeric_limits
           <real_type>::infinity() ;

            real_type _qt[3] ;

            _ht =  null_hits ;

            __projface3d(
                    _pb, _pd, _pc) ;
            __projface3d(
                    _pc, _pd, _pa) ;
            __projface3d(
                    _pa, _pd, _pb) ;
            __projface3d(
                    _pa, _pb, _pc) ;

            return ( true )  ;
            }
            else
            {

            _ht =  tria_hits ;

            copy_node_3d(_qq, _pp) ;

            return ( true )  ;
            }
        }

    #   undef  __projface3d
    }

    /*
    --------------------------------------------------------
     * linear intersection kernels
    --------------------------------------------------------
     */

    __normal_call double cleave2d (
    __const_ptr  (double) _pa,
    __const_ptr  (double) _pb,
    __const_ptr  (double) _pp
        )
    {
        double _ab[2] ;
        _ab[0] = _pb[0]-_pa[0] ;
        _ab[1] = _pb[1]-_pa[1] ;

        double _ap[2] ;
        _ap[0] = _pp[0]-_pa[0] ;
        _ap[1] = _pp[1]-_pa[1] ;

        double _dp =
            _ab[0] * _ap[0] +
            _ab[1] * _ap[1] ;

        return _dp ;
    }

    __normal_call double cleave3d (
    __const_ptr  (double) _pa,
    __const_ptr  (double) _pb,
    __const_ptr  (double) _pp
        )
    {
        double _ab[3] ;
        _ab[0] = _pb[0]-_pa[0] ;
        _ab[1] = _pb[1]-_pa[1] ;
        _ab[2] = _pb[2]-_pa[2] ;

        double _ap[3] ;
        _ap[0] = _pp[0]-_pa[0] ;
        _ap[1] = _pp[1]-_pa[1] ;
        _ap[2] = _pp[2]-_pa[2] ;

        double _dp =
            _ab[0] * _ap[0] +
            _ab[1] * _ap[1] +
            _ab[2] * _ap[2] ;

        return _dp ;
    }


    /*
    --------------------------------------------------------
     * NODE-NODE-KD: node-node intersections
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __inline_call hits_type node_node_2d (
    __const_ptr  (real_type) _pa, // node
    __const_ptr  (real_type) _pb, // node
    __write_ptr  (real_type) _qq  // intersection
        )
    {
        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] )
        {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pa[0] ;
            _qq[1] =  _pa[1] ;

            return node_hits ;
        }
        else
        {
    /*----------------------- no intersections: null hits */
            return null_hits ;
        }
    }

    template <
        typename  real_type
             >
    __inline_call hits_type node_node_3d (
    __const_ptr  (real_type) _pa, // node
    __const_ptr  (real_type) _pb, // node
    __write_ptr  (real_type) _qq  // intersection
        )
    {
        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] &&
            _pa[2] == _pb[2] )
        {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pa[0] ;
            _qq[1] =  _pa[1] ;
            _qq[2] =  _pa[2] ;

            return node_hits ;
        }
        else
        {
    /*----------------------- no intersections: null hits */
            return null_hits ;
        }
    }

    /*
    --------------------------------------------------------
     * NODE-LINE-KD: node-line intersections
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __normal_call hits_type node_line_2d (
    __const_ptr  (real_type) _pp, // node
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true
        )
    {
        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] )
        {
    /*----------------------- test node-node intersection */
            return
            node_node_2d(_pa, _pp, _qq);
        }
        else
        {
    /*----------------------- test node-line intersection */
            double _PA[2] ;
            _PA[0] =  _pa[0] ;
            _PA[1] =  _pa[1] ;

            double _PB[2] ;
            _PB[0] =  _pb[0] ;
            _PB[1] =  _pb[1] ;

            double _PP[2] ;
            _PP[0] =  _pp[0] ;
            _PP[1] =  _pp[1] ;

            double _ss =
                geompred::orient2d (
                (double*)_PA ,
                (double*)_PB ,
                (double*)_PP )   ;

            if (_ss == +0.0)
            {
    /*----------------------- orient w.r.t. line endpoint */
            double _sa = cleave2d (
                (double*)_PA ,
                (double*)_PB ,
                (double*)_PP )   ;

            double _sb = cleave2d (
                (double*)_PB ,
                (double*)_PA ,
                (double*)_PP )   ;

            if (_sa*_sb > +0.0)
            {
    /*----------------------- have node-line intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return edge_hits ;
            }
            else
            if (_sa == +0.0)
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pa[0] ;
            _qq[1] =  _pa[1] ;

            return node_hits ;
            }
            else
            if (_sb == +0.0)
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pb[0] ;
            _qq[1] =  _pb[1] ;

            return node_hits ;
            }
            else
            if ( !_bind )
            {
    /*----------------------- is unconstrained: edge hits */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return edge_hits ;
            }
            else
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }

            }
            else
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }
        }
    }

    template <
        typename  real_type
             >
    __normal_call hits_type node_line_3d (
    __const_ptr  (real_type) _pp, // node
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true
        )
    {
        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] &&
            _pa[2] == _pb[2] )
        {
    /*----------------------- test node-node intersection */
            return
            node_node_3d(_pa, _pp, _qq);
        }
        else
        {
    /*----------------------- test node-line intersection */
            double _PA[3] ;
            _PA[0] =  _pa[0] ;
            _PA[1] =  _pa[1] ;
            _PA[2] =  _pa[2] ;

            double _PB[3] ;
            _PB[0] =  _pb[0] ;
            _PB[1] =  _pb[1] ;
            _PB[2] =  _pb[2] ;

            double _PP[3] ;
            _PP[0] =  _pp[0] ;
            _PP[1] =  _pp[1] ;
            _PP[2] =  _pp[2] ;

    /*----------------------- get orientation in xy-plane */
            double _A1[2] ;
            _A1[0] =  _PA[0] ;
            _A1[1] =  _PA[1] ;

            double _B1[2] ;
            _B1[0] =  _PB[0] ;
            _B1[1] =  _PB[1] ;

            double _P1[2] ;
            _P1[0] =  _PP[0] ;
            _P1[1] =  _PP[1] ;

            double _s1 =
                geompred::orient2d (
                (double*)_A1 ,
                (double*)_B1 ,
                (double*)_P1 )   ;

    /*----------------------- get orientation in xz-plane */
            double _A2[2] ;
            _A2[0] =  _PA[0] ;
            _A2[1] =  _PA[2] ;

            double _B2[2] ;
            _B2[0] =  _PB[0] ;
            _B2[1] =  _PB[2] ;

            double _P2[2] ;
            _P2[0] =  _PP[0] ;
            _P2[1] =  _PP[2] ;

            double _s2 =
                geompred::orient2d (
                (double*)_A2 ,
                (double*)_B2 ,
                (double*)_P2 )   ;

    /*----------------------- get orientation in yz-plane */
            double _A3[2] ;
            _A3[0] =  _PA[1] ;
            _A3[1] =  _PA[2] ;

            double _B3[2] ;
            _B3[0] =  _PB[1] ;
            _B3[1] =  _PB[2] ;

            double _P3[2] ;
            _P3[0] =  _PP[1] ;
            _P3[1] =  _PP[2] ;

            double _s3 =
                geompred::orient2d (
                (double*)_A3 ,
                (double*)_B3 ,
                (double*)_P3 )   ;

    /*----------------------- test intersection hierarchy */
            if (_s1 == +0.0 &&
                _s2 == +0.0 &&
                _s3 == +0.0 )
            {
    /*----------------------- orient w.r.t. line endpoint */
            double _sa = cleave3d (
                (double*)_PA ,
                (double*)_PB ,
                (double*)_PP )   ;

            double _sb = cleave3d (
                (double*)_PB ,
                (double*)_PA ,
                (double*)_PP )   ;

            if (_sa*_sb > +0.0)
            {
    /*----------------------- have node-line intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return edge_hits ;
            }
            else
            if (_sa == +0.0)
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pa[0] ;
            _qq[1] =  _pa[1] ;
            _qq[2] =  _pa[2] ;

            return node_hits ;
            }
            else
            if (_sb == +0.0)
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pb[0] ;
            _qq[1] =  _pb[1] ;
            _qq[2] =  _pb[2] ;

            return node_hits ;
            }
            else
            if ( !_bind )
            {
    /*----------------------- is unconstrained: edge hits */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return edge_hits ;
            }
            else
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }

            }
            else
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * LINE-LINE-KD: line-line intersections
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __normal_call hits_type line_line_2d (
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc, // line
    __const_ptr  (real_type) _pd,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true ,
                  char_type _part = +1
        )
    {
        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] )
        {
    /*----------------------- test node-line intersection */
            return node_line_2d (
                _pa, _pc, _pd, _qq, _bind) ;
        }
        else
        if (_pc[0] == _pd[0] &&
            _pc[1] == _pd[1] )
        {
    /*----------------------- test node-line intersection */
            return node_line_2d (
                _pc, _pa, _pb, _qq, _bind) ;
        }
        else
        {
    /*----------------------- test line-line intersection */
            double _PA[2] ;
            _PA[0] =  _pa[0] ;
            _PA[1] =  _pa[1] ;

            double _PB[2] ;
            _PB[0] =  _pb[0] ;
            _PB[1] =  _pb[1] ;

            double _PC[2] ;
            _PC[0] =  _pc[0] ;
            _PC[1] =  _pc[1] ;

            double _PD[2] ;
            _PD[0] =  _pd[0] ;
            _PD[1] =  _pd[1] ;

    /*----------------------- orient w.r.t. line endpoint */
            double  _sa =
                geompred::orient2d (
                (double*)_PA ,
                (double*)_PC ,
                (double*)_PD )   ;

            double  _sb =
                geompred::orient2d (
                (double*)_PB ,
                (double*)_PD ,
                (double*)_PC )   ;

            if (_bind)
            {
    /*----------------------- no intersections: null hits */
            if (_sa * _sb < +0.0 )
            return null_hits;
            }

            double  _sc =
                geompred::orient2d (
                (double*)_PC ,
                (double*)_PA ,
                (double*)_PB )   ;

            double  _sd =
                geompred::orient2d (
                (double*)_PD ,
                (double*)_PB ,
                (double*)_PA )   ;

            if (_bind)
            {
    /*----------------------- no intersections: null hits */
            if (_sc * _sd < +0.0 )
            return null_hits;
            }

            if (_sa == +0.0 &&
                _sb == +0.0 )
            {
    /*----------------------- colinear: test for overlaps */

          //std::cout << "line-line-2d" << std::endl;

            return null_hits ;

            }
            else
            if (_sa == +0.0 )
            {
    /*----------------------- have node-line intersection */
            _qq[0] =  _pa[0] ;
            _qq[1] =  _pa[1] ;

            return node_hits ;
            }
            else
            if (_sb == +0.0 )
            {
    /*----------------------- have node-line intersection */
            _qq[0] =  _pb[0] ;
            _qq[1] =  _pb[1] ;

            return node_hits ;
            }
            else
            if (_sc == +0.0 )
            {
    /*----------------------- have node-line intersection */
            _qq[0] =  _pc[0] ;
            _qq[1] =  _pc[1] ;

            return node_hits ;
            }
            else
            if (_sd == +0.0 )
            {
    /*----------------------- have node-line intersection */
            _qq[0] =  _pd[0] ;
            _qq[1] =  _pd[1] ;

            return node_hits ;
            }
            else
            {
    /*----------------------- have line-line intersection */
            double _mm [2*2] ;
            _mm[__ij(0,0,2)] = _PB[0]-_PA[0] ;
            _mm[__ij(1,0,2)] = _PB[1]-_PA[1] ;
            _mm[__ij(0,1,2)] = _PC[0]-_PD[0] ;
            _mm[__ij(1,1,2)] = _PC[1]-_PD[1] ;

            double _im [2*2] , _dm;
            math::inv_2x2(
                +2, _mm, +2, _im, _dm) ;

            double _rv [2*1] ;
            _rv[__ij(0,0,2)] = _PC[0]-_PA[0] ;
            _rv[__ij(1,0,2)] = _PC[1]-_PA[1] ;

            if (_dm == +0.0)
            return null_hits ;

            double _tu =
            _im[__ij(0,0,2)] * _rv [0] +
            _im[__ij(0,1,2)] * _rv [1] ;
            double _tv =
            _im[__ij(1,0,2)] * _rv [0] +
            _im[__ij(1,1,2)] * _rv [1] ;

            _tu       /= _dm ;
            _tv       /= _dm ;

            if (_bind)
            {
                _tu =
            std::min(+1.,std::max(+0., _tu)) ;
                _tv =
            std::min(+1.,std::max(+0., _tv)) ;
            }

            if (_part == +1)
            {
    /*----------------------- calc. intersection on [a,b] */
                dd_flt _WA = 1. - _tu ;
                dd_flt _WB = 0. + _tu ;

                dd_flt _FA[2] ;
                _FA[0] =  _pa[0] ;
                _FA[1] =  _pa[1] ;

                dd_flt _FB[2] ;
                _FB[0] =  _pb[0] ;
                _FB[1] =  _pb[1] ;

                dd_flt _QQ[2] ;
                _QQ[0] = _FA[0] * _WA +
                         _FB[0] * _WB ;
                _QQ[1] = _FA[1] * _WA +
                         _FB[1] * _WB ;

                _qq[0] =  _QQ[0] ;
                _qq[1] =  _QQ[1] ;
            }
            else
            if (_part == +2)
            {
    /*----------------------- calc. intersection on [c,d] */
                dd_flt _WC = 1. - _tv ;
                dd_flt _WD = 0. + _tv ;

                dd_flt _FC[2] ;
                _FC[0] =  _pc[0] ;
                _FC[1] =  _pc[1] ;

                dd_flt _FD[2] ;
                _FD[0] =  _pd[0] ;
                _FD[1] =  _pd[1] ;

                dd_flt _QQ[2] ;
                _QQ[0] = _FC[0] * _WC +
                         _FD[0] * _WD ;
                _QQ[1] = _FC[1] * _WC +
                         _FD[1] * _WD ;

                _qq[0] =  _QQ[0] ;
                _qq[1] =  _QQ[1] ;
            }
            else
            {
            __assert( false &&
            "line_line_2d: invalid part!!") ;
            }

            return edge_hits ;
            }

        }
    }

    template <
        typename  real_type
             >
    __normal_call hits_type line_line_3d (
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _pc, // line
    __const_ptr  (real_type) _pd,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true ,
                  char_type _part = +1
        )
    {

        __unreferenced (_part);

        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] &&
            _pa[2] == _pb[2] )
        {
    /*----------------------- test node-line intersection */
            return node_line_3d (
                _pa, _pc, _pd, _qq, _bind) ;
        }
        else
        if (_pc[0] == _pd[0] &&
            _pc[1] == _pd[1] &&
            _pc[2] == _pd[2] )
        {
    /*----------------------- test node-line intersection */
            return node_line_3d (
                _pc, _pa, _pb, _qq, _bind) ;
        }
        else
        {
    /*----------------------- test line-line intersection */
            /*
            double _PA[3] ;
            _PA[0] =  _pa[0] ;
            _PA[1] =  _pa[1] ;
            _PA[2] =  _pa[2] ;

            double _PB[3] ;
            _PB[0] =  _pb[0] ;
            _PB[1] =  _pb[1] ;
            _PB[2] =  _pb[2] ;

            double _PC[3] ;
            _PC[0] =  _pc[0] ;
            _PC[1] =  _pc[1] ;
            _PC[2] =  _pc[2] ;

            double _PD[3] ;
            _PD[0] =  _pd[0] ;
            _PD[1] =  _pd[1] ;
            _PD[2] =  _pd[2] ;
            */

        //!! actual line intersection...

        //!!std::cout << "line-line-3d" << std::endl;

            return null_hits ;

        }
    }

    /*
    --------------------------------------------------------
     * NODE-TRIA-KD: node-tria intersections
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __normal_call hits_type node_tria_2d (
    __const_ptr  (real_type) _pp, // node
    __const_ptr  (real_type) _p1, // tria
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true
        )
    {
        if (_p1[0] == _p2[0] &&
            _p1[1] == _p2[1] )
        {
    /*----------------------- test node-line intersection */
            return node_line_2d (
                _pp, _p2, _p3, _qq, _bind) ;
        }
        else
        if (_p2[0] == _p3[0] &&
            _p2[1] == _p3[1] )
        {
    /*----------------------- test node-line intersection */
            return node_line_2d (
                _pp, _p3, _p1, _qq, _bind) ;
        }
        else
        if (_p3[0] == _p1[0] &&
            _p3[1] == _p1[1] )
        {
    /*----------------------- test node-line intersection */
            return node_line_2d (
                _pp, _p1, _p2, _qq, _bind) ;
        }
        else
        {
    /*----------------------- test node-tria intersection */
            double _PP[2] ;
            _PP[0] =  _pp[0] ;
            _PP[1] =  _pp[1] ;

            double _P1[2] ;
            _P1[0] =  _p1[0] ;
            _P1[1] =  _p1[1] ;

            double _P2[2] ;
            _P2[0] =  _p2[0] ;
            _P2[1] =  _p2[1] ;

            double _P3[2] ;
            _P3[0] =  _p3[0] ;
            _P3[1] =  _p3[1] ;

    /*----------------------- orient w.r.t. tria vertices */
            double _s3 =
                geompred::orient2d (
                (double*)_P1 ,
                (double*)_P2 ,
                (double*)_PP )   ;

            double _s1 =
                geompred::orient2d (
                (double*)_P2 ,
                (double*)_P3 ,
                (double*)_PP )   ;

            double _s2 =
                geompred::orient2d (
                (double*)_P3 ,
                (double*)_P1 ,
                (double*)_PP )   ;

    /*----------------------- test intersection hierarchy */
            if (_s1*_s2<0.0 ||
                _s1*_s3<0.0 ||
                _s2*_s3<0.0 )
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }
            else
            if (_s1 == +0.0 &&
                _s2 == +0.0 &&
                _s3 == +0.0 )
            {
    /*----------------------- degenerate tria: check line */

        //!!std::cout << "node-tria-2d" << std::endl;

            return null_hits ;

            }
            else
            if (_s1 == +0.0 &&
                _s2 == +0.0 )
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return node_hits ;
            }
            else
            if (_s2 == +0.0 &&
                _s3 == +0.0 )
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return node_hits ;
            }
            else
            if (_s3 == +0.0 &&
                _s1 == +0.0 )
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return node_hits ;
            }
            else
            if (_s1 == +0.0 )
            {
    /*----------------------- have node-edge intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return edge_hits ;
            }
            else
            if (_s2 == +0.0 )
            {
    /*----------------------- have node-edge intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return edge_hits ;
            }
            else
            if (_s3 == +0.0 )
            {
    /*----------------------- have node-edge intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return edge_hits ;
            }
            else
            {
    /*----------------------- have node-tria intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;

            return face_hits ;
            }
        }
    }

    template <
        typename  real_type
             >
    __normal_call hits_type node_tria_3d (
    __const_ptr  (real_type) _pp, // node
    __const_ptr  (real_type) _p1, // tria
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true
        )
    {
        if (_p1[0] == _p2[0] &&
            _p1[1] == _p2[1] &&
            _p1[2] == _p2[2] )
        {
    /*----------------------- test node-line intersection */
            return node_line_3d (
                _pp, _p2, _p3, _qq, _bind) ;
        }
        else
        if (_p2[0] == _p3[0] &&
            _p2[1] == _p3[1] &&
            _p2[2] == _p3[2] )
        {
    /*----------------------- test node-line intersection */
            return node_line_3d (
                _pp, _p3, _p1, _qq, _bind) ;
        }
        else
        if (_p3[0] == _p1[0] &&
            _p3[1] == _p1[1] &&
            _p3[2] == _p1[2] )
        {
    /*----------------------- test node-line intersection */
            return node_line_3d (
                _pp, _p1, _p2, _qq, _bind) ;
        }
        else
        {
    /*----------------------- test node-tria intersection */
            double _PP[3] ;
            _PP[0] =  _pp[0] ;
            _PP[1] =  _pp[1] ;
            _PP[2] =  _pp[2] ;

            double _P1[3] ;
            _P1[0] =  _p1[0] ;
            _P1[1] =  _p1[1] ;
            _P1[2] =  _p1[2] ;

            double _P2[3] ;
            _P2[0] =  _p2[0] ;
            _P2[1] =  _p2[1] ;
            _P2[2] =  _p2[2] ;

            double _P3[3] ;
            _P3[0] =  _p3[0] ;
            _P3[1] =  _p3[1] ;
            _P3[2] =  _p3[2] ;

            double _ss =
                geompred::orient3d (
                (double*)_P1 ,
                (double*)_P2 ,
                (double*)_P3 ,
                (double*)_PP )  ;

            if (_ss == +0.0)
            {
    /*----------------------- get orientation in xy-plane */
            double _TP[2] ;
            _TP[0] =  _PP[0] ;
            _TP[1] =  _PP[1] ;

            double _T1[2] ;
            _T1[0] =  _P1[0] ;
            _T1[1] =  _P1[1] ;

            double _T2[2] ;
            _T2[0] =  _P2[0] ;
            _T2[1] =  _P2[1] ;

            double _T3[2] ;
            _T3[0] =  _P3[0] ;
            _T3[1] =  _P3[1] ;

            double _s3[3] ;
            _s3[0] = geompred::orient2d (
                (double*)_T1 ,
                (double*)_T2 ,
                (double*)_TP )  ;

            double _s1[3] ;
            _s1[0] = geompred::orient2d (
                (double*)_T2 ,
                (double*)_T3 ,
                (double*)_TP )  ;

            double _s2[3] ;
            _s2[0] = geompred::orient2d (
                (double*)_T3 ,
                (double*)_T1 ,
                (double*)_TP )  ;

            if (_s1[0]*_s2[0] < +0.0 ||
                _s2[0]*_s3[0] < +0.0 ||
                _s3[0]*_s1[0] < +0.0 )
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }

    /*----------------------- get orientation in xz-plane */
            _TP[0] =  _PP[0] ;
            _TP[1] =  _PP[2] ;

            _T1[0] =  _P1[0] ;
            _T1[1] =  _P1[2] ;

            _T2[0] =  _P2[0] ;
            _T2[1] =  _P2[2] ;

            _T3[0] =  _P3[0] ;
            _T3[1] =  _P3[2] ;

            _s3[1] = geompred::orient2d (
                (double*)_T1 ,
                (double*)_T2 ,
                (double*)_TP )  ;

            _s1[1] = geompred::orient2d (
                (double*)_T2 ,
                (double*)_T3 ,
                (double*)_TP )  ;

            _s2[1] = geompred::orient2d (
                (double*)_T3 ,
                (double*)_T1 ,
                (double*)_TP )  ;

            if (_s1[1]*_s2[1] < +0.0 ||
                _s2[1]*_s3[1] < +0.0 ||
                _s3[1]*_s1[1] < +0.0 )
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }

    /*----------------------- get orientation in yz-plane */
            _TP[0] =  _PP[1] ;
            _TP[1] =  _PP[2] ;

            _T1[0] =  _P1[1] ;
            _T1[1] =  _P1[2] ;

            _T2[0] =  _P2[1] ;
            _T2[1] =  _P2[2] ;

            _T3[0] =  _P3[1] ;
            _T3[1] =  _P3[2] ;

            _s3[2] = geompred::orient2d (
                (double*)_T1 ,
                (double*)_T2 ,
                (double*)_TP )  ;

            _s1[2] = geompred::orient2d (
                (double*)_T2 ,
                (double*)_T3 ,
                (double*)_TP )  ;

            _s2[2] = geompred::orient2d (
                (double*)_T3 ,
                (double*)_T1 ,
                (double*)_TP )  ;

            if (_s1[2]*_s2[2] < +0.0 ||
                _s2[2]*_s3[2] < +0.0 ||
                _s3[2]*_s1[2] < +0.0 )
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }

    /*----------------------- test intersection hierarchy */
            bool_type _z1 =
                _s1[0] == +0.0 &&
                _s1[1] == +0.0 &&
                _s1[2] == +0.0  ;

            bool_type _z2 =
                _s2[0] == +0.0 &&
                _s2[1] == +0.0 &&
                _s2[2] == +0.0  ;

            bool_type _z3 =
                _s3[0] == +0.0 &&
                _s3[1] == +0.0 &&
                _s3[2] == +0.0  ;

            if (_z1 && _z2 && _z3 )
            {
    /*----------------------- degenerate tria: check line */

        //!!std::cout << "node-tria-3d" << std::endl;

            return null_hits ;

            }
            else
            if (_z1 && _z2 )
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return edge_hits ;
            }
            else
            if (_z2 && _z3 )
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return node_hits ;
            }
            else
            if (_z3 && _z1 )
            {
    /*----------------------- have node-node intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return node_hits ;
            }
            else
            if (_z1 )
            {
    /*----------------------- have node-edge intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return edge_hits ;
            }
            else
            if (_z2 )
            {
    /*----------------------- have node-edge intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return edge_hits ;
            }
            else
            if (_z3 )
            {
    /*----------------------- have node-edge intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return edge_hits ;
            }
            else
            {
    /*----------------------- have node-tria intersection */
            _qq[0] =  _pp[0] ;
            _qq[1] =  _pp[1] ;
            _qq[2] =  _pp[2] ;

            return face_hits ;
            }

            }
            else
            {
    /*----------------------- no intersections: null hits */
            return null_hits ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * LINE-TRIA-KD: line-tria intersections
    --------------------------------------------------------
     */

    template <
        typename  real_type
             >
    __normal_call hits_type line_tria_3d (
    __const_ptr  (real_type) _pa, // line
    __const_ptr  (real_type) _pb,
    __const_ptr  (real_type) _p1, // tria
    __const_ptr  (real_type) _p2,
    __const_ptr  (real_type) _p3,
    __write_ptr  (real_type) _qq, // intersection
                  bool_type _bind = true ,
                  char_type _part = +1
        )
    {

        __unreferenced (_part);

        if (_pa[0] == _pb[0] &&
            _pa[1] == _pb[1] &&
            _pa[2] == _pb[2] )
        {
    /*----------------------- test node-tria intersection */
            return node_tria_3d (
                _pa, _p1, _p2, _p3, _qq, _bind) ;
        }
        else
        if (_p1[0] == _p2[0] &&
            _p1[1] == _p2[1] &&
            _p1[2] == _p2[2] )
        {
    /*----------------------- test line-line intersection */
            return line_line_3d (
            _pa, _pb, _p2, _p3, _qq, _bind) ;
        }
        else
        if (_p2[0] == _p3[0] &&
            _p2[1] == _p3[1] &&
            _p2[2] == _p3[2] )
        {
    /*----------------------- test line-line intersection */
            return line_line_3d (
            _pa, _pb, _p3, _p1, _qq, _bind) ;
        }
        else
        if (_p3[0] == _p1[0] &&
            _p3[1] == _p1[1] &&
            _p3[2] == _p1[2] )
        {
    /*----------------------- test line-line intersection */
            return line_line_3d (
            _pa, _pb, _p1, _p2, _qq, _bind) ;
        }
        else
        {
    /*----------------------- test line-tria intersection */
            double _PA[3] ;
            _PA[0] =  _pa[0] ;
            _PA[1] =  _pa[1] ;
            _PA[2] =  _pa[2] ;

            double _PB[3] ;
            _PB[0] =  _pb[0] ;
            _PB[1] =  _pb[1] ;
            _PB[2] =  _pb[2] ;

            double _P1[3] ;
            _P1[0] =  _p1[0] ;
            _P1[1] =  _p1[1] ;
            _P1[2] =  _p1[2] ;

            double _P2[3] ;
            _P2[0] =  _p2[0] ;
            _P2[1] =  _p2[1] ;
            _P2[2] =  _p2[2] ;

            double _P3[3] ;
            _P3[0] =  _p3[0] ;
            _P3[1] =  _p3[1] ;
            _P3[2] =  _p3[2] ;

    /*----------------------- test if line straddles tria */
            double _sa =
                geompred::orient3d (
                (double*)_P1 ,
                (double*)_P2 ,
                (double*)_P3 ,
                (double*)_PA )   ;

            double _sb =
                geompred::orient3d (
                (double*)_P1 ,
                (double*)_P3 ,
                (double*)_P2 ,
                (double*)_PB )   ;

            if (_bind)
            {
    /*----------------------- no intersections: null hits */
            if (_sa*_sb < 0.0)
            return null_hits ;
            }

    /*----------------------- test if tria straddles line */
            double _s1 =
                geompred::orient3d (
                (double*)_P1 ,
                (double*)_P2 ,
                (double*)_PA ,
                (double*)_PB )   ;

            double _s2 =
                geompred::orient3d (
                (double*)_P2 ,
                (double*)_P3 ,
                (double*)_PA ,
                (double*)_PB )   ;

            double _s3 =
                geompred::orient3d (
                (double*)_P3 ,
                (double*)_P1 ,
                (double*)_PA ,
                (double*)_PB )   ;

            if (_bind)
            {
    /*----------------------- no intersections: null hits */
            if (_s1*_s2 < 0.0)
            return null_hits ;

            if (_s2*_s3 < 0.0)
            return null_hits ;

            if (_s3*_s1 < 0.0)
            return null_hits ;
            }

            if (_sa == +0.0 &&
                _sb == +0.0 )
            {
        // line + tria in same plane

        //!!std::cout << "line-tria-3d" << std::endl;

            return null_hits ;

            }
            else
            if (_s1 == +0.0 &&
                _s2 == +0.0 )
            {
    /*----------------------- have line-node intersection */
            _qq[0] =  _p2[0] ;
            _qq[1] =  _p2[1] ;
            _qq[2] =  _p2[2] ;

            return node_hits ;
            }
            else
            if (_s2 == +0.0 &&
                _s3 == +0.0 )
            {
    /*----------------------- have line-node intersection */
            _qq[0] =  _p3[0] ;
            _qq[1] =  _p3[1] ;
            _qq[2] =  _p3[2] ;

            return node_hits ;
            }
            else
            if (_s3 == +0.0 &&
                _s1 == +0.0 )
            {
    /*----------------------- have line-node intersection */
            _qq[0] =  _p1[0] ;
            _qq[1] =  _p1[1] ;
            _qq[2] =  _p1[2] ;

            return node_hits ;
            }
            else
            if (_s1 == +0.0 )
            {
    /*----------------------- have line-edge intersection */
            double _WS = _s2 + _s3 ;
            double _W1 = _s2 / _WS ;
            double _W2 = _s3 / _WS ;

            dd_flt _F1[3] ;
            _F1[0] =  _p1[0] ;
            _F1[1] =  _p1[1] ;
            _F1[2] =  _p1[2] ;

            dd_flt _F2[3] ;
            _F2[0] =  _p2[0] ;
            _F2[1] =  _p2[1] ;
            _F2[2] =  _p2[2] ;

            dd_flt _QQ[3] ;
            _QQ[0] =  _F1[0] * _W1 +
                      _F2[0] * _W2 ;
            _QQ[1] =  _F1[1] * _W1 +
                      _F2[1] * _W2 ;
            _QQ[2] =  _F1[2] * _W1 +
                      _F2[2] * _W2 ;

            _qq[0] =  _QQ[0] ;
            _qq[1] =  _QQ[1] ;
            _qq[2] =  _QQ[2] ;

            return edge_hits ;
            }
            else
            if (_s2 == +0.0 )
            {
    /*----------------------- have line-edge intersection */
            double _WS = _s1 + _s3 ;
            double _W2 = _s3 / _WS ;
            double _W3 = _s1 / _WS ;

            dd_flt _F2[3] ;
            _F2[0] =  _p2[0] ;
            _F2[1] =  _p2[1] ;
            _F2[2] =  _p2[2] ;

            dd_flt _F3[3] ;
            _F3[0] =  _p3[0] ;
            _F3[1] =  _p3[1] ;
            _F3[2] =  _p3[2] ;

            dd_flt _QQ[3] ;
            _QQ[0] =  _F2[0] * _W2 +
                      _F3[0] * _W3 ;
            _QQ[1] =  _F2[1] * _W2 +
                      _F3[1] * _W3 ;
            _QQ[2] =  _F2[2] * _W2 +
                      _F3[2] * _W3 ;

            _qq[0] =  _QQ[0] ;
            _qq[1] =  _QQ[1] ;
            _qq[2] =  _QQ[2] ;

            return edge_hits ;
            }
            else
            if (_s3 == +0.0 )
            {
    /*----------------------- have line-edge intersection */
            double _WS = _s1 + _s2 ;
            double _W3 = _s1 / _WS ;
            double _W1 = _s2 / _WS ;

            dd_flt _F3[3] ;
            _F3[0] =  _p3[0] ;
            _F3[1] =  _p3[1] ;
            _F3[2] =  _p3[2] ;

            dd_flt _F1[3] ;
            _F1[0] =  _p1[0] ;
            _F1[1] =  _p1[1] ;
            _F1[2] =  _p1[2] ;

            dd_flt _QQ[3] ;
            _QQ[0] =  _F3[0] * _W3 +
                      _F1[0] * _W1 ;
            _QQ[1] =  _F3[1] * _W3 +
                      _F1[1] * _W1 ;
            _QQ[2] =  _F3[2] * _W3 +
                      _F1[2] * _W1 ;

            _qq[0] =  _QQ[0] ;
            _qq[1] =  _QQ[1] ;
            _qq[2] =  _QQ[2] ;

            return edge_hits ;
            }
            else
            {
    /*----------------------- have line-face intersection */
            double _WS = _s1 + _s2 + _s3 ;
            double _W1 = _s2 / _WS ;
            double _W2 = _s3 / _WS ;
            double _W3 = _s1 / _WS ;

            dd_flt _F1[3] ;
            _F1[0] =  _p1[0] ;
            _F1[1] =  _p1[1] ;
            _F1[2] =  _p1[2] ;

            dd_flt _F2[3] ;
            _F2[0] =  _p2[0] ;
            _F2[1] =  _p2[1] ;
            _F2[2] =  _p2[2] ;

            dd_flt _F3[3] ;
            _F3[0] =  _p3[0] ;
            _F3[1] =  _p3[1] ;
            _F3[2] =  _p3[2] ;

            dd_flt _QQ[3] ;
            _QQ[0] =  _F1[0] * _W1 +
                      _F2[0] * _W2 +
                      _F3[0] * _W3 ;
            _QQ[1] =  _F1[1] * _W1 +
                      _F2[1] * _W2 +
                      _F3[1] * _W3 ;
            _QQ[2] =  _F1[2] * _W1 +
                      _F2[2] * _W2 +
                      _F3[2] * _W3 ;

            _qq[0] =  _QQ[0] ;
            _qq[1] =  _QQ[1] ;
            _qq[2] =  _QQ[2] ;

            return face_hits ;
            }
        }
    }


    }

#   endif   //__INTERSECT_K__



