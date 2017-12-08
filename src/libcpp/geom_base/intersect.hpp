

#pragma once

#ifndef __INTERSECT__
#define __INTERSECT__

//!! TODO:
//!! should use dd_flt wherever an intersection
//!!    is interpolated from existing points
//!!
//!! should apply a careful interpolation calc.
//!!    a'la line_tria -- any existing routines
//!!    that compute a "midpoint" need updating
//!!

namespace geometry
{

/*
------------------------------------------------------------
 intersection predicates in R^2/R^3.
------------------------------------------------------------
 */
    template <
    typename      data_type
             >
    __normal_call bool_type polyroots (
        data_type _aa,      // aa*xx^2+bb*xx+cc=0
        data_type _bb,
        data_type _cc,
    __write_ptr  (data_type) _xx
        )
    {
        bool_type _real = false ;

        data_type _sq = _bb * _bb -
       (data_type)+4. * _aa * _cc ;

		if (_sq >= (data_type)+0.)  // real roots
		{
        	_sq  = std::sqrt(_sq) ;
		
            _real = true;

            _xx[0] = (-_bb + _sq) ;
			_xx[1] = (-_bb - _sq) ;

            data_type _xm = std::max (
                std::abs(_xx[0]), 
                std::abs(_xx[1])) ;

            data_type _rt = 
        +std::numeric_limits<data_type>::epsilon();

            if (_aa >= _xm * _rt)
            {
                _aa *=(data_type)+2.;

                _xx[0] /= _aa ;
                _xx[1] /= _aa ;
            }
            else
            {
                _xx[0] = -_cc / _bb ;
                _xx[1] = -_cc / _bb ;
            }
        }

        return _real ;
    }

    template <
    typename      data_type
             >
	__normal_call size_t ball_line_2d (
    __const_ptr  (data_type) _pc, // ball BB(pc,rc)
		data_type            _rc,
	__const_ptr  (data_type) _pa, // line
	__const_ptr  (data_type) _pb,
	__write_ptr  (data_type) _qa, // intersections
	__write_ptr  (data_type) _qb
		)
	{
		size_t _nn = +0;
		
		data_type _pm[2] = {
	   (data_type)+.5 * (_pa[0]+_pb[0]),
	   (data_type)+.5 * (_pa[1]+_pb[1])
		    } ;
		data_type _pd[2] = {
	   (data_type)+.5 * (_pb[0]-_pa[0]),
	   (data_type)+.5 * (_pb[1]-_pa[1])
		    } ;
		data_type _mc[2] = {
	   (data_type)+1. * (_pm[0]-_pc[0]),
	   (data_type)+1. * (_pm[1]-_pc[1])
		    } ;
		
        data_type _aa = dot_2d(_pd, _pd) ;
        data_type _bb = dot_2d(_pd, _mc) *
                       (data_type) +2. ;	
        data_type _cc = dot_2d(_mc, _mc) ;
        _cc -= _rc * _rc ;
		
        data_type _tt[2] ;
        if (polyroots(_aa, _bb, _cc, _tt))
        {
        if (_tt[0] >= (data_type)-1. &&
			_tt[0] <= (data_type)+1. )
		{
		    data_type *_qq = _nn++ == +0 
		                   ? _qa 
		                   : _qb ;
		    
		    dd_flt _WB = _tt[0] ;
		    dd_flt _WA = _tt[0] ;
		    
		    _WB = (dd_flt)+1. + _WB;
		    _WB = (dd_flt)+.5 * _WB;
		    
		    _WA = (dd_flt)+1. - _WA;
		    _WA = (dd_flt)+.5 * _WA;
		  
		    if (_WA > _WB)
		    _WA = (dd_flt)+1. - _WB;
		    else
		    _WB = (dd_flt)+1. - _WA;
		    
		    dd_flt _PA[2] ;
		    _PA[0]=_pa[0] ;
		    _PA[1]=_pa[1] ;
		    
		    dd_flt _PB[2] ;
		    _PB[0]=_pb[0] ;
		    _PB[1]=_pb[1] ;
		    
		    dd_flt _QQ[2] ;
		    _QQ[0]=_PA[0] * _WA +
		           _PB[0] * _WB ;
		    _QQ[1]=_PA[1] * _WA +
		           _PB[1] * _WB ;
		           
		    _qq[0]=_QQ[0] ;
		    _qq[1]=_QQ[1] ;
		}

		if (_tt[1] >= (data_type)-1. &&
			_tt[1] <= (data_type)+1. )
		{
		    data_type *_qq = _nn++ == +0  
		                   ? _qa 
		                   : _qb ;
		                   
		    dd_flt _WB = _tt[1] ;
		    dd_flt _WA = _tt[1] ;
		    
		    _WB = (dd_flt)+1. + _WB;
		    _WB = (dd_flt)+.5 * _WB;
		    
		    _WA = (dd_flt)+1. - _WA;
		    _WA = (dd_flt)+.5 * _WA;
		  
		    if (_WA > _WB)
		    _WA = (dd_flt)+1. - _WB;
		    else
		    _WB = (dd_flt)+1. - _WA;
		    
		    dd_flt _PA[2] ;
		    _PA[0]=_pa[0] ;
		    _PA[1]=_pa[1] ;
		    
		    dd_flt _PB[2] ;
		    _PB[0]=_pb[0] ;
		    _PB[1]=_pb[1] ;
		    
		    dd_flt _QQ[2] ;
		    _QQ[0]=_PA[0] * _WA +
		           _PB[0] * _WB ;
		    _QQ[1]=_PA[1] * _WA +
		           _PB[1] * _WB ;
		           
		    _qq[0]=_QQ[0] ;
		    _qq[1]=_QQ[1] ;
		}
        }

        return ( _nn ) ; // return num roots
	}

    template <
    typename      data_type
             >
	__normal_call size_t ball_line_3d (
    __const_ptr  (data_type) _pc, // ball BB(pc,rc)
		data_type            _rc,
	__const_ptr  (data_type) _pa, // line
	__const_ptr  (data_type) _pb,
	__write_ptr  (data_type) _qa, // intersections
	__write_ptr  (data_type) _qb
		)
	{
        size_t _nn = +0;
		
		data_type _pm[3] = {
	   (data_type)+.5 * (_pa[0]+_pb[0]),
	   (data_type)+.5 * (_pa[1]+_pb[1]),
	   (data_type)+.5 * (_pa[2]+_pb[2])
		    } ;
		data_type _pd[3] = {
	   (data_type)+.5 * (_pb[0]-_pa[0]),
	   (data_type)+.5 * (_pb[1]-_pa[1]),
	   (data_type)+.5 * (_pb[2]-_pa[2])
		    } ;
		data_type _mc[3] = {
	   (data_type)+1. * (_pm[0]-_pc[0]),
	   (data_type)+1. * (_pm[1]-_pc[1]),
	   (data_type)+1. * (_pm[2]-_pc[2])
		    } ;
		
        data_type _aa = dot_3d(_pd, _pd) ;
        data_type _bb = dot_3d(_pd, _mc) *
                       (data_type) +2. ;	
        data_type _cc = dot_3d(_mc, _mc) ;
        _cc -= _rc * _rc ;
		
        data_type _tt[2] ;
        if (polyroots(_aa, _bb, _cc, _tt))
        {
        if (_tt[0] >= (data_type)-1. &&
			_tt[0] <= (data_type)+1. )
		{
		    data_type *_qq = _nn++ == +0 
		                   ? _qa 
		                   : _qb ;
		     
		    dd_flt _WB = _tt[0] ;
		    dd_flt _WA = _tt[0] ;
		    
		    _WB = (dd_flt)+1. + _WB;
		    _WB = (dd_flt)+.5 * _WB;
		    
		    _WA = (dd_flt)+1. - _WA;
		    _WA = (dd_flt)+.5 * _WA;
		  
		    if (_WA > _WB)
		    _WA = (dd_flt)+1. - _WB;
		    else
		    _WB = (dd_flt)+1. - _WA;
		    
		    dd_flt _PA[3] ;
		    _PA[0]=_pa[0] ;
		    _PA[1]=_pa[1] ;
		    _PA[2]=_pa[2] ;
		    
		    dd_flt _PB[3] ;
		    _PB[0]=_pb[0] ;
		    _PB[1]=_pb[1] ;
		    _PB[2]=_pb[2] ;
		    
		    dd_flt _QQ[3] ;
		    _QQ[0]=_PA[0] * _WA +
		           _PB[0] * _WB ;
		    _QQ[1]=_PA[1] * _WA +
		           _PB[1] * _WB ;
		    _QQ[2]=_PA[2] * _WA +
		           _PB[2] * _WB ; 
		           
		    _qq[0]=_QQ[0] ;
		    _qq[1]=_QQ[1] ;
		    _qq[2]=_QQ[2] ;
		}

		if (_tt[1] >= (data_type)-1. &&
			_tt[1] <= (data_type)+1. )
		{
		    data_type *_qq = _nn++ == +0  
		                   ? _qa 
		                   : _qb ;
		   
		    dd_flt _WB = _tt[1] ;
		    dd_flt _WA = _tt[1] ;
		    
		    _WB = (dd_flt)+1. + _WB;
		    _WB = (dd_flt)+.5 * _WB;
		    
		    _WA = (dd_flt)+1. - _WA;
		    _WA = (dd_flt)+.5 * _WA;
		  
		    if (_WA > _WB)
		    _WA = (dd_flt)+1. - _WB;
		    else
		    _WB = (dd_flt)+1. - _WA;
		    
		    dd_flt _PA[3] ;
		    _PA[0]=_pa[0] ;
		    _PA[1]=_pa[1] ;
		    _PA[2]=_pa[2] ;
		    
		    dd_flt _PB[3] ;
		    _PB[0]=_pb[0] ;
		    _PB[1]=_pb[1] ;
		    _PB[2]=_pb[2] ;
		    
		    dd_flt _QQ[3] ;
		    _QQ[0]=_PA[0] * _WA +
		           _PB[0] * _WB ;
		    _QQ[1]=_PA[1] * _WA +
		           _PB[1] * _WB ;
		    _QQ[2]=_PA[2] * _WA +
		           _PB[2] * _WB ; 
		           
		    _qq[0]=_QQ[0] ;
		    _qq[1]=_QQ[1] ;
		    _qq[2]=_QQ[2] ;
		}
        }

        return ( _nn ) ; // return num roots
	}

    template <
    typename      data_type
             >
	__normal_call bool line_flat_3d (
	__const_ptr  (data_type) _pp, // (xx-pp).nv=0
	__const_ptr  (data_type) _nv,
	__const_ptr  (data_type) _pa, // line
	__const_ptr  (data_type) _pb,
	__write_ptr  (data_type) _qq, // intersection
	    bool _bind = true
		)
	{
		data_type _pm[3] = {
	   (data_type)+.5 * (_pa[0]+_pb[0]),
	   (data_type)+.5 * (_pa[1]+_pb[1]),
	   (data_type)+.5 * (_pa[2]+_pb[2])
		    } ;
		data_type _pd[3] = {
	   (data_type)+.5 * (_pb[0]-_pa[0]),
	   (data_type)+.5 * (_pb[1]-_pa[1]),
	   (data_type)+.5 * (_pb[2]-_pa[2])
		    } ;
		data_type _pv[3] = {
	   (data_type)+1. * (_pp[0]-_pm[0]),
	   (data_type)+1. * (_pp[1]-_pm[1]),
	   (data_type)+1. * (_pp[2]-_pm[2])
		    } ;

        data_type _ep = 
        std::numeric_limits<data_type>::epsilon() ;

		data_type _d1 = 
		geometry::dot_3d(_pv, _nv) ; 
		data_type _d2 = 
		geometry::dot_3d(_pd, _nv) ;

        if (std::abs(_d2) <= _ep * std::abs(_d1))
            return ( false ) ;
		
		data_type _tt =  _d1 / _d2 ;
		
		if (_bind)
		{
		if (_tt  < (data_type)-1.)
			return ( false ) ;
		if (_tt  > (data_type)+1.)
			return ( false ) ;
        }

		_qq[0] = _pm[0] + _tt * _pd[0] ;
		_qq[1] = _pm[1] + _tt * _pd[1] ;
		_qq[2] = _pm[2] + _tt * _pd[2] ;	

		return ( true ) ;
	}

    template <
    typename      data_type
             >
	__normal_call bool line_tria_3d__ (
	__const_ptr  (data_type) _p1, // tria
	__const_ptr  (data_type) _p2,
	__const_ptr  (data_type) _p3,
	__const_ptr  (data_type) _pa, // line
	__const_ptr  (data_type) _pb,
	    data_type& _tt,
	__write_ptr  (data_type) _qq,
	    bool _bind = false
		)
	{
		data_type _nv[3] ;
		tria_norm_3d (_p1,_p2,_p3,_nv) ;

        data_type _pm[3] = {
	   (data_type)+.5 * (_pa[0]+_pb[0]),
	   (data_type)+.5 * (_pa[1]+_pb[1]),
	   (data_type)+.5 * (_pa[2]+_pb[2])
		    } ;
		data_type _pd[3] = {
	   (data_type)+.5 * (_pb[0]-_pa[0]),
	   (data_type)+.5 * (_pb[1]-_pa[1]),
	   (data_type)+.5 * (_pb[2]-_pa[2])
		    } ;
		data_type _pv[3] = {
	   (data_type)+1. * (_p1[0]-_pm[0]),
	   (data_type)+1. * (_p1[1]-_pm[1]),
	   (data_type)+1. * (_p1[2]-_pm[2])
		    } ;

		data_type _ta = 
		    geometry::dot_3d(_pv, _nv) ; 
		data_type _tb = 
		    geometry::dot_3d(_pd, _nv) ;

        _tt = _ta/_tb ;

   //!! std::numeric_limits<data_type>::epsilon();
        data_type _rt = (data_type)+1.0E-14;

        bool_type _ok = std::abs(_tb) > std::abs(_ta) * _rt ;

        if (_bind)
        {
        _tt = std::min((data_type)+1., _tt);
        _tt = std::max((data_type)-1., _tt);
        }

        _qq[0] = _pm[0] + _tt * _pd[0] ;
		_qq[1] = _pm[1] + _tt * _pd[1] ;
		_qq[2] = _pm[2] + _tt * _pd[2] ;
		
		return ( _ok ) ;
		
    //!!return ( std::isfinite(_tt) ) ;
	}

    template <
    typename      data_type
             >
	__normal_call size_t tria_flat_3d (
	__const_ptr  (data_type) _pp, // (xx-pp).nv=0
	__const_ptr  (data_type) _nv,
	__const_ptr  (data_type) _pa, // tria
	__const_ptr  (data_type) _pb,
	__const_ptr  (data_type) _pc,
	__write_ptr  (data_type) _qa, // intersection
	__write_ptr  (data_type) _qb
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









template <typename _Iptr>
	__inline_call _Iptr coin_flip (
		_Iptr  _inum
		)
	{
		double _rval = 
			((double) std::rand()) 
				/ (RAND_MAX + 1.0)  ;
		double _ival = +1.0 / _inum ;
		return(_Iptr)(_rval / _ival);
	}



template <typename _Type> 
	__inline_call bool line_line_2d (
		__const_ptr(_Type) _p1,
		__const_ptr(_Type) _p2,
		__const_ptr(_Type) _q1,
		__const_ptr(_Type) _q2,
		_Type &_tp,
		_Type &_tq
		)
    {
        _Type _vp[2], _vq[2], _vv[2];
        geometry::vector_2d(_p1, _p2, _vp);
        geometry::vector_2d(_q1, _q2, _vq);
        geometry::vector_2d(_q1, _p1, _vv);

        _Type _dd = _vq[1] * _vp[0] - 
                    _vq[0] * _vp[1] ;

        _tp = (_vq[0] * _vv[1] - 
               _vq[1] * _vv[0]) / _dd;
        _tq = (_vp[0] * _vv[1] - 
               _vp[1] * _vv[0]) / _dd;		

        return ( std::isfinite(_tp) &&
                 std::isfinite(_tq) );
    } 

template <typename _Type> 
	__inline_call void point_line_2d (
		__const_ptr(_Type) _pp,
		__const_ptr(_Type) _p1,
		__const_ptr(_Type) _p2,
		_Type &_tt
		)
	{
		_Type _P1[2], _P2[2], _VV[2];
        geometry::vector_2d(_pp, _p1, _P1); // translate to origin
        geometry::vector_2d(_pp, _p2, _P2);
        geometry::vector_2d(_P1, _P2, _VV);

        _Type _LL = geometry::lensqr_2d(_VV);
        
        _tt = -(_P1[0] * _VV[0] + 
                _P1[1] * _VV[1] ) / _LL ;
    }
template <typename _Type> 
	__inline_call void point_line_2d (
		__const_ptr(_Type) _pp,
		__const_ptr(_Type) _p1,
		__const_ptr(_Type) _p2,
		__write_ptr(_Type) _qq,
        bool _bind = true
		)
	{
        _Type _tt;
        geometry::point_line_2d(_pp,_p1,_p2,_tt) ;

        if (_tt < (_Type)+.0)
            _tt = (_Type)+.0;
        else
        if (_tt > (_Type)+1.)
            _tt = (_Type)+1.;

        _qq[0] = _p1[0] + _tt * (_p2[0] - _p1[0]);
        _qq[1] = _p1[1] + _tt * (_p2[1] - _p1[1]);
    }
template <typename _Type> 
	__inline_call void point_line_3d (
		__const_ptr(_Type) _pp,
		__const_ptr(_Type) _p1,
		__const_ptr(_Type) _p2,
		_Type &_tt
		)
	{
		_Type _P1[3], _P2[3], _VV[3];
        geometry::vector_3d(_pp, _p1, _P1); // translate to origin
        geometry::vector_3d(_pp, _p2, _P2);
        geometry::vector_3d(_P1, _P2, _VV);

        _Type _LL = geometry::lensqr_3d(_VV);
        
        _tt = -(_P1[0] * _VV[0] +
                _P1[1] * _VV[1] +
                _P1[2] * _VV[2] ) / _LL ;
    }
template <typename _Type> 
	__inline_call void point_line_3d (
		__const_ptr(_Type) _pp,
		__const_ptr(_Type) _p1,
		__const_ptr(_Type) _p2,
		__write_ptr(_Type) _qq,
        bool _bind = true
		)
	{
        _Type _tt;
        geometry::point_line_3d(_pp,_p1,_p2,_tt) ;

        if (_tt < (_Type)+.0)
            _tt = (_Type)+.0;
        else
        if (_tt > (_Type)+1.)
            _tt = (_Type)+1.;

        _qq[0] = _p1[0] + _tt * (_p2[0] - _p1[0]);
        _qq[1] = _p1[1] + _tt * (_p2[1] - _p1[1]);
        _qq[2] = _p1[2] + _tt * (_p2[2] - _p1[2]);
    }

template <typename _Type>
	__normal_call void point_plane_3d (
		__const_ptr(_Type) _pp,
		__write_ptr(_Type) _qq,
		__const_ptr(_Type) _pi,
		__const_ptr(_Type) _nv
		)
	{
		_Type _vv[3];
		_vv[0] = _pi[0] - _pp[0] ;
		_vv[1] = _pi[1] - _pp[1] ;
		_vv[2] = _pi[2] - _pp[2] ;
		
		_Type _tt = 
		geometry::dot_3d(_vv, _nv) / 
		geometry::dot_3d(_nv, _nv) ;

		_qq[0] = _pp[0]+_tt*_nv[0] ;
		_qq[1] = _pp[1]+_tt*_nv[1] ;
		_qq[2] = _pp[2]+_tt*_nv[2] ;
	}

template <typename _Type>
	__normal_call void point_tria_3d (
		__const_ptr(_Type) _pp,
		__write_ptr(_Type) _qq,
		__const_ptr(_Type) _p1,
		__const_ptr(_Type) _p2,
		__const_ptr(_Type) _p3
		)
	{
		_Type _vv[3];
		_vv[0] = _p1[0] - _pp[0] ;
		_vv[1] = _p1[1] - _pp[1] ;
		_vv[2] = _p1[2] - _pp[2] ;
		
		_Type _nv[3];
		geometry::tria_norm_3d(_p1, _p2, _p3, _nv);

		_Type _tt = 
		geometry::dot_3d(_vv, _nv) / 
		geometry::dot_3d(_nv, _nv) ;

		_qq[0] = _pp[0]+_tt*_nv[0] ;
		_qq[1] = _pp[1]+_tt*_nv[1] ;
		_qq[2] = _pp[2]+_tt*_nv[2] ;
	}





    typedef char_type hits_type;

    hits_type null_hits = +0;
    hits_type node_hits = +1;
    hits_type edge_hits = +2;
    hits_type face_hits = +3;

    __normal_call hits_type line_line_2d (
	__const_ptr  (double) _pa,
	__const_ptr  (double) _pb,
	__const_ptr  (double) _pc,
	__const_ptr  (double) _pd,
	__write_ptr  (double) _qq,
		bool _bind = true ,
		char _part = +1
		)
	{
        hits_type _hits = edge_hits;

        double  _sa = 
            geompred::orient2d (
       (double*)_pa, (double*)_pc, 
       (double*)_pd) ;
        double  _sb = 
            geompred::orient2d (
       (double*)_pb, (double*)_pd, 
       (double*)_pc) ;

    //!! all zero cases??

        if (_bind)
        {
        if (_sa * _sb < (double)+.0) 
            return ( null_hits )  ;
        }

        double  _sc = 
            geompred::orient2d (
       (double*)_pc, (double*)_pa, 
       (double*)_pb) ;
        double  _sd = 
            geompred::orient2d (
       (double*)_pd, (double*)_pb, 
       (double*)_pa) ;

        if (_bind)
        {
        if (_sc * _sd < (double)+.0) 
            return ( null_hits )  ;
        }

        if (_sa == (double) +0.)
        {
            _qq[0] = _pa[0] ;
            _qq[1] = _pa[1] ;
            
            return ( node_hits ) ;
        }
        if (_sb == (double) +0.)
        {
            _qq[0] = _pb[0] ;
            _qq[1] = _pb[1] ;
            
            return ( node_hits ) ;
        }     
        if (_sc == (double) +0.)
        {
            _qq[0] = _pc[0] ;
            _qq[1] = _pc[1] ;
            
            return ( node_hits ) ;
        }
        if (_sd == (double) +0.)
        {
            _qq[0] = _pd[0] ;
            _qq[1] = _pd[1] ;
            
            return ( node_hits ) ;
        }

        if (_part == +1) //!! outside bounds?
        {
            dd_flt _SA = _sb ;
            dd_flt _SB = _sa ;
            dd_flt _SS = _SA + _SB ;

            _SA /= _SS ;
            _SB /= _SS ;
            
            if (_SA > _SB)
            _SA = (dd_flt)+1. -_SB ;
            else
            _SB = (dd_flt)+1. -_SA ;
            
            dd_flt  _PA[2] ;
            _PA[0] =_pa[0] ;
            _PA[1] =_pa[1] ;

            dd_flt  _PB[2] ;
            _PB[0] =_pb[0] ;
            _PB[1] =_pb[1] ;

            dd_flt  _QQ[2] ;
            _QQ[0] =_SA * _PA[0] + 
                    _SB * _PB[0] ;
            _QQ[1] =_SA * _PA[1] + 
                    _SB * _PB[1] ;
                    
            _qq[0] =_QQ[0] ;
            _qq[1] =_QQ[1] ;
        }
        else
        if (_part == +2)
        {
            dd_flt _SC = _sd ;
            dd_flt _SD = _sc ;
            dd_flt _SS = _SC + _SD ;

            _SC /= _SS ;
            _SD /= _SS ;
            
            if (_SC > _SD)
            _SC = (dd_flt)+1. -_SD ;
            else
            _SD = (dd_flt)+1. -_SC ;
            
            dd_flt  _PC[2] ;
            _PC[0] =_pc[0] ;
            _PC[1] =_pc[1] ;

            dd_flt  _PD[2] ;
            _PD[0] =_pd[0] ;
            _PD[1] =_pd[1] ;

            dd_flt  _QQ[2] ;
            _QQ[0] =_SC * _PC[0] + 
                    _SD * _PD[0] ;
            _QQ[1] =_SC * _PC[1] + 
                    _SD * _PD[1] ;
                    
            _qq[0] =_QQ[0] ;
            _qq[1] =_QQ[1] ;     
        }
        else
        {
        __assert( false && "line_line_2d: invalid part!") ;
        }

        return ( _hits ) ;
    }
          
    __normal_call hits_type line_tria_3d (
	__const_ptr  (double) _pa,  // line
	__const_ptr  (double) _pb,
	__const_ptr  (double) _p1,  // tria
	__const_ptr  (double) _p2,
	__const_ptr  (double) _p3,
	__write_ptr  (double) _qq,  // intersection
		bool _bind = true ,
		char _part =   +1
		)
	{
        hits_type _hits ;

    // does the line straddle the tria??
        double  _sa = 
            geompred::orient3d (
       (double*)_p1, (double*)_p2, 
       (double*)_p3, (double*)_pa) ;
        
        double  _sb = 
            geompred::orient3d (
       (double*)_p1, (double*)_p3, 
       (double*)_p2, (double*)_pb) ;

        if (_bind)
        {
        if (_sa*_sb<(double)+.0)
            return null_hits ;
            
        if (_sa == (double)+.0 && // totally degenerate??
            _sb == (double)+.0 )
        {     
        if (_pa [0] == _pb [0] &&
            _pa [1] == _pb [1] &&
            _pa [2] == _pb [2] )
        {
            ;  // allow line to degenerate into a point??
        }
        else
            return null_hits ;
        }
        }

    // does the tria straddle the line??
        double  _s1 = 
            geompred::orient3d (
       (double*)_p1, (double*)_p2,
       (double*)_pa, (double*)_pb) ;

        double  _s2 = 
            geompred::orient3d (
       (double*)_p2, (double*)_p3,
       (double*)_pa, (double*)_pb) ;
        
        double  _s3 = 
            geompred::orient3d (
       (double*)_p3, (double*)_p1,
       (double*)_pa, (double*)_pb) ;

        // genuine outside
        if (_s1*_s2<(double)+.0)
            return null_hits ;

        if (_s1*_s3<(double)+.0)
            return null_hits ;
        
        if (_s1 == (double)+.0 && // totally degenerate??
            _s2 == (double)+.0 &&
            _s3 == (double)+.0 )
            return null_hits ;
        
        // on node
        if (_sa == (double)+.0 )
        {
            _qq[0] = _pa[0];
            _qq[1] = _pa[1];
            _qq[2] = _pa[2];
            
            return node_hits ;
        }
        if (_sb == (double)+.0 )
        {
            _qq[0] = _pb[0];
            _qq[1] = _pb[1];
            _qq[2] = _pb[2];
            
            return node_hits ;
        }
        
        // on node
        if (_s1 == (double)+.0 &&
            _s2 == (double)+.0 )
        {
            _qq[0] = _p2[0];
            _qq[1] = _p2[1];
            _qq[2] = _p2[2];
            
            return node_hits ;
        }
        if (_s2 == (double)+.0 &&
            _s3 == (double)+.0 )
        {
            _qq[0] = _p3[0];
            _qq[1] = _p3[1];
            _qq[2] = _p3[2];
            
            return node_hits ;
        }   
        if (_s3 == (double)+.0 &&
            _s1 == (double)+.0 )
        {
            _qq[0] = _p1[0];
            _qq[1] = _p1[1];
            _qq[2] = _p1[2];
            
            return node_hits ;
        }
        
        // on edge
        if (_s1 == (double)+.0 ||
            _s2 == (double)+.0 ||
            _s3 == (double)+.0 )
            _hits = edge_hits ;
        else
            _hits = face_hits ;



        if (_part == +1)
        {
        
        
        
        /*
            double _tt ;
            if (!line_tria_3d__ (
	       (double*) _p1 ,
           (double*) _p2 , 
           (double*) _p3 ,
	       (double*) _pa ,
	       (double*) _pb ,
	        _tt, _qq ) )
	        return null_hits ;
	        
	        
	        
	        if (_tt < 0.)
	        _tt += 1.e-8;
	        else
	        _tt -= 1.e-8;
	        
	        
	        
	        dd_flt _WB = _tt ;
	        dd_flt _WA = _tt ;
	        
	        _WB = (dd_flt)+1. + _WB;
	        _WB = (dd_flt)+.5 * _WB;
	        
	        _WA = (dd_flt)+1. - _WA;
	        _WA = (dd_flt)+.5 * _WA;
	      
	        if (_WA > _WB)
	        _WA = (dd_flt)+1. - _WB;
	        else
	        _WB = (dd_flt)+1. - _WA;
	        
	        dd_flt _PA[3] ;
	        _PA[0]=_pa[0] ;
	        _PA[1]=_pa[1] ;
	        _PA[2]=_pa[2] ;
	        
	        dd_flt _PB[3] ;
	        _PB[0]=_pb[0] ;
	        _PB[1]=_pb[1] ;
	        _PB[2]=_pb[2] ;
	        
	        dd_flt _QQ[3] ;
	        _QQ[0]=_PA[0] * _WA +
	               _PB[0] * _WB ;
	        _QQ[1]=_PA[1] * _WA +
	               _PB[1] * _WB ;
	        _QQ[2]=_PA[2] * _WA +
	               _PB[2] * _WB ; 
	               
	        _qq[0]=_QQ[0] ;
	        _qq[1]=_QQ[1] ;
	        _qq[2]=_QQ[2] ;
	            
	    */
	        
	    
        
        
        
        
        
            dd_flt _SA = _sa ;
            dd_flt _SB = _sb ;
            dd_flt _SS = _SA + _SB ;

            dd_flt _WW[2];
            _WW[0] = _SB / _SS ;
            _WW[1] = _SA / _SS ;
            
            
            
            /*
            char _ii, _id = +0;
            for (_ii = +2; _ii-- != +1; )
            {
            if (_WW[_ii] > _WW[_id])
                _id = _ii;
            }
            
            if (_id == +0)
            _WW[0] = (dd_flt)+1.-_WW[1] ;
            else
            if (_id == +1)
            _WW[1] = (dd_flt)+1.-_WW[0] ;
             */
             
             
            
              
            if (_WW[0] < (dd_flt)+0.)
                _WW[0] = (dd_flt)+0.;
            else
            if (_WW[0] > (dd_flt)+1.)
                _WW[0] = (dd_flt)+1.;
       
            _WW[1] = (dd_flt)+1. - _WW[0] ; 
            
            
            
            dd_flt  _PA[3] ;
            _PA[0] =_pa[0] ;
            _PA[1] =_pa[1] ;
            _PA[2] =_pa[2] ;

            dd_flt  _PB[3] ;
            _PB[0] =_pb[0] ;
            _PB[1] =_pb[1] ;
            _PB[2] =_pb[2] ;

            dd_flt  _QQ[3] ;
            _QQ[0] =_WW[0] * _PA[0] + 
                    _WW[1] * _PB[0] ;
            _QQ[1] =_WW[0] * _PA[1] + 
                    _WW[1] * _PB[1] ;
            _QQ[2] =_WW[0] * _PA[2] + 
                    _WW[1] * _PB[2] ;
                    
            _qq[0] =_QQ[0] ;
            _qq[1] =_QQ[1] ;
            _qq[2] =_QQ[2] ;
           
            
            
        }
        else
        if (_part == +2)
        {   
            dd_flt _S1 = _s1 ;
            dd_flt _S2 = _s2 ;
            dd_flt _S3 = _s3 ;
            
            dd_flt _WW[4];
            _WW[3] = _S1 + _S2 + _S3 ;
            _WW[0] = _S2 / _WW[3] ;
            _WW[1] = _S3 / _WW[3] ;
            _WW[2] = _S1 / _WW[3] ;
   
   
            /*
            char _ii, _id = +0;
            for (_ii = +3; _ii-- != +1; )
            {
            if (_WW[_ii] > _WW[_id])
                _id = _ii;
            }
            
            if (_id == +0)
            _WW[0] = (dd_flt)+1.-_WW[1] -
                                 _WW[2] ;
            else
            if (_id == +1)
            _WW[1] = (dd_flt)+1.-_WW[0] -
                                 _WW[2] ;
            else
            if (_id == +2)
            _WW[2] = (dd_flt)+1.-_WW[0] -
                                 _WW[1] ;
              
             */
                 
                 
            if (_WW[1] < (dd_flt)+0.)
                _WW[1] = (dd_flt)+0.;
            else
            if (_WW[1] > (dd_flt)+1.)
                _WW[1] = (dd_flt)+1.;
                
            if (_WW[2] < (dd_flt)+0.)
                _WW[2] = (dd_flt)+0.;
            else
            if (_WW[2] > (dd_flt)+1.)
                _WW[2] = (dd_flt)+1.;
       
            _WW[0] = (dd_flt)+1. - _WW[1] - _WW[2] ; 
                 
                  
            dd_flt  _P1[3] ;
            _P1[0] =_p1[0] ;
            _P1[1] =_p1[1] ;
            _P1[2] =_p1[2] ;

            dd_flt  _P2[3] ;
            _P2[0] =_p2[0] ;
            _P2[1] =_p2[1] ;
            _P2[2] =_p2[2] ;

            dd_flt  _P3[3] ;
            _P3[0] =_p3[0] ;
            _P3[1] =_p3[1] ;
            _P3[2] =_p3[2] ;

            dd_flt  _QQ[3] ;   
            _QQ[0] =_WW[0] * _P1[0] + 
                    _WW[1] * _P2[0] +
                    _WW[2] * _P3[0] ;
            _QQ[1] =_WW[0] * _P1[1] + 
                    _WW[1] * _P2[1] +
                    _WW[2] * _P3[1] ;
            _QQ[2] =_WW[0] * _P1[2] + 
                    _WW[1] * _P2[2] +
                    _WW[2] * _P3[2] ;
                    
            _qq[0] =_QQ[0] ;
            _qq[1] =_QQ[1] ;
            _qq[2] =_QQ[2] ;
        }
        else
        {
        __assert( false && "line_tria_3d: invalid part!") ;
        }
        
        return _hits ;
	}


    template <
    typename      poly_list
             >
    __normal_call hits_type line_poly_3d (
	__const_ptr  (double) _pa,  // line
	__const_ptr  (double) _pb,
      poly_list         & _pi,
	__write_ptr  (double) _qq   // intersection
		)
	{
        hits_type _hits = null_hits ;

    /*--------------------------------- line within poly? */

    //!! this is exact if the points are

        double _sj = +0.;

        for(size_t _ii  = +0 ;
                   _ii != _pi.count() ;
                 ++_ii  )
        {
            size_t _i1 = 
            (_ii + 0) % _pi.count() ;
            size_t _i2 = 
            (_ii + 1) % _pi.count() ;

            double _si = 
            geompred::orient3d (
           (double*)&_pi[_i1][0] , 
           (double*)&_pi[_i2][0] ,
           (double*) _pa, 
           (double*) _pb)    ;

            if (_si*_sj < +0.) 
            return null_hits ;

            if (_si    == +0.)
            _hits= edge_hits ;

            if (_si    != +0.)
            _sj = _si   ;
        }

    /*--------------------------------- poly within line? */

        size_t _nn = +0;

        _qq[0] = +0. ;
        _qq[1] = +0. ;
        _qq[2] = +0. ;

        for(size_t _ii  = +1 ;
                   _ii != _pi.count() ;
                 ++_ii  )
        {
            size_t _i1 = 
            (_ii + 0) % _pi.count() ;
            size_t _i2 = 
            (_ii + 1) % _pi.count() ;

            double _qi[ +3];
            if (line_tria_3d (
	            _pa, _pb,
	           &_pi[ +0][0],  // tria
	           &_pi[_i1][0],
	           &_pi[_i2][0],
	            _qi, true, +1) 
	            != null_hits )
	        {
	            _hits = face_hits ;
	            
	            _nn    += +1 ;
	            
	            _qq[0] += _qi[0];
	            _qq[1] += _qi[1];
	            _qq[2] += _qi[2];
	        }

        }
        
        _qq[0] /=  _nn ;
        _qq[1] /=  _nn ;
        _qq[2] /=  _nn ;
        
        return ( _hits ) ;
    }

}

#endif//__INTERSECT__


