
    /*
    --------------------------------------------------------
     * TRIA-BALL-K: various circumscribing ball calc.'s. 
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
     * Last updated: 13 November, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
#   pragma once

#   ifndef __TRIA_BALL_K__
#   define __TRIA_BALL_K__

    namespace geometry {

    /*
    --------------------------------------------------------
     * small matrix utilities. 
    --------------------------------------------------------
     */
      
#   define __ij(__ir , __ic , __nr) \
              ((__ir)+(__ic)*(__nr))

    template <
    typename      real_type, 
    typename      size_type
             >
    __inline_call real_type det_2x2 (
	    size_type  _la,
    __const_ptr  (real_type) _aa
	    )
    {   return 
        _aa[__ij(0,0,_la)] *
		_aa[__ij(1,1,_la)] -
		_aa[__ij(0,1,_la)] *
		_aa[__ij(1,0,_la)] ;
    }

    template <
    typename      real_type, 
    typename      size_type
             >
    __inline_call real_type det_3x3 (
	    size_type  _la,
    __const_ptr  (real_type) _aa
	    )
    {   return 
        _aa[__ij(0,0,_la)] * (
		_aa[__ij(1,1,_la)] * 
		_aa[__ij(2,2,_la)] - 
		_aa[__ij(1,2,_la)] * 
		_aa[__ij(2,1,_la)] ) -

		_aa[__ij(0,1,_la)] * (
		_aa[__ij(1,0,_la)] * 
		_aa[__ij(2,2,_la)] - 
		_aa[__ij(1,2,_la)] * 
		_aa[__ij(2,0,_la)] ) +

		_aa[__ij(0,2,_la)] * (
		_aa[__ij(1,0,_la)] * 
		_aa[__ij(2,1,_la)] - 
		_aa[__ij(1,1,_la)] * 
		_aa[__ij(2,0,_la)] ) ;
    }

    template <
    typename      real_type,
    typename      size_type 
             >
    __inline_call void_type inv_2x2 (
	    size_type  _la, 
    __const_ptr  (real_type) _aa,
	    size_type  _lx,
    __write_ptr  (real_type) _xx,
        real_type &_da
	    )
    {
	    _da = det_2x2(_la, _aa) ;
    
        _xx[__ij(0,0,_lx)] = 
        _aa[__ij(1,1,_la)] ;
        _xx[__ij(1,1,_lx)] = 
        _aa[__ij(0,0,_la)] ;
        _xx[__ij(0,1,_lx)] = 
       -_aa[__ij(0,1,_la)] ;
        _xx[__ij(1,0,_lx)] = 
       -_aa[__ij(1,0,_la)] ;
    }

    template <
    typename      real_type,
    typename      size_type 
             >
    __inline_call void_type inv_3x3 (
	    size_type  _la, 
    __const_ptr  (real_type) _aa,
	    size_type  _lx,
    __write_ptr  (real_type) _xx,
        real_type &_da
	    )
    {
	    _da = det_3x3(_la, _aa) ;
    
	    _xx[__ij(0,0,_lx)] =
	    _aa[__ij(2,2,_la)] * 
	    _aa[__ij(1,1,_la)] - 
	    _aa[__ij(2,1,_la)] * 
	    _aa[__ij(1,2,_la)] ;
	
	    _xx[__ij(0,1,_lx)] =
	    _aa[__ij(2,1,_la)] * 
	    _aa[__ij(0,2,_la)] - 
	    _aa[__ij(2,2,_la)] * 
	    _aa[__ij(0,1,_la)] ;

	    _xx[__ij(0,2,_lx)] =
	    _aa[__ij(1,2,_la)] * 
	    _aa[__ij(0,1,_la)] - 
	    _aa[__ij(1,1,_la)] * 
	    _aa[__ij(0,2,_la)] ;

	    _xx[__ij(1,0,_lx)] =
	    _aa[__ij(2,0,_la)] * 
	    _aa[__ij(1,2,_la)] - 
	    _aa[__ij(2,2,_la)] * 
	    _aa[__ij(1,0,_la)] ;
	
	    _xx[__ij(1,1,_lx)] =
	    _aa[__ij(2,2,_la)] * 
	    _aa[__ij(0,0,_la)] - 
	    _aa[__ij(2,0,_la)] * 
	    _aa[__ij(0,2,_la)] ;
	
	    _xx[__ij(1,2,_lx)] =
	    _aa[__ij(1,0,_la)] * 
	    _aa[__ij(0,2,_la)] - 
	    _aa[__ij(1,2,_la)] * 
	    _aa[__ij(0,0,_la)] ;

	    _xx[__ij(2,0,_lx)] =
	    _aa[__ij(2,1,_la)] * 
	    _aa[__ij(1,0,_la)] - 
	    _aa[__ij(2,0,_la)] * 
	    _aa[__ij(1,1,_la)] ;

	    _xx[__ij(2,1,_lx)] =
	    _aa[__ij(2,0,_la)] * 
	    _aa[__ij(0,1,_la)] - 
	    _aa[__ij(2,1,_la)] * 
	    _aa[__ij(0,0,_la)] ;
	
	    _xx[__ij(2,2,_lx)] =
	    _aa[__ij(1,1,_la)] * 
	    _aa[__ij(0,0,_la)] - 
	    _aa[__ij(1,0,_la)] * 
	    _aa[__ij(0,1,_la)] ;
    }
     
    /*
    --------------------------------------------------------
     * TRIA-BALL: circumscribing balls. 
    --------------------------------------------------------
     */
     
    template <
    typename      data_type
             >
	__inline_call void_type tria_norm_3d (  // forward dec.
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
	__const_ptr  (data_type) _p3,
	__write_ptr  (data_type) _nv
		 ) ;
     
    template <
    typename      real_type
             >
	__normal_call void_type tria_ball_2d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2,
	__const_ptr  (real_type) _p3
		)
	{
		real_type _xm[2*2] ;
		real_type _xi[2*2] ;
		_xm[__ij(0,0,2)] = _p2[0]-_p1[0] ;
		_xm[__ij(0,1,2)] = _p2[1]-_p1[1] ;

		_xm[__ij(1,0,2)] = _p3[0]-_p1[0] ;
		_xm[__ij(1,1,2)] = _p3[1]-_p1[1] ;

        real_type _xr[2*1] ;
        _xr[0] = (real_type)+.5 * (
        _xm[__ij(0,0,2)] * 
        _xm[__ij(0,0,2)] +
        _xm[__ij(0,1,2)] *
        _xm[__ij(0,1,2)] ) ;
        
        _xr[1] = (real_type)+.5 * (
        _xm[__ij(1,0,2)] * 
        _xm[__ij(1,0,2)] +
        _xm[__ij(1,1,2)] *
        _xm[__ij(1,1,2)] ) ;

        real_type  _dd;
		inv_2x2(2, _xm, 2, _xi, _dd) ;
        _bb[0] = (
        _xi[__ij(0,0,2)] * _xr[0] +
		_xi[__ij(0,1,2)] * _xr[1] )  ;

        _bb[1] = (
        _xi[__ij(1,0,2)] * _xr[0] +
		_xi[__ij(1,1,2)] * _xr[1] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        
        real_type _ee[2*1];
        real_type _db[2*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,2)] * _bb[0] +
        _xm[__ij(0,1,2)] * _bb[1] )  ;
        
        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,2)] * _bb[0] +
        _xm[__ij(1,1,2)] * _bb[1] )  ;
        
        _db[0] = (
        _xi[__ij(0,0,2)] * _ee[0] +
		_xi[__ij(0,1,2)] * _ee[1] )  ;
		
		_db[1] = (
        _xi[__ij(1,0,2)] * _ee[0] +
		_xi[__ij(1,1,2)] * _ee[1] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        }

        _bb[ 2]  = lensqr_2d(_bb) ;

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
	}

    template <
    typename      real_type
             >
	__normal_call void_type tria_ball_3d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2,
	__const_ptr  (real_type) _p3
		)
	{
        real_type _nv[4*1] ;
		tria_norm_3d(_p1, _p2, _p3, _nv) ;

		real_type _xm[3*3] ;
		real_type _xi[3*3] ;
		_xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
		_xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;
        
        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
		_xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;
        
        _xm[__ij(2,0,3)] = _nv[0] ;
		_xm[__ij(2,1,3)] = _nv[1] ;
        _xm[__ij(2,2,3)] = _nv[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        _xm[__ij(0,0,3)] * 
        _xm[__ij(0,0,3)] +
        _xm[__ij(0,1,3)] *
        _xm[__ij(0,1,3)] +
        _xm[__ij(0,2,3)] *
        _xm[__ij(0,2,3)] ) ;

        _xr[1] = (real_type)+.5 * (
        _xm[__ij(1,0,3)] * 
        _xm[__ij(1,0,3)] +
        _xm[__ij(1,1,3)] *
        _xm[__ij(1,1,3)] +
        _xm[__ij(1,2,3)] *
        _xm[__ij(1,2,3)] ) ;
        
        _xr[2] = (real_type)+.0 ;

        real_type  _dd;
		inv_3x3(3, _xm, 3, _xi, _dd) ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
		_xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;
		
        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
		_xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
	    _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;
        
        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;
        
        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;
        
        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
		_xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;
		
        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
		_xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
	    _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;
        
        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;     
        }

        _bb[ 3]  = lensqr_3d(_bb) ;

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;
	}

    template <
    typename      real_type
             >
	__normal_call void_type tria_ball_3d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2,
	__const_ptr  (real_type) _p3,
    __const_ptr  (real_type) _p4
		)
	{
	    real_type _xm[3*3] ;
		real_type _xi[3*3] ;
		_xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
		_xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;
        
        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
		_xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;
        
        _xm[__ij(2,0,3)] = _p4[0]-_p1[0] ;
		_xm[__ij(2,1,3)] = _p4[1]-_p1[1] ;
        _xm[__ij(2,2,3)] = _p4[2]-_p1[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        _xm[__ij(0,0,3)] * 
        _xm[__ij(0,0,3)] +
        _xm[__ij(0,1,3)] *
        _xm[__ij(0,1,3)] +
        _xm[__ij(0,2,3)] *
        _xm[__ij(0,2,3)] ) ;

        _xr[1] = (real_type)+.5 * (
        _xm[__ij(1,0,3)] * 
        _xm[__ij(1,0,3)] +
        _xm[__ij(1,1,3)] *
        _xm[__ij(1,1,3)] +
        _xm[__ij(1,2,3)] *
        _xm[__ij(1,2,3)] ) ;

        _xr[2] = (real_type)+.5 * (
        _xm[__ij(2,0,3)] * 
        _xm[__ij(2,0,3)] +
        _xm[__ij(2,1,3)] *
        _xm[__ij(2,1,3)] +
        _xm[__ij(2,2,3)] *
        _xm[__ij(2,2,3)] ) ;

    	real_type  _dd;
		inv_3x3(3, _xm, 3, _xi, _dd) ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
		_xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;
		
        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
		_xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
	    _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;
        
        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;
        
        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;
        
        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
		_xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;
		
        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
		_xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
	    _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ; 
        }

        _bb[ 3]  = lensqr_3d(_bb) ;

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;
        
	}
	
	/*
    --------------------------------------------------------
     * ORTHOBALL: (scalar) weighted orthogonal balls. 
    --------------------------------------------------------
     */
    
    template <
    typename      real_type
             >
	__normal_call void_type orthoball_2d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2
		)
	{
	    real_type _dd[3*1] ;
	    _dd[0] = _p1[0] - _p2[0] ;
	    _dd[1] = _p1[1] - _p2[1] ;
	    _dd[2] = _p1[2] - _p2[2] ;
	    
	    real_type _dp = 
	        geometry::lensqr_2d(_dd) ;

        real_type _tt = 
            (real_type).5 * 
                (_dd[2] + _dp) / _dp ;
	    
	    _bb[0] = 
	    _p1[0] -  _tt * _dd[0] ;
	    _bb[1] = 
	    _p1[1] -  _tt * _dd[1] ;
	    
	    real_type _r1 = 
	    geometry::lensqr_2d(_p1, _bb);
	    real_type _r2 = 
	    geometry::lensqr_2d(_p2, _bb);
	    
	    _r1 -= _p1[2] ;
	    _r2 -= _p2[2] ;
	    
	    _bb[2] = 
	    (_r1 + _r2) / (real_type)+2. ;
	}
	
	template <
    typename      real_type
             >
	__normal_call void_type orthoball_3d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2
		)
	{
	    real_type _dd[4*1] ;
	    _dd[0] = _p1[0] - _p2[0] ;
	    _dd[1] = _p1[1] - _p2[1] ;
	    _dd[2] = _p1[2] - _p2[2] ;
	    _dd[3] = _p1[3] - _p2[3] ;
	    
	    real_type _dp = 
	        geometry::lensqr_3d(_dd) ;

        real_type _tt = 
            (real_type).5 * 
                (_dd[3] + _dp) / _dp ;
	    
	    _bb[0] = 
	    _p1[0] -  _tt * _dd[0] ;
	    _bb[1] = 
	    _p1[1] -  _tt * _dd[1] ;
	    _bb[2] = 
	    _p1[2] -  _tt * _dd[2] ;
	    
	    real_type _r1 = 
	    geometry::lensqr_3d(_p1, _bb);
	    real_type _r2 = 
	    geometry::lensqr_3d(_p2, _bb);
	    
	    _r1 -= _p1[3] ;
	    _r2 -= _p2[3] ;
	    
	    _bb[3] = 
	    (_r1 + _r2) / (real_type)+2. ;
	}
     
    template <
    typename      real_type
             >
	__normal_call void_type orthoball_2d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2,
	__const_ptr  (real_type) _p3
		)
	{
		real_type _xm[2*2] ;
		real_type _xi[2*2] ;
		_xm[__ij(0,0,2)] = _p2[0]-_p1[0] ;
		_xm[__ij(0,1,2)] = _p2[1]-_p1[1] ;

		_xm[__ij(1,0,2)] = _p3[0]-_p1[0] ;
		_xm[__ij(1,1,2)] = _p3[1]-_p1[1] ;

        real_type _xr[2*1] ;
        _xr[0] = (real_type)+.5 * (
        _xm[__ij(0,0,2)] * 
        _xm[__ij(0,0,2)] +
        _xm[__ij(0,1,2)] *
        _xm[__ij(0,1,2)] ) ;
        
        _xr[1] = (real_type)+.5 * (
        _xm[__ij(1,0,2)] * 
        _xm[__ij(1,0,2)] +
        _xm[__ij(1,1,2)] *
        _xm[__ij(1,1,2)] ) ;
    
        real_type _w21 = _p2[2]-_p1[2] ;
        real_type _w31 = _p3[2]-_p1[2] ;
    
        _xr[0]-= (real_type)+.5 * _w21 ;
        _xr[1]-= (real_type)+.5 * _w31 ;   
        
        real_type  _dd;
		inv_2x2(2, _xm, 2, _xi, _dd) ;
        _bb[0] = (
        _xi[__ij(0,0,2)] * _xr[0] +
		_xi[__ij(0,1,2)] * _xr[1] )  ;

        _bb[1] = (
        _xi[__ij(1,0,2)] * _xr[0] +
		_xi[__ij(1,1,2)] * _xr[1] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        
        real_type _ee[2*1];
        real_type _db[2*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,2)] * _bb[0] +
        _xm[__ij(0,1,2)] * _bb[1] )  ;
        
        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,2)] * _bb[0] +
        _xm[__ij(1,1,2)] * _bb[1] )  ;
        
        _db[0] = (
        _xi[__ij(0,0,2)] * _ee[0] +
		_xi[__ij(0,1,2)] * _ee[1] )  ;
		
		_db[1] = (
        _xi[__ij(1,0,2)] * _ee[0] +
		_xi[__ij(1,1,2)] * _ee[1] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        }

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
	
	    real_type _r1 = 
	    geometry::lensqr_2d(_p1, _bb) ;
	    real_type _r2 = 
	    geometry::lensqr_2d(_p2, _bb) ;
	    real_type _r3 = 
	    geometry::lensqr_2d(_p3, _bb) ;
	
	    _r1     -= _p1[2]  ;
	    _r2     -= _p2[2]  ;
	    _r3     -= _p3[2]  ;
	
	    _bb[2] = (_r1+_r2+_r3) / (real_type)3. ;
	}

    template <
    typename      real_type
             >
	__normal_call void_type orthoball_3d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2,
	__const_ptr  (real_type) _p3
		)
	{
        real_type _nv[4*1] ;
		tria_norm_3d(_p1, _p2, _p3, _nv) ;

		real_type _xm[3*3] ;
		real_type _xi[3*3] ;
		_xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
		_xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;
        
        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
		_xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;
        
        _xm[__ij(2,0,3)] = _nv[0] ;
		_xm[__ij(2,1,3)] = _nv[1] ;
        _xm[__ij(2,2,3)] = _nv[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        _xm[__ij(0,0,3)] * 
        _xm[__ij(0,0,3)] +
        _xm[__ij(0,1,3)] *
        _xm[__ij(0,1,3)] +
        _xm[__ij(0,2,3)] *
        _xm[__ij(0,2,3)] ) ;

        _xr[1] = (real_type)+.5 * (
        _xm[__ij(1,0,3)] * 
        _xm[__ij(1,0,3)] +
        _xm[__ij(1,1,3)] *
        _xm[__ij(1,1,3)] +
        _xm[__ij(1,2,3)] *
        _xm[__ij(1,2,3)] ) ;
        
        real_type _w21 = _p2[3]-_p1[3] ;
        real_type _w31 = _p3[3]-_p1[3] ;
    
        _xr[0]-= (real_type)+.5 * _w21 ;
        _xr[1]-= (real_type)+.5 * _w31 ;

        _xr[2] = (real_type)+.0 ;

        real_type  _dd;
		inv_3x3(3, _xm, 3, _xi, _dd) ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
		_xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;
		
        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
		_xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
	    _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;

        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;
        
        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;
        
        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;
        
        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
		_xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;
		
        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
		_xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
	    _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;
        
        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;
        }

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;
	
	    real_type _r1 = 
	    geometry::lensqr_3d(_p1, _bb) ;
	    real_type _r2 = 
	    geometry::lensqr_3d(_p2, _bb) ;
	    real_type _r3 = 
	    geometry::lensqr_3d(_p3, _bb) ;
	
	    _r1     -= _p1[3]  ;
	    _r2     -= _p2[3]  ;
	    _r3     -= _p3[3]  ;
	
	    _bb[3] = (_r1+_r2+_r3) / (real_type)3. ;
	}

    template <
    typename      real_type
             >
	__normal_call void_type orthoball_3d (
	__write_ptr  (real_type) _bb,
	__const_ptr  (real_type) _p1,
	__const_ptr  (real_type) _p2,
	__const_ptr  (real_type) _p3,
    __const_ptr  (real_type) _p4
		)
	{
	    real_type _xm[3*3] ;
		real_type _xi[3*3] ;
		_xm[__ij(0,0,3)] = _p2[0]-_p1[0] ;
		_xm[__ij(0,1,3)] = _p2[1]-_p1[1] ;
        _xm[__ij(0,2,3)] = _p2[2]-_p1[2] ;
        
        _xm[__ij(1,0,3)] = _p3[0]-_p1[0] ;
		_xm[__ij(1,1,3)] = _p3[1]-_p1[1] ;
        _xm[__ij(1,2,3)] = _p3[2]-_p1[2] ;
        
        _xm[__ij(2,0,3)] = _p4[0]-_p1[0] ;
		_xm[__ij(2,1,3)] = _p4[1]-_p1[1] ;
        _xm[__ij(2,2,3)] = _p4[2]-_p1[2] ;

        real_type _xr[3*1] ;
        _xr[0] = (real_type)+.5 * (
        _xm[__ij(0,0,3)] * 
        _xm[__ij(0,0,3)] +
        _xm[__ij(0,1,3)] *
        _xm[__ij(0,1,3)] +
        _xm[__ij(0,2,3)] *
        _xm[__ij(0,2,3)] ) ;

        _xr[1] = (real_type)+.5 * (
        _xm[__ij(1,0,3)] * 
        _xm[__ij(1,0,3)] +
        _xm[__ij(1,1,3)] *
        _xm[__ij(1,1,3)] +
        _xm[__ij(1,2,3)] *
        _xm[__ij(1,2,3)] ) ;

        _xr[2] = (real_type)+.5 * (
        _xm[__ij(2,0,3)] * 
        _xm[__ij(2,0,3)] +
        _xm[__ij(2,1,3)] *
        _xm[__ij(2,1,3)] +
        _xm[__ij(2,2,3)] *
        _xm[__ij(2,2,3)] ) ;

    	real_type  _dd;
		inv_3x3(3, _xm, 3, _xi, _dd) ;
        _bb[0] = (
        _xi[__ij(0,0,3)] * _xr[0] +
		_xi[__ij(0,1,3)] * _xr[1] +
        _xi[__ij(0,2,3)] * _xr[2] )  ;
		
        _bb[1] = (
        _xi[__ij(1,0,3)] * _xr[0] +
		_xi[__ij(1,1,3)] * _xr[1] +
        _xi[__ij(1,2,3)] * _xr[2] )  ;

        _bb[2] = (
        _xi[__ij(2,0,3)] * _xr[0] +
	    _xi[__ij(2,1,3)] * _xr[1] +
        _xi[__ij(2,2,3)] * _xr[2] )  ;

        _bb[0] /= _dd ;
        _bb[1] /= _dd ;
        _bb[2] /= _dd ;


        real_type _ee[3*1];
        real_type _db[3*1];
        for(int _ii = +1; _ii-- != +0; )
        {
        _ee[0] = _xr[0] - (
        _xm[__ij(0,0,3)] * _bb[0] +
        _xm[__ij(0,1,3)] * _bb[1] +
        _xm[__ij(0,2,3)] * _bb[2] )  ;
        
        _ee[1] = _xr[1] - (
        _xm[__ij(1,0,3)] * _bb[0] +
        _xm[__ij(1,1,3)] * _bb[1] +
        _xm[__ij(1,2,3)] * _bb[2] )  ;
        
        _ee[2] = _xr[2] - (
        _xm[__ij(2,0,3)] * _bb[0] +
        _xm[__ij(2,1,3)] * _bb[1] +
        _xm[__ij(2,2,3)] * _bb[2] )  ;
        
        _db[0] = (
        _xi[__ij(0,0,3)] * _ee[0] +
		_xi[__ij(0,1,3)] * _ee[1] +
        _xi[__ij(0,2,3)] * _ee[2] )  ;
		
        _db[1] = (
        _xi[__ij(1,0,3)] * _ee[0] +
		_xi[__ij(1,1,3)] * _ee[1] +
        _xi[__ij(1,2,3)] * _ee[2] )  ;

        _db[2] = (
        _xi[__ij(2,0,3)] * _ee[0] +
	    _xi[__ij(2,1,3)] * _ee[1] +
        _xi[__ij(2,2,3)] * _ee[2] )  ;

        _bb[0]+= _db[0] / _dd ;
        _bb[1]+= _db[1] / _dd ;
        _bb[2]+= _db[2] / _dd ;  
        }

        _bb[ 0] += _p1[0]  ;
        _bb[ 1] += _p1[1]  ;
        _bb[ 2] += _p1[2]  ;
	
	    real_type _r1 = 
	    geometry::lensqr_3d(_p1, _bb) ;
	    real_type _r2 = 
	    geometry::lensqr_3d(_p2, _bb) ;
	    real_type _r3 = 
	    geometry::lensqr_3d(_p3, _bb) ;
	    real_type _r4 = 
	    geometry::lensqr_3d(_p4, _bb) ;
	
	    _r1     -= _p1[3]  ;
	    _r2     -= _p2[3]  ;
	    _r3     -= _p3[3]  ;
	    _r4     -= _p4[3]  ;
	
	    _bb[3] = (_r1+_r2+_r3+_r4)/(real_type)4. ;
	}

    #undef  __ij
    
    
    /*
    --------------------------------------------------------
     * MINI-BALL: min. enclosing balls. 
    --------------------------------------------------------
     */
    
    #define __miniface12(__cc, __pi, __pj, __pq) \
        do {                                \
    /*---------- calc. face circumball */   \
        mini_ball_2d(__bb, __pi, __pj) ;    \
    /*---------- keep if min-enclosing */   \
        if (__bb[2] >= lensqr_2d(__bb, __pq))    \
        {                                   \
            if (__bb[2] < __cc[2])          \
            {                               \
                __cc[0] = __bb[0];          \
                __cc[1] = __bb[1];          \
                __cc[2] = __bb[2];          \
            }                               \
        }                                   \
        } while (false) ;

    #define __miniface13(__cc, __pi, __pj, __pq) \
        do {                                \
    /*---------- calc. face circumball */   \
        mini_ball_3d(__bb, __pi, __pj) ;    \
    /*---------- keep if min-enclosing */   \
        if (__bb[3] >= lensqr_3d(__bb, __pq))    \
        {                                   \
            if (__bb[3] < __cc[3])          \
            {                               \
                __cc[0] = __bb[0];          \
                __cc[1] = __bb[1];          \
                __cc[2] = __bb[2];          \
                __cc[3] = __bb[3];          \
            }                               \
        }                                   \
        } while (false) ;

    #define __miniface23(__cc, __pi, __pj, __pk, __pq) \
        do {                                \
    /*---------- calc. face circumball */   \
        mini_ball_3d(__bb, __pi, __pj, __pk);    \
    /*---------- keep if min-enclosing */   \
        if (__bb[3] >= lensqr_3d(__bb, __pq))    \
        {                                   \
            if (__bb[3] < __cc[3])          \
            {                               \
                __cc[0] = __bb[0];          \
                __cc[1] = __bb[1];          \
                __cc[2] = __bb[2];          \
                __cc[3] = __bb[3];          \
            }                               \
        }                                   \
        } while (false) ;

/*------------------------------------ min-ball 1dim face */
    template <
    typename      data_type
             >
    __normal_call void_type mini_ball_2d (
    __write_ptr  (data_type) _cc,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2
        )
    {
        _cc[0] = 
        (data_type).5*(_p1[0]+_p2[0]);
        _cc[1] = 
        (data_type).5*(_p1[1]+_p2[1]);

        _cc[2] = (data_type).0 ;
		_cc[2] = std::max (
		_cc[2] , lensqr_2d(_cc, _p1));
		_cc[2] = std::max (
		_cc[2] , lensqr_2d(_cc, _p2));
    }

    template <
    typename      data_type
             >
    __normal_call void_type mini_ball_3d (
    __write_ptr  (data_type) _cc,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2
        )
    {
        _cc[0] = 
        (data_type).5*(_p1[0]+_p2[0]);
        _cc[1] = 
        (data_type).5*(_p1[1]+_p2[1]);
        _cc[2] = 
        (data_type).5*(_p1[2]+_p2[2]);

        _cc[3] = (data_type).0 ;
		_cc[3] = std::max (
		_cc[3] , lensqr_3d(_cc, _p1));
		_cc[3] = std::max (
		_cc[3] , lensqr_3d(_cc, _p2));
    }

/*------------------------------------ min-ball 2dim face */
    template <
    typename      data_type
             >
    __normal_call void_type mini_ball_2d (
    __write_ptr  (data_type) _cc,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        tria_ball_2d(_cc, _p1, _p2, _p3);

        data_type __bb[3] ;
        __miniface12(_cc, _p1, _p2, _p3);
        __miniface12(_cc, _p2, _p3, _p1);
        __miniface12(_cc, _p3, _p1, _p2);
    }

    template <
    typename      data_type
             >
    __normal_call void_type mini_ball_3d (
    __write_ptr  (data_type) _cc,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3
        )
    {
        tria_ball_3d(_cc, _p1, _p2, _p3);

        data_type __bb[4] ;
        __miniface13(_cc, _p1, _p2, _p3);
        __miniface13(_cc, _p2, _p3, _p1);
        __miniface13(_cc, _p3, _p1, _p2);
    }

/*------------------------------------ min-ball 3dim face */
    template <
    typename      data_type
             >
    __normal_call void_type mini_ball_3d (
    __write_ptr  (data_type) _cc,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
    __const_ptr  (data_type) _p3,
    __const_ptr  (data_type) _p4
        )
    {
        tria_ball_3d( _cc, 
            _p1, _p2, _p3, _p4) ;

        data_type __bb[4];
        __miniface23( _cc, 
            _p1, _p2, _p3, _p4) ;
        __miniface23( _cc, 
            _p1, _p2, _p4, _p3) ;
        __miniface23( _cc, 
            _p2, _p3, _p4, _p1) ;
        __miniface23( _cc, 
            _p3, _p1, _p4, _p2) ;
    }

    #undef  __miniface12
    #undef  __miniface13
    #undef  __miniface23
    
    /*
    --------------------------------------------------------
     * MASS-BALL: centre-of-mass balls. 
    --------------------------------------------------------
     */
    
    template <
    typename      data_type
             >
    __normal_call void_type mass_ball_2d (
    __write_ptr  (data_type) _bb,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]/= (data_type) +2. ;
        
        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]/= (data_type) +2. ;
        
        data_type _r1 = 
        geometry::lensqr_2d(_bb, _p1);
        data_type _r2 = 
        geometry::lensqr_2d(_bb, _p2);

        _bb[2] = 
        (_r1+_r2)/(data_type)+2. ;
    }
    
    template <
    typename      data_type
             >
    __normal_call void_type mass_ball_3d (
    __write_ptr  (data_type) _bb,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]/= (data_type) +2. ;
        
        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]/= (data_type) +2. ;
        
        _bb[2] = _p1[2] ;
        _bb[2]+= _p2[2] ;
        _bb[2]/= (data_type) +2. ;

        data_type _r1 = 
        geometry::lensqr_3d(_bb, _p1);
        data_type _r2 = 
        geometry::lensqr_3d(_bb, _p2);

        _bb[3] = 
        (_r1+_r2)/(data_type)+2. ;
    }
    
    template <
    typename      data_type
             >
    __normal_call void_type mass_ball_2d (
    __write_ptr  (data_type) _bb,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
	__const_ptr  (data_type) _p3
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]+= _p3[0] ;
        _bb[0]/= (data_type) +3. ;
        
        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]+= _p3[1] ;
        _bb[1]/= (data_type) +3. ;
        
        data_type _r1 = 
        geometry::lensqr_2d(_bb, _p1);
        data_type _r2 = 
        geometry::lensqr_2d(_bb, _p2);
        data_type _r3 = 
        geometry::lensqr_2d(_bb, _p3);

        _bb[2] = std::max ( _r3, 
                 std::max ( _r1, _r2)
                          ) ;
    }
    
    template <
    typename      data_type
             >
    __normal_call void_type mass_ball_3d (
    __write_ptr  (data_type) _bb,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
	__const_ptr  (data_type) _p3
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]+= _p3[0] ;
        _bb[0]/= (data_type) +3. ;
        
        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]+= _p3[1] ;
        _bb[1]/= (data_type) +3. ;
        
        _bb[2] = _p1[2] ;
        _bb[2]+= _p2[2] ;
        _bb[2]+= _p3[2] ;
        _bb[2]/= (data_type) +3. ;

        data_type _r1 = 
        geometry::lensqr_3d(_bb, _p1);
        data_type _r2 = 
        geometry::lensqr_3d(_bb, _p2);
        data_type _r3 = 
        geometry::lensqr_3d(_bb, _p3);

        _bb[3] = std::max ( _r3, 
                 std::max ( _r1, _r2)
                          ) ;
    }
    
    template <
    typename      data_type
             >
    __normal_call void_type mass_ball_3d (
    __write_ptr  (data_type) _bb,
	__const_ptr  (data_type) _p1,
	__const_ptr  (data_type) _p2,
	__const_ptr  (data_type) _p3,
	__const_ptr  (data_type) _p4
        )
    {
        _bb[0] = _p1[0] ;
        _bb[0]+= _p2[0] ;
        _bb[0]+= _p3[0] ;
        _bb[0]+= _p4[0] ;
        _bb[0]/= (data_type) +4. ;
        
        _bb[1] = _p1[1] ;
        _bb[1]+= _p2[1] ;
        _bb[1]+= _p3[1] ;
        _bb[1]+= _p4[1] ;
        _bb[1]/= (data_type) +4. ;
        
        _bb[2] = _p1[2] ;
        _bb[2]+= _p2[2] ;
        _bb[2]+= _p3[2] ;
        _bb[2]+= _p4[2] ;
        _bb[2]/= (data_type) +4. ;
    
        data_type _r1 = 
        geometry::lensqr_3d(_bb, _p1);
        data_type _r2 = 
        geometry::lensqr_3d(_bb, _p2);
        data_type _r3 = 
        geometry::lensqr_3d(_bb, _p3);
        data_type _r4 = 
        geometry::lensqr_3d(_bb, _p4);

        _bb[3] = std::max ( _r4,
                 std::max ( _r3,
                 std::max ( _r1, _r2)
                         )) ;
    }
          
    }
    
#   endif //__TRIA_BALL_K__



