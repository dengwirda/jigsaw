
    /*
    --------------------------------------------------------
     * WALK-MESH-NODE: find enclosing tria. via traversal.
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
     * Last updated: 12 May, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from delaunay_tri_k.hpp
    

    /*
    --------------------------------------------------------
     * WALK-MESH-NODE: find enclosing tria. via traversal.
    --------------------------------------------------------
     */
 
	__normal_call bool_type walk_mesh_node (
		real_type *_ppos,
		iptr_type &_elem,
		iptr_type  _hint =  null_flag ()
		)
	{
	/*------------- walk mesh until we find enclosing tri */
		bool_type _done = false;
		bool_type _scan = false;
        iptr_type _ntri = 
       (iptr_type)this->_tset.count () ;

        if (this->_nset.count() == +0 ||
            this->_tset.count() == +0 )
        {
    /*----------------------------- return when empty */
            return  _done ;
        }
        if (_hint < +0 || _hint>=_ntri)
        {
    /*----------------------------- not a valid index */
            _elem = +0 ;
            _scan =  true ;
        }
        else
        if (tria(_hint)->mark() <= -1 )
        {
    /*----------------------------- not a valid _tria */
            _elem = +0 ;
            _scan =  true ;
        }
        else
        {
            _elem = _hint ;
        }

	/*------------- randomised scan: locate initial _elem */
		if (_scan)
		{
		
	    real_type _best = 
       +std::numeric_limits<real_type>::infinity() ;
	
    /*--------------------- scale no. iter with dims. */
        iptr_type _imax = (iptr_type)
            std::pow ((double) _ntri, 
        +1.0 / tria_pred::_dims) / 4 + 1 ;

		for (auto _iter = _imax; _iter-- != +0; )
		{
        /*----------------------- randomise selection */
            iptr_type _ipos =
	            (iptr_type)std::rand() 
	                % this->_nset.count();

        /*----------------------- reject "null" nodes */
            if (node(_ipos)->mark() < +0)
            { 
                continue  ;
            }
		/*----------------------- calc. dist. to node */
			real_type _dist = 
                tria_pred::lensqr_kd(_ppos, 
                    &node(_ipos)->pval(0));

			if (_best > _dist)
			{
        /*----------------------- keep closest so far */
				_best = _dist;
				
				_elem = node(_ipos)->next() ;
			}
		}
		
		}

	/*------------- ensure that we have a non--null _elem */
		for ( ; _elem != _ntri; ++_elem)
		{
		    if (tria(_elem)->mark()>=+0) break  ;
		}

	/*------------- walk elements toward _ppos from _elem */
	
	    iptr_type _next = +0 ;
	    iptr_type _face = +0 ;
	    
		for (auto _iter = _ntri; 
		          _iter-- != +0; _elem = _next)
		{
		/*--------------- test node-in-tria condition */
            _done = typename tria_pred::
                template walk_pred<self_type>() (
                    *this, _ppos,
                    _elem, _face) ;
		
			if(!_done )
            {
        /*--------------- push next element onto walk */
			    _next = __unflip (
                    tria(_elem)->next(_face));

                if(_next == this->null_flag())
                {
        /*--------------- stop at external boundaries */
                    break ;
                }
            }
        /*--------------- stop when _ppos is enclosed */
			else break ;
		}

		return ( _done )  ;
	}
	
	
	
