
    /*
    --------------------------------------------------------
     * STAR-TRIA-VOID: triangulate a star-shaped cavity.
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
     * STAR-TRIA-VOID: re-tria. star-shaped cavity.
    --------------------------------------------------------
     */
    
    template <
    typename      list_type 
			 >
	__normal_call void_type star_tria_void (
		iptr_list &_tset ,
		iptr_type  _npos ,
		char_type  _flag ,
	    list_type *_tnew = nullptr ,
        list_type *_told = nullptr
		)
    {
    /*---------------------- mark list of trias in cavity */
        for(auto _iter  = _tset.head(); 
                 _iter != _tset.tend(); 
               ++_iter  )
		{
			__mark(*_iter) += _flag ;
		}

    /*---------------------- scan list of trias in cavity */
        iptr_type static constexpr
            FACE_SIZE = tria_pred::_dims+1;
        
        this->_fset._lptr.set_count (
            _tset.count()*tria_dims ,
        containers::loose_alloc, nullptr) ;

        for(auto _iter  = _tset.head() ; 
                    _iter != _tset.tend() ; 
                        ++_iter  )
		{
            if (__mark(*_iter)<_flag) continue ;

        /*-------------------------- scan for cavity face */
			for(iptr_type _fpos = tria_dims +1 ; 
                          _fpos-- != +0 ; )
			{
				iptr_type  _tadj, _fadj, _fmrk ;
				find_pair(*_iter, _tadj, 
					_fpos, _fadj, _fmrk) ;

            /*--------------- continue if cavity face */
				if (_tadj == this->null_flag() || 
					__mark(_tadj) <    _flag )
				{

                iptr_type _tria = _get_tria () ;

            /*--------------- push new tria onto list */
				if (_tnew != nullptr)
                {
				    _tnew->push_tail(_tria) ;
                }
            /*------------------------- face indexing */
                iptr_type _tnod[FACE_SIZE];
                tria_type::face_node(
                    _tnod,_fpos ,
                tria_dims-0, tria_dims-1) ;        
                iptr_type  _inod;
                for(_inod = tria_dims + 0 ; 
                    _inod-- != +0 ; )
                {
                    _tnod[ _inod] = 
                        tria( *_iter)->
                            node(_tnod[ _inod]);
                }
				_tnod[ tria_dims] = _npos ;

            /*------------------------- node indexing */
                for(_inod = tria_dims + 1 ;
				    _inod-- != +0 ; )
				{
					node(
                   _tnod[_inod])->next()=_tria ;
				}
                for(_inod = tria_dims + 1 ; 
                    _inod-- != +0 ; )
                {
					tria(_tria)->
                    node(_inod) =  _tnod[_inod];
                }
                
            /*---------------- push about cavity boundary */
				push_pair(_tria, _tadj, 0 , 
					      _fadj, _fmrk) ;

            /*---------------- face indexing for new tria */
				for(iptr_type _posf = tria_dims + 1; 
                              _posf-- != +1 ; )
				{
                    iptr_type _fnod[FACE_SIZE];
                    tria_type::face_node(
                        _fnod,_posf ,
                    tria_dims-0, tria_dims-1) ;
                    
                    face_ptrs _fdat, _same ;
                    for(_inod = tria_dims + 0 ; 
                        _inod-- != +0 ; )
                    {
                        _fdat._node[_inod] = 
                            tria(_tria)-> 
                                node(_fnod[_inod]) ;
                    }
                    _fdat._tadj = _tria ;
                    _fdat._fadj = _posf ;

                    algorithms::isort (
                        _fdat._node.head(),
                            _fdat._node.tend(),
                               std::less<iptr_type>()) ;

                    if(!this->_fset._pop(_fdat,_same))
                    {
                        this->_fset.push(_fdat) ;
                    }
                    else
                    {
                        push_pair(_fdat._tadj ,
                            _same._tadj, 
                            _fdat._fadj, 
					        _same._fadj, +1 ) ;
                    }
				}
				}
			}
		}

    /*------------------------------ un-mark cavity trias */
		for(auto _iter  = _tset.head() ; 
                 _iter != _tset.tend() ; 
               ++_iter  )
		{
            if (__mark(*_iter) >= _flag)
            {
	            __mark(*_iter) -= _flag;
    
                if (_told != nullptr)
                {
            /*------------------ push "old" tria list */
			        _told->push_tail(*_iter) ;
                }
                else
                {
            /*------------------ recycle cavity trias */
                    _put_tria(*_iter)  ;
                }
            }
		}

    /*------------------------------ mesh is now delaunay */
	}
	
	
	
