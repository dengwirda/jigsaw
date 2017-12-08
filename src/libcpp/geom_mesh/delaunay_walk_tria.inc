
    /*
    --------------------------------------------------------
     * WALK-TRIA-LIST: breadth-first-search about tria.
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
     * WALK-TRIA-LIST: breadth-first-search about tria.
    --------------------------------------------------------
     */
     
	template <
    typename      func_type
			 >
	__normal_call void_type walk_tria_list (
		iptr_type  _elem,
		fptr_type  _flag,
	    func_type  _pred,
		iptr_list &_list
		)
	{
        _list.push_tail(_elem);
        
		__mark(_elem)+= _flag ;
		
	/*--------- traversal about _elem while _pred is true */
		for (typename iptr_list::
                size_type _iter = +0 ; 
                    _iter != _list.count(); ++_iter)
		{
		/*--------- scan all faces within current element */
			for(iptr_type _fpos = tria_pred::_dims + 1; 
                          _fpos-- != +0 ; )
			{
			/*--------- find neighbour about current edge */
				iptr_type  _tadj, _fadj , _fmrk;
				find_pair( _list[ _iter], _tadj,
					_fpos, _fadj, _fmrk);

            /*---------------------- skip if null or seen */
				if (_tadj == this->null_flag())
					continue ;
				if (__mark(_tadj) >=   _flag )
					continue ;

			/*--------- push neighbour when _pred is true */
				if (_pred(*this, _tadj, _fadj))
				{
					__mark(_tadj)+= _flag ;
					_list.push_tail(_tadj);
				}
			}
		}
		
    /*------------------------------ unmark list of trias */
		for (typename iptr_list::
                size_type _iter = +0 ; 
                    _iter != _list.count(); ++_iter)
		{
            __mark( _list[_iter]) -=_flag ;
        }
	}



