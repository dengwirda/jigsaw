
    /*
    --------------------------------------------------------
     * ITER-DIVS-2: optim. schemes to split edges.
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
     * Last updated: 12 August, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
    
    // from iter_mesh_2.hpp
    

    /*
    --------------------------------------------------------
     * _DIV-EDGE: try edge split to improve adj. cost-fun.
    --------------------------------------------------------
     */
 
    __static_call
    __normal_call void_type _div_edge_3 (
        geom_type &_geom ,
        mesh_type &_mesh , 
        size_type &_hfun ,
        pred_type &_pred ,
        real_list &_hval ,
        iter_opts &_opts ,
        iptr_type  _enum ,
        bool_type &_okay ,
        iptr_list &_tset ,
        real_list &_csrc ,
        real_list &_cdst ,
        real_type  _lmin
            = (real_type) +0.67E+00,
        real_type  _good 
            = (real_type) +9.25E-01,
        real_type  _qinc 
            = (real_type) +1.00E-08
        )
    {    
        __unreferenced(_pred) ;         // for MSVC...
        __unreferenced(_opts) ;

    /*--------------------------------- get edge indexing */
         auto _edge = 
        _mesh._set2.head() + _enum ;
        
        iptr_type _enod[2] ;
        _enod[ 0] = _edge->node(0) ;
        _enod[ 1] = _edge->node(1) ;
        
         auto _iptr = _mesh.
        _set1.head()+ _edge->node(0) ;        
         auto _jptr = _mesh.
        _set1.head()+ _edge->node(1) ;
        
        _tset.set_count(+0);
        _mesh.edge_tri4(_enum, _tset);
    
        _okay = false ;
    
    /*--------------------------------- get edge h-sizing */
        real_type _ipos[_dims] ;
        real_type _jpos[_dims] ;
        iptr_type _idim = +0;
        for (_idim = _dims+0; _idim-- != 0; )
        {
            _ipos[_idim] =
                _iptr->pval(_idim) ;
            _jpos[_idim] =
                _jptr->pval(_idim) ;
        }

        real_type  _isiz = 
        _hfun.eval(_ipos, _iptr->hidx ()) ;
        real_type  _jsiz = 
        _hfun.eval(_jpos, _jptr->hidx ()) ;
        
        real_type  _lsqr =
            _pred.length_sq(_ipos, _jpos) ;

        real_type  _hbar = 
            std::min(_isiz , _jsiz);

    /*--------------------------------- exit if too small */
        if (_lsqr <= _hbar * _lmin *
                     _hbar * _lmin )
            return  ;
         
    /*--------------------------------- get adjacent face */
         
         
         
    }
    
    
    
