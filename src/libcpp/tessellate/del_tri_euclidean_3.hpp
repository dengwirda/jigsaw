
/* 
------------------------------------------------------------
 * DEL-TRI-EUCLIDEAN-3: kernel for 3-dim. delaunay tria.
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
 * Last updated: 17 March, 2018
 *
 * Copyright 2013-2018
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __DELAUNAY_TRI_EUCLIDEAN_3__
#   define __DELAUNAY_TRI_EUCLIDEAN_3__

    namespace mesh {

/*
------------------------------------------------------------
 * DELAUNAY-TRI-EUCLIDEAN-3: for del.-tri. in R^3.
------------------------------------------------------------
 * IPTR-TYPE - signed-integer typedef.
 * REAL-TYPE - floating-point typedef.
------------------------------------------------------------
 */
     
    template <
    typename I,
    typename R
             >
    class delaunay_tri_euclidean_3
    {
/*----------- predicate for delaunay triangulation in R^2 */
    public  :
    typedef R               real_type ;
    typedef I               iptr_type ;

    iptr_type static constexpr _dims = +3 ;

    public  :
/*----------------------- initialise geometric predicates */
    __static_call
    __inline_call void_type exactinit (
        )
    {   geompred::exactinit() ;
    }

/*----------------------- (squared) node-to-node distance */
    __static_call
    __inline_call real_type lensqr_kd (
    __const_ptr ( real_type) _ipos,
    __const_ptr ( real_type) _jpos
        ) 
    {   return  ( geometry::lensqr_3d (
                    _ipos, 
                    _jpos) ) ; 
    }
    
    template <typename mesh_type> class walk_pred
        { 
/*------------------------------------ walk--simplex test */
        public  :
        __inline_call bool_type operator() (
            mesh_type    &_mesh,
        __const_ptr ( real_type) _ppos,
            iptr_type     _tpos,
            iptr_type    &_fpos
            ) const
        {   
            double     _xpos[3];
            _xpos[0] = _ppos[0];
            _xpos[1] = _ppos[1];
            _xpos[2] = _ppos[2];
        
            for(_fpos = 4; _fpos-- != 0; )
            {
        /*--------------- test orientation wrt. k-th face */            
            iptr_type  _fnod[4];
            mesh_type::tria_type::
            face_node(_fnod, _fpos, 3, 2) ;
            _fnod[0] = _mesh.
             tria(_tpos)->node(_fnod[0]);
            _fnod[1] = _mesh.
             tria(_tpos)->node(_fnod[1]);
            _fnod[2] = _mesh.
             tria(_tpos)->node(_fnod[2]);
             
            double _ipos[3] = {
            _mesh.node(_fnod[0])->pval(0) ,
            _mesh.node(_fnod[0])->pval(1) ,
            _mesh.node(_fnod[0])->pval(2)
                } ;
            double _jpos[3] = {
            _mesh.node(_fnod[1])->pval(0) ,
            _mesh.node(_fnod[1])->pval(1) ,
            _mesh.node(_fnod[1])->pval(2)
                } ;
            double _kpos[3] = {
            _mesh.node(_fnod[2])->pval(0) ,
            _mesh.node(_fnod[2])->pval(1) ,
            _mesh.node(_fnod[2])->pval(2)
                } ;
             
            double _sign; 
            _sign = geompred::orient3d (
                &_ipos[ 0] ,
                &_jpos[ 0] ,
                &_kpos[ 0] ,
                &_xpos[ 0] ) ;
             
            if (_sign > (double) +0.00 )
                return false ;            
            }
            
            return  true ;
        }
        } ;
        
    template <typename mesh_type> class circ_pred
        { 
/*------------------------------------ in-circumball test */
        public  :
        __write_ptr ( real_type)    _ppos ;

        public  :
        __inline_call circ_pred (
            real_type *_psrc
            ) : _ppos( _psrc){}
            
        __inline_call bool_type operator()(
            mesh_type &_mesh,
            iptr_type  _tpos,
            iptr_type//_fpos
            ) const
        {
            iptr_type _tnod[4] = {
            _mesh.tria(_tpos   )->node(0) ,
            _mesh.tria(_tpos   )->node(1) ,
            _mesh.tria(_tpos   )->node(2) ,
            _mesh.tria(_tpos   )->node(3)
                } ;
                
            double _xpos[3] ;
            _xpos[0] = this->_ppos[0] ;
            _xpos[1] = this->_ppos[1] ;
            _xpos[2] = this->_ppos[2] ;
        
            double _ipos[3] = {
            _mesh.node(_tnod[0])->pval(0) ,
            _mesh.node(_tnod[0])->pval(1) ,
            _mesh.node(_tnod[0])->pval(2)
                } ;
            double _jpos[3] = {
            _mesh.node(_tnod[1])->pval(0) ,
            _mesh.node(_tnod[1])->pval(1) ,
            _mesh.node(_tnod[1])->pval(2)
                } ;
            double _kpos[3] = {
            _mesh.node(_tnod[2])->pval(0) ,
            _mesh.node(_tnod[2])->pval(1) ,
            _mesh.node(_tnod[2])->pval(2)
                } ;
            double _lpos[3] = {
            _mesh.node(_tnod[3])->pval(0) ,
            _mesh.node(_tnod[3])->pval(1) ,
            _mesh.node(_tnod[3])->pval(2)
                } ;
        
            double _sign; 
            _sign = geompred::insphere  (
                &_ipos[ 0] ,
                &_jpos[ 0] ,
                &_kpos[ 0] ,
                &_lpos[ 0] ,
                &_xpos[ 0] ) ;
                
            return (_sign < (double)+0.0) ;
        }
        } ; 

    class face_ptrs
        {
/*------------------------------ face adjacency structure */
        public :
            containers::fixed_array <
            iptr_type, +3 > _node ;
            
            iptr_type       _tadj ;
            iptr_type       _fadj ;
        } ;
        
    class face_same
        {
/*------------------------------ face indexing comparison */
        public  :
        __inline_call bool_type operator() (
            face_ptrs const&_fi,
            face_ptrs const&_fj
            ) const
        {
            return  _fi._node[0] ==
                    _fj._node[0] &&
                    _fi._node[1] ==
                    _fj._node[1] &&
                    _fi._node[2] ==
                    _fj._node[2]  ;
        }
        } ;
        
    class face_hash
        {
/*------------------------------ face indexing hash-value */
        public  :
        __inline_call iptr_type operator() (
            face_ptrs const&_ff
            ) const
        {
            return hash::hashword (
            (uint32_t*)&_ff._node[ 0], 
            +3* sizeof(iptr_type) 
              / sizeof( uint32_t),13);
        }
        } ; 
    
    } ;
 
 
}
    
#   endif  //__DELAUNAY_TRI_EUCLIDEAN_3__   
    
    
    
