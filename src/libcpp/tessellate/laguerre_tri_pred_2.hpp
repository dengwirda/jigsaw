
/*
------------------------------------------------------------
 * LAG-TRI-EUCLIDEAN-2: kernel for 2-dim. laguerre tri.
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
 * Last updated: 01 March, 2020
 *
 * Copyright 2013-2020
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __LAGUERRE_TRI_EUCLIDEAN_2__
#   define __LAGUERRE_TRI_EUCLIDEAN_2__

    namespace mesh {

/*
------------------------------------------------------------
 * LAGUERRE-TRI-EUCLIDEAN-2: for lag.-tri. in R^2.
------------------------------------------------------------
 * IPTR-TYPE - signed-integer typedef.
 * REAL-TYPE - floating-point typedef.
------------------------------------------------------------
 */

    template <
    typename I,
    typename R
             >
    class laguerre_tri_euclidean_2
    {
/*----------- predicate for laguerre triangulation in R^2 */
    public  :
    typedef R               real_type ;
    typedef I               iptr_type ;

    iptr_type static constexpr geom_dims = +2 ;
    iptr_type static constexpr real_dims = +3 ;
    iptr_type static constexpr topo_dims = +2 ;

    public  :
/*----------------------- (squared) node-to-node distance */
    __static_call
    __inline_call real_type lensqr_kd (
    __const_ptr ( real_type) _ipos,
    __const_ptr ( real_type) _jpos
        )
    {   return  ( geometry::lensqr_2d (
                    _ipos,
                    _jpos) ) ;
    }

    template <typename mesh_type> class near_pred
        {
/*------------------------------------ walk--simplex test */
        public  :
        __inline_call bool_type operator() (
            mesh_type    &_mesh,
        __const_ptr ( real_type) _ppos,
            real_type    &_dmin,
            iptr_type    &_nmin,
            iptr_type     _tpos,
            iptr_type    &_fmin
            ) const
        {
            bool_type _done = true  ;
            iptr_type _fpos;
            iptr_type _fadj;
            iptr_type _tadj;
            iptr_type _fnod[3] = {0};

            for(_fpos = 3; _fpos-- != 0; )
            {
        /*--------------- check dist. wrt. k-th face apex */
            _tadj = _mesh.
                tria(_tpos)->next(_fpos) ;
            _fadj = _mesh.
                tria(_tpos)->fpos(_fpos) ;

            _tadj = __unflip(_tadj) ;

            if(_tadj == _mesh.null_flag())
                continue ;

            mesh_type::tria_type::
            face_node(_fnod, _fadj, 2, 1) ;
            /*
            _fnod[0] = _mesh.
             tria(_tadj)->node(_fnod[0]);
            _fnod[1] = _mesh.
             tria(_tadj)->node(_fnod[1]);
             */
            _fnod[2] = _mesh.
             tria(_tadj)->node(_fnod[2]);

            iptr_type _apex =  _fnod[2] ;

            real_type _dsqr = lensqr_kd (
                _ppos,
           &_mesh.node(_fnod[2])->pval(0));

            if (_dsqr < _dmin)
            {
                _done = false;
                _dmin = _dsqr;
                _nmin = _apex;
                _fmin = _fpos;
            }
            }

            return _done ;
        }
        } ;

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
            double     _xpos[2];
            _xpos[0] = _ppos[0];
            _xpos[1] = _ppos[1];

            for(_fpos = 3; _fpos-- != 0; )
            {
        /*--------------- test orientation wrt. k-th face */
            iptr_type  _fnod[ 3] ;
            mesh_type::tria_type::
            face_node(_fnod, _fpos, 2, 1) ;
            _fnod[0] = _mesh.
             tria(_tpos)->node(_fnod[0]);
            _fnod[1] = _mesh.
             tria(_tpos)->node(_fnod[1]);

            double _ipos[2] = {
            _mesh.node(_fnod[0])->pval(0) ,
            _mesh.node(_fnod[0])->pval(1)
                } ;
            double _jpos[2] = {
            _mesh.node(_fnod[1])->pval(0) ,
            _mesh.node(_fnod[1])->pval(1)
                } ;

            double _sign;
            _sign = geompred::orient2d (
                &_ipos[ 0] ,
                &_jpos[ 0] ,
                &_xpos[ 0] ) ;

            if (_sign < (double) +0.00 )
                return false ;
            }

            return  true ;
        }
        } ;

    template <typename mesh_type> class circ_pred
        {
/*------------------------------------ in-circumball test */
        public  :
        __const_ptr ( real_type)    _ppos ;

        public  :
        __inline_call circ_pred  (
        __const_ptr ( real_type) _psrc
            ) : _ppos( _psrc ) { ; }

        __inline_call bool_type operator()(
            mesh_type &_mesh,
            iptr_type  _tpos,
            iptr_type//_fpos
            ) const
        {
            iptr_type _tnod[3] = {
            _mesh.tria(_tpos   )->node(0) ,
            _mesh.tria(_tpos   )->node(1) ,
            _mesh.tria(_tpos   )->node(2)
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

            double _sign;
            _sign = geompred::inball2w  (
                &_ipos[ 0] ,
                &_jpos[ 0] ,
                &_kpos[ 0] ,
                &_xpos[ 0] ) ;

            return (_sign >= (double)+0.) ;
        }
        } ;

    class face_ptrs
        {
/*------------------------------ face adjacency structure */
        public :
            containers::fixed_array <
            iptr_type, +2 > _node ;

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
                    _fj._node[1]  ;
        }
        } ;

    class face_hash
        {
/*------------------------------ face indexing hash-value */
        public  :
        __inline_call uint32_t  operator() (
            face_ptrs const&_ff
            ) const
        {
            return hash::hashword (
            (uint32_t*)&_ff._node[ 0],
            +2* sizeof(iptr_type)
              / sizeof( uint32_t),13);
        }
        } ;

    } ;


    }

#   endif  //__LAGUERRE_TRI_EUCLIDEAN_2__



