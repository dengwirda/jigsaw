
/*
------------------------------------------------------------
 * LAG-TRI-EUCLIDEAN-4: kernel for 4-dim. laguerre tri.
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
 * Disclaimer:  Neither I nor THE CONTRIBUTORS warrant
 * this code in any way whatsoever.  This code is
 * provided "as-is" to be used at your own risk.
 *
 * THE CONTRIBUTORS include:
 * (a) The University of Sydney
 * (b) The Massachusetts Institute of Technology
 * (c) Columbia University
 * (d) The National Aeronautics & Space Administration
 * (e) Los Alamos National Laboratory
 *
------------------------------------------------------------
 *
 * Last updated: 30 Apr., 2020
 *
 * Copyright 2013-2020
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __LAGUERRE_TRI_EUCLIDEAN_4__
#   define __LAGUERRE_TRI_EUCLIDEAN_4__

    namespace mesh {

/*
------------------------------------------------------------
 * LAGUERRE-TRI-EUCLIDEAN-4: for lag.-tri. in R^4.
------------------------------------------------------------
 * IPTR-TYPE - signed-integer typedef.
 * REAL-TYPE - floating-point typedef.
------------------------------------------------------------
 */

    template <
    typename I,
    typename R
             >
    class laguerre_tri_euclidean_4
    {
/*----------- predicate for laguerre triangulation in R^4 */
    public  :
    typedef R               real_type ;
    typedef I               iptr_type ;

    iptr_type static constexpr geom_dims = +4 ;
    iptr_type static constexpr real_dims = +5 ;
    iptr_type static constexpr topo_dims = +4 ;

    public  :
/*----------------------- (squared) node-to-node distance */
    __static_call
    __inline_call real_type lensqr_kd (
    __const_ptr ( real_type) _ipos,
    __const_ptr ( real_type) _jpos
        )
    {   return  ( geometry::lensqr_4d (
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
            iptr_type _fnod[5] = {0};

            for(_fpos = 5; _fpos-- != 0; )
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
            face_node(_fnod, _fadj, 4, 3) ;
            /*
            _fnod[0] = _mesh.
             tria(_tpos)->node(_fnod[0]);
            _fnod[1] = _mesh.
             tria(_tpos)->node(_fnod[1]);
            _fnod[2] = _mesh.
             tria(_tpos)->node(_fnod[2]);
            _fnod[3] = _mesh.
             tria(_tpos)->node(_fnod[3]);
             */
            _fnod[4] = _mesh.
             tria(_tadj)->node(_fnod[4]);

            iptr_type _apex =  _fnod[4] ;

            real_type _dsqr = lensqr_kd (
                _ppos,
           &_mesh.node(_fnod[4])->pval(0));

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
            double     _xpos[4];
            _xpos[0] = _ppos[0];
            _xpos[1] = _ppos[1];
            _xpos[2] = _ppos[2];
            _xpos[3] = _ppos[3];

            for(_fpos = 5; _fpos-- != 0; )
            {
        /*--------------- test orientation wrt. k-th face */
            iptr_type  _fnod[ 5] ;
            mesh_type::tria_type::
            face_node(_fnod, _fpos, 4, 3) ;
            _fnod[0] = _mesh.
             tria(_tpos)->node(_fnod[0]);
            _fnod[1] = _mesh.
             tria(_tpos)->node(_fnod[1]);
            _fnod[2] = _mesh.
             tria(_tpos)->node(_fnod[2]);
            _fnod[3] = _mesh.
             tria(_tpos)->node(_fnod[3]);

            double _ipos[4] = {
            _mesh.node(_fnod[0])->pval(0) ,
            _mesh.node(_fnod[0])->pval(1) ,
            _mesh.node(_fnod[0])->pval(2) ,
            _mesh.node(_fnod[0])->pval(3)
                } ;
            double _jpos[4] = {
            _mesh.node(_fnod[1])->pval(0) ,
            _mesh.node(_fnod[1])->pval(1) ,
            _mesh.node(_fnod[1])->pval(2) ,
            _mesh.node(_fnod[1])->pval(3)
                } ;
            double _kpos[4] = {
            _mesh.node(_fnod[2])->pval(0) ,
            _mesh.node(_fnod[2])->pval(1) ,
            _mesh.node(_fnod[2])->pval(2) ,
            _mesh.node(_fnod[2])->pval(3)
                } ;
            double _lpos[4] = {
            _mesh.node(_fnod[3])->pval(0) ,
            _mesh.node(_fnod[3])->pval(1) ,
            _mesh.node(_fnod[3])->pval(2) ,
            _mesh.node(_fnod[3])->pval(3)
                } ;

            double _sign;
            _sign = geompred::orient4d (
                &_ipos[ 0] ,
                &_jpos[ 0] ,
                &_kpos[ 0] ,
                &_lpos[ 0] ,
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
            iptr_type _tnod[5] = {
            _mesh.tria(_tpos   )->node(0) ,
            _mesh.tria(_tpos   )->node(1) ,
            _mesh.tria(_tpos   )->node(2) ,
            _mesh.tria(_tpos   )->node(3) ,
            _mesh.tria(_tpos   )->node(4)
                } ;

            double _xpos[4] ;
            _xpos[0] = this->_ppos[0] ;
            _xpos[1] = this->_ppos[1] ;
            _xpos[2] = this->_ppos[2] ;
            _xpos[3] = this->_ppos[3] ;

            double _ipos[4] = {
            _mesh.node(_tnod[0])->pval(0) ,
            _mesh.node(_tnod[0])->pval(1) ,
            _mesh.node(_tnod[0])->pval(2) ,
            _mesh.node(_tnod[0])->pval(3)
                } ;
            double _jpos[4] = {
            _mesh.node(_tnod[1])->pval(0) ,
            _mesh.node(_tnod[1])->pval(1) ,
            _mesh.node(_tnod[1])->pval(2) ,
            _mesh.node(_tnod[1])->pval(3)
                } ;
            double _kpos[4] = {
            _mesh.node(_tnod[2])->pval(0) ,
            _mesh.node(_tnod[2])->pval(1) ,
            _mesh.node(_tnod[2])->pval(2) ,
            _mesh.node(_tnod[2])->pval(3)
                } ;
            double _lpos[4] = {
            _mesh.node(_tnod[3])->pval(0) ,
            _mesh.node(_tnod[3])->pval(1) ,
            _mesh.node(_tnod[3])->pval(2) ,
            _mesh.node(_tnod[3])->pval(3)
                } ;
            double _mpos[4] = {
            _mesh.node(_tnod[4])->pval(0) ,
            _mesh.node(_tnod[4])->pval(1) ,
            _mesh.node(_tnod[4])->pval(2) ,
            _mesh.node(_tnod[4])->pval(3)
                } ;

            __unreferenced(_xpos);
            __unreferenced(_ipos);
            __unreferenced(_jpos);
            __unreferenced(_kpos);
            __unreferenced(_lpos);
            __unreferenced(_mpos);

            double _sign = 0.;
            /*
            _sign = geompred::inball4w  (
                &_ipos[ 0] ,
                &_jpos[ 0] ,
                &_kpos[ 0] ,
                &_lpos[ 0] ,
                &_mpos[ 0] ,
                &_xpos[ 0] ) ;
            */
            return (_sign <= (double)0.0) ;
        }
        } ;

    class face_ptrs
        {
/*------------------------------ face adjacency structure */
        public :
            containers::fixed_array <
            iptr_type, +4 > _node ;

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
                    _fj._node[2] &&
                    _fi._node[3] ==
                    _fj._node[3]  ;
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
            +4* sizeof(iptr_type)
              / sizeof( uint32_t),13);
        }
        } ;

    } ;


}

#   endif  //__LAGUERRE_TRI_EUCLIDEAN_4__



