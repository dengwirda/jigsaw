
/*
------------------------------------------------------------
 * AABB-MESH: setup AABB-tree for k-dim. faces.
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
 * Last updated: 11 Sept., 2017
 *
 * Copyright 2013-2017
 * Darren Engwirda
 * d.engwirda@gmail.com
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __AABB_MESH_K__
#   define __AABB_MESH_K__

    namespace mesh {

    /*
    --------------------------------------------------------
     * AABB-MESH: setup AABB-tree for k-dim. faces.
    --------------------------------------------------------
     */

#   define REAL_TYPE typename tree_type \
                            ::item_type \
                            ::real_type

#   define IPTR_TYPE typename face_list \
                            ::data_type \
                            ::iptr_type

#   define FACE_TYPE typename face_list \
                            ::data_type

#   define TREE_ITEM typename tree_type \
                            ::item_type

    class push_aabb
        {
    /*--------------------- default aabb "push" predicate */
        public  :
        template <
            typename  face_type
             >
        __inline_call
            bool_type operator() (
            face_type const&_fdat
            ) const
        {   __unreferenced (_fdat) ;

            return ( true ) ;
        }
        } ;


    /*
    --------------------------------------------------------
     * AABB-MESH: compute AABB tree for k-face list.
    --------------------------------------------------------
     */

    template <
        typename  node_list ,
        typename  face_list ,
        typename  tree_type ,
        typename  push_pred = push_aabb ,
        typename  allocator =
                allocators::basic_alloc
             >
    __normal_call void_type aabb_mesh (
        node_list& _nset ,  // array of mesh nodes
        face_list& _fset ,  // array of mesh cells
        tree_type& _tree ,  // tree to load into
        REAL_TYPE* _btol ,  // expand tree bbox over axes
        IPTR_TYPE  _nbox =  + 32,  // cells-per-bbox lim.
        push_pred  _push =  push_aabb()  // filter scheme
        )
    {
        IPTR_TYPE constexpr _nnod =
            face_list::data_type::_dims + 1;
        IPTR_TYPE constexpr _ndim =
            node_list::data_type::_dims + 0;

        containers::block_array <
                TREE_ITEM,
                allocator> _bbox;

    /*------------------------- form aabb's for all faces */
        IPTR_TYPE _fpos  = + 0 ;
        for (auto _iter  = _fset.head() ;
                  _iter != _fset.tend() ;
                ++_iter, ++_fpos)
        {
            if ( !_push(*_iter) ) continue ;

        /*------------------------- calc. current subtree */
            REAL_TYPE _xmin[_ndim];
            REAL_TYPE _xmax[_ndim];

        /*------------------------- init. aabb at -+ inf. */
            for(IPTR_TYPE _idim = _ndim; _idim-- != +0 ; )
            {
                _xmin[_idim] =
            +std::numeric_limits<REAL_TYPE>::infinity();
                _xmax[_idim] =
            -std::numeric_limits<REAL_TYPE>::infinity();
            }

        /*------------------------- calc. face aabb geom. */
            for(IPTR_TYPE _inod = _nnod; _inod-- != +0 ; )
            for(IPTR_TYPE _idim = _ndim; _idim-- != +0 ; )
            {
                IPTR_TYPE _node = _iter->node(_inod) ;

                if (_xmin[_idim] > (REAL_TYPE)
                    _nset[_node] .pval (_idim))
                {
                    _xmin[_idim] = (REAL_TYPE)
                    _nset[_node] .pval (_idim);
                }
                if (_xmax[_idim] < (REAL_TYPE)
                    _nset[_node] .pval (_idim))
                {
                    _xmax[_idim] = (REAL_TYPE)
                    _nset[_node] .pval (_idim);
                }
            }

        /*------------------------- inflate aabb's by tol */
            for(IPTR_TYPE _idim = _ndim; _idim-- != +0 ; )
            {
                _xmin [_idim]  -=   _btol[_idim] ;
                _xmax [_idim]  +=   _btol[_idim] ;
            }

        /*------------------------- push aabb's onto list */
            _bbox.push_tail();
            for(IPTR_TYPE _idim = _ndim; _idim-- != +0 ; )
            {
                _bbox.tail()
                    ->pmin(_idim) = _xmin[_idim] ;
                _bbox.tail()
                    ->pmax(_idim) = _xmax[_idim] ;

                _bbox.tail()
                    ->ipos()      = _fpos ;
            }
        }

    /*----------------------------- form tree from aabb's */
        _tree.load (
        _bbox.head(), _bbox.tend(), _nbox ) ;
    }

#   undef  IPTR_TYPE
#   undef  REAL_TYPE

#   undef  TREE_ITEM

#   undef  FACE_TYPE


    }

#   endif //__AABB_MESH_K__



