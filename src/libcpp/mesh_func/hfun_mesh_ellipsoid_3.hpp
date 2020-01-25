
    /*
    --------------------------------------------------------
     * HFUN-MESH-ELLIPSOID-3D: ellipsoidal H(x) in R^3.
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
     * Last updated: 12 August, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_MESH_ELLIPSOID_3__
#   define __HFUN_MESH_ELLIPSOID_3__

    namespace mesh {

    template <
    typename R ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class hfun_mesh_ellipsoid_3d
        : public hfun_base_kd <I, R>
    {
    public  :

    /*------------------------- euclidean size-fun in R^2 */

    typedef R                       real_type ;
    typedef I                       iptr_type ;
    typedef A                       allocator ;

    typedef hfun_mesh_ellipsoid_3d  <
            real_type ,
            iptr_type >             hfun_type ;

    typedef typename  hfun_base_kd  <
            iptr_type ,
            real_type >::hint_type  hint_type ;

    typedef containers::array   <
            real_type ,
            allocator >             real_list ;


    typedef mesh_complex_node_3<I, R>
                                    node_type ;

    typedef mesh_complex_edge_2<I>  edge_type ;
    typedef mesh_complex_tria_3<I>  tri3_type ;

    typedef mesh::tria_complex_2<
                node_type,
                edge_type,
                tri3_type,
                allocator       >   mesh_type ;

    typedef geom_tree::aabb_node_base_k
                                    tree_node ;

    typedef geom_tree::aabb_item_rect_k <
                    float,
                iptr_type,
                + 3             >   tree_item ;

    typedef geom_tree::aabb_tree<
                tree_item,
                + 3 ,
                tree_node,
                allocator       >   tree_type ;

    typedef geom_tree::aabb_pred_node_3 <
                    float,
                iptr_type       >   tree_pred ;

    public  :

    /*--------------- (x/a)**2 + (y/b)**2 + (z/c)**2 = 1. */

    real_type                      _radA = 1. ;
    real_type                      _radB = 1. ;
    real_type                      _radC = 1. ;

    containers::
        fixed_array<real_type,3>   _bmin ;
    containers::
        fixed_array<real_type,3>   _bmax ;

    mesh_type                      _mesh ;
    tree_type                      _tree ;

    containers::array<
        real_type, allocator >     _hval ;

    containers::array<
        real_type, allocator >     _dhdx ;

    public  :

    __inline_call void_type toR3 (
      __const_ptr(real_type) _apos ,
      __write_ptr(real_type) _ppos
        ) const
    {
    /*------------ helper: convert from S^2 to R^3 coord. */
        _ppos[0] = this->_radA *
            std::cos( _apos[0] ) *
            std::cos( _apos[1] ) ;

        _ppos[1] = this->_radB *
            std::sin( _apos[0] ) *
            std::cos( _apos[1] ) ;

        _ppos[2] = this->_radC *
            std::sin( _apos[1] ) ;
    }

    __inline_call void_type toS2 (
      __const_ptr(real_type) _ppos ,
      __write_ptr(real_type) _apos
        ) const
    {
    /*------------ helper: convert from R^3 to S^2 coord. */
        real_type _xmul =
            _ppos[0] * this->_radB ;
        real_type _ymul =
            _ppos[1] * this->_radA ;
        real_type _zrat =
            _ppos[2] / this->_radC ;

        _zrat = std::min(+1.,_zrat);
        _zrat = std::max(-1.,_zrat);

        _apos[0]= std::atan2(_ymul,
                             _xmul);
        _apos[1]= std::asin (_zrat);
    }

    /*
    --------------------------------------------------------
     * INIT-HFUN: make tree - check signs.
    --------------------------------------------------------
     */

    __normal_call void_type init (
        )
    {
        class tria_pred
            {
    /*-------------------- TRUE if tree should index tria */
            public  :
            __inline_call
                bool_type operator () (
                tri3_type const& _tdat
                ) const
            {   return _tdat.mark() >= 0 ;
            }
            } ;

    /*----------------------------- init. aabb at -+ inf. */
        for(auto _idim = 3; _idim-- != 0 ; )
        {
            this->_bmin[_idim] =
                +std::numeric_limits<
                    real_type>::infinity() ;

            this->_bmax[_idim] =
                -std::numeric_limits<
                    real_type>::infinity() ;
        }

    /*----------------------------- convert to R^3 coord. */
        for (auto _iter  =
             this->_mesh._set1.head() ;
                  _iter !=
             this->_mesh._set1.tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
                real_type _apos[2];
                real_type _ppos[3];

                _apos[0] = _iter->pval(0);
                _apos[1] = _iter->pval(1);

                toR3(_apos, _ppos);

                _iter->pval(0) = _ppos[0];
                _iter->pval(1) = _ppos[1];
                _iter->pval(2) = _ppos[2];
            }
        }

    /*----------------------------- calc. aabb for inputs */
        for (auto  _iter  =
             this->_mesh._set1.head() ;
                   _iter !=
             this->_mesh._set1.tend() ;
                 ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
            this->_bmin[0] = std::min (
            this->_bmin[0] ,
            _iter->pval(0) ) ;
            this->_bmin[1] = std::min (
            this->_bmin[1] ,
            _iter->pval(1) ) ;
            this->_bmin[2] = std::min (
            this->_bmin[2] ,
            _iter->pval(2) ) ;

            this->_bmax[0] = std::max (
            this->_bmax[0] ,
            _iter->pval(0) ) ;
            this->_bmax[1] = std::max (
            this->_bmax[1] ,
            _iter->pval(1) ) ;
            this->_bmax[2] = std::max (
            this->_bmax[2] ,
            _iter->pval(2) ) ;
            }
        }

        float     static const _RTOL =
            std::pow (
            std::numeric_limits<float>
            ::epsilon(), (float) +0.9) ;

        iptr_type static
        constexpr _NBOX=(iptr_type)+8  ;

        float      _BTOL[3] ;
        _BTOL[0] =
       (float)    ( this->_bmax[ 0] -
                    this->_bmin[ 0] )
                 * _RTOL ;
        _BTOL[1] =
       (float)    ( this->_bmax[ 1] -
                    this->_bmin[ 1] )
                 * _RTOL ;
        _BTOL[2] =
       (float)    ( this->_bmax[ 2] -
                    this->_bmin[ 2] )
                 * _RTOL ;

    /*-------------------- make aabb-tree and init. bbox. */
        aabb_mesh( this->_mesh._set1 ,
                   this->_mesh._set3 ,
                   this->_tree,_BTOL ,
                  _NBOX , tria_pred()) ;
    }

    /*
    --------------------------------------------------------
     * CLIP-HFUN: impose |dh/dx| limits.
    --------------------------------------------------------
     */

    __normal_call void_type clip (
        )
    {
        class  less_than
        {
    /*-------------------- "LESS-THAN" operator for queue */
        public  :
            typename
            real_list::_write_it _hptr;

        public  :
        __inline_call less_than  (
            typename
            real_list::_write_it _hsrc
            ) : _hptr(_hsrc) {}

        __inline_call
            bool_type operator() (
            iptr_type _ipos,
            iptr_type _jpos
            )
        {   return *(this->_hptr+_ipos) <
                   *(this->_hptr+_jpos) ;
        }
        } ;

        typedef typename
            allocator:: size_type   uint_type ;

        uint_type static constexpr
            _null =
        std::numeric_limits<uint_type>::max() ;

        containers::prioritymap <
            iptr_type ,
            less_than ,
            allocator >
        _sort((less_than(this->_hval.head())));

        containers:: array      <
            typename
            allocator:: size_type,
            allocator >     _keys;

        containers:: array      <
            iptr_type ,
            allocator >     _tset;

    /*-------------------- push nodes onto priority queue */
        _keys.set_count (
            _mesh._set1.count(),
        containers::tight_alloc, _null) ;

        iptr_type _inum  = +0;
        for (auto _iter  =
             this->_mesh._set1.head() ;
                  _iter !=
             this->_mesh._set1.tend() ;
                ++_iter , ++_inum)
        {
            if (_iter->mark() >= +0 )
            {
                _keys[_inum] =
                    _sort.push(_inum) ;
            }
        }

    /*-------------------- compute h(x) via fast-marching */
        for ( ; !_sort.empty() ; )
        {
            iptr_type _base ;
            _sort._pop_root(_base) ;

            _keys[_base] = _null ;

            _tset.set_count( +0) ;
             this->
            _mesh.node_tri3(_base, _tset);

            for (auto _next  = _tset.head();
                      _next != _tset.tend();
                    ++_next )
            {
                 auto _inod = this->
                _mesh._set3[*_next].node(0);
                 auto _jnod = this->
                _mesh._set3[*_next].node(1);
                 auto _knod = this->
                _mesh._set3[*_next].node(2);

    /*-------------------- skip any cells with null nodes */
                if (_keys[_inod] == _null &&
                    _inod != _base) continue ;
                if (_keys[_jnod] == _null &&
                    _jnod != _base) continue ;
                if (_keys[_knod] == _null &&
                    _knod != _base) continue ;

                if (this->_dhdx.count() >1)
                {
    /*-------------------- update adj. set, g = g(x) case */
                if (eikonal_tria_3d (
                   &this->
                _mesh._set1[ _inod].pval(0),
                   &this->
                _mesh._set1[ _jnod].pval(0),
                   &this->
                _mesh._set1[ _knod].pval(0),
                    this->_hval[_inod],
                    this->_hval[_jnod],
                    this->_hval[_knod],
                    this->_dhdx[_inod],
                    this->_dhdx[_jnod],
                    this->_dhdx[_knod]) )
                {

                if (_keys[_inod] != _null)
                    _sort.update(
                    _keys[_inod] ,  _inod) ;

                if (_keys[_jnod] != _null)
                    _sort.update(
                    _keys[_jnod] ,  _jnod) ;

                if (_keys[_knod] != _null)
                    _sort.update(
                    _keys[_knod] ,  _knod) ;

                }
                }
                else
                if (this->_dhdx.count()==1)
                {
    /*-------------------- update adj. set, const. g case */
                if (eikonal_tria_3d (
                   &this->
                _mesh._set1[ _inod].pval(0),
                   &this->
                _mesh._set1[ _jnod].pval(0),
                   &this->
                _mesh._set1[ _knod].pval(0),
                    this->_hval[_inod],
                    this->_hval[_jnod],
                    this->_hval[_knod],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ],
                    this->_dhdx[  +0 ]) )
                {

                if (_keys[_inod] != _null)
                    _sort.update(
                    _keys[_inod] ,  _inod) ;

                if (_keys[_jnod] != _null)
                    _sort.update(
                    _keys[_jnod] ,  _jnod) ;

                if (_keys[_knod] != _null)
                    _sort.update(
                    _keys[_knod] ,  _knod) ;

                }
                }
            }
        }

    }

    /*
    --------------------------------------------------------
     * FIND-TRIA: scan for nearest tria.
    --------------------------------------------------------
     */

    class find_tria
        {
        public  :
        real_type              *_ppos ;

        mesh_type              *_mesh ;

        bool_type               _find ;
        iptr_type               _tpos ;

        public  :

    /*------------------------ make a tree-tria predicate */
        __inline_call find_tria (
            real_type*_psrc = nullptr ,
            mesh_type*_msrc = nullptr
            ) : _ppos(_psrc) ,
                _mesh(_msrc) ,
                _find(false) ,
                _tpos(   -1)   {}

    /*------------------------ call pred. on tree matches */
        __inline_call
            void_type operator () (
                typename
            tree_type::item_data  *_iptr
            )
        {
            if (this->_find) return ;

            for ( ; _iptr != nullptr;
                    _iptr = _iptr->_next)
            {
                real_type  _qtmp[+3];
                iptr_type  _TPOS =
                    _iptr->_data.ipos() ;

                if (near_pred ( _ppos ,
                        _qtmp ,*_mesh ,
                        _TPOS ) )
                {
    /*------------------------ is fully inside: finished! */
                this->_find =  true ;
                this->_tpos = _TPOS ;

                break ;
                }
            }
        }

        } ;

    class near_tria
        {
        public  :
        real_type              *_ppos ;
        real_type              *_qpos ;

        real_type               _dsqr ;

        mesh_type              *_mesh ;

        bool_type               _find ;
        iptr_type               _tpos ;

        public  :

    /*------------------------ make a tree-tria predicate */
        __inline_call near_tria (
            real_type*_psrc = nullptr ,
            real_type*_qsrc = nullptr ,
            mesh_type*_msrc = nullptr
            ) : _ppos(_psrc) ,
                _qpos(_qsrc) ,
           _dsqr(+std::numeric_limits
           <real_type>::infinity()) ,
                _mesh(_msrc) ,
                _find(false) ,
                _tpos(   -1)   {}

    /*------------------------ call pred. on tree matches */
        __inline_call float operator () (
                typename
            tree_type::item_data  *_iptr
            )
        {
            if (this->_find) return +0. ;

            for ( ; _iptr != nullptr;
                    _iptr = _iptr->_next)
            {
                real_type  _qtmp[+3];
                iptr_type  _TPOS =
                    _iptr->_data.ipos() ;

                if (near_pred ( _ppos ,
                        _qtmp ,*_mesh ,
                        _TPOS ) )
                {
    /*------------------------ projected match: keep best */
                real_type _dtmp =
            geometry::lensqr_3d(_ppos, _qtmp);

                if (_dtmp<_dsqr )   // n.b. proj. on tria
                {
                _qpos[0] = _qtmp[0] ;
                _qpos[1] = _qtmp[1] ;
                _qpos[2] = _qtmp[2] ;

                this->_dsqr = _dtmp ;
                this->_tpos = _TPOS ;
                }

                }
                else
                {
    /*------------------------ projected match: keep best */
                real_type _dtmp =
            geometry::lensqr_3d(_ppos, _qtmp);

                if (_dtmp<_dsqr )   // n.b. proj. to rest
                {
                _qpos[0] = _qtmp[0] ;
                _qpos[1] = _qtmp[1] ;
                _qpos[2] = _qtmp[2] ;

                this->_dsqr = _dtmp ;
                this->_tpos = _TPOS ;
                }

                }
            }

            return (float) this->_dsqr ;
        }

        } ;

    /*
    --------------------------------------------------------
     * NEAR-PRED: TRUE if PPOS is in TPOS.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type near_pred (
        real_type*_ppos ,
        real_type*_qpos ,
        mesh_type&_mesh ,
        iptr_type _tpos
        )
    {
        geometry::hits_type _hits ;
        if (geometry::proj_tria_3d(_ppos,
           &_mesh._set1[
            _mesh._set3[
            _tpos].node(0)].pval(0) ,
           &_mesh._set1[
            _mesh._set3[
            _tpos].node(1)].pval(0) ,
           &_mesh._set1[
            _mesh._set3[
            _tpos].node(2)].pval(0) ,
            _qpos,_hits) )
        {
            return (_hits ==
                 geometry::face_hits) ;
        }
        else
        {
            return ( false ) ;
        }
    }

    /*
    --------------------------------------------------------
     * HINT: check for valid index.
    --------------------------------------------------------
     */

    __inline_call bool_type hint_okay (
        hint_type _hint
        )
    {
    /*------------------------ test whether hint is valid */
        return _hint >= (iptr_type)0
         && _hint < (iptr_type)
            this->_mesh._set3.count()
         && this->_mesh.
           _set3 [_hint].mark() >= 0 ;
    }

    /*
    --------------------------------------------------------
     * EVAL: eval. size-fun. value.
    --------------------------------------------------------
     */

    __normal_call real_type eval (
        real_type *_ppos ,
        hint_type &_hint
        )
    /*------------------------ find tria + linear interp. */
    {
        real_type  _QPOS[3] = {
       (real_type) _ppos[0] ,
       (real_type) _ppos[1] ,
       (real_type) _ppos[2] } ;

        float      _PPOS[3] = {
       (float)     _ppos[0] ,
       (float)     _ppos[1] ,
       (float)     _ppos[2] } ;

        real_type _hOUT =
    +std::numeric_limits<real_type>::infinity() ;

        if (hint_okay(_hint))
        {
    /*------------------------ test whether hint is valid */
            if(!near_pred( _ppos,
                    _QPOS, _mesh,
                    _hint)  )
            {
            _hint =  this->null_hint();
            }
        }
        else
        {
    /*------------------------ hint is definitely invalid */
            _hint =  this->null_hint();
        }

        if (_hint == this->null_hint())
        {
    /*------------------------ outside: find nearest tria */
            near_tria _func (_ppos,
                      _QPOS,&_mesh) ;

            this->
           _tree.near(_PPOS, _func) ;

           _hint =  _func._tpos ;
        }

        if (_hint != this->null_hint())
        {
    /*------------------------ linear interp. within tria */
        real_type _hsum = (real_type)+.0 ;
        real_type _vsum = (real_type)+.0 ;

        for(auto _fpos = 3; _fpos-- != 0; )
        {
            iptr_type  _fnod [3] ;
            tri3_type::
            face_node(_fnod, _fpos, 2, 1) ;

            _fnod[0] = this->_mesh.
            _set3[_hint].node(_fnod[0]);
            _fnod[1] = this->_mesh.
            _set3[_hint].node(_fnod[1]);
            _fnod[2] = this->_mesh.
            _set3[_hint].node(_fnod[2]);

            real_type _tvol =
                geometry::tria_area_3d (
               &this->_mesh.
               _set1[_fnod[0]].pval(0) ,
               &this->_mesh.
               _set1[_fnod[1]].pval(0) ,
               _QPOS) ;

            _hsum += _tvol *
                    this->_hval[_fnod[2]] ;

            _vsum += _tvol ;
        }


        _hOUT = _hsum / _vsum ;

        }

    /*------------------------- size-fun interp. to ppos. */
        return  _hOUT ;
    }

    } ;


    }

#   endif   //__HFUN_MESH_ELLIPSOID_3__



