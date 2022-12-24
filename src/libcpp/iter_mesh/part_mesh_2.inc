
    /*
    --------------------------------------------------------
     * PART-MESH-2: form dynamic parallel partitions.
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
    --------------------------------------------------------
     *
     * Last updated: 15 Dec., 2022
     *
     * Copyright 2013-2022
     * Marc Tunnell, 
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from iter_mesh_2.hpp


    __static_call
    __normal_call void_type part_mesh (
        mesh_type &_mesh,   // full mesh obj.
        part_sets &_part,   // partition data-structure
        iptr_type  _nprt,   // number of partitions to make
        iptr_list  _aset,   // list of active nodes
        iptr_type  _isub,   // inner iteration marker
        iter_stat &_tcpu    // CPU timers/info
        )
    {
    #   define PUSH_BNDS( _PNUM )   \
        __isflip(_PNUM) ? _PNUM : __doflip(_PNUM)

    #   define PART_LIMS( _PNUM )   \
        if (_PNUM != -1)            \
        {                           \
            _have += +1 ;           \
            _pmin = std::min(       \
                    _pmin, __unflip(_PNUM));    \
            _pmax = std::max(       \
                    _pmax, __unflip(_PNUM));    \
        }

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

        __unreferenced(_isub) ;

        if (_nprt > +1) 
        {
        /*------ define AABB-tree for recursive partition */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            typedef geom_tree::aabb_item_node_k <
                    real_type,
                    iptr_type,
                    pred_type::geom_dims>   tree_item ;

            typedef geom_tree::aabb_tree<
                    tree_item,
                    pred_type::geom_dims>   tree_type ;

            tree_type   _tree;

            auto _npop =  // no. nodes per part for nprts
                (8 * _aset.count()) / _nprt / 3 ;

        /*-------------------- populate tree: xyz + index */
            containers::array<tree_item>  _item ;           
            for (auto _apos  = _aset.tail() ;
                      _apos != _aset.hend() ;
                    --_apos  ) // backward maintains sort?
            {
                 auto _aptr = 
                _mesh.node().head() + *_apos;

                if (_aptr->mark() >= 0 && 
                   (_aptr->ndeg(TRIA3_tag) > 0 ||
                    _aptr->ndeg(QUAD4_tag) > 0 ))
                {
                    _item.push_tail();
                    _item.tail()->ipos() = *_apos;

                    for (auto _idim = geom_dims;
                              _idim-- != 0; )
                        _item.tail()->pval(_idim) =
                              _aptr ->pval(_idim);
                }
            }

        /*-------------------- split: npop nodes per leaf */
            _tree.load(
                _item.head(), _item.tend(),
                _npop ,  // pop. threshold, for leaf div.
                _tree.median_split ,
                 +1.0 ,  // set LONG=1, only div. by pop.
                 +0.0) ; // set VTOL=0, only div. by pop.

        /*-------------------- form sparse lists for part */

        //  LIST[LPTR[K]:LPTR[K+1]-1]: list of nodes that 
        //  comprise the K-th partition

            _part._lptr.clear() ;
            _part._list.clear() ;
            _part._part.clear() ;
            _part._seqs.clear() ;

            _tree.dcmp(_part._lptr,_part._list) ;
            
    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._tree_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
            _tcpu._full_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
    #       endif//__use_timers

        /*-------------------- and build reverse indexing */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            _part._part.set_count(
                _mesh.node().count(), 
                    containers::loose_alloc, -1) ;

        //  PART[J] = K: the partition for J-th mesh node

            for (auto 
                _lnum = _part._lptr.count() - 1  ;
                _lnum-- != 0 ;  )
            {
                for (auto _iptr = _part._lptr[_lnum+0];
                          _iptr < _part._lptr[_lnum+1];
                        ++_iptr )
                _part._part[_part._list[_iptr]] = _lnum ;
            }

        //  PART[J] =-1: the J-th node is on an interface

            for (auto _tria  = _mesh.tri3().head() ;
                      _tria != _mesh.tri3().tend() ;
                    ++_tria  )
            {
                if (_tria->mark() >= +0)
                {
        /*-------------- interface if cell straddles part */
                    iptr_type _inod, _jnod, _knod;
                    _inod = _tria->node( 0);
                    _jnod = _tria->node( 1);
                    _knod = _tria->node( 2);

                    iptr_type _have = 0;
                    iptr_type _pmin = 
                std::numeric_limits<iptr_type>::max() ;
                    iptr_type _pmax = 
                std::numeric_limits<iptr_type>::min() ;

                    PART_LIMS(_part._part[_inod]);
                    PART_LIMS(_part._part[_jnod]);
                    PART_LIMS(_part._part[_knod]);

                    if (_have==0 || _pmin == _pmax) 
                        continue ;

                    _part._part[_inod] = PUSH_BNDS (
                        _part._part[_inod]);
                    _part._part[_jnod] = PUSH_BNDS (
                        _part._part[_jnod]);
                    _part._part[_knod] = PUSH_BNDS (
                        _part._part[_knod]);
                }
            }

            for (auto _quad  = _mesh.quad().head() ;
                      _quad != _mesh.quad().tend() ;
                    ++_quad  )
            {
                if (_quad->mark() >= +0)
                {
        /*-------------- interface if cell straddles part */
                    iptr_type _inod, _jnod, _knod,
                              _lnod;
                    _inod = _quad->node( 0);
                    _jnod = _quad->node( 1);
                    _knod = _quad->node( 2);
                    _lnod = _quad->node( 3);

                    iptr_type _have = 0;
                    iptr_type _pmin = 
                std::numeric_limits<iptr_type>::max() ;
                    iptr_type _pmax = 
                std::numeric_limits<iptr_type>::min() ;

                    PART_LIMS(_part._part[_inod]);
                    PART_LIMS(_part._part[_jnod]);
                    PART_LIMS(_part._part[_knod]);
                    PART_LIMS(_part._part[_lnod]);

                    if (_have==0 || _pmin == _pmax) 
                        continue ;

                    _part._part[_inod] = PUSH_BNDS (
                        _part._part[_inod]);
                    _part._part[_jnod] = PUSH_BNDS (
                        _part._part[_jnod]);
                    _part._part[_knod] = PUSH_BNDS (
                        _part._part[_knod]);
                    _part._part[_lnod] = PUSH_BNDS (
                        _part._part[_lnod]);
                }
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._part_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
            _tcpu._full_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
    #       endif//__use_timers

        /*-------------- add interface to sequential list */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            for (auto _apos  = _aset.head() ;
                      _apos != _aset.tend() ;
                    ++_apos  ) // aset sweep ensures sort
            {
                if (_part._part[*_apos] < -1)
                {
                    _part._seqs.push_tail(*_apos);
                }
            }

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._seqs_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
            _tcpu._full_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
    #       endif//__use_timers
        } 
        else 
        {
        /*-------------- dummy part lists for single part */
    #       ifdef  __use_timers
            _ttic = _time.now() ;
    #       endif//__use_timers

            _part._lptr.set_count(+2);
            _part._list.set_alloc(_aset.count()) ;
            _part._part.set_count(
                _mesh.node().count(), 
                    containers::loose_alloc, +0) ;

            _part._lptr[0] = +0;
            _part._lptr[1] = _aset.count() ;
            _part._list.push_tail(
                _aset.head(), _aset.tend());

    #       ifdef  __use_timers
            _ttoc = _time.now() ;
            _tcpu._seqs_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
            _tcpu._full_part += 
                  _tcpu.time_span(_ttic, _ttoc) ;
    #       endif//__use_timers
        }

        /*
        iptr_type _ipos = +0;
        for (auto _iter  = _mesh.node().head();
                  _iter != _mesh.node().tend();
                ++_iter, ++_ipos)
        {
        //  just for debug: output part idx. to node itag 
            if (_iter->mark() >= 0)
            {
                _iter->itag() = _part._part[_ipos];
            }
        }
    
        iptr_type _pmin = 
            std::numeric_limits<iptr_type>::max();
        iptr_type _pmax = 
            std::numeric_limits<iptr_type>::min();
        for (auto _rank  = 0u ;
            _rank<_part._lptr.count()-1; ++_rank)
        {
        // just for debug: output min/max sizes of part's
            iptr_type _size = 
                _part._lptr[_rank+1] - 
                _part._lptr[_rank+0] ;
            _pmin = std::min(_pmin, _size);
            _pmax = std::max(_pmax, _size);
        }
        std::cout << _pmin << ", " << _pmax << std::endl;
        */

    #   undef   PART_LIMS
    #   undef   PUSH_BNDS
    }
    
