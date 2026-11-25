
    /*
    --------------------------------------------------------
     * HFUN-SEED-kD: initialise the MARCHE heap.
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
     * Last updated: 24 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HFUN_SEED_K__
#   define __HFUN_SEED_K__

    namespace mesh {

    template <    size_t N ,
        typename  vals_list, 
        typename  seed_list,
        typename  node_func
             >
    __inline_call void_type hfun_mark_seed (
        vals_list const& _vals,
        seed_list      & _seed,
        node_func const& _node
        )
    {
        auto const _data =_vals[_node( 0 )] ;
        size_t _ii = 1;
        while (_ii < N && _vals[_node(_ii)] == _data) ++_ii;

        if (_ii == N) return;  // cell == const.

        for (_ii = 0; _ii < N; ++_ii) _seed[_node(_ii)] = 1;
    }

    template <
        typename  seed_list, 
        typename  heap_type
             >
    __inline_call void_type hfun_push_seed (
        seed_list const& _seed,
        heap_type      & _sort
        )
    {
        using iptr_type = typename heap_type::iptr_type;

        for (size_t _ii = 0; _ii < _seed.count(); ++_ii)
            if (_seed[_ii])
                _sort.push(
                   (iptr_type)_ii, (iptr_type)_ii) ;
    }

    template <    size_t N , 
        typename  vals_list, 
        typename  node_list,
        typename  cell_list, 
        typename  heap_type
             >
    __normal_call void_type hfun_seed_mesh (
        vals_list const& _vals,
        node_list const& _node,
        cell_list const& _cell,
        heap_type      & _sort
        )
    {
        containers::array<unsigned char,
            typename heap_type::allocator> _seed ;
        _seed.set_count(
        _node.count(), containers::tight_alloc, 0) ;

        for (auto _icel  = _cell.head(); 
                  _icel != _cell.tend(); ++_icel)
        {
            if (_icel->mark() < 0) continue;

            hfun_mark_seed<N>(_vals, _seed ,
                [&_icel](size_t _ipos) { 
                    return _icel->node(_ipos); }) ;
        }

        hfun_push_seed(_seed, _sort) ;
    }

    template <
        typename  vals_list, 
        typename  heap_type
             >
    __normal_call void_type hfun_seed_grid (
        vals_list const& _vals,
        size_t _nx, size_t _ny,
        heap_type      & _sort
        )
    {
        if (_nx < 2 || _ny < 2) return ;

        containers::array<unsigned char,
            typename heap_type::allocator> _seed ;
        _seed.set_count(
        _vals.count(), containers::tight_alloc, 0) ;

        for (size_t _jj = 0; _jj + 1 < _nx; ++_jj)
        for (size_t _ii = 0; _ii + 1 < _ny; ++_ii)
        {
            auto const _kk = _jj * _ny + _ii;
            size_t const _node[] = {
                _kk, _kk+1, _kk+_ny+1, _kk+_ny
                } ;
            hfun_mark_seed<4>(_vals, _seed,
                [&_node](size_t _ipos) { return _node[_ipos]; }) ;
        }

        hfun_push_seed(_seed, _sort) ;
    }

    template <
        typename  vals_list, 
        typename  heap_type
             >
    __normal_call void_type hfun_seed_grid (
        vals_list const& _vals,
        size_t _nx, size_t _ny, size_t _nz,
        heap_type      & _sort
        )
    {
        if (_nx < 2 || _ny < 2 || _nz < 2) return;

        containers::array<unsigned char,
            typename heap_type::allocator> _seed ;
        _seed.set_count(
        _vals.count(), containers::tight_alloc, 0) ;

        auto const _xy = _nx * _ny;

        for (size_t _kk = 0; _kk + 1 < _nz; ++_kk)
        for (size_t _jj = 0; _jj + 1 < _nx; ++_jj)
        for (size_t _ii = 0; _ii + 1 < _ny; ++_ii)
        {
            auto const _aa = 
                _kk * _xy + _jj * _ny + _ii;
            auto const _bb = _aa + _xy ;
            size_t const _node[] = {
                _aa, _aa+1, _aa+_ny+1, _aa+_ny,
                _bb, _bb+1, _bb+_ny+1, _bb+_ny
                } ;
            hfun_mark_seed<8>(_vals, _seed,
                [&_node](size_t _ipos) { return _node[_ipos]; }) ;
        }

        hfun_push_seed(_seed, _sort) ;
    }

    }

#   endif//__HFUN_SEED_K__
