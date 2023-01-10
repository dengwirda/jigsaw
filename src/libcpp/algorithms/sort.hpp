
    /*
    --------------------------------------------------------
     * SORT: sorting algorithms.
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
     * Last updated: 24 Jul., 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __SORT__
#   define __SORT__

#   include "../containers/iter_base.hpp"

    namespace algorithms {

    /*
    --------------------------------------------------------
     * IS-SORTED: returns TRUE if sequence is sorted.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __normal_call bool_type is_sorted (
        iter_type _head ,  // head of sequence to test
        iter_type _tail ,  // last of sequence to test
        pred_type _less    // comparison predicate
        )
    {
        for ( ; _head+1 != _tail; ++_head )
        { /* check pair ordering */
            if (_less(*(_head+1),  *_head))
            {
                return ( false )  ;
            }
        }

        return ( true ) ;
    }

    /*
    --------------------------------------------------------
     * SORT-N: explicit sort of n-tuples for small n.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __inline_call void_type sort_2 ( // sort 2-tuple
        iter_type _aa,
        iter_type _bb,
        pred_type const&_less
        )
    {
        if (_less(*_bb, *_aa))
        std::swap(*_aa, *_bb);
    }

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __inline_call void_type sort_3 ( // sort 3-tuple
        iter_type _aa,
        iter_type _bb,
        iter_type _cc,
        pred_type const&_less
        )
    {
        if (_less(*_bb, *_aa))
        std::swap(*_aa, *_bb);

        if (_less(*_cc, *_bb))
        {
        std::swap(*_bb, *_cc);

        if (_less(*_bb, *_aa))
        std::swap(*_aa, *_bb);
        }
    }

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __inline_call void_type sort_4 ( // sort 4-tuple
        iter_type _aa,
        iter_type _bb,
        iter_type _cc,
        iter_type _dd,
        pred_type const&_less
        )
    {
        if (_less(*_bb, *_aa))
        std::swap(*_aa, *_bb);

        if (_less(*_dd, *_cc))
        std::swap(*_cc, *_dd);

        if (_less(*_cc, *_aa))
        std::swap(*_aa, *_cc);

        if (_less(*_dd, *_bb))
        std::swap(*_bb, *_dd);

        if (_less(*_cc, *_bb))
        std::swap(*_bb, *_cc);
    }

    /*
    --------------------------------------------------------
     * I-SORT: insertion sort.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __normal_call void_type isort (
        iter_type _head,  // head of sequence to sort
        iter_type _tend,  // last of sequence to sort
        pred_type _less   // comparison predicate
        )
    {
    /*--------------------------- sort small input ranges */
        switch (_tend - _head)
        {
        /*-------------------------------- quick exit */
        case 0:
        case 1: return ;
        case 2:
        /*-------------------------------- length = 2 */
            {
            sort_2(_head+0, _head+1, _less);
            return ;
            }
        case 3:
        /*-------------------------------- length = 3 */
            {
            sort_3(_head+0,
                   _head+1, _head+2, _less);
            return ;
            }
        case 4:
        /*-------------------------------- length = 4 */
            {
            sort_4(_head+0, _head+1,
                   _head+2, _head+3, _less);
            return ;
            }
        }

    /*-------- scan sequence, find min., swap as sentinel */
        iter_type _mark  = _head + 0 ;
        for (auto _item  = _head + 1 ;
                  _item != _tend ; ++ _item )
        {
        /*-------- _item becomes _mark if new minimum */
            if(_less( *_item,*_mark)) _mark = _item;
        }

        if (_mark != _head)
        {
        /*-------- swap _mark to _head if new minimum */
            std::swap(*_mark, *_head);
        }

    /* insertion sort: sift down toward head until sorted */
        for(_head += 2; _head != _tend; ++_head)
        {
            if (_less( *_head, *(_head-1)))
            {
                typename containers::
                    iterator_traits<iter_type>::
                        data_type _temp =
                            std::move (*_head) ;

                iter_type _item = _head ;

            /*---------- sift pivot into sorted order */
                do
                {
                    *_item = std::move(*(_item-1)) ;
                   --_item ;
                }
                while ( _less( _temp,  *(_item-1)));

            /*---------- swap pivot into new position */
               *_item = std::move( _temp) ;
            }
        }
    }

    /*
    --------------------------------------------------------
     * S-SORT: shell sort.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __normal_call void_type ssort (
        iter_type _head,  // head of sequence to sort
        iter_type _tend,  // head of sequence to sort
        pred_type _less   // comparison predicate
        )
    {
    /*--------------------------- sort small input ranges */
        switch (_tend - _head)
        {
        /*-------------------------------- quick exit */
        case 0:
        case 1: return ;
        case 2:
        /*-------------------------------- length = 2 */
            {
            sort_2(_head+0, _head+1, _less);
            return ;
            }
        case 3:
        /*-------------------------------- length = 3 */
            {
            sort_3(_head+0,
                   _head+1, _head+2, _less);
            return ;
            }
        case 4:
        /*-------------------------------- length = 4 */
            {
            sort_4(_head+0, _head+1,
                   _head+2, _head+3, _less);
            return ;
            }
        }

    /* Niels Pardons, A154393: "Empirically good sequence of
     * increments for the shell sort algorithm", 2009, OEIS.
     */
        typename containers::
            iterator_traits<iter_type>::
                size_type _finc = _tend - _head;
        typename containers::
            iterator_traits<iter_type>::
                size_type _incn = +12 ;
        typename containers::
            iterator_traits<iter_type>::
            size_type static constexpr _incv[12]
            = {   1,        9,        34,
                182,      836,      4025,
              19001,    90358,    428481,
            2034035,  9651787,  45806244}  ;

    /* insertion sort: sift down toward head until sorted */
        for ( ; _incn != 0 ;  --_incn)
        {
            if (_finc <= _incv[_incn-1]) continue;

        /* just a pass of insertion sort with stride ival */
            typename containers::
                iterator_traits<iter_type>::
                    data_type _temp ;
            typename containers::
                iterator_traits<iter_type>::
                    size_type _ival ;

            _ival = _incv[_incn - 1];

            iter_type _ipos =_head + _ival ;
            iter_type _iend =_head + _ival ;
            iter_type _item ;
            for ( ; _ipos < _tend; ++_ipos )
            {
            if (_less(*_ipos,*(_ipos-_ival)))
            {
                _temp = std::move(*_ipos) ;
                _item = _ipos ;
            /*---------- sift pivot into sorted order */
                do
                {
                    *_item = std::move(*(_item-_ival)) ;
                }
                while ((_item-=_ival) >= _iend &&
                        _less( _temp, *( _item-_ival)));

            /*---------- swap pivot into new position */
               *_item = std::move(_temp) ;
            }
            }
        }

    }

    /*
    --------------------------------------------------------
     * Q-SORT: quick sort.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __inline_call iter_type pivot (     // pivot selection
        iter_type _head,
        iter_type _tail,
        pred_type const&_less
        )
    {
        //typename containers::
        //    iterator_traits<iter_type>::
        //diff_type static constexpr _LONG = 256;

        iter_type _imid =
            _head + (_tail - _head) / 2 ;

        if (false) //(_tail - _head >= _LONG)
        {
    /*-------------- median-of-5 choice for pivot element */
            typename containers::
                iterator_traits<iter_type>::
                    size_type _step =
                        (_tail - _head) / 4;

            iter_type _pos1 =
                _head + _step * 1 ;
            iter_type _pos3 =
                _head + _step * 3 ;

            sort_3(_head, _pos1, _imid, _less);
            sort_3(_imid, _pos3, _tail, _less);
            sort_3(_pos1, _imid, _pos3, _less);
        }
        else
        {
    /*-------------- median-of-3 choice for pivot element */
            sort_3(_head, _imid, _tail, _less);
        }

        return ( _imid ) ;
    }

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __normal_call void_type qsort (  // unrolled quick sort
        iter_type _head,
        iter_type _tend,  // last of sequence to sort
        pred_type _less   // comparison predicate
        )
    {
        typedef containers::
            iterator_traits<iter_type> iter_base ;

    /*------- tiny-list threshold, push to insertion sort */
        typename iter_base::
            diff_type static constexpr _LONG = 16;

    /*------- stack max depth = log2(maxint/sizeof(data)) */
        typename iter_base::
            size_type static constexpr _DEEP =
        8 * sizeof(typename iter_base::size_type);

        class node_type
            {
            public  :
            iter_type         _obj[2] ;
            public  :
            __inline_call void push (
                iter_type  _ll ,
                iter_type  _rr
                )
            {   _obj[0] =  _ll ;
                _obj[1] =  _rr ;
            }
            __inline_call void _pop (
                iter_type &_ll,
                iter_type &_rr
                )
            {   _ll = _obj [0] ;
                _rr = _obj [1] ;
            }
            } ;

        node_type _node[_DEEP], *_nptr =_node ;

        _nptr->push( _head, _tend-1); ++_nptr ;

    /*------------ non-recursive quicksort implementation */
        do
        {
            iter_type _hh, _tt, _mm;
        /* _pop next partition bounds from stack */
            ( --_nptr)->_pop(_hh, _tt);

            if (_tt - _hh + 1 < _LONG)
            { /* leave small partitions */
                continue ;
            }

        /* find pivot item for current partition */
            _mm = pivot(_hh, _tt, _less);

        /* reduce partition around pivot element */
            iter_type _ll = _hh + 1;
            iter_type _rr = _tt - 1;
            for ( ; _ll < _rr; )
            {
                for(;_less(*_ll, *_mm);) ++_ll;
                for(;_less(*_mm, *_rr);) --_rr;

                if (_ll <  _rr)
                { /* swap elements */
                if (_ll == _mm)
                    _mm =  _rr;
                else
                if (_rr == _mm)
                    _mm =  _ll;

                std::swap(*_ll,
                          *_rr) ;
                }
                if (_ll <= _rr)
                { /* shrink bounds */
                  ++_ll; --_rr  ;
                }
            }

        /* push partitions onto stack, big first */
            if (_rr - _hh > _tt - _ll)
            {
                _nptr->push(_hh, _rr); ++_nptr;
                _nptr->push(_ll, _tt); ++_nptr;
            }
            else
            {
                _nptr->push(_ll, _tt); ++_nptr;
                _nptr->push(_hh, _rr); ++_nptr;
            }
        }
        while (_nptr !=_node) ;

    /*------- sort remaining sequences via insertion sort */
        isort (_head , _tend, _less) ;
    }

    /*
    --------------------------------------------------------
     * Q-SIFT: quick select; sort kth element
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  pred_type
             >
    __normal_call void_type qsift (  // quick select
        iter_type _head,
        iter_type _kk,    // pos. in sequence to find
        iter_type _tend,  // last of sequence to parition
        pred_type _less   // comparison predicate
        )
    {
        iter_type  _hh, _tt, _mm;
        _hh = _head + 0;
        _tt = _tend - 1;

        do {
        /* deal with "spec.-case" tiny partition */
            if (_tt - _hh + 1 <= 1)
                break ;

            if (_tt - _hh + 1 == 2)
            {
                if (_less(*_tt, *_hh))
                std::swap(*_hh, *_tt);
                break ;
            }

            /*
            if (_tt - _hh + 1 <= 4)
            {
                isort ( _hh, _tt, _less);
                break ;
            }
            */

        /* find pivot item for current partition */
            _mm = pivot(_hh, _tt, _less);

        /* reduce partition around pivot element */
            iter_type _ll = _hh + 1;
            iter_type _rr = _tt - 1;
            for ( ; _ll < _rr; )
            {
                for(;_less(*_ll, *_mm);) ++_ll;
                for(;_less(*_mm, *_rr);) --_rr;

                if (_ll <  _rr)
                { /* swap elements */
                if (_ll == _mm)
                    _mm =  _rr;
                else
                if (_rr == _mm)
                    _mm =  _ll;

                std::swap(*_ll,
                          *_rr) ;
                }
                if (_ll <= _rr)
                { /* shrink bounds */
                  ++_ll; --_rr  ;
                }
            }

          //if (_mm == _kk) break ;

            if (_kk <  _rr) _tt = _rr;
            else
            if (_kk >  _ll) _hh = _ll;
            else
            {
                _hh =  _rr+1; _tt = _ll-1;
            }
        }
        while (true) ;
    }


    }

#   endif   //__SORT__



