
    /*
    --------------------------------------------------------
     * FIND: search-within-sequence algorithms.
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
     * Last updated: 01 Sept., 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __FIND__
#   define __FIND__

#   include "../containers/iter_base.hpp"

    namespace algorithms {

    /*
    --------------------------------------------------------
     * LOWER-BOUND: returns lowermost index "not <" XVAL.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  data_type ,
        typename  pred_type
             >
    __normal_call iter_type lower_bound (
        iter_type _head ,  // head of sequence to search
        iter_type _tend ,  // last of sequence to search
        data_type _xval ,  // value to search for
        pred_type _less    // comparison predicate
        )
    {
        iter_type _iter ;
        typename containers::
            iterator_traits<iter_type>::
        diff_type _span  = _tend-_head;
        typename containers::
            iterator_traits<iter_type>::
        diff_type _step = +0;

        while (_span > +0)
        {
            _step = _span / +2;

            _iter = _head + _step;

            if (_less (*_iter, _xval))
            {
                _head = _iter + 1;
                _span-= _step + 1;
            }
            else
            {
                _span = _step + 0;
            }

        }

        return _head ;
    }

    /*
    --------------------------------------------------------
     * UPPER-BOUND: returns uppermost index ">" XVAL.
    --------------------------------------------------------
     */

    template <
        typename  iter_type ,
        typename  data_type ,
        typename  pred_type
             >
    __normal_call iter_type upper_bound (
        iter_type _head ,  // head of sequence to search
        iter_type _tend ,  // last of sequence to search
        data_type _xval ,  // value to search for
        pred_type _less    // comparison predicate
        )
    {
        iter_type _iter ;
        typename containers::
            iterator_traits<iter_type>::
        diff_type _span  = _tend-_head;
        typename containers::
            iterator_traits<iter_type>::
        diff_type _step = +0;

        while (_span > +0)
        {
            _step = _span / +2;

            _iter = _head + _step;

            if(!_less ( _xval,*_iter))
            {
                _head = _iter + 1;
                _span-= _step + 1;
            }
            else
            {
                _span = _step + 0;
            }

        }

        return _head ;
    }


    }

#   endif   //__FIND__



