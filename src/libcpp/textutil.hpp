
    /*
    --------------------------------------------------------
     * TEXTUTIL: utilities for string/file parsing.
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
     * Last updated: 07 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __TEXTUTIL__
#   define __TEXTUTIL__

    /*
    --------------------------------------------------------
     * TO-STRING: to_string, but controls precision.
    --------------------------------------------------------
     */

    template <
    typename      real_type
             >
    __static_call
    __inline_call std::string to_string_prec (
        real_type const&_real,
        std::size_t     _prec
        )
    {
        std::stringstream  _sstr ;
        _sstr << std::scientific
              << std::setprecision(_prec)
              << _real ;
        return   _sstr.str() ;
    }

    /*
    --------------------------------------------------------
     * TRIM: trim leading and trailing "whitespace".
    --------------------------------------------------------
     */

    __static_call
    __inline_call std::string trim (
        std::string const&_line
        )
    {
        typename std::string::size_type _pos1 =
            _line.find_first_not_of(" \n\r\t") ;

        typename std::string::size_type _pos2 =
            _line. find_last_not_of(" \n\r\t") ;

        if (_pos1 != std::string::npos &&
            _pos2 != std::string::npos )
        {
            return _line.substr(
                _pos1, _pos2-_pos1 + 1 );
        }
        else
        {
            return ( std::string() ) ;
        }
    }

    /*
    --------------------------------------------------------
     * CENTRED: "centre" a string within a string.
    --------------------------------------------------------
     */

    __static_call
    __inline_call std::string centred (
        std::size_t _wide,
        const std::string& _ssrc)
    {
        std::size_t _slen = _ssrc.length() ;
        if (_wide < _slen)
        {
            return  _ssrc;
        }

        std::size_t _diff = _wide - _slen ;
        std::size_t _pad1 = _diff / +2 ;
        std::size_t _pad2 = _diff - _pad1 ;

        return std::string(_pad1, ' ')
            + _ssrc
            +  std::string(_pad2, ' ') ;
    }

    /*
    --------------------------------------------------------
     * FIND-TOKS: push delimited string onto tokens.
    --------------------------------------------------------
     */

    __normal_call bool_type test_char (
        char               _char,
        std::string const& _dlim
        )
    {
        typename std::
            string::const_iterator
        _dpos = _dlim.begin(),
        _tend = _dlim.end  ();
        for ( ; _dpos != _tend; ++_dpos)
        {
        if (_char ==*_dpos) return true;
        }
        return ( false )  ;
    }

    template <
    typename      list_type
             >
    __normal_call void_type find_toks (
        std::string const& _istr,
        std::string const& _dlim,
        list_type        & _list
        )
    {
        typename std::
            string::const_iterator
        _hpos = _istr.begin() ,
        _tpos = _istr.begin() ,
        _tend = _istr.end  () ;
        for ( ; _tpos != _tend; )
        {
        /*--------------------------- find next delimiter */
            for (; _tpos != _tend; ++_tpos )
            {
                if (test_char(*_tpos,_dlim))
                    break ;
            }
            if (_tpos == _hpos)
            {
        /*--------------------------- push zero-len token */
                _list.push_tail("") ;
            }
            else
            {
        /*-------------------------- push non-empty token */
                _list.push_tail  (
                    std::move(std::string (
                        _hpos, _tpos))) ;
        /*-------------------------- reached end-of-input */
                if (_tpos == _tend) break ;
            }
            _hpos = _tpos + 1;
            _tpos = _tpos + 1;
        }
    }

    /*
    --------------------------------------------------------
     * FILE-PART: split a file name into path-name-fext.
    --------------------------------------------------------
     */

    __normal_call void_type file_part (
        std::string const& _fsrc,
        std::string      & _path,
        std::string      & _name,
        std::string      & _fext
        )
    {
        typename std::string::size_type
            _spos = _fsrc.find_last_of("\\/");

        typename std::string::size_type
            _dpos = _fsrc.find_last_of("."  );

        typename std::string::const_iterator
            _pos0, _pos1, _pos2,
            _pos3, _pos4, _pos5;

        if (_spos != std::string::npos )
        {
            _pos0 = _fsrc.begin();
            _pos1 = _fsrc.begin()+_spos-0;
            _pos2 = _fsrc.begin()+_spos+1;
        }
        else
        {
            _pos0 = _fsrc.begin();
            _pos1 = _fsrc.begin();
            _pos2 = _fsrc.begin();
        }

        if (_dpos != std::string::npos &&
           (_spos == std::string::npos ||
            _dpos >= _spos) )
        {
            _pos3 = _fsrc.begin()+_dpos-0;
            _pos4 = _fsrc.begin()+_dpos+1;
            _pos5 = _fsrc.end  ();
        }
        else
        {
            _pos3 = _fsrc.end  ();
            _pos4 = _fsrc.end  ();
            _pos5 = _fsrc.end  ();
        }

        _path = std::string(_pos0, _pos1);
        _name = std::string(_pos2, _pos3);
        _fext = std::string(_pos4, _pos5);
    }


#   endif   //__TEXTUTIL__



