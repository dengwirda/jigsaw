
    /*
    --------------------------------------------------------
     * HASH-COMPLEX-k: utilities for hashing cells in T^k.
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
     * Last updated: 27 Apr., 2020
     *
     * Copyright 2020--
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __HASH_COMPLEX_K__
#   define __HASH_COMPLEX_K__

    namespace mesh {

    /*
    --------------------------------------------------------
     * HASH-COMPLEX-k hashes for k-cell complexes.
    --------------------------------------------------------
     */

    template <
        typename  node_list
             >
    class node_hash
        {
        public  :
        __const_ptr(node_list) _nset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call node_hash  (
            node_list *_nsrc
            ) : _nset( _nsrc) {}
    /*----------------------- hash node indexing for node */
        __inline_call uint32_t  operator() (
            uint32_t  _npos
            ) const
        {
            return (uint32_t)this->_nset->
                operator[](_npos).node( 0) ;
        }
        } ;

    template <
        typename  edge_list
             >
    class edge_hash
        {
        public  :
        __const_ptr(edge_list) _eset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call edge_hash  (
            edge_list *_esrc
            ) : _eset( _esrc) {}
    /*----------------------- hash node indexing for edge */
        __inline_call uint32_t  operator() (
            uint32_t  _epos
            ) const
        {
            uint32_t  _enod[2];
            _enod[0] = this->_eset->
                operator[](_epos).node( 0) ;
            _enod[1] = this->_eset->
                operator[](_epos).node( 1) ;

            algorithms::isort (
                &_enod[0], &_enod[2] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_enod[0], +2, +137);
        }
        } ;

    template <
        typename  tri3_list
             >
    class tri3_hash
        {
        public  :
        __const_ptr(tri3_list) _tset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call tri3_hash  (
            tri3_list *_tsrc
            ) : _tset( _tsrc) {}
    /*----------------------- hash node indexing for face */
        __inline_call uint32_t  operator() (
            uint32_t  _tpos
            ) const
        {
            uint32_t  _tnod[3];
            _tnod[0] = this->_tset->
                operator[](_tpos).node( 0) ;
            _tnod[1] = this->_tset->
                operator[](_tpos).node( 1) ;
            _tnod[2] = this->_tset->
                operator[](_tpos).node( 2) ;

            algorithms::isort (
                &_tnod[0], &_tnod[3] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_tnod[0], +3, +137);
        }
        } ;

    template <
        typename  quad_list
             >
    class quad_hash
        {
        public  :
        __const_ptr(quad_list) _qset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call quad_hash  (
            quad_list *_qsrc
            ) : _qset( _qsrc) {}
    /*----------------------- hash node indexing for face */
        __inline_call uint32_t  operator() (
            uint32_t  _qpos
            ) const
        {
            uint32_t  _qnod[4];
            _qnod[0] = this->_qset->
                operator[](_qpos).node( 0) ;
            _qnod[1] = this->_qset->
                operator[](_qpos).node( 1) ;
            _qnod[2] = this->_qset->
                operator[](_qpos).node( 2) ;
            _qnod[3] = this->_qset->
                operator[](_qpos).node( 3) ;

            algorithms::isort (
                &_qnod[0], &_qnod[4] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_qnod[0], +4, +137);
        }
        } ;

    template <
        typename  tri4_list
             >
    class tri4_hash
        {
        public  :
        __const_ptr(tri4_list) _tset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call tri4_hash  (
            tri4_list *_tsrc
            ) : _tset( _tsrc) {}
    /*----------------------- hash node indexing for face */
        __inline_call uint32_t  operator() (
            uint32_t  _tpos
            ) const
        {
            uint32_t  _tnod[4];
            _tnod[0] = this->_tset->
                operator[](_tpos).node( 0) ;
            _tnod[1] = this->_tset->
                operator[](_tpos).node( 1) ;
            _tnod[2] = this->_tset->
                operator[](_tpos).node( 2) ;
            _tnod[3] = this->_tset->
                operator[](_tpos).node( 3) ;

            algorithms::isort (
                &_tnod[0], &_tnod[4] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_tnod[0], +4, +137);
        }
        } ;

    template <
        typename  hexa_list
             >
    class hexa_hash
        {
        public  :
        __const_ptr(hexa_list) _hset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call hexa_hash  (
            hexa_list *_hsrc
            ) : _hset( _hsrc) {}
    /*----------------------- hash node indexing for face */
        __inline_call uint32_t  operator() (
            uint32_t  _hpos
            ) const
        {
            uint32_t  _hnod[8];
            _hnod[0] = this->_hset->
                operator[](_hpos).node( 0) ;
            _hnod[1] = this->_hset->
                operator[](_hpos).node( 1) ;
            _hnod[2] = this->_hset->
                operator[](_hpos).node( 2) ;
            _hnod[3] = this->_hset->
                operator[](_hpos).node( 3) ;
            _hnod[4] = this->_hset->
                operator[](_hpos).node( 4) ;
            _hnod[5] = this->_hset->
                operator[](_hpos).node( 5) ;
            _hnod[6] = this->_hset->
                operator[](_hpos).node( 6) ;
            _hnod[7] = this->_hset->
                operator[](_hpos).node( 7) ;

            algorithms::isort (
                &_hnod[0], &_hnod[8] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_hnod[0], +8, +137);
        }
        } ;

    template <
        typename  wedg_list
             >
    class wedg_hash
        {
        public  :
        __const_ptr(wedg_list) _wset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call wedg_hash  (
            wedg_list *_wsrc
            ) : _wset( _wsrc) {}
    /*----------------------- hash node indexing for face */
        __inline_call uint32_t  operator() (
            uint32_t  _wpos
            ) const
        {
            uint32_t _wnod[6];
            _wnod[0] = this->_wset->
                operator[](_wpos).node( 0) ;
            _wnod[1] = this->_wset->
                operator[](_wpos).node( 1) ;
            _wnod[2] = this->_wset->
                operator[](_wpos).node( 2) ;
            _wnod[3] = this->_wset->
                operator[](_wpos).node( 3) ;
            _wnod[4] = this->_wset->
                operator[](_wpos).node( 4) ;
            _wnod[5] = this->_wset->
                operator[](_wpos).node( 5) ;

            algorithms::isort (
                &_wnod[0], &_wnod[6] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_wnod[0], +6, +137);
        }
        } ;

    template <
        typename  pyra_list
             >
    class pyra_hash
        {
        public  :
        __const_ptr(pyra_list) _pset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call pyra_hash  (
            pyra_list *_psrc
            ) : _pset( _psrc) {}
    /*----------------------- hash node indexing for face */
        __inline_call uint32_t  operator() (
            uint32_t  _ppos
            ) const
        {
            uint32_t  _pnod[5];
            _pnod[0] = this->_pset->
                operator[](_ppos).node( 0) ;
            _pnod[1] = this->_pset->
                operator[](_ppos).node( 1) ;
            _pnod[2] = this->_pset->
                operator[](_ppos).node( 2) ;
            _pnod[3] = this->_pset->
                operator[](_ppos).node( 3) ;
            _pnod[4] = this->_pset->
                operator[](_ppos).node( 4) ;

            algorithms::isort (
                &_pnod[0], &_pnod[5] ,
                    std::less<uint32_t >());

            return hash::hashword (
           (uint32_t *)&_pnod[0], +5, +137);
        }
        } ;

    template <
        typename  node_list
             >
    class node_pred
        {
        public  :
        __const_ptr(node_list) _nset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call node_pred  (
            node_list *_nsrc
            ) : _nset( _nsrc) {}
    /*----------------------- compute "equal-to" for node */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[1];
            _inod[0] = this->_nset->
                operator[](_ipos).node( 0) ;

            uint32_t  _jnod[1];
            _jnod[0] = this->_nset->
                operator[](_jpos).node( 0) ;

            return _inod[0] == _jnod[0]  ;
        }
        } ;

    template <
        typename  edge_list
             >
    class edge_pred
        {
        public  :
        __const_ptr(edge_list) _eset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call edge_pred  (
            edge_list *_esrc
            ) : _eset( _esrc) {}
    /*----------------------- compute "equal-to" for edge */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[2];
            _inod[0] = this->_eset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_eset->
                operator[](_ipos).node( 1) ;

            uint32_t  _jnod[2];
            _jnod[0] = this->_eset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_eset->
                operator[](_jpos).node( 1) ;

            algorithms::isort (
                &_inod[0], &_inod[2],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[2],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1]  ;
        }
        } ;

    template <
        typename  tri3_list
             >
    class tri3_pred
        {
        public  :
        __const_ptr(tri3_list) _tset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call tri3_pred  (
            tri3_list *_tsrc
            ) : _tset( _tsrc) {}
    /*----------------------- compute "equal-to" for face */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[3];
            _inod[0] = this->_tset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_tset->
                operator[](_ipos).node( 1) ;
            _inod[2] = this->_tset->
                operator[](_ipos).node( 2) ;

            uint32_t  _jnod[3];
            _jnod[0] = this->_tset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_tset->
                operator[](_jpos).node( 1) ;
            _jnod[2] = this->_tset->
                operator[](_jpos).node( 2) ;

            algorithms::isort (
                &_inod[0], &_inod[3],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[3],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1] &&
                   _inod[2] == _jnod[2]  ;
        }
        } ;

    template <
        typename  quad_list
             >
    class quad_pred
        {
        public  :
        __const_ptr(quad_list) _qset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call quad_pred  (
            quad_list *_qsrc
            ) : _qset( _qsrc) {}
    /*----------------------- compute "equal-to" for face */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[4];
            _inod[0] = this->_qset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_qset->
                operator[](_ipos).node( 1) ;
            _inod[2] = this->_qset->
                operator[](_ipos).node( 2) ;
            _inod[3] = this->_qset->
                operator[](_ipos).node( 3) ;

            uint32_t  _jnod[4];
            _jnod[0] = this->_qset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_qset->
                operator[](_jpos).node( 1) ;
            _jnod[2] = this->_qset->
                operator[](_jpos).node( 2) ;
            _jnod[3] = this->_qset->
                operator[](_jpos).node( 3) ;

            algorithms::isort (
                &_inod[0], &_inod[4],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[4],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1] &&
                   _inod[2] == _jnod[2] &&
                   _inod[3] == _jnod[3]  ;
        }
        } ;

    template <
        typename  tri4_list
             >
    class tri4_pred
        {
        public  :
        __const_ptr(tri4_list) _tset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call tri4_pred  (
            tri4_list *_tsrc
            ) : _tset( _tsrc) {}
    /*----------------------- compute "equal-to" for face */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[4];
            _inod[0] = this->_tset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_tset->
                operator[](_ipos).node( 1) ;
            _inod[2] = this->_tset->
                operator[](_ipos).node( 2) ;
            _inod[3] = this->_tset->
                operator[](_ipos).node( 3) ;

            uint32_t  _jnod[4];
            _jnod[0] = this->_tset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_tset->
                operator[](_jpos).node( 1) ;
            _jnod[2] = this->_tset->
                operator[](_jpos).node( 2) ;
            _jnod[3] = this->_tset->
                operator[](_jpos).node( 3) ;

            algorithms::isort (
                &_inod[0], &_inod[4],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[4],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1] &&
                   _inod[2] == _jnod[2] &&
                   _inod[3] == _jnod[3]  ;
        }
        } ;

    template <
        typename  hexa_list
             >
    class hexa_pred
        {
        public  :
        __const_ptr(hexa_list) _hset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call hexa_pred  (
            hexa_list *_hsrc
            ) : _hset( _hsrc) {}
    /*----------------------- compute "equal-to" for face */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[8];
            _inod[0] = this->_hset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_hset->
                operator[](_ipos).node( 1) ;
            _inod[2] = this->_hset->
                operator[](_ipos).node( 2) ;
            _inod[3] = this->_hset->
                operator[](_ipos).node( 3) ;
            _inod[4] = this->_hset->
                operator[](_ipos).node( 4) ;
            _inod[5] = this->_hset->
                operator[](_ipos).node( 5) ;
            _inod[6] = this->_hset->
                operator[](_ipos).node( 6) ;
            _inod[7] = this->_hset->
                operator[](_ipos).node( 7) ;

            uint32_t  _jnod[8];
            _jnod[0] = this->_hset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_hset->
                operator[](_jpos).node( 1) ;
            _jnod[2] = this->_hset->
                operator[](_jpos).node( 2) ;
            _jnod[3] = this->_hset->
                operator[](_jpos).node( 3) ;
            _jnod[4] = this->_hset->
                operator[](_jpos).node( 4) ;
            _jnod[5] = this->_hset->
                operator[](_jpos).node( 5) ;
            _jnod[6] = this->_hset->
                operator[](_jpos).node( 6) ;
            _jnod[7] = this->_hset->
                operator[](_jpos).node( 7) ;

            algorithms::isort (
                &_inod[0], &_inod[8],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[8],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1] &&
                   _inod[2] == _jnod[2] &&
                   _inod[3] == _jnod[3] &&
                   _inod[4] == _jnod[4] &&
                   _inod[5] == _jnod[5] &&
                   _inod[6] == _jnod[6] &&
                   _inod[7] == _jnod[7]  ;
        }
        } ;

    template <
        typename  wedg_list
             >
    class wedg_pred
        {
        public  :
        __const_ptr(wedg_list) _wset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call wedg_pred  (
            wedg_list *_wsrc
            ) : _wset( _wsrc) {}
    /*----------------------- compute "equal-to" for face */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[6];
            _inod[0] = this->_wset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_wset->
                operator[](_ipos).node( 1) ;
            _inod[2] = this->_wset->
                operator[](_ipos).node( 2) ;
            _inod[3] = this->_wset->
                operator[](_ipos).node( 3) ;
            _inod[4] = this->_wset->
                operator[](_ipos).node( 4) ;
            _inod[5] = this->_wset->
                operator[](_ipos).node( 5) ;

            uint32_t  _jnod[6];
            _jnod[0] = this->_wset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_wset->
                operator[](_jpos).node( 1) ;
            _jnod[2] = this->_wset->
                operator[](_jpos).node( 2) ;
            _jnod[3] = this->_wset->
                operator[](_jpos).node( 3) ;
            _jnod[4] = this->_wset->
                operator[](_jpos).node( 4) ;
            _jnod[5] = this->_wset->
                operator[](_jpos).node( 5) ;

            algorithms::isort (
                &_inod[0], &_inod[6],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[6],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1] &&
                   _inod[2] == _jnod[2] &&
                   _inod[3] == _jnod[3] &&
                   _inod[4] == _jnod[4] &&
                   _inod[5] == _jnod[5]  ;
        }
        } ;

    template <
        typename  pyra_list
             >
    class pyra_pred
        {
        public  :
        __const_ptr(pyra_list) _pset ;
        public  :
    /*----------------------- build hash _obj. from _src. */
        __inline_call pyra_pred  (
            pyra_list *_psrc
            ) : _pset( _psrc) {}
    /*----------------------- compute "equal-to" for face */
        __inline_call bool_type operator() (
            uint32_t  _ipos,
            uint32_t  _jpos
            ) const
        {
            uint32_t  _inod[5];
            _inod[0] = this->_pset->
                operator[](_ipos).node( 0) ;
            _inod[1] = this->_pset->
                operator[](_ipos).node( 1) ;
            _inod[2] = this->_pset->
                operator[](_ipos).node( 2) ;
            _inod[3] = this->_pset->
                operator[](_ipos).node( 3) ;
            _inod[4] = this->_pset->
                operator[](_ipos).node( 4) ;

            uint32_t  _jnod[5];
            _jnod[0] = this->_pset->
                operator[](_jpos).node( 0) ;
            _jnod[1] = this->_pset->
                operator[](_jpos).node( 1) ;
            _jnod[2] = this->_pset->
                operator[](_jpos).node( 2) ;
            _jnod[3] = this->_pset->
                operator[](_jpos).node( 3) ;
            _jnod[4] = this->_pset->
                operator[](_jpos).node( 4) ;

            algorithms::isort (
                &_inod[0], &_inod[5],
                    std::less<uint32_t >());

            algorithms::isort (
                &_jnod[0], &_jnod[5],
                    std::less<uint32_t >());

            return _inod[0] == _jnod[0] &&
                   _inod[1] == _jnod[1] &&
                   _inod[2] == _jnod[2] &&
                   _inod[3] == _jnod[3] &&
                   _inod[4] == _jnod[4]  ;
        }
        } ;


    }

#   endif  //__HASH_COMPLEX_k__



