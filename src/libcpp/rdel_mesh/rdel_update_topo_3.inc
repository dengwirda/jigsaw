
    /*
    --------------------------------------------------------
     * RDEL-UPDATE-3: update nodes for topo. checks.
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
     * Last updated: 03 July, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_3.hpp


    /*
    --------------------------------------------------------
     * FILL-TOPO: enqueue node for topo.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type fill_topo (
        mesh_type &_mesh,
        iptr_type  _pass,
        node_heap &_etpq,
        iptr_list &_emrk,
        node_heap &_ftpq,
        iptr_list &_fmrk,
        edat_list &_enew,
        edat_list &_eold,
        fdat_list &_fnew,
        fdat_list &_fold,
        rdel_opts &_opts
        )
    {
        if (_opts.top1())
        {
    /*----------------------------- check nodes adj. to edges */
        if (_eold.count() == +0)
        {
            _emrk.set_count (
                _mesh._tria._nset.count(),
                    containers::loose_alloc , -1) ;

        /*------------------------- add all: everything "new" */
            node_data _ndat ;
            for (auto _iter  = _enew.head();
                      _iter != _enew.tend();
                    ++_iter  )
            {
                _ndat._pass = _pass ;

                if (_emrk[_iter->_node[0]] == -1)
                {
                    _emrk[_iter->_node[0]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[0] ;

                    _etpq.push(_ndat) ;
                }

                if (_emrk[_iter->_node[1]] == -1)
                {
                    _emrk[_iter->_node[1]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[1] ;

                    _etpq.push(_ndat) ;
                }
            }
        }
        else
        {
            _emrk.set_count (
                _mesh._tria._nset.count(),
                    containers::loose_alloc , -1) ;

            {
        /*------------------------- init. for local hash obj. */
            typename
                mesh_type::edge_list _eset (
            typename mesh_type::edge_hash(),
            typename mesh_type::edge_pred(),
               +.8, _mesh._eset.get_alloc()) ;

        /*------------------------- push alloc. for hash obj. */
            _eset.set_slots (
                _enew.count() * 2 ,
                    containers::tight_alloc) ;

        /*------------------------- add "old" absent in "new" */
            node_data _ndat ;
            for (auto _iter  = _enew.head();
                      _iter != _enew.tend();
                    ++_iter  )
            {
                _eset.push(*_iter) ;
            }
            for (auto _iter  = _eold.head();
                      _iter != _eold.tend();
                    ++_iter  )
            {
                _ndat._pass = _pass;
                typename mesh_type::
                         edge_list::
                     item_type *_mptr = nullptr ;
                if(!_eset.find(*_iter, _mptr))
                {

                if (_emrk[_iter->_node[0]] == -1)
                {
                    _emrk[_iter->_node[0]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[0] ;

                    _etpq.push(_ndat) ;
                }

                if (_emrk[_iter->_node[1]] == -1)
                {
                    _emrk[_iter->_node[1]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[1] ;

                    _etpq.push(_ndat) ;
                }

                }
            }
            }

            {
        /*------------------------- init. for local hash obj. */
            typename
                mesh_type::edge_list _eset (
            typename mesh_type::edge_hash(),
            typename mesh_type::edge_pred(),
               +.8, _mesh._eset.get_alloc()) ;

        /*------------------------- push alloc. for hash obj. */
            _eset.set_slots (
                _eold.count() * 2 ,
                    containers::tight_alloc) ;

        /*------------------------- add "new" absent in "old" */
            node_data _ndat ;
            for (auto _iter  = _eold.head();
                      _iter != _eold.tend();
                    ++_iter  )
            {
                _eset.push(*_iter) ;
            }
            for (auto _iter  = _enew.head();
                      _iter != _enew.tend();
                    ++_iter  )
            {
                _ndat._pass = _pass;
                typename mesh_type::
                         edge_list::
                     item_type *_mptr = nullptr ;
                if(!_eset.find(*_iter, _mptr))
                {

                if (_emrk[_iter->_node[0]] == -1)
                {
                    _emrk[_iter->_node[0]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[0] ;

                    _etpq.push(_ndat) ;
                }

                if (_emrk[_iter->_node[1]] == -1)
                {
                    _emrk[_iter->_node[1]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[1] ;

                    _etpq.push(_ndat) ;
                }

                }
            }
            }
        }
        }

        if (_opts.top2())
        {
    /*----------------------------- check nodes adj. to faces */
        if (_fold.count() == +0)
        {
            _fmrk.set_count (
                _mesh._tria._nset.count(),
                    containers::loose_alloc , -1) ;

        /*------------------------- add all: everything "new" */
            node_data _ndat ;
            for (auto _iter  = _fnew.head();
                      _iter != _fnew.tend();
                    ++_iter  )
            {
                _ndat._pass = _pass ;

                if (_fmrk[_iter->_node[0]] == -1)
                {
                    _fmrk[_iter->_node[0]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[0] ;

                    _ftpq.push(_ndat) ;
                }

                if (_fmrk[_iter->_node[1]] == -1)
                {
                    _fmrk[_iter->_node[1]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[1] ;

                    _ftpq.push(_ndat) ;
                }

                if (_fmrk[_iter->_node[2]] == -1)
                {
                    _fmrk[_iter->_node[2]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[2] ;

                    _ftpq.push(_ndat) ;
                }
            }
        }
        else
        {
            _fmrk.set_count (
                _mesh._tria._nset.count(),
                    containers::loose_alloc , -1) ;

            {
        /*------------------------- init. for local hash obj. */
            typename
                mesh_type::face_list _fset (
            typename mesh_type::face_hash(),
            typename mesh_type::face_pred(),
               +.8, _mesh._fset.get_alloc()) ;

        /*------------------------- push alloc. for hash obj. */
            _fset.set_slots (
                _fnew.count() * 2 ,
                    containers::tight_alloc) ;

        /*------------------------- add "old" absent in "new" */
            node_data _ndat ;
            for (auto _iter  = _fnew.head();
                      _iter != _fnew.tend();
                    ++_iter  )
            {
                _fset.push(*_iter) ;
            }
            for (auto _iter  = _fold.head();
                      _iter != _fold.tend();
                    ++_iter  )
            {
                _ndat._pass = _pass;
                typename mesh_type::
                         face_list::
                     item_type *_mptr = nullptr ;
                if(!_fset.find(*_iter, _mptr))
                {

                if (_fmrk[_iter->_node[0]] == -1)
                {
                    _fmrk[_iter->_node[0]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[0] ;

                    _ftpq.push(_ndat) ;
                }

                if (_fmrk[_iter->_node[1]] == -1)
                {
                    _fmrk[_iter->_node[1]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[1] ;

                    _ftpq.push(_ndat) ;
                }

                if (_fmrk[_iter->_node[2]] == -1)
                {
                    _fmrk[_iter->_node[2]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[2] ;

                    _ftpq.push(_ndat) ;
                }

                }
            }
            }

            {
        /*------------------------- init. for local hash obj. */
            typename
                mesh_type::face_list _fset (
            typename mesh_type::face_hash(),
            typename mesh_type::face_pred(),
               +.8, _mesh._fset.get_alloc()) ;

        /*------------------------- push alloc. for hash obj. */
            _fset.set_slots (
                _fold.count() * 2 ,
                    containers::tight_alloc) ;

        /*------------------------- add "new" absent in "old" */
            node_data _ndat ;
            for (auto _iter  = _fold.head();
                      _iter != _fold.tend();
                    ++_iter  )
            {
                _fset.push(*_iter) ;
            }
            for (auto _iter  = _fnew.head();
                      _iter != _fnew.tend();
                    ++_iter  )
            {
                _ndat._pass = _pass;
                typename mesh_type::
                         face_list::
                     item_type *_mptr = nullptr ;
                if(!_fset.find(*_iter, _mptr))
                {

                if (_fmrk[_iter->_node[0]] == -1)
                {
                    _fmrk[_iter->_node[0]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[0] ;

                    _ftpq.push(_ndat) ;
                }

                if (_fmrk[_iter->_node[1]] == -1)
                {
                    _fmrk[_iter->_node[1]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[1] ;

                    _ftpq.push(_ndat) ;
                }

                if (_fmrk[_iter->_node[2]] == -1)
                {
                    _fmrk[_iter->_node[2]] =  +1;

                    _ndat. _node[0] =
                    _iter->_node[2] ;

                    _ftpq.push(_ndat) ;
                }

                }
            }
            }
        }
        }
    }



