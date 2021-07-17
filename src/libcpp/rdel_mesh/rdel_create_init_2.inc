
    /*
    --------------------------------------------------------
     * RDEL-INITIAL-MESH-2: set initial conditions in R^2.
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
     * Last updated: 15 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from rdel_mesh_2.hpp


    /*
    --------------------------------------------------------
     * INIT-SORT: create an insertion order.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_sort (
        init_type &_init,
        iptr_list &_iset
        )
    {
        typedef geom_tree::aabb_node_base_k
                           tree_node ;

        typedef geom_tree::
                    aabb_item_node_k <
            real_type,
            iptr_type, 2>  tree_item ;

        typedef geom_tree::aabb_tree <
            tree_item, 2,
            tree_node,
            allocator   >  tree_type ;

        containers::array<tree_item> _bbox;

    /*------------------------------ initialise aabb-tree */
        iptr_type _npos  = 0 ;
        tree_type _tree  ;
        for (auto _node  =
            _init._mesh.node().head();
                  _node !=
            _init._mesh.node().tend();
                ++_node, ++_npos)
        {
            if (_node->mark() >= +0)
            {
            _bbox.push_tail() ;
            _bbox.tail()->
               pval(0) = _node->pval(0) ;
            _bbox.tail()->
               pval(1) = _node->pval(1) ;

            _bbox.tail()->
               ipos () = _npos ;
            }
        }

        iptr_type _NBOX =
            (iptr_type) std::pow (8, 2) ;       // 8^ndim

        _tree.load(_bbox.head(),
                   _bbox.tend(), _NBOX) ;

    /*------------------------------ randomised tree sort */
        _tree.brio(_iset) ;
    }

    /*
    --------------------------------------------------------
     * HAVE-EDGE: TRUE if exists in DT.
    --------------------------------------------------------
     */

    __static_call
    __normal_call bool_type have_edge (
        mesh_type &_mesh,
        iptr_list &_imap,
        char_type &_eadj,
        iptr_type &_tadj,
        iptr_type *_enod
        )
    {
        class node_pred
        {
    /*--------------------- find adj. set of tria-to-node */
        public  :

        typedef typename
            mesh_type::tria_type tria_type ;

        public  :
            iptr_type _npos;

        public  :
        __inline_call node_pred (
            iptr_type _nsrc
            ) : _npos(_nsrc) {}
    /*--------------------- find adj. set of tria-to-node */
        __inline_call bool_type operator() (
            tria_type&_tria,
            iptr_type _tpos,
            iptr_type _fpos
            )
        {
            iptr_type _tnod[3] = {
            _tria.tria(_tpos)->node(0) ,
            _tria.tria(_tpos)->node(1) ,
            _tria.tria(_tpos)->node(2)
                } ;

            __unreferenced(_fpos) ;

            if (_tnod[0] == this->_npos)
            return  true ;
            if (_tnod[1] == this->_npos)
            return  true ;
            if (_tnod[2] == this->_npos)
            return  true ;

            return false ;
        }
        } ;

        if (_imap[_enod[0]] != -1 &&
            _imap[_enod[1]] != -1 )
        {
    /*--------------------- TRUE if we find edge in DTRI. */
        iptr_type _emap[2];
        _emap[0] = _imap[_enod[ 0]] ;
        _emap[1] = _imap[_enod[ 1]] ;

    /*--------------------- find all tria's about ENOD[0] */
        iptr_list _tset;
        char_type _epos;
        _tset.set_alloc( +32) ;

        _mesh._tria.walk_node (
            _emap[0],
         node_pred(_emap[0]), _tset);

    /*--------------------- can find exact match to edge? */
        for (auto _tpos  = _tset.head();
                  _tpos != _tset.tend();
                ++_tpos  )
        {
        for (_epos = 3; _epos-- != 0; )
        {
            iptr_type _ENOD [ +3];
            mesh_type::tria_type::
            tria_type::face_node(
                _ENOD, _epos, +2, +1 );

            _ENOD[0] = _mesh.
                _tria  .tria(
               *_tpos)->node(_ENOD[0]);
            _ENOD[1] = _mesh.
                _tria  .tria(
               *_tpos)->node(_ENOD[1]);

            iptr_type _same = +0 ;
            if (_emap[0] == _ENOD[0] ||
                _emap[0] == _ENOD[1] )
                _same += +1 ;
            if (_emap[1] == _ENOD[0] ||
                _emap[1] == _ENOD[1] )
                _same += +1 ;

            if (_same == +2 )
            {
    /*--------------------- return TRUE for edge in DTRI. */
                _eadj = _epos ;

                _tadj =*_tpos ;

                return true ;
            }
        }
        }

        }

        return false  ;
    }

    /*
    --------------------------------------------------------
     * INIT-INIT: push IC nodes + faces.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_init (
        init_type &_init ,
        geom_type &_geom ,
        mesh_type &_mesh ,
        typename
    mesh_type::edge_list & _epro,
        typename
    mesh_type::edge_list & _ebad,
        real_type  _NEAR
        )
    {
    /*------------------------------ form insertion order */
        iptr_type _hint  = -1 ;
        iptr_type _null  = -1 ;
        iptr_list _iset, _imap  ;
        init_sort(_init, _iset) ;

    /*------------------------------ find "central" point */
        iptr_type _imid  = -1 ;
        real_type _dmin  =
            std::numeric_limits
                <real_type>::infinity();

        _imap.set_count (
            _iset.count(),
        containers::tight_alloc, _null);

        real_type _pmid[2] ;
        _pmid[0] = (real_type) +0. ;
        _pmid[1] = (real_type) +0. ;

        for (auto _iter  = _iset.head();
                  _iter != _iset.tend();
                ++_iter  )
        {
             auto _node =
           &_init._mesh. node (*_iter) ;

            _pmid[0] +=
                _node->pval(0);
            _pmid[1] +=
                _node->pval(1);
        }

        _pmid[0] /= _iset.count () ;
        _pmid[1] /= _iset.count () ;

        for (auto _iter  = _iset.head();
                  _iter != _iset.tend();
                ++_iter  )
        {
             auto _node =
           &_init._mesh. node (*_iter) ;

            real_type _dsqr =
            geometry::lensqr_2d(
               &_node->pval(0), _pmid) ;

            if (_dsqr < _dmin)
            {
                _dmin = _dsqr ;

                _imid =*_iter ;
            }
        }

    /*------------------------------ seed node from init. */
        iptr_type _base =
            (iptr_type)
            _mesh._tria._nset.count () ;

        if (_imid > -1)
        {
             auto _node =
           &_init._mesh. node ( _imid) ;

            iptr_type _npos =
                _base+_imid -  0 ;
            iptr_type _near = -1 ;

            if (_node->itag () >= +0 )
                _geom.projector(
               &_node->pval(0) ,
                _node->itag () ,
               &_node->pval(0) ) ;

            if (_mesh._tria.find_node(
               &_node->pval(0) ,
                _near, _hint ) )
            {
    /*------------------------------ find len. to nearest */
             auto _nptr =
            _mesh._tria. node ( _near) ;

            _hint = _mesh._tria.
                 node(_near)->next() ;

            real_type _dsqr =
                 geometry::lensqr_2d (
               &_node->pval(0) ,
               &_nptr->pval(0) ) ;

            if (_dsqr >= _NEAR )
            {
            if (_mesh._tria.push_node(
               &_node->pval(0) ,
                _npos, _hint ) )
            {
    /*------------------------------ is "new" if distinct */
            _imap[_imid] = _npos ;

            _mesh._tria.node
                (_npos)->fdim()
                    = _node->fdim() ;

            _mesh._tria.node
                (_npos)->feat()
                    = _node->feat() ;

            _mesh._tria.node
                (_npos)->topo() = 2 ;

            _mesh._tria.node
                (_npos)->part()
                    = _node->itag() ;

            _hint        = _npos ;

            }
            }
            else
            {
    /*------------------------------ else move to nearest */
            _imap[_imid] = _near ;

            _mesh._tria.node
                (_near)->part()
                    = _node->itag() ;

            _hint        = _near ;

            }

            }
        }

    /*------------------------------ seed node from init. */
        for (auto _iter  = _iset.head();
                  _iter != _iset.tend();
                ++_iter  )
        {
            if (*_iter == _imid)
                continue ;

             auto _inum =*_iter;
             auto _node =
           &_init._mesh . node(*_iter) ;

            iptr_type _npos =
               *_iter+_base -  0 ;
            iptr_type _near = -1 ;

            if (_node->itag () >= +0 )
                _geom.projector(
               &_node->pval(0) ,
                _node->itag () ,
               &_node->pval(0) ) ;

            if (_mesh._tria.find_node(
               &_node->pval(0) ,
                _near, _hint ) )
            {
    /*------------------------------ find len. to nearest */
             auto _nptr =
            _mesh._tria. node ( _near) ;

            _hint = _mesh._tria.
                 node(_near)->next() ;

            real_type _dsqr =
                 geometry::lensqr_2d (
               &_node->pval(0) ,
               &_nptr->pval(0) ) ;

            if (_dsqr >= _NEAR )
            {
            if (_mesh._tria.push_node(
               &_node->pval(0) ,
                _npos, _hint ) )
            {
    /*------------------------------ is "new" if distinct */
            _imap[_inum] = _npos ;

            _mesh._tria.node
                (_npos)->fdim()
                    = _node->fdim() ;

            _mesh._tria.node
                (_npos)->feat()
                    = _node->feat() ;

            _mesh._tria.node
                (_npos)->topo() = 2 ;

            _mesh._tria.node
                (_npos)->part()
                    = _node->itag() ;

            _hint        = _npos ;

            }
            }
            else
            {
    /*------------------------------ else move to nearest */
            _imap[_inum] = _near ;

            _mesh._tria.node
                (_near)->part()
                    = _node->itag() ;

            _hint        = _near ;

            }

            }
        }

    /*------------------------------ seed edge from init. */
        for (auto _iter  =
            _init._mesh.edge().head() ;
                  _iter !=
            _init._mesh.edge().tend() ;
                ++_iter  )
        {
            if (_iter->mark() >= +0)
            if (_iter->itag() <= -1)
            {
    /*------------------------------ find edge "in" DTRI. */
            char_type _eadj ;
            iptr_type _tadj ;
            if (have_edge( _mesh,
               _imap,_eadj,_tadj,
              &_iter->node( +0 )) )
            {
    /*------------------------------ mark as un-refinable */
            edge_data  _edat;
            _edat._node[0] =
            _imap[_iter->node( 0 )] ;
            _edat._node[1] =
            _imap[_iter->node( 1 )] ;

            _edat._pass    =     +0 ;
            _edat._kind    = 
                     mesh::good_item;
            _edat._part    =
                    _iter->itag() ;

            _edat._eadj    = _eadj  ;
            _edat._tadj    = _tadj  ;

            algorithms::isort (
            _edat._node.head(),
            _edat._node.tend(),
            std::less<iptr_type>()) ;

            _epro.push     ( _edat) ;

            _mesh.push_edge( _edat) ;

            }
            else
            {
    /*------------------------------ edge is non-Delaunay */
            edge_data  _edat;
            _edat._node[0] =
                    _iter->node( 0) ;
            _edat._node[1] =
                    _iter->node( 1) ;

            _edat._part    =
                    _iter->itag() ;

            _ebad.push     ( _edat) ;

            }

            }
        }

    }



