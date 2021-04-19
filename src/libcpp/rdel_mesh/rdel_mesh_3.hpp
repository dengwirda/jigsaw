
    /*
    --------------------------------------------------------
     * RDEL-MESH-3: restricted delaunay mesh-gen. in R^3.
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
     * Last updated: 02 Feb., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     *
     * This class defines the "restricted" delaunay mesh
     * generation scheme for domains in R^3. This is the
     * top-level class, implementing the main method loops:
     * "sampling" the geometry and incrementally refining
     * to convergence. In brief: a set of priority queues
     * are maintained for the nodes, edges, faces and cells
     * in the rDT, with new refinement points inserted to
     * eliminate any remaining "bad" simplexes.
     *
     * The algorithm is parameterised by various templated
     * types: MESH-TYPE which holds the rDT, MESH-PRED that
     * defines the refinement strategy, GEOM-TYPE which
     * represents the domain geometry, and HFUN-TYPE which
     * represents the "mesh-spacing" function h(x). These
     * predicates derive from base classes defining common
     * functionality.
     *
     * My implementation is described in:
     *
     * D. Engwirda and D. Ivers, (2016): Off-centre Steiner
     * points for Delaunay-refinement on curved surfaces,
     * Computer-Aided Design, 72, pp. 157-171,
     * http://dx.doi.org/10.1016/j.cad.2015.10.007
     *
     * D. Engwirda, (2016): "Conforming restricted Delaunay
     * mesh generation for piecewise smooth complexes",
     * Procedia Engineering, 163, pp. 84-96,
     * http://dx.doi.org/10.1016/j.proeng.2016.11.024
     *
     * D. Engwirda, (2014): "Locally-optimal Delaunay-
     * refinement and optimisation-based mesh generation",
     * Ph.D. Thesis, School of Mathematics and Statistics,
     * Univ. of Sydney.
     * http://hdl.handle.net/2123/13148
     *
     * building on various previous works on rDT methods,
     * including (primarily):
     *
     * J.D. Boissonnat, S. Oudot, (2005): "Provably Good
     * Sampling and Meshing of Surfaces", Graphical Models,
     * 67, pp. 405-451,
     * https://doi.org/10.1016/j.gmod.2005.01.004
     *
     * C. Jamin, P. Alliez, M. Yvinec, and J.D. Boissonnat,
     * (2015): "CGALmesh: a generic framework for Delaunay
     * mesh generation", ACM Transactions on Mathematical
     * Software (TOMS), 41, pp. 23
     * https://doi.org/10.1145/2699463
     *
     * L. Rineau, M. Yvinec, (2008): "Meshing 3D Domains
     * Bounded by Piecewise Smooth Surfaces", Proc. of the
     * 16th International Meshing Roundtable, pp. 443-460,
     * https://doi.org/10.1007/978-3-540-75103-8_25
     *
     * S.W. Cheng, T.K. Dey, E.A. Ramos, (2010): "Delaunay
     * Refinement for Piecewise Smooth Complexes",
     * Discrete & Computational Geometry, 43, pp. 121-166,
     * https://doi.org/10.1007/s00454-008-9109-3
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_MESH_3__
#   define __RDEL_MESH_3__

    namespace mesh {

    template <
    typename M ,
    typename P ,
    typename G ,
    typename H ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class rdel_mesh_3d
    {
    public  :

    /*-------- restricted delaunay mesh-generation in R^3 */

    typedef M                           mesh_type ;
    typedef P                           mesh_pred ;
    typedef G                           geom_type ;
    typedef H                           hfun_type ;
    typedef I                           init_type ;
    typedef A                           allocator ;

    typedef typename
            mesh_type::real_type        real_type ;
    typedef typename
            mesh_type::iptr_type        iptr_type ;

    typedef typename
            allocator::size_type        uint_type ;

    typedef char_type  mode_type ;

    char_type static constexpr null_mode = +0 ;
    char_type static constexpr node_mode = +1 ;
    char_type static constexpr edge_mode = +2 ;
    char_type static constexpr etop_mode = +3 ;
    char_type static constexpr face_mode = +4 ;
    char_type static constexpr ftop_mode = +5 ;
    char_type static constexpr tria_mode = +6 ;

    class node_pred ; class ball_pred ;
    class edge_pred ; class face_pred ;
    class tria_pred ;

    class edge_cost : public mesh_pred::edge_data
        {
        public  :
        iptr_type           _node[ +2] ;
        iptr_type           _pass;
        } ;
    class face_cost : public mesh_pred::face_data
        {
        public  :
        iptr_type           _node[ +3] ;
        iptr_type           _pass;
        } ;
    class tria_cost : public mesh_pred::tria_data
        {
        public  :
        iptr_type           _node[ +4] ;
        iptr_type           _pass;
        } ;

    typedef typename
            mesh_type::node_data        node_data ;
    typedef typename
            mesh_type::ball_data        ball_data ;
    typedef typename
            mesh_type::edge_data        edge_data ;
    typedef typename
            mesh_type::face_data        face_data ;
    typedef typename
            mesh_type::tria_data        tria_data ;

    typedef mesh::mesh_params       <
                real_type,
                iptr_type           >   rdel_opts ;

    typedef mesh::rdel_timers       <
                real_type ,
                iptr_type           >   rdel_stat ;

    typedef containers::array       <
                iptr_type           >   iptr_list ;

/*------------------------------------------ cavity lists */
    typedef containers::array       <
                edge_data           >   edat_list ;

    typedef containers::array       <
                edge_cost           >   escr_list ;

    typedef containers::array       <
                face_data           >   fdat_list ;

    typedef containers::array       <
                face_cost           >   fscr_list ;

    typedef containers::array       <
                tria_data           >   tdat_list ;

    typedef containers::array       <
                tria_cost           >   tscr_list ;

/*------------------------------------------ refine lists */
    typedef containers::priorityset <
                node_data,
                node_pred           >   node_heap ;

    typedef containers::priorityset <
                edge_cost,
                edge_pred           >   edge_heap ;

    typedef containers::priorityset <
                face_cost,
                face_pred           >   face_heap ;

    typedef containers::priorityset <
                tria_cost,
                tria_pred           >   tria_heap ;

/*------------------------------------------ collar lists */
    typedef containers::array       <
                ball_data           >   ball_list ;

    typedef containers::priorityset <
                ball_data,
                ball_pred           >   ball_heap ;

    class node_pred
        {
    /*---------------------- "less-than" for node objects */
        public  :
        __inline_call bool_type operator() (
            node_data const& _idat,
            node_data const& _jdat
            ) const
        {   return _idat._node[0] <
                   _jdat._node[0] ;
        }
        } ;
    class ball_pred
        {
    /*---------------------- "less-than" for ball objects */
        public  :
        __inline_call bool_type operator() (
            ball_data const& _idat,
            ball_data const& _jdat
            ) const
        {   return _idat._ball[3] >
                   _jdat._ball[3] ;
        }
        } ;
    class edge_pred
        {
    /*---------------------- "less-than" for edge objects */
        public  :
        __inline_call bool_type operator() (
            edge_cost const& _idat,
            edge_cost const& _jdat
            ) const
        {   return mesh_pred
                 ::edge_pred(_idat, _jdat) ;
        }
        } ;
    class face_pred
        {
    /*---------------------- "less-than" for face objects */
        public  :
        __inline_call bool_type operator() (
            face_cost const& _idat,
            face_cost const& _jdat
            ) const
        {   return mesh_pred
                 ::face_pred(_idat, _jdat) ;
        }
        } ;
    class tria_pred
        {
    /*---------------------- "less-than" for tria objects */
        public  :
        __inline_call bool_type operator() (
            tria_cost const& _idat,
            tria_cost const& _jdat
            ) const
        {   return mesh_pred
                 ::tria_pred(_idat, _jdat) ;
        }
        } ;


    /*
    --------------------------------------------------------
     * RDEL-CREATE: set-up initial nodes/faces.
    --------------------------------------------------------
     */

    #include "rdel_create_init_3.inc"


    /*
    --------------------------------------------------------
     * RDEL-UPDATE: update the restricted-tria.
    --------------------------------------------------------
     */

    #include "rdel_update_face_3.inc"
    #include "rdel_update_topo_3.inc"


    /*
    --------------------------------------------------------
     * RDEL-BOUNDS: test restricted boundaries.
    --------------------------------------------------------
     */

    #include "rdel_test_bounds_3.inc"


    /*
    --------------------------------------------------------
     * RDEL-REFINE: refine restricted subfaces.
    --------------------------------------------------------
     */

    #include "rdel_refine_base_3.inc"
    #include "rdel_refine_ball_3.inc"
    #include "rdel_refine_face_3.inc"
    #include "rdel_refine_topo_3.inc"


    /*
    --------------------------------------------------------
     * TRIM-LIST: prune null faces from queues.
    --------------------------------------------------------
     */

    template <
    typename      list_type
             >
    __static_call
    __inline_call void_type trim_list (
        list_type &_list
        )
    {
        typedef typename
            list_type::
                size_type size_type;

        size_type _amin = +512;
        size_type _alim = +256;
        size_type _amax =
       (size_type)+3 * _list.count() ;
        size_type _anew =
       (size_type)+2 * _list.count() ;

        _anew =
            std::max(_alim, _anew) ;

        if (_list.alloc() > _amin)
        if (_list.alloc() > _amax)
            _list.set_alloc(_anew) ;
    }

    /*
    --------------------------------------------------------
     * TRIM-EEPQ: prune null edges from queues.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type trim_eepq (
        mesh_type &_mesh,
        edge_heap &_eepq
        )
    {
        if(!_eepq.empty())
        {

        iptr_type _dead = +0 ;
        iptr_type _okay = +0 ;

        for (auto _hpos = _eepq.count() - 1 ;
                  _hpos > +0 ;
                --_hpos )
        {
            if (_okay > +1024 &&
                    _dead < _okay / 4) break;

            iptr_type _pass;
            edge_data _edat;
            _edat._node[0] =
            _eepq. peek(_hpos)._node[0];
            _edat._node[1] =
            _eepq. peek(_hpos)._node[1];

            _pass =
            _eepq. peek(_hpos)._pass;

            typename mesh_type::
                     edge_list::
                 item_type *_eptr = nullptr ;
            if(!_mesh.find_edge(_edat,_eptr))
            {
                _eepq._pop(_hpos);

                _dead += +1 ;
            }
            else
            if (_eptr->_data._pass != _pass )
            {
                _eepq._pop(_hpos);

                _dead += +1 ;
            }
            else
            {
                _okay += +1 ;
            }
        }
        }

        trim_list ( _eepq ) ;
    }

    /*
    --------------------------------------------------------
     * TRIM-FFPQ: prune null faces from queues.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type trim_ffpq (
        mesh_type &_mesh,
        face_heap &_ffpq
        )
    {
        if(!_ffpq.empty())
        {

        iptr_type _dead = +0 ;
        iptr_type _okay = +0 ;

        for (auto _hpos = _ffpq.count() - 1 ;
                  _hpos > +0 ;
                --_hpos )
        {
            if (_okay > +1024 &&
                    _dead < _okay / 4) break;

            iptr_type _pass ;
            face_data _fdat;
            _fdat._node[0] =
            _ffpq. peek(_hpos)._node[0];
            _fdat._node[1] =
            _ffpq. peek(_hpos)._node[1];
            _fdat._node[2] =
            _ffpq. peek(_hpos)._node[2];

            _pass =
            _ffpq. peek(_hpos)._pass;

            typename mesh_type::
                     face_list::
                 item_type *_fptr = nullptr ;
            if(!_mesh.find_face(_fdat,_fptr))
            {
                _ffpq._pop(_hpos);

                _dead += +1 ;
            }
            else
            if (_fptr->_data._pass != _pass )
            {
                _ffpq._pop(_hpos);

                _dead += +1 ;
            }
            else
            {
                _okay += +1 ;
            }
        }
        }

        trim_list ( _ffpq ) ;
    }

    /*
    --------------------------------------------------------
     * TRIM-TTPQ: prune null tria. from queues.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type trim_ttpq (
        mesh_type &_mesh,
        tria_heap &_ttpq
        )
    {
        if(!_ttpq.empty())
        {

        iptr_type _dead = +0 ;
        iptr_type _okay = +0 ;

        for (auto _hpos = _ttpq.count() - 1 ;
                  _hpos > +0 ;
                --_hpos )
        {
            if (_okay > +1024 &&
                    _dead < _okay / 4) break;

            iptr_type _pass;
            tria_data _tdat;
            _tdat._node[0] =
            _ttpq. peek(_hpos)._node[0];
            _tdat._node[1] =
            _ttpq. peek(_hpos)._node[1];
            _tdat._node[2] =
            _ttpq. peek(_hpos)._node[2];
            _tdat._node[3] =
            _ttpq. peek(_hpos)._node[3];

            _pass =
            _ttpq. peek(_hpos)._pass;

            typename mesh_type::
                     tria_list::
                 item_type *_tptr = nullptr ;
            if(!_mesh.find_tria(_tdat,_tptr))
            {
                _ttpq._pop(_hpos);

                _dead += +1 ;
            }
            else
            if (_tptr->_data._pass != _pass )
            {
                _ttpq._pop(_hpos);

                _dead += +1 ;
            }
            else
            {
                _okay += +1 ;
            }
        }
        }

        trim_list ( _ttpq ) ;
    }

    /*
    --------------------------------------------------------
     * INIT-RDEL: init. face-set in rDT.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_rdel (
        geom_type &_geom ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        bool_type  _init ,
        iptr_list &_nnew ,
        iptr_list &_tnew ,
        edat_list &_ecav ,
        escr_list &_escr ,
        fdat_list &_fcav ,
        fscr_list &_fscr ,
        tdat_list &_tcav ,
        tscr_list &_tscr ,
        ball_list &_bcav ,
        ball_list &_bscr ,
        iptr_type  _pass ,
        mode_type  _fdim ,
        rdel_opts &_args
        )
    {
    /*-------------------- mark all existing elem. as new */
        iptr_type _npos = 0, _tpos = 0 ;
        for (auto _iter  =
            _mesh._tria._nset.head() ;
                  _iter !=
            _mesh._tria._nset.tend() ;
                ++_iter, ++_npos)
        {
            if (_iter->mark() >= +0)
            {
            _nnew. push_tail( _npos) ;
            }
        }
        for (auto _iter  =
            _mesh._tria._tset.head() ;
                  _iter !=
            _mesh._tria._tset.tend() ;
                ++_iter, ++_tpos)
        {
            if (_iter->mark() >= +0)
            {
            _tnew. push_tail( _tpos) ;
            }
        }
    /*-------------------- init. restricted triangulation */
        push_rdel( _geom, _hfun,
            _mesh, _init,
            _nnew, _tnew,
            _escr, _ecav,
            _fscr, _fcav,
            _tscr, _tcav,
            _bscr, _bcav,
               -1, _pass,
            _fdim, _fdim, _args) ;
    }

    /*
    --------------------------------------------------------
     * INIT-MESH: init. the bounding DT.
    --------------------------------------------------------
     */

    __static_call
    __normal_call void_type init_mesh (
        geom_type &_geom,
        init_type &_init,
        hfun_type &_hfun,
        mesh_type &_mesh,
        typename
    mesh_type::edge_list & _epro,
        typename
    mesh_type::edge_list & _ebad,
        typename
    mesh_type::face_list & _fpro,
        typename
    mesh_type::face_list & _fbad,
        rdel_opts &_opts
        )
    {
        __unreferenced(_hfun) ;

    /*------------------------------ initialise mesh bbox */
        real_type _pmin[ 3] ;
        real_type _pmax[ 3] ;
        _pmin[ 0] = _geom._bmin[ 0] ;
        _pmin[ 1] = _geom._bmin[ 1] ;
        _pmin[ 2] = _geom._bmin[ 2] ;
        _pmax[ 0] = _geom._bmax[ 0] ;
        _pmax[ 1] = _geom._bmax[ 1] ;
        _pmax[ 2] = _geom._bmax[ 2] ;

        for (auto _node  =
            _init._mesh.node().head() ;
                  _node !=
            _init._mesh.node().tend() ;
                ++_node  )
        {
        if (_node->mark() >= + 0 )
        {
        _pmin[ 0] = std::min(
        _pmin[ 0], _node->pval(0)) ;
        _pmax[ 0] = std::max(
        _pmax[ 0], _node->pval(0)) ;

        _pmin[ 1] = std::min(
        _pmin[ 1], _node->pval(1)) ;
        _pmax[ 1] = std::max(
        _pmax[ 1], _node->pval(1)) ;

        _pmin[ 2] = std::min(
        _pmin[ 2], _node->pval(2)) ;
        _pmax[ 2] = std::max(
        _pmax[ 2], _node->pval(2)) ;
        }
        }

        real_type _plen[ 3] = {
        _pmax[ 0] - _pmin[ 0] ,
        _pmax[ 1] - _pmin[ 1] ,
        _pmax[ 2] - _pmin[ 2] }  ;

        real_type _scal =
                    (real_type)+0.0 ;
        _scal =  std::max(
            _scal , _plen[ 0]);
        _scal =  std::max(
            _scal , _plen[ 1]);
        _scal =  std::max(
            _scal , _plen[ 2]);

        _plen[ 0]*= (real_type)+4.0 ;
        _plen[ 1]*= (real_type)+4.0 ;
        _plen[ 2]*= (real_type)+4.0 ;

        _pmin[ 0]-= _plen[ 0] ;
        _pmin[ 1]-= _plen[ 1] ;
        _pmin[ 2]-= _plen[ 2] ;
        _pmax[ 0]+= _plen[ 0] ;
        _pmax[ 1]+= _plen[ 1] ;
        _pmax[ 2]+= _plen[ 2] ;

        _mesh.
        _tria.push_root(_pmin, _pmax) ;

    /*------------------------------ initialise mesh root */
        _mesh.
        _tria.node(+0)->fdim() = +4 ;
        _mesh.
        _tria.node(+1)->fdim() = +4 ;
        _mesh.
        _tria.node(+2)->fdim() = +4 ;
        _mesh.
        _tria.node(+3)->fdim() = +4 ;

        _mesh.
        _tria.node(+0)->feat() = +0 ;
        _mesh.
        _tria.node(+1)->feat() = +0 ;
        _mesh.
        _tria.node(+2)->feat() = +0 ;
        _mesh.
        _tria.node(+3)->feat() = +0 ;

        _mesh.
        _tria.node(+0)->topo() = +0 ;
        _mesh.
        _tria.node(+1)->topo() = +0 ;
        _mesh.
        _tria.node(+2)->topo() = +0 ;
        _mesh.
        _tria.node(+3)->topo() = +0 ;

    /*------------------------------ seed feat from geom. */
        _geom.seed_root(_mesh, _opts) ;
        _geom.seed_feat(_mesh, _opts) ;

    /*------------------------------ seed node from init. */
         real_type _NEAR =
            _scal *_opts.near() *
            _scal *_opts.near() ;

         init_init(_init, _geom ,
                   _mesh,
                   _epro, _ebad ,
                   _fpro, _fbad ,
                   _NEAR) ;

    /*------------------------------ seed mesh from geom. */
        _geom.seed_mesh(_mesh, _opts) ;
    }




    __static_call
    __normal_call void_type init_dual (
        mesh_type &_mesh
        )
    {

      //return;

        iptr_list _nsum;
        _nsum.set_count(
            _mesh._tria._nset.count(),
                containers::tight_alloc, +0) ;

        for (auto _iter  =
                  _mesh._eset._lptr.head() ;
                  _iter !=
                  _mesh._eset._lptr.tend() ;
                ++_iter  )
        {
            if ( *_iter == nullptr) continue ;

            for (auto _item  = *_iter ;
                _item != nullptr;
                _item  = _item->_next )
            {
            if (_item->_data._feat == mesh::hard_feat)
            {
             auto _inod = _item->_data._node[0] ;
             auto _jnod = _item->_data._node[1] ;

             auto _iptr =
                & _mesh._tria._nset[_inod] ;
             auto _jptr =
                & _mesh._tria._nset[_jnod] ;

             auto _lsqr = geometry::lensqr_3d (
                & _iptr->pval( 0),
                & _jptr->pval( 0)) ;

            _lsqr *= (real_type) std::pow(2./3., +2) ;

            _iptr->pval(3) += _lsqr ;
            _nsum[_inod]   += +1 ;

            _jptr->pval(3) += _lsqr ;
            _nsum[_jnod]   += +1 ;
            }
            }
        }

        for (auto _iter  = _nsum.head() ,
                  _npos  = +0 ;
                  _iter != _nsum.tend() ;
                ++_iter, ++_npos)
        {
            if ( *_iter != +0 )
            _mesh._tria._nset[_npos].pval(3) /= *_iter ;
        }

    }




    /*
    --------------------------------------------------------
     * RDEL-MESH: build an rDT mesh in R^3.
    --------------------------------------------------------
     */

    template <
    typename      jlog_file
             >
    __static_call
    __normal_call void_type rdel_mesh (
        geom_type &_geom ,
        init_type &_init ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        rdel_opts &_args ,
        jlog_file &_dump
        )
    {
        mode_type _mode = null_mode ;

    /*------------------------------ push log-file header */
        if (_args.verb() >= 0 )
        {
            _dump.push(
    "#------------------------------------------------------------\n"
    "#    |ITER.|      |DEL-1|      |DEL-2|      |DEL-3| \n"
    "#------------------------------------------------------------\n"
                    ) ;
        }

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::time_point  _ttic ;
        typename std ::chrono::
        high_resolution_clock::time_point  _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*------------------------------ ensure deterministic */
        std::srand( +1 ) ;

        rdel_stat _tcpu  ;

    /*------------------------------ init. list workspace */
        iptr_list _nnew, _nold ;
        iptr_list _tnew, _told ;

        iptr_list _emrk, _fmrk ;

        escr_list _escr  ;
        fscr_list _fscr  ;
        tscr_list _tscr  ;
        ball_list _bscr  ;

        edat_list _edat, _eprv ;
        fdat_list _fdat, _fprv ;
        tdat_list _tdat  ;
        ball_list _bdat  ;

    /*------------------------------ refinement p.-queues */
        edge_heap _eepq  ;
        face_heap _ffpq  ;
        tria_heap _ttpq  ;
        ball_heap _nbpq  ;

        node_heap _etpq, _ftpq ;

    /*------------------------------ alloc. for hash obj. */
        _mesh._eset.set_slots (
            _mesh._tria._tset.count()*+6 ,
                containers::tight_alloc) ;

        _mesh._fset.set_slots (
            _mesh._tria._tset.count()*+4 ,
                containers::tight_alloc) ;

        _mesh._tset.set_slots (
            _mesh._tria._tset.count()*+1 ,
                containers::tight_alloc) ;

    /*------------------------------ init. topo hash obj. */
        typename
            mesh_type::edge_list _epro (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(),
            +.8,_mesh._eset.get_alloc()) ;

        typename
            mesh_type::edge_list _ebad (
        typename mesh_type::edge_hash(),
        typename mesh_type::edge_pred(),
            +.8,_mesh._eset.get_alloc()) ;

        typename
            mesh_type::face_list _fpro (
        typename mesh_type::face_hash(),
        typename mesh_type::face_pred(),
            +.8,_mesh._fset.get_alloc()) ;

        typename
            mesh_type::face_list _fbad (
        typename mesh_type::face_hash(),
        typename mesh_type::face_pred(),
            +.8,_mesh._fset.get_alloc()) ;

    /*------------------------------ init. point counters */
        containers:: fixed_array<
            iptr_type,
            rdel_opts::last_kind> _enod;
        _enod.fill( +0 ) ;

        containers:: fixed_array<
            iptr_type,
            rdel_opts::last_kind> _fnod;
        _fnod.fill( +0 ) ;

        containers:: fixed_array<
            iptr_type,
            rdel_opts::last_kind> _tnod;
        _tnod.fill( +0 ) ;

    /*------------------------------ initialise mesh obj. */
    #   ifdef  __use_timers
        _ttic = _time.now() ;
    #   endif//__use_timers

        init_mesh( _geom , _init, _hfun,
            _mesh, _epro ,
            _ebad, _fpro ,
            _fbad, _args )  ;

        if (_ebad.count() > +0)
        {
            _dump.push("\n");
            _dump.push(
    "**input warning: initial edge(s) not recovered! \n") ;
            _dump.push("\n");
        }
        if (_fbad.count() > +0)
        {
            _dump.push("\n");
            _dump.push(
    "**input warning: initial face(s) not recovered! \n") ;
            _dump.push("\n");
        }

    #   ifdef  __use_timers
        _ttoc = _time.now() ;
        _tcpu._mesh_seed +=
            _tcpu.time_span(_ttic,_ttoc) ;
    #   endif//__use_timers

    /*------------------------------ calc. hfun. at seeds */
        for (auto _node  =
            _mesh._tria._nset.head() ;
                  _node !=
            _mesh._tria._nset.tend() ;
                ++_node  )
        {
            if (_node->mark() >= +0)
            {
                _node->idxh()  =
                 hfun_type::null_hint () ;
            }
        }

    /*-------------------- main: refine edges/faces/trias */
        iptr_type _pass =   +0 ;
        uint_type _Nbpq =   +0 ;
        uint_type _Nepq =   +0 ;
        uint_type _Nfpq =   +0 ;
        uint_type _Ntpq =   +0 ;

        for(bool_type _done=false; !_done ; )
        {
            iptr_type _trim_freq = +10000 ;
    #       ifdef _DEBUG
            iptr_type _jlog_freq = +250 ;
    #       else
            iptr_type _jlog_freq = +12500 ;
    #       endif

            if(++_pass>_args.iter()) break;

        /*------------------------- init. array workspace */

            _Nbpq = std::max(
                _Nbpq, _nbpq.alloc()) ;
            _Nepq = std::max(
                _Nepq, _eepq.alloc()) ;
            _Nfpq = std::max(
                _Nfpq, _ffpq.alloc()) ;
            _Ntpq = std::max(
                _Ntpq, _ttpq.alloc()) ;

            _nnew.set_count(  +0 ) ; // del-tri idx lists
            _nold.set_count(  +0 ) ;
            _tnew.set_count(  +0 ) ;
            _told.set_count(  +0 ) ;

            _escr.set_count(  +0 ) ; // face "cost" lists
            _fscr.set_count(  +0 ) ;
            _tscr.set_count(  +0 ) ;
            _bscr.set_count(  +0 ) ;

            _eprv.set_count(  +0 ) ; // "old" edge in cav
            _edat.set_count(  +0 ) ; // "new" edge in cav
            _fprv.set_count(  +0 ) ; // "old" face in cav
            _fdat.set_count(  +0 ) ; // "new" face in cav
            _tdat.set_count(  +0 ) ; // "new" cell in cav
            _bdat.set_count(  +0 ) ; // "new" ball in cav

        /*--------- calc. "restricted-ness" incrementally */

            bool_type _irDT = false ;

            if (_mode == null_mode )
            {
        /*------------------------- init. protecting ball */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _mode  = node_mode;

                _irDT  = true;   // init. new face in rDT

                init_rdel( _geom, _hfun,
                    _mesh, false,
                    _nnew, _tnew,
                    _edat, _escr,
                    _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr, _pass,
                    _mode, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._node_init +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }

            if (_mode == node_mode &&
                    _nbpq. empty() &&
                    _bscr. empty() &&
                    _bdat. empty() )
            {
        /*------------------------- init. restricted edge */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _mode  = edge_mode;

                if(_geom.have_feat(+1))
                {
                _irDT  = true;   // init. new face in rDT

                init_ball( _geom, _hfun,
                    _mesh,
                    _epro, _fpro, _pass,
                    _mode, _args) ;

                init_rdel( _geom, _hfun,
                    _mesh,  true,// init. circum. for rDT
                    _nnew, _tnew,
                    _edat, _escr,
                    _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr, _pass,
                    _mode, _args) ;
                }

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._edge_init +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            if (_mode == edge_mode &&
                    _eepq. empty() &&
                    _escr. empty() &&
                    _edat. empty() )
            {
        /*------------------------- init. restricted topo */
                _mode  = etop_mode ;
            }

            if (_mode == etop_mode &&
                    _etpq. empty() &&
                    _escr. empty() &&
                    _edat. empty() )
            {
        /*------------------------- init. restricted face */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _mode  = face_mode;

                if(_geom.have_feat(+2))
                {
                _irDT  = true;   // init. new face in rDT

                /*
                init_disc( _geom, _hfun,
                    _mesh,
                    _epro, _fpro, _pass,
                    _mode, _args) ;
                */

              //init_dual( _mesh) ;

                init_rdel( _geom, _hfun,
                    _mesh,  true,
                    _nnew, _tnew,
                    _edat, _escr,
                    _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr, _pass,
                    _mode, _args) ;
                }

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._face_init +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            if (_mode == face_mode &&
                    _ffpq. empty() &&
                    _fscr. empty() &&
                    _fdat. empty() )
            {
        /*------------------------- init. restricted topo */
                _mode  = ftop_mode ;

            }

            if (_mode == ftop_mode &&
                    _ftpq. empty() &&
                    _fscr. empty() &&
                    _fdat. empty() )
            {
        /*------------------------- init. restricted tria */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _mode  = tria_mode;

                if(_geom.have_feat(+2))
                {
                _irDT  = true;   // init. new face in rDT

                init_rdel( _geom, _hfun,
                    _mesh, false,
                    _nnew, _tnew,
                    _edat, _escr,
                    _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr, _pass,
                    _mode, _args) ;
                }

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._tria_init +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }

        /*------------- refine "bad" sub-faces until done */

            if (_irDT == false )
            {

            char_type _tdim = -1;
            typename rdel_opts::node_kind
            _kind =  rdel_opts::null_kind ;

            if (!_nbpq.empty() )
            {
        /*----------------------------- refine "bad" ball */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _kind =_bad_ball( _geom,
                    _hfun, _mesh, _mode,
                    _nnew, _nold,
                    _tnew, _told, _nbpq,
                    _eprv, _edat, _escr,
                    _fprv, _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr,
                    _tdim, _pass, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._node_rule +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            else
            if (!_eepq.empty() )
            {
        /*----------------------------- refine "bad" edge */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _kind =_bad_edge( _geom,
                    _hfun, _mesh, _mode,
                    _epro, _fpro,
                    _nnew, _nold,
                    _tnew, _told, _eepq,
                    _eprv, _edat, _escr,
                    _fprv, _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr,
                    _tdim, _pass, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._edge_rule +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            else
            if (!_etpq.empty() )
            {
        /*----------------------------- refine "bad" topo */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _kind =_bad_etop( _geom,
                    _hfun, _mesh, _mode,
                    _epro, _fpro,
                    _nnew, _nold,
                    _tnew, _told,
                    _etpq, _emrk,
                    _eprv, _edat, _escr,
                    _fprv, _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr,
                    _tdim, _pass, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._edge_rule +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            else
            if (!_ffpq.empty() )
            {
        /*----------------------------- refine "bad" face */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _kind =_bad_face( _geom,
                    _hfun, _mesh, _mode,
                    _epro, _fpro,
                    _nnew, _nold,
                    _tnew, _told, _ffpq,
                    _eprv, _edat, _escr,
                    _fprv, _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr,
                    _tdim, _pass, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._face_rule +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            else
            if (!_ftpq.empty() )
            {
        /*----------------------------- refine "bad" topo */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _kind =_bad_ftop( _geom,
                    _hfun, _mesh, _mode,
                    _epro, _fpro,
                    _nnew, _nold,
                    _tnew, _told,
                    _ftpq, _fmrk,
                    _eprv, _edat, _escr,
                    _fprv, _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr,
                    _tdim, _pass, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._face_rule +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
            else
            if (!_ttpq.empty() )
            {
        /*----------------------------- refine "bad" tria */
    #           ifdef  __use_timers
                _ttic = _time.now() ;
    #           endif//__use_timers

                _kind =_bad_tria( _geom,
                    _hfun, _mesh, _mode,
                    _epro, _fpro,
                    _nnew, _nold,
                    _tnew, _told, _ttpq,
                    _eprv, _edat, _escr,
                    _fprv, _fdat, _fscr,
                    _tdat, _tscr,
                    _bdat, _bscr,
                    _tdim, _pass, _args) ;

    #           ifdef  __use_timers
                _ttoc = _time.now() ;
                _tcpu._tria_rule +=
                    _tcpu.time_span(_ttic,_ttoc) ;
    #           endif//__use_timers
            }
        /*----------------------------- meshing converged */
            else { _done = true ; }

            if (_pass%_jlog_freq==+0 || _done)
            {
        /*----------------------------- output to logfile */
                if (_args.verb() >= +0)
                {
                std::stringstream _sstr;
                _sstr << std::setw(+11) <<
                _pass << std::setw(+13) <<
                _mesh._eset.count ()
                      << std::setw(+13) <<
                _mesh._fset.count ()
                      << std::setw(+13) <<
                _mesh._tset.count ()
                      <<   "\n"    ;
                _dump.push(_sstr.str());
                }
            }

            if (_kind != rdel_opts::null_kind)
            {
        /*----------------------------- inc. steiner info */
                if (_tdim == +1)
                {
                    _enod[_kind] += +1 ;
                }
                else
                if (_tdim == +2)
                {
                    _fnod[_kind] += +1 ;
                }
                else
                if (_tdim == +3)
                {
                    _tnod[_kind] += +1 ;
                }
            }

            }

            if (_pass%_trim_freq == +0 )
            {
        /*--------------- trim null PQ items "on-the-fly" */
                trim_list( _nbpq ) ;
                trim_eepq( _mesh ,
                           _eepq ) ;
                trim_ffpq( _mesh ,
                           _ffpq ) ;
                trim_ttpq( _mesh ,
                           _ttpq ) ;
                trim_list( _etpq ) ;
                trim_list( _ftpq ) ;

                trim_list( _nnew ) ;
                trim_list( _nold ) ;
                trim_list( _tnew ) ;
                trim_list( _told ) ;

                trim_list( _eprv ) ;
                trim_list( _edat ) ;
                trim_list( _escr ) ;
                trim_list( _fprv ) ;
                trim_list( _fdat ) ;
                trim_list( _fscr ) ;
                trim_list( _tdat ) ;
                trim_list( _tscr ) ;
                trim_list( _bdat ) ;
                trim_list( _bscr ) ;
            }

        /*--------------- enqueue nodes for topol. checks */

            fill_topo( _mesh, _pass,
                _etpq, _emrk,
                _ftpq, _fmrk,
                _edat, _eprv,
                _fdat, _fprv, _args)   ;

        /*--------------- update restricted triangulation */

            for (auto _npos  = _nold.head() ;
                      _npos != _nold.tend() ;
                    ++_npos  )
            {
                ball_data  _ball, _same;
                _ball._node[0] = *_npos;
                _ball._kind = feat_ball;
                _mesh.
                _pop_ball( _ball, _same) ;
            }

            for (auto _tpos  = _told.head() ;
                      _tpos != _told.tend() ;
                    ++_tpos  )
            {
                _pop_edge(_mesh, *_tpos) ;
                _pop_face(_mesh, *_tpos) ;
                _pop_tria(_mesh, *_tpos) ;
            }
            for (auto _tpos  = _told.head() ;
                      _tpos != _told.tend() ;
                    ++_tpos  )
            {
                _mesh.
                _tria._put_tria( *_tpos) ;
            }

            for (auto _iter  = _bscr.head() ;
                      _iter != _bscr.tend() ;
                    ++_iter  )
            {
                _nbpq .push( *_iter ) ;
            }
            for (auto _iter  = _escr.head() ;
                      _iter != _escr.tend() ;
                    ++_iter  )
            {
                _eepq .push( *_iter ) ;
            }
            for (auto _iter  = _fscr.head() ;
                      _iter != _fscr.tend() ;
                    ++_iter  )
            {
                _ffpq .push( *_iter ) ;
            }
            for (auto _iter  = _tscr.head() ;
                      _iter != _tscr.tend() ;
                    ++_iter  )
            {
                _ttpq .push( *_iter ) ;
            }

            for (auto _iter  = _bdat.head() ;
                      _iter != _bdat.tend() ;
                    ++_iter  )
            {
                _mesh.push_ball( *_iter) ;
            }
            for (auto _iter  = _edat.head() ;
                      _iter != _edat.tend() ;
                    ++_iter  )
            {
                _mesh.push_edge( *_iter) ;
            }
            for (auto _iter  = _fdat.head() ;
                      _iter != _fdat.tend() ;
                    ++_iter  )
            {
                _mesh.push_face( *_iter) ;
            }
            for (auto _iter  = _tdat.head() ;
                      _iter != _tdat.tend() ;
                    ++_iter  )
            {
                _mesh.push_tria( *_iter) ;
            }

        }

        if (_args.verb() >= +2 )
        {
    /*-------------------- push refinement scheme metrics */

        _dump.push("\n")  ;
        _dump.push("**REFINE statistics... \n") ;
        _dump.push("\n")  ;

        _dump.push("**FUNCTION timing: ") ;
        _dump.push("\n")  ;

        _dump.push("  MESH-SEED = ") ;
        _dump.push(
        std::to_string (_tcpu._mesh_seed));
        _dump.push("\n")  ;

        _dump.push("  NODE-INIT = ") ;
        _dump.push(
        std::to_string (_tcpu._node_init));
        _dump.push("\n")  ;
        _dump.push("  NODE-RULE = ") ;
        _dump.push(
        std::to_string (_tcpu._node_rule));
        _dump.push("\n")  ;

        _dump.push("  EDGE-INIT = ") ;
        _dump.push(
        std::to_string (_tcpu._edge_init));
        _dump.push("\n")  ;
        _dump.push("  EDGE-RULE = ") ;
        _dump.push(
        std::to_string (_tcpu._edge_rule));
        _dump.push("\n")  ;

        _dump.push("  FACE-INIT = ") ;
        _dump.push(
        std::to_string (_tcpu._face_init));
        _dump.push("\n")  ;
        _dump.push("  FACE-RULE = ") ;
        _dump.push(
        std::to_string (_tcpu._face_rule));
        _dump.push("\n")  ;

        _dump.push("  TRIA-INIT = ") ;
        _dump.push(
        std::to_string (_tcpu._tria_init));
        _dump.push("\n")  ;
        _dump.push("  TRIA-RULE = ") ;
        _dump.push(
        std::to_string (_tcpu._tria_rule));
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        _dump.push("**INSERTION rules: ") ;
        _dump.push("\n")  ;

        _dump.push("  EDGE-CIRC = ") ;
        _dump.push(std::to_string(
             _enod[rdel_opts::circ_kind]));
        _dump.push("\n")  ;
        _dump.push("  EDGE-OFFH = ") ;
        _dump.push(std::to_string(
             _enod[rdel_opts::offH_kind]));
        _dump.push("\n")  ;
        _dump.push("  EDGE-OFFT = ") ;
        _dump.push(std::to_string(
             _enod[rdel_opts::offT_kind]));
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        _dump.push("  FACE-CIRC = ") ;
        _dump.push(std::to_string(
             _fnod[rdel_opts::circ_kind]));
        _dump.push("\n")  ;
        _dump.push("  FACE-SINK = ") ;
        _dump.push(std::to_string(
             _fnod[rdel_opts::sink_kind]));
        _dump.push("\n")  ;
        _dump.push("  FACE-OFFH = ") ;
        _dump.push(std::to_string(
             _fnod[rdel_opts::offH_kind]));
        _dump.push("\n")  ;
        _dump.push("  FACE-OFFC = ") ;
        _dump.push(std::to_string(
             _fnod[rdel_opts::offC_kind]));
        _dump.push("\n")  ;
        _dump.push("  FACE-OFFT = ") ;
        _dump.push(std::to_string(
             _fnod[rdel_opts::offT_kind]));
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        _dump.push("  TRIA-CIRC = ") ;
        _dump.push(std::to_string(
             _tnod[rdel_opts::circ_kind]));
        _dump.push("\n")  ;
        _dump.push("  TRIA-SINK = ") ;
        _dump.push(std::to_string(
             _tnod[rdel_opts::sink_kind]));
        _dump.push("\n")  ;
        _dump.push("  TRIA-OFFH = ") ;
        _dump.push(std::to_string(
             _tnod[rdel_opts::offH_kind]));
        _dump.push("\n")  ;
        _dump.push("  TRIA-OFFC = ") ;
        _dump.push(std::to_string(
             _tnod[rdel_opts::offC_kind]));
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        }

        if (_args.verb() >= +3 )
        {
    /*-------------------- push refinement memory metrics */

        _dump.push("\n")  ;
        _dump.push("**MEMORY statistics... \n") ;
        _dump.push("\n")  ;

        _dump.push("**DELAUNAY-OBJECT: ") ;
        _dump.push("\n")  ;

        _dump.push("  NODE-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(typename mesh_type::
                tria_type:: node_type)) ) ;
        _dump.push("\n")  ;
        _dump.push("  NODE-LIST = ") ;
        _dump.push(std::to_string(
            _mesh._tria._nset.alloc())) ;
        _dump.push("\n")  ;

        _dump.push("  TRIA-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(typename mesh_type::
                tria_type:: tria_type)) ) ;
        _dump.push("\n")  ;
        _dump.push("  TRIA-LIST = ") ;
        _dump.push(std::to_string(
            _mesh._tria._tset.alloc())) ;
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        _dump.push("**RESTRICTED-TRIA: ") ;
        _dump.push("\n")  ;

        _dump.push("  BALL-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(
        typename mesh_type::ball_item)) ) ;
        _dump.push("\n")  ;
        _dump.push("  BALL-HASH = ") ;
        _dump.push(std::to_string(
            _mesh._bset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push("  POOL-BYTE = ") ;
        _dump.push(std::to_string(
            _mesh._bpol.bytes () ) ) ;
        _dump.push("\n")  ;

        _dump.push("  NODE-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(
        typename mesh_type::node_item)) ) ;
        _dump.push("\n")  ;
        _dump.push("  NODE-HASH = ") ;
        _dump.push(std::to_string(
            _mesh._nset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push("  POOL-BYTE = ") ;
        _dump.push(std::to_string(
            _mesh._npol.bytes () ) ) ;
        _dump.push("\n")  ;

        _dump.push("  EDGE-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(
        typename mesh_type::edge_item)) ) ;
        _dump.push("\n")  ;
        _dump.push("  EDGE-HASH = ") ;
        _dump.push(std::to_string(
            _mesh._eset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push("  POOL-BYTE = ") ;
        _dump.push(std::to_string(
            _mesh._epol.bytes () ) ) ;
        _dump.push("\n")  ;

        _dump.push("  FACE-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(
        typename mesh_type::face_item)) ) ;
        _dump.push("\n")  ;
        _dump.push("  FACE-HASH = ") ;
        _dump.push(std::to_string(
            _mesh._fset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push("  POOL-BYTE = ") ;
        _dump.push(std::to_string(
            _mesh._fpol.bytes () ) ) ;
        _dump.push("\n")  ;

        _dump.push("  TRIA-BYTE = ") ;
        _dump.push(std::to_string(
            sizeof(
        typename mesh_type::tria_item)) ) ;
        _dump.push("\n")  ;
        _dump.push("  TRIA-HASH = ") ;
        _dump.push(std::to_string(
            _mesh._tset._lptr.alloc())) ;
        _dump.push("\n")  ;
        _dump.push("  POOL-BYTE = ") ;
        _dump.push(std::to_string(
            _mesh._tpol.bytes () ) ) ;
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        _dump.push("**PRIORITY-QUEUES: ") ;
        _dump.push("\n")  ;

        _dump.push("  BSCR-BYTE = ") ;
        _dump.push(
        std::to_string(sizeof(ball_data)));
        _dump.push("\n")  ;
        _dump.push("  BBPQ-ITEM = ") ;
        _dump.push(std::to_string(_Nbpq)) ;
        _dump.push("\n")  ;

        _dump.push("  ESCR-BYTE = ") ;
        _dump.push(
        std::to_string(sizeof(edge_cost)));
        _dump.push("\n")  ;
        _dump.push("  EEPQ-ITEM = ") ;
        _dump.push(std::to_string(_Nepq)) ;
        _dump.push("\n")  ;

        _dump.push("  FSCR-BYTE = ") ;
        _dump.push(
        std::to_string(sizeof(face_cost)));
        _dump.push("\n")  ;
        _dump.push("  FFPQ-ITEM = ") ;
        _dump.push(std::to_string(_Nfpq)) ;
        _dump.push("\n")  ;

        _dump.push("  TSCR-BYTE = ") ;
        _dump.push(
        std::to_string(sizeof(tria_cost)));
        _dump.push("\n")  ;
        _dump.push("  TTPQ-ITEM = ") ;
        _dump.push(std::to_string(_Ntpq)) ;
        _dump.push("\n")  ;
        _dump.push("\n")  ;

        }

        _dump.push("\n")  ;
    }

    } ;


    }

#   endif   //__RDEL_MESH_3__



