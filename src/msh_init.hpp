
    /*
    --------------------------------------------------------
     * MSH-INIT: init MESH data-structures.
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
     * Last updated: 27 Jul., 2020
     *
     * Copyright 2013-2020
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_INIT__
#   define __MSH_INIT__

     /*
    --------------------------------------------------------
     * INIT-MESH: initialise MESH data.
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type init_mesh (
        jcfg_data &_jcfg ,
        jlog_data &_jlog ,
        geom_data &_geom ,
        mesh_data &_mesh
        )
    {
        iptr_type _errv  = __no_error ;

        __unreferenced(_jcfg) ;
        __unreferenced(_jlog) ;

        if (_mesh._ndim == +2 &&
            _mesh._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            for (auto _iter  = _mesh.
            _euclidean_mesh_2d._mesh.node().head() ;
                      _iter != _mesh.
            _euclidean_mesh_2d._mesh.node().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                if (_iter->itag () > +0)
                if (_geom._kind ==
                jmsh_kind::euclidean_mesh)
                    _geom.
                _euclidean_mesh_2d.projector (
                   &_iter->pval(0) ,
                    _iter->itag () ,
                   &_iter->pval(0) ) ;

                if (_iter->feat()
                        == mesh::null_feat)
                    _iter->fdim() = 2  ;
                else
                if (_iter->feat()
                        == mesh::user_feat)
                    _iter->fdim() = 2  ;
                else
                if (_iter->feat()
                        == mesh::soft_feat)
                    _iter->fdim() = 0  ;
                else
                if (_iter->feat()
                        == mesh::hard_feat)
                    _iter->fdim() = 0  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_2d._mesh.tri3().head() ;
                      _iter != _mesh.
            _euclidean_mesh_2d._mesh.tri3().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_2d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;
                 auto _knod  =
                _head+_iter->node( +2) ;

                if (_inod->fdim() > 2)
                    _inod->fdim() = 2  ;
                if (_jnod->fdim() > 2)
                    _jnod->fdim() = 2  ;
                if (_knod->fdim() > 2)
                    _knod->fdim() = 2  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_2d._mesh.quad().head() ;
                      _iter != _mesh.
            _euclidean_mesh_2d._mesh.quad().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_2d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;
                 auto _knod  =
                _head+_iter->node( +2) ;
                 auto _lnod  =
                _head+_iter->node( +3) ;

                if (_inod->fdim() > 2)
                    _inod->fdim() = 2  ;
                if (_jnod->fdim() > 2)
                    _jnod->fdim() = 2  ;
                if (_knod->fdim() > 2)
                    _knod->fdim() = 2  ;
                if (_lnod->fdim() > 2)
                    _lnod->fdim() = 2  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_2d._mesh.edge().head() ;
                      _iter != _mesh.
            _euclidean_mesh_2d._mesh.edge().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_2d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;

                if (_inod->fdim() > 1)
                    _inod->fdim() = 1  ;
                if (_jnod->fdim() > 1)
                    _jnod->fdim() = 1  ;
            }
        }
        else
        if (_mesh._ndim == +3 &&
            _mesh._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            for (auto _iter  = _mesh.
            _euclidean_mesh_3d._mesh.node().head() ;
                      _iter != _mesh.
            _euclidean_mesh_3d._mesh.node().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                if (_iter->itag () > +0)
                if (_geom._kind ==
                jmsh_kind::euclidean_mesh)
                    _geom.
                _euclidean_mesh_3d.projector (
                   &_iter->pval(0) ,
                    _iter->itag () ,
                   &_iter->pval(0) ) ;

                if (_iter->itag () > +0)
                if (_geom._kind ==
                jmsh_kind::ellipsoid_mesh)
                    _geom.
                _ellipsoid_mesh_3d.projector (
                   &_iter->pval(0) ,
                    _iter->itag () ,
                   &_iter->pval(0) ) ;

                if (_iter->feat()
                        == mesh::null_feat)
                    _iter->fdim() = 3  ;
                else
                if (_iter->feat()
                        == mesh::user_feat)
                    _iter->fdim() = 3  ;
                else
                if (_iter->feat()
                        == mesh::soft_feat)
                    _iter->fdim() = 0  ;
                else
                if (_iter->feat()
                        == mesh::hard_feat)
                    _iter->fdim() = 0  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_3d._mesh.tri4().head() ;
                      _iter != _mesh.
            _euclidean_mesh_3d._mesh.tri4().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_3d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;
                 auto _knod  =
                _head+_iter->node( +2) ;
                 auto _lnod  =
                _head+_iter->node( +3) ;

                if (_inod->fdim() > 3)
                    _inod->fdim() = 3  ;
                if (_jnod->fdim() > 3)
                    _jnod->fdim() = 3  ;
                if (_knod->fdim() > 3)
                    _knod->fdim() = 3  ;
                if (_lnod->fdim() > 3)
                    _lnod->fdim() = 3  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_3d._mesh.tri3().head() ;
                      _iter != _mesh.
            _euclidean_mesh_3d._mesh.tri3().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_3d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;
                 auto _knod  =
                _head+_iter->node( +2) ;

                if (_inod->fdim() > 2)
                    _inod->fdim() = 2  ;
                if (_jnod->fdim() > 2)
                    _jnod->fdim() = 2  ;
                if (_knod->fdim() > 2)
                    _knod->fdim() = 2  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_3d._mesh.quad().head() ;
                      _iter != _mesh.
            _euclidean_mesh_3d._mesh.quad().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_3d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;
                 auto _knod  =
                _head+_iter->node( +2) ;
                 auto _lnod  =
                _head+_iter->node( +3) ;

                if (_inod->fdim() > 2)
                    _inod->fdim() = 2  ;
                if (_jnod->fdim() > 2)
                    _jnod->fdim() = 2  ;
                if (_knod->fdim() > 2)
                    _knod->fdim() = 2  ;
                if (_lnod->fdim() > 2)
                    _lnod->fdim() = 2  ;
            }

            for (auto _iter  = _mesh.
            _euclidean_mesh_3d._mesh.edge().head() ;
                      _iter != _mesh.
            _euclidean_mesh_3d._mesh.edge().tend() ;
                    ++_iter )
            {
                if (_iter->mark() < 0) continue ;

                 auto _head  = _mesh.
            _euclidean_mesh_3d._mesh.node().head() ;

                 auto _inod  =
                _head+_iter->node( +0) ;
                 auto _jnod  =
                _head+_iter->node( +1) ;

                if (_inod->fdim() > 1)
                    _inod->fdim() = 1  ;
                if (_jnod->fdim() > 1)
                    _jnod->fdim() = 1  ;
            }
        }

        return (  _errv ) ;
    }

#   endif   //__MSH_INIT__



