
    /*
    --------------------------------------------------------
     * OFFSET: return domain separators.
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
     * Last updated: 18 Aug., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __OFFSET__
#   define __OFFSET__

    /*
    --------------------------------------------------------
     * SEPARATOR: return offset to domain "centre"
    --------------------------------------------------------
     */

    template <
    typename      jlog_data
             >
    __normal_call iptr_type separator (
        jcfg_data &_jcfg,
        jlog_data &_jlog,
        geom_data &_geom,
        float     *_xoff
        )
    {
        iptr_type _errv = __no_error ;

        __unreferenced(_jcfg) ;
        __unreferenced(_jlog) ;

    /*--------------------------------- find GEOM offsets */
        double _xmid[3] = { +0. } ;

        if (_geom._ndim == +2 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            for (auto _iter  = _geom.
            _euclidean_mesh_2d._tria.node().head();
                      _iter != _geom.
            _euclidean_mesh_2d._tria.node().tend();
                    ++_iter  )
            {
                _xmid[0] += _iter->pval(0);
                _xmid[1] += _iter->pval(1);
            }

            _xoff[0] = (float)(_xmid[0] /
            _geom._euclidean_mesh_2d._tria.node().count());
            _xoff[1] = (float)(_xmid[1] /
            _geom._euclidean_mesh_2d._tria.node().count());
        }
        else
        if (_geom._ndim == +3 &&
            _geom._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            for (auto _iter  = _geom.
            _euclidean_mesh_3d._tria.node().head();
                      _iter != _geom.
            _euclidean_mesh_3d._tria.node().tend();
                    ++_iter  )
            {
                _xmid[0] += _iter->pval(0);
                _xmid[1] += _iter->pval(1);
                _xmid[2] += _iter->pval(2);
            }

            _xoff[0] = (float)(_xmid[0] /
            _geom._euclidean_mesh_3d._tria.node().count());
            _xoff[1] = (float)(_xmid[1] /
            _geom._euclidean_mesh_3d._tria.node().count());
            _xoff[2] = (float)(_xmid[2] /
            _geom._euclidean_mesh_3d._tria.node().count());
        }

        return (  _errv ) ;
    }

    template <
    typename      jlog_data
             >
    __normal_call iptr_type separator (
        jcfg_data &_jcfg,
        jlog_data &_jlog,
        hfun_data &_hfun,
        float     *_xoff
        )
    {
        iptr_type _errv = __no_error ;

        __unreferenced(_jcfg) ;
        __unreferenced(_jlog) ;

    /*--------------------------------- find GEOM offsets */
        double _xmid[3] = { +0. } ;

        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-2d */
            for (auto _iter  = _hfun.
            _euclidean_mesh_2d._mesh.node().head();
                      _iter != _hfun.
            _euclidean_mesh_2d._mesh.node().tend();
                    ++_iter  )
            {
                _xmid[0] += _iter->pval(0);
                _xmid[1] += _iter->pval(1);
            }

            _xoff[0] = (float)(_xmid[0] /
            _hfun._euclidean_mesh_2d._mesh.node().count());
            _xoff[1] = (float)(_xmid[1] /
            _hfun._euclidean_mesh_2d._mesh.node().count());
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_mesh)
        {
    /*--------------------------------- euclidean-mesh-3d */
            for (auto _iter  = _hfun.
            _euclidean_mesh_3d._mesh.node().head();
                      _iter != _hfun.
            _euclidean_mesh_3d._mesh.node().tend();
                    ++_iter  )
            {
                _xmid[0] += _iter->pval(0);
                _xmid[1] += _iter->pval(1);
                _xmid[2] += _iter->pval(2);
            }

            _xoff[0] = (float)(_xmid[0] /
            _hfun._euclidean_mesh_3d._mesh.node().count());
            _xoff[1] = (float)(_xmid[1] /
            _hfun._euclidean_mesh_3d._mesh.node().count());
            _xoff[2] = (float)(_xmid[2] /
            _hfun._euclidean_mesh_3d._mesh.node().count());
        }
        else
        if (_hfun._ndim == +2 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-2d */
            for (auto _iter  = _hfun.
            _euclidean_grid_2d._xpos.head();
                      _iter != _hfun.
            _euclidean_grid_2d._xpos.tend();
                    ++_iter  )
            {
                _xmid[0] +=  *_iter;
            }

            _xoff[0] = (float)(_xmid[0] /
            _hfun._euclidean_grid_2d._xpos.count()) ;

            for (auto _iter  = _hfun.
            _euclidean_grid_2d._ypos.head();
                      _iter != _hfun.
            _euclidean_grid_2d._ypos.tend();
                    ++_iter  )
            {
                _xmid[1] +=  *_iter;
            }

            _xoff[1] = (float)(_xmid[1] /
            _hfun._euclidean_grid_2d._ypos.count()) ;
        }
        else
        if (_hfun._ndim == +3 &&
            _hfun._kind ==
             jmsh_kind::euclidean_grid)
        {
    /*--------------------------------- euclidean-grid-3d */
            for (auto _iter  = _hfun.
            _euclidean_grid_3d._xpos.head();
                      _iter != _hfun.
            _euclidean_grid_3d._xpos.tend();
                    ++_iter  )
            {
                _xmid[0] +=  *_iter;
            }

            _xoff[0] = (float)(_xmid[0] /
            _hfun._euclidean_grid_3d._xpos.count()) ;

            for (auto _iter  = _hfun.
            _euclidean_grid_3d._ypos.head();
                      _iter != _hfun.
            _euclidean_grid_3d._ypos.tend();
                    ++_iter  )
            {
                _xmid[1] +=  *_iter;
            }

            _xoff[1] = (float)(_xmid[1] /
            _hfun._euclidean_grid_3d._ypos.count()) ;

            for (auto _iter  = _hfun.
            _euclidean_grid_3d._zpos.head();
                      _iter != _hfun.
            _euclidean_grid_3d._zpos.tend();
                    ++_iter  )
            {
                _xmid[2] +=  *_iter;
            }

            _xoff[2] = (float)(_xmid[2] /
            _hfun._euclidean_grid_3d._zpos.count()) ;
        }

        return (  _errv ) ;
    }


#   endif   //__OFFSET__



