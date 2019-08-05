
    /*
    --------------------------------------------------------
     * RDEL-PRED-BASE-2: base-class for refinement in R^2. 
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
     * Last updated: 10 July, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     *
     * This class defines the basic "restricted" delaunay 
     * predicates for meshing in R^2, as described in:
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
     * S.W. Cheng, T.K. Dey, E.A. Ramos, (2010): "Delaunay 
     * Refinement for Piecewise Smooth Complexes", 
     * Discrete & Computational Geometry, 43, pp. 121-166,
     * https://doi.org/10.1007/s00454-008-9109-3
     *
     * Testing "restricted-ness" by incrementally clipping 
     * against dual half-planes is an elaboration on ideas
     * introduced in:
     *
     * B. Levy, (2016): "Robustness and efficiency of 
     * geometric programs: The Predicate Construction Kit 
     * (PCK)", 72, pp. 3-12,
     * https://doi.org/10.1016/j.cad.2015.10.004
     *
     * though I make use of an "epsilon"-based strategy 
     * designed to support geometry predicates that return 
     * intersections that may be "imprecise".
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_PRED_BASE_2__
#   define __RDEL_PRED_BASE_2__

    namespace mesh {
    
    template <
    typename G ,
    typename M
             >
    class rdel_pred_base_2
    {
    public  :
    
    /*------------ base-class for refinement kern. in R^2 */
    
    typedef G                           geom_type ;
    typedef M                           mesh_type ;

    typedef typename 
            mesh_type::real_type        real_type ;
    typedef typename 
            mesh_type::iptr_type        iptr_type ;
            
    typedef containers::
            fixed_array<iptr_type, 2>   half_type ;

    typedef typename 
            geom_type::line_type        line_type ;
    typedef typename 
            geom_type::ball_type        ball_type ;
    
    /*
    --------------------------------------------------------
     * CLIP-DUAL: test pt. wrt dual halfplanes.
    --------------------------------------------------------
     */

    __static_call
    __normal_call double half_sign (
    __const_ptr  (double) _pp,
    __const_ptr  (double) _pa,
    __const_ptr  (double) _pb
        )
    {
    /*-------- helper: eval. sign w.r.t. half-plane [a,b] */
        double _pm[2];
        _pm[0] = _pa[0] * 0.5 ;
        _pm[1] = _pa[1] * 0.5 ;
        _pm[0]+= _pb[0] * 0.5 ;
        _pm[1]+= _pb[1] * 0.5 ;
        
        double _ab[2];
        _ab[0] = _pb[0] ;
        _ab[1] = _pb[1] ;
        _ab[0]-= _pa[0] ;
        _ab[1]-= _pa[1] ;
        
        double _mp[2];
        _mp[0] = _pp[0] ;
        _mp[1] = _pp[1] ;
        _mp[0]-= _pm[0] ;
        _mp[1]-= _pm[1] ;
        
        double _dp = 
        _ab[0] * _mp[0] + 
        _ab[1] * _mp[1] ;
        
        return ((double) _dp );
    }
    
    template <
        typename  half_list
             >
    __static_call
    __normal_call bool_type clip_dual (
        mesh_type &_mesh,
        half_list &_hset,
        real_type *_ppos,
        bool_type &_safe,
        real_type  _rtol
        )
    {
    /*-------- test PPOS against adjacent dual halfplanes */
        double _PPOS[2] ;
        _PPOS[0] =  _ppos[0] ;
        _PPOS[1] =  _ppos[1] ;

        _safe = true ;

        for(auto _hpos = _hset.count() ;
                 _hpos-- != +0 ; )
        {
            auto _anod = _hset[_hpos][0] ;
            auto _bnod = _hset[_hpos][1] ;

            double _APOS[2] ;
            _APOS[0] = _mesh.
            _tria.node(_anod)->pval(0) ;
            _APOS[1] = _mesh.
            _tria.node(_anod)->pval(1) ;
            
            double _BPOS[2] ;
            _BPOS[0] = _mesh.
            _tria.node(_bnod)->pval(0) ;
            _BPOS[1] = _mesh.
            _tria.node(_bnod)->pval(1) ;
            
            double _sign = +0.0 ;
            if (_bnod > _anod)
            _sign = + half_sign (
                (double*) _PPOS,
                (double*) _APOS,
                (double*) _BPOS) ;
            else
            _sign = - half_sign (
                (double*) _PPOS,
                (double*) _BPOS,
                (double*) _APOS) ;
            
    /*-------- "fatten" dual cavity to filter imprecision */
            if (_sign >= -_rtol &&
                _sign <= +_rtol)
                _safe  =  false;

            if (_sign <  -_rtol) return false ;
        }

        return true ;
    }

    /*
    --------------------------------------------------------
     * EDGE-BALL: calc. edge-based circumballs.
    --------------------------------------------------------
     */
     
    __static_call 
    __normal_call bool_type edge_ball (
        geom_type &_geom,
        mesh_type &_mesh,
        iptr_type  _tadj,
        iptr_type  _eadj,
        real_type *_ebal,
        real_type *_sbal,
        char_type &_feat,
        char_type &_topo,
        iptr_type &_part
        )
    {
        real_type static const _rEPS = 
            std::pow(std::numeric_limits
                <real_type>::epsilon(),+.67);

        real_type _blen = 
        _geom._bmax[0]-_geom._bmin[0]
      + _geom._bmax[1]-_geom._bmin[1] ;
        _blen   /= (real_type) +2. ;

    /*--------------------------- init. output balls = 0. */
        _ebal[0] = (real_type) +0. ;
        _ebal[1] = (real_type) +0. ;
        _ebal[2] = (real_type) +0. ;
        
        _sbal[0] = (real_type) +0. ;
        _sbal[1] = (real_type) +0. ;
        _sbal[2] = (real_type) +0. ;
        
    /*--------------------------- get local neighbourhood */
        iptr_type _topp = +0 ;
        iptr_type _eopp = +0 ;
        iptr_type _tmrk = +0 ;
        _mesh. _tria.find_pair (
        _tadj, _topp,
        _eadj, _eopp, _tmrk) ;
    
    /*--------------------------- skip faces adj. to hull */
        if (_topp == _mesh._tria.null_flag())
            return ( false ) ;
    
    /*--------------------------- assemble local indexing */
        if (_topp < _tadj)
        {        
        std::swap(_tadj, _topp) ;
        std::swap(_eadj, _eopp) ;
        }
    
        iptr_type _enod[ +3] ;
        mesh_type::tria_type::tria_type::
        face_node(_enod, _eadj, 2, 1);
        _enod[0] =_mesh._tria.
         tria(_tadj)->node(_enod[ 0]);
        _enod[1] =_mesh._tria.
         tria(_tadj)->node(_enod[ 1]);
        _enod[2] =_mesh._tria.
         tria(_tadj)->node(_enod[ 2]);
        
        iptr_type _onod[ +3] ;
        mesh_type::tria_type::tria_type::
        face_node(_onod, _eopp, 2, 1);
        _onod[0] =_mesh._tria.
         tria(_topp)->node(_onod[ 0]);
        _onod[1] =_mesh._tria.
         tria(_topp)->node(_onod[ 1]);
        _onod[2] =_mesh._tria.
         tria(_topp)->node(_onod[ 2]);

    /*--------------------------- calc. faceted face-ball */
        iptr_type  _ENOD[ 2] ;
        _ENOD[0] = _enod[ 0] ;
        _ENOD[1] = _enod[ 1] ;

        algorithms::isort(
            &_ENOD[0], &_ENOD[2] , 
                std::less<iptr_type>()) ;

        geometry::
            circ_ball_2d( _ebal, 
            &_mesh._tria.
        node(_ENOD[0])->pval(0),
            &_mesh._tria.
        node(_ENOD[1])->pval(0)) ;
    
    /*--------------------------- est. dual voronoi patch */
        line_type _line;
        _line._ipos[0] =  
             _mesh._tria.
        tria(_tadj   )->circ(0) ;
        _line._ipos[1] =  
             _mesh._tria.
        tria(_tadj   )->circ(1) ;
        
        _line._jpos[0] =  
             _mesh._tria.
        tria(_topp   )->circ(0) ;
        _line._jpos[1] =  
             _mesh._tria.
        tria(_topp   )->circ(1) ;

        real_type _nvec[3];
        geometry::line_norm_2d( 
            &_mesh._tria.
        node(_enod[0])->pval(0) ,
            &_mesh._tria.
        node(_enod[1])->pval(0) ,
             _nvec) ;
          
        _nvec[2] = 
        geometry::length_2d(_nvec) ;
        _nvec[0]/= _nvec[2] ;
        _nvec[1]/= _nvec[2] ;

        real_type  _ival, _jval;
        geometry::proj_line_2d(
            _line. _ipos,
            _ebal, _nvec, _ival) ;
        
        geometry::proj_line_2d(
            _line. _jpos,
            _ebal, _nvec, _jval) ;

        real_type _diff ;
        _diff =   _jval - _ival;

        if (_diff > (real_type)0.)
        {
            _diff = 
        std::max( _diff, +_blen) ;
        }
        else
        {
            _diff = 
        std::min( _diff, -_blen) ;
        }

        _ival -=  _rEPS * _diff;
      
        _line._ipos[0] = 
        _ebal[0] + _ival*_nvec [0] ;
        _line._ipos[1] = 
        _ebal[1] + _ival*_nvec [1] ;
        
        _jval +=  _rEPS * _diff;
      
        _line._jpos[0] = 
        _ebal[0] + _jval*_nvec [0] ;
        _line._jpos[1] = 
        _ebal[1] + _jval*_nvec [1] ;
        
    /*--------------------------- find loc. intersections */
        mesh::keep_all_2d <
            real_type ,
            iptr_type > _pred ;

        if(!_geom.intersect (
            _line, _pred) )
    /*--------------------------- face cant be restricted */
            return false  ; 

    /*--------------------------- form list of halfplanes */
        containers::
        fixed_array<half_type, +4> _hset ;
        
        _hset[0][0] = _enod[2] ;
        _hset[0][1] = _enod[0] ;
        _hset[1][0] = _enod[2] ;
        _hset[1][1] = _enod[1] ;
        
        _hset[2][0] = _onod[2] ;
        _hset[2][1] = _enod[0] ;
        _hset[3][0] = _onod[2] ;
        _hset[3][1] = _enod[1] ;
        
    /*--------------------------- test loc. intersections */
        auto _iful = _pred._list.tend() ;
        auto _imax = _pred._list.tend() ;
        
        real_type _RTOL = _rEPS*_ebal[2];

        real_type _dmax  =
            -std::numeric_limits
                <real_type>::infinity() ;
        real_type _dful  =
            -std::numeric_limits
                <real_type>::infinity() ;
        
        bool_type _safe  ;
        for (auto _iter  = 
                  _pred._list.head() ;
                  _iter != 
                  _pred._list.tend() ;
                ++_iter  )
        {
            if (clip_dual( _mesh, _hset , 
                   &_iter->pval( 0), 
                    _safe, _RTOL) )
            {
    /*--------------------------- prune near-degeneracies */
                if (!_safe)
                {                
                typename mesh_type::
                         tria_type::
                         tria_pred::
                template circ_pred<
                typename mesh_type::tria_type> 
                    _isDT (&_iter-> pval( 0));
          
                bool_type _okay =
                _isDT(_mesh._tria, _tadj, 0)||
                _isDT(_mesh._tria, _topp, 0) ;

                if (!_okay) continue ;   
                }

    /*--------------------------- dist to face circumball */
                real_type _dsqr = 
                geometry::lensqr_2d(
                    _ebal, 
                   &_iter->pval( 0)) ;  

    /*--------------------------- keep furthest from ball */
                if (_dsqr > _dmax )
                {
                    _dmax = _dsqr ;
                    _imax = _iter ;
                }
                if (_dsqr > _dful &&
                    _safe )
                {
                    _dful = _dsqr ;
                    _iful = _iter ;
                }
            }
        }
        
        if (_iful !=
                _pred._list.tend() )
        {
    /*--------------------------- keep best intersections */
        _sbal[ 0] = _iful->pval(0);
        _sbal[ 1] = _iful->pval(1);
            
        _part     = _iful->itag ();
        _feat     = _iful->feat ();
        _topo     = _iful->topo ();
        
    /*--------------------------- eval. surf. ball radius */
        _sbal[ 2]+= 
        geometry::lensqr_2d(_sbal , 
       &_mesh._tria.
            node(_enod[0])->pval(0)) ;
        _sbal[ 2]+= 
        geometry::lensqr_2d(_sbal , 
       &_mesh._tria.
            node(_enod[1])->pval(0)) ;
             
        _sbal[ 2]/= (real_type) +2.;   
    
    /*--------------------------- return restricted state */
        return (  true ) ;
        }
        else
        if (_imax !=
                _pred._list.tend() )
        {
    /*--------------------------- keep best intersections */
        _sbal[ 0] = _imax->pval(0);
        _sbal[ 1] = _imax->pval(1);
            
        _part     = _imax->itag ();
        _feat     = _imax->feat ();
        _topo     = _imax->topo ();
        
    /*--------------------------- eval. surf. ball radius */
        _sbal[ 2]+= 
        geometry::lensqr_2d(_sbal , 
       &_mesh._tria.
            node(_enod[0])->pval(0)) ;
        _sbal[ 2]+= 
        geometry::lensqr_2d(_sbal , 
       &_mesh._tria.
            node(_enod[1])->pval(0)) ;
             
        _sbal[ 2]/= (real_type) +2.;   
    
    /*--------------------------- return restricted state */
        return (  true ) ;
        }
    
        return ( false ) ;
    }

    /*
    --------------------------------------------------------
     * TRIA-BALL: calc. tria-based circumballs.
    --------------------------------------------------------
     */
     
    __static_call 
    __inline_call bool_type tria_ball (
        geom_type &_geom,
        mesh_type &_mesh,
        iptr_type  _tpos,
        real_type *_tbal,
        iptr_type &_part
        )
    {
    /*--------------------------- assemble local indexing */
        iptr_type  _tnod[3] = {
            _mesh.
        _tria.tria(_tpos)->node(+0) ,
            _mesh.
        _tria.tria(_tpos)->node(+1) ,
            _mesh.
        _tria.tria(_tpos)->node(+2)
            } ;

    /*--------------------------- init. local output ball */
        _tbal[0] = _mesh.
            _tria.tria(_tpos)->circ(+0) ;
        _tbal[1] = _mesh.
            _tria.tria(_tpos)->circ(+1) ;
  
        _tbal[2] = (real_type)+.0 ; 
        _tbal[2]+= 
            geometry::lensqr_2d (_tbal, 
               &_mesh._tria.node(
                    _tnod[0])->pval(0)) ;
        _tbal[2]+= 
            geometry::lensqr_2d (_tbal, 
               &_mesh._tria.node(
                    _tnod[1])->pval(0)) ;
        _tbal[2]+= 
            geometry::lensqr_2d (_tbal, 
               &_mesh._tria.node(
                    _tnod[2])->pval(0)) ;
                    
        _tbal[2]/= (real_type)+3. ;
        
    /*------------------------- evaluate "in--out" status */
        if (_part <= -1 && (_part = 
            _geom.is_inside(_tbal)) < +0) 
        {
    /*------------------------- is not a restricted facet */
            return false ;
        }
        else
        {
            return  true ;
        }    
    }

    } ;

    
    }
    
#   endif   //__RDEL_PRED_BASE_2__



