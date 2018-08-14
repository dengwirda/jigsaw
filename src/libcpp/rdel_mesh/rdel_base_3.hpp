
    /*
    --------------------------------------------------------
     * RDEL-PRED-BASE-3: base-class for refinement in R^3. 
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
     * Last updated: 10 December, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_PRED_BASE_3__
#   define __RDEL_PRED_BASE_3__

    namespace mesh {
    
    template <
    typename G ,
    typename M
             >
    class rdel_pred_base_3
    {
    public  :
    
    /*------------ base-class for refinement kern. in R^3 */
    
    typedef G                           geom_type ;
    typedef M                           mesh_type ;

    typedef typename 
            mesh_type::real_type        real_type ;
    typedef typename 
            mesh_type::iptr_type        iptr_type ;

    /*
    --------------------------------------------------------
     * EDGE-LOOP: assemble tria.'s adj. to edge.
    --------------------------------------------------------
     */
     
    template <
    typename      list_type
             >
    __static_call
    __normal_call void_type edge_loop (
        mesh_type &_mesh,
        iptr_type *_enod,
        iptr_type  _tadj,
        iptr_type  _fadj,
        list_type &_loop
        )
    {
        iptr_type _tcur = _tadj ;
        iptr_type _fcur = _fadj ;

        iptr_type _null = 
            _mesh._tria.null_flag ();

        while (true)
        {   
            _loop.push_tail ( _tcur);

            iptr_type _tpos = _tcur ;
            iptr_type _fpos ;
            iptr_type _tmrk ;
            for(_fpos = +4; _fpos-- != +0; )
            {
                if(_fpos == _fcur) continue;

                iptr_type _fnod[  +4] ;
                mesh_type::tria_type::
                    tria_type::
                face_node(_fnod, _fpos, 3, 2);
                _fnod[ +0] = _mesh._tria.
                 tria(_tpos)->node(_fnod[ 0]);
                _fnod[ +1] = _mesh._tria.
                 tria(_tpos)->node(_fnod[ 1]);
                _fnod[ +2] = _mesh._tria.
                 tria(_tpos)->node(_fnod[ 2]);

                iptr_type _same = +0;
                if (_fnod[0] == _enod[0] ||
                    _fnod[0] == _enod[1] )
                    _same += +1 ;
                if (_fnod[1] == _enod[0] ||
                    _fnod[1] == _enod[1] )
                    _same += +1 ;
                if (_fnod[2] == _enod[0] ||
                    _fnod[2] == _enod[1] )
                    _same += +1 ;

                if (_same == +2) { break ; }
            }

            _mesh._tria.find_pair(
            _tpos, _tcur, 
            _fpos, _fcur, _tmrk) ;

            if (_tcur == _tadj)  { break ; }
            if (_tcur == _null)  { break ; }
        }
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
        real_type *_pmax,
        char_type &_hits,
        char_type &_feat,
        char_type &_topo,
        iptr_type &_part
        )
    {
        typedef containers::
        fixed_array<real_type, +3> ball_type;

    /*--------------------------- assemble local indexing */
        iptr_type _enod[ +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_enod, _eadj, 3, 1);
        _enod[0] =_mesh._tria.
         tria(_tadj)->node(_enod[ 0]);
        _enod[1] =_mesh._tria.
         tria(_tadj)->node(_enod[ 1]);
         
        if (_enod[1] < _enod[0])
            std::swap( _enod[0],_enod[1]);
            
    /*--------------------------- init. output balls = 0. */
        _ebal[0] = (real_type) +0. ;
        _ebal[1] = (real_type) +0. ;
        _ebal[2] = (real_type) +0. ;
        _ebal[3] = (real_type) +0. ;
        
        _pmax[0] = (real_type) +0. ;
        _pmax[1] = (real_type) +0. ;
        _pmax[2] = (real_type) +0. ;
        _pmax[3] = (real_type) +0. ;
        
    /*--------------------------- get local neighbourhood */        
        iptr_type _fadj;
        for (_fadj = +4; _fadj-- != +0; )
        {
            iptr_type _fnod[ +4] ;
            mesh_type::tria_type::
                tria_type::
            face_node(_fnod, _fadj, 3, 2);
            _fnod[0] =_mesh._tria.
             tria(_tadj)->node(_fnod[ 0]);
            _fnod[1] =_mesh._tria.
             tria(_tadj)->node(_fnod[ 1]);
            _fnod[2] =_mesh._tria.
             tria(_tadj)->node(_fnod[ 2]);

            iptr_type _same = +0;
            if (_fnod[0] == _enod[0] ||
                _fnod[0] == _enod[1] )
                _same += +1 ;
            if (_fnod[1] == _enod[0] ||
                _fnod[1] == _enod[1] )
                _same += +1 ;
            if (_fnod[2] == _enod[0] ||
                _fnod[2] == _enod[1] )
                _same += +1 ;

            if (_same == +2) { break ; }
        }

        containers::array<iptr_type> _tset ;
        containers::array<ball_type> _tbal ;
        _tset. set_alloc(+32);
        _tbal. set_alloc(+32);
        edge_loop( _mesh, _enod, 
            _tadj, _fadj, _tset)  ;
        
        if (_tset.count()<+3) return false ;
           
        _ebal[0] =
        (real_type)+.50*_mesh._tria.
         node(_enod[0])->pval(0) +
        (real_type)+.50*_mesh._tria.
         node(_enod[1])->pval(0) ;
        _ebal[1] = 
        (real_type)+.50*_mesh._tria.
         node(_enod[0])->pval(1) +
        (real_type)+.50*_mesh._tria.
         node(_enod[1])->pval(1) ;
        _ebal[2] = 
        (real_type)+.50*_mesh._tria.
         node(_enod[0])->pval(2) +
        (real_type)+.50*_mesh._tria.
         node(_enod[1])->pval(2) ;
         
         real_type _nvec[ 3] = {
            _mesh._tria.
         node(_enod[1])->pval(0) -
            _mesh._tria.
         node(_enod[0])->pval(0) ,
            _mesh._tria.
         node(_enod[1])->pval(1) -
            _mesh._tria.
         node(_enod[0])->pval(1) ,
            _mesh._tria.
         node(_enod[1])->pval(2) -
            _mesh._tria.
         node(_enod[0])->pval(2) } ;
         
        _ebal[3] = (real_type)+.25 * 
            geometry::lensqr_3d (
       &_mesh._tria.
         node(_enod[0])->pval(0) ,
       &_mesh._tria.
         node(_enod[1])->pval(0) ) ;

        mesh::keep_minmax_length_3d<
            real_type, 
            iptr_type> _pred(_ebal ) ;

    #   ifdef __testdual
        real_type _tvec[3], _last[3] ;
        _last[0] = (real_type)+0. ;
        _last[1] = (real_type)+0. ;
        _last[2] = (real_type)+0. ;

        for (auto _tpos  = _tset.head()+2;
                  _tpos != _tset.tend()+0;
                ++_tpos  )
        {
            geometry::tria_norm_3d (
           &_mesh._tria.
             tria(  _tset[0])->circ(0),
           &_mesh._tria.
             tria(*(_tpos-1))->circ(0),
           &_mesh._tria.
             tria(*(_tpos+0))->circ(0),
            _tvec) ;
           
            real_type _vdot = 
            geometry::dot_3d(_tvec, _last) ;
           
            if (_vdot < (real_type)0.) 
            {
                std::cout << 
                    "bad voronoi face! \n" ;
            }
           
            _last[0] = _tvec[0] ;
            _last[1] = _tvec[1] ;
            _last[2] = _tvec[2] ;
        }
    #   endif
    
        _tbal.set_count(_tset.count ());

        iptr_type _tnum  = +0;
        for (auto _tpos  = _tset.head();
                  _tpos != _tset.tend();
                ++_tpos, ++_tnum)
        {
            _tbal[_tnum][0] = 
            _mesh._tria.
                tria(*_tpos)->circ( 0) ;
            _tbal[_tnum][1] = 
            _mesh._tria.
                tria(*_tpos)->circ( 1) ;
            _tbal[_tnum][2] = 
            _mesh._tria.
                tria(*_tpos)->circ( 2) ;
        }

        _geom.intersect (_ebal, _nvec, _tbal, _pred) ;

        if ( _pred._find )
        {
    /*------------------------- store surface ball param. */    
        _part=_pred._pmin.itag () ;
        
        _hits=_pred._pmin.hits () ;
        _feat=_pred._pmin.feat () ;
        _topo=_pred._pmin.topo () ;
        
        _pmax[ 0] 
            = _pred._pmin.pval(0) ;
        _pmax[ 1] 
            = _pred._pmin.pval(1) ;
        _pmax[ 2] 
            = _pred._pmin.pval(2) ;
        
    /*------------------------- eval. surface ball radius */
        _pmax[ 3]+= 
        geometry::lensqr_3d(_pmax , 
       &_mesh._tria.
            node(_enod[0])->pval(0)) ;
        _pmax[ 3]+= 
        geometry::lensqr_3d(_pmax , 
       &_mesh._tria.
            node(_enod[1])->pval(0)) ;
             
        _pmax[ 3]/= (real_type)+2.;
        
        }
        
    /*------------------------- return TRUE if restricted */        
        return ( _pred._find ) ;
    }

    /*
    --------------------------------------------------------
     * FACE-BALL: calc. face-based circumballs.
    --------------------------------------------------------
     */
     
    __static_call 
    __normal_call bool_type face_ball (
        geom_type &_geom,
        mesh_type &_mesh,
        iptr_type  _tadj,
        iptr_type  _fadj,
        real_type *_fbal,
        real_type *_sbal,
        char_type &_feat,
        char_type &_topo,
        iptr_type &_part
        )
    {
    /*--------------------------- assemble local indexing */
        iptr_type _fnod[ +4] ;
        mesh_type::tria_type::tria_type::
        face_node(_fnod, _fadj, 3, 2);
        _fnod[0] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 0]);
        _fnod[1] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 1]);
        _fnod[2] =_mesh._tria.
         tria(_tadj)->node(_fnod[ 2]);

    /*--------------------------- init. output balls = 0. */
        _fbal[0] = (real_type) +0. ;
        _fbal[1] = (real_type) +0. ;
        _fbal[2] = (real_type) +0. ;
        _fbal[3] = (real_type) +0. ;
        
        _sbal[0] = (real_type) +0. ;
        _sbal[1] = (real_type) +0. ;
        _sbal[2] = (real_type) +0. ;
        _sbal[3] = (real_type) +0. ;

    /*--------------------------- get local neighbourhood */
        iptr_type _topp = +0 ;
        iptr_type _fopp = +0 ;
        iptr_type _tmrk = +0 ;
        _mesh. _tria.find_pair (
        _tadj, _topp,
        _fadj, _fopp, _tmrk) ;

    /*--------------------------- skip faces adj. to hull */
        if (_topp == 
            _mesh._tria.null_flag())
            return ( false ) ;

        if (_topp < _tadj)
            std::swap(_tadj, _topp);

    #   ifdef __testdual
        real_type _nvec[ +3] ;
        geometry::tria_norm_3d(
       &_mesh._tria.
         node(_fnod[ 0])->pval(0) ,
       &_mesh._tria.
         node(_fnod[ 1])->pval(0) ,
       &_mesh._tria.
         node(_fnod[ 2])->pval(0) ,
        _nvec ) ;
        
        real_type _cvec[3] = {
        _mesh._tria.
            tria( _tadj)->circ(0) - 
        _mesh._tria.
            tria( _topp)->circ(0) ,
        _mesh._tria.
            tria( _tadj)->circ(1) - 
        _mesh._tria.
            tria( _topp)->circ(1) ,
        _mesh._tria.
            tria( _tadj)->circ(2) - 
        _mesh._tria.
            tria( _topp)->circ(2)
            } ;

        real_type _vdot = 
        geometry::dot_3d(_nvec, _cvec) ;
        
        if (_vdot < (real_type)0.)
        {
            std::cout << 
                "bad voronoi edge! \n" ;
        }
    #   endif

    /*------------------------------ bipolar voronoi edge */
        real_type  _ibal[ 3] = {
        _mesh.
        _tria.tria(_tadj)->circ( 0),
        _mesh.
        _tria.tria(_tadj)->circ( 1),
        _mesh.
        _tria.tria(_tadj)->circ( 2)
            } ;
        
        real_type  _jbal[ 3] = {
        _mesh.
        _tria.tria(_topp)->circ( 0),
        _mesh.
        _tria.tria(_topp)->circ( 1),
        _mesh.
        _tria.tria(_topp)->circ( 2)
            } ;
 
    /*------------------------------ calc. diametric ball */
        geometry::circ_ball_3d(_fbal , 
       &_mesh._tria.
            node(_fnod[0])->pval(0),
       &_mesh._tria.
            node(_fnod[1])->pval(0),
       &_mesh._tria.
            node(_fnod[2])->pval(0)) ;

    /*------------------------- find surface intersection */
        mesh::keep_minmax_length_3d<
            real_type, 
            iptr_type> _pred(_fbal ) ;
        if (!_geom.intersect (
               &_ibal[+0] ,
               &_jbal[+0] ,  _pred))
    /*------------------------- is not a restricted facet */
        return ( _pred._find ) ;
   
    /*------------------------- store surface ball param. */
        _sbal[ 0] = _pred._pmin.pval(0);
        _sbal[ 1] = _pred._pmin.pval(1);
        _sbal[ 2] = _pred._pmin.pval(2);
               
        _feat     = _pred._pmin.feat ();
        _topo     = _pred._pmin.topo ();
   
        _part     = _pred._pmin.itag ();
   
    /*------------------------- eval. surface ball radius */
        _sbal[ 3]+= 
        geometry::lensqr_3d (_sbal, 
            &_mesh._tria.
              node(_fnod[0])->pval(0)) ;
        _sbal[ 3]+= 
        geometry::lensqr_3d (_sbal, 
            &_mesh._tria.
              node(_fnod[1])->pval(0)) ;
        _sbal[ 3]+= 
        geometry::lensqr_3d (_sbal, 
            &_mesh._tria.
              node(_fnod[2])->pval(0)) ;

        _sbal[ 3]/= (real_type)+3. ;

    /*----------------------- return TRUE if "restricted" */
        return ( _pred._find ) ;
    }
    
    } ;
    
    
    }
    
#   endif   //__RDEL_PRED_BASE_3__
    
    
    
