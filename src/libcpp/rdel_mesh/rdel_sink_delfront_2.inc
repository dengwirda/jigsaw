
    /*
    --------------------------------------------------------
     * RDEL-PRED-DELFRONT-2: "frontal-DEL" kernel in R^2. 
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
     * Last updated: 04 October, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
    // from rdel_pred_delfront_2.hpp
    
  
    /*
    --------------------------------------------------------
     * TRIA-SINK: "sink"-based off-centre.
    --------------------------------------------------------
     */

    __static_call 
    __normal_call 
    typename rdel_opts::node_kind tria_sink (
        geom_type &_geom,
        hfun_type &_size,
        mesh_type &_mesh,
        iptr_type  _tpos,
        real_type *_tbal,
        real_type *_ppos,
        rdel_opts &_args
        )
    {
        class sink_pred
            {
    /*-------------------- helper: find "sink" candidates */
            public  :
            
            typedef typename 
                mesh_type::tria_type tria_type ;
            
            public  :    
                mesh_type*_mesh ;
                real_type _safe ;
                
                real_type _tbal[ 3] ;

            public  :
            __inline_call sink_pred (
                mesh_type*_msrc ,
                real_type*_tbal ,
                real_type _stol
                ) : _mesh(_msrc),
                    _safe(_stol)
            {   
                this->_tbal[0] = _tbal[ 0] ;
                this->_tbal[1] = _tbal[ 1] ;
                this->_tbal[2] = _tbal[ 2] ;
            }
            
            __inline_call bool_type operator() (
                tria_type&_tria,
                iptr_type _tpos,
                iptr_type _fpos
                )
            {
                real_type  _ball[2] = {
                _tria.tria(_tpos)->circ(+0),
                _tria.tria(_tpos)->circ(+1)
                    } ;

                __unreferenced(_fpos) ;

                real_type _blen = 
                    geometry::lensqr_2d (
                        this->_tbal, _ball);

                return _blen < this->_safe ;
            }
            } ;

        typename rdel_opts::node_kind 
        _kind  = rdel_opts::null_kind ;

        __unreferenced(_geom);
        __unreferenced(_size);
        
    /*-------------------------------- find list of sinks */
        real_type _safe = _args.snk2()*
                          _args.snk2();
        _safe *= _tbal[2] ;

        containers::
            array <iptr_type>  _tset;
        _tset. set_alloc(+32);
        _mesh._tria.walk_tria( _tpos, 
            sink_pred(&_mesh, 
                _tbal, _safe), _tset) ;

        real_type static const _ftol = 
            std::pow(
        std::numeric_limits<real_type>
        ::epsilon(),(real_type) +.80) ;

        real_type _best = _tbal[2] ;
        real_type _btol = 
            (real_type)+1. + _ftol ;
        _best  *= _btol ;

    /*-------------------------------- keep the best sink */
        for (auto _iter  = _tset.head();
                  _iter != _tset.tend();
                ++_iter )
        {
            if ( *_iter == _tpos) continue;

            iptr_type _tnod[ 3] ;
            _tnod[0] = _mesh.
            _tria.tria(*_iter)->node(0);
            _tnod[1] = _mesh.
            _tria.tria(*_iter)->node(1);
            _tnod[2] = _mesh.
            _tria.tria(*_iter)->node(2);

            typename mesh_type::
                     tria_data _tdat ;
            _tdat._node[ 0] = _tnod[ 0];
            _tdat._node[ 1] = _tnod[ 1];
            _tdat._node[ 2] = _tnod[ 2];
            
            typename mesh_type::
                     tria_list::
                 item_type *_tptr = nullptr ;
            if (_mesh.find_tria(_tdat,_tptr)) 
            {
                if (_tptr->_data._kind 
                        ==  mesh::good_item )
                continue ;// skip good tria
            }
            else
            {   continue ;// not restricted
            }

            real_type  _ball[3];
            _ball[0] = _mesh.
               _tria.tria(*_iter)->circ(0);
            _ball[1] = _mesh.
               _tria.tria(*_iter)->circ(1);
      
            _ball[ 2] = (real_type)+0. ; 
            _ball[ 2]+= 
                geometry::lensqr_2d (_ball, 
                   &_mesh._tria.node(
                        _tnod[0])->pval(0)) ;
            _ball[ 2]+= 
                geometry::lensqr_2d (_ball, 
                   &_mesh._tria.node(
                        _tnod[1])->pval(0)) ;
            _ball[ 2]+= 
                geometry::lensqr_2d (_ball, 
                   &_mesh._tria.node(
                        _tnod[2])->pval(0)) ;
                        
            _ball[ 2]/= (real_type)+3. ;            

            if (_best    < _ball[2])
            {
                _best    = _ball[2];

                _kind    =  
                 rdel_opts::sink_kind  ;

                _ppos[0] = _ball[0];
                _ppos[1] = _ball[1];
            }
        }

        return ( _kind ) ;
    }
    
    
    
