
    /*
    --------------------------------------------------------
     * TREE-MESH-2: tree-like bisection for 2-manifolds. 
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
     * Last updated: 22 November, 2019
     *
     * Copyright 2013-2019
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    #   pragma once

#   ifndef __TREE_MESH_2__
#   define __TREE_MESH_2__

    namespace mesh {
    
    template <
    typename M ,
    typename P ,
    typename G ,
    typename H ,
    typename I ,
    typename A = allocators::basic_alloc
             >
    class tree_mesh_2d
    { 
    public  : 
    
    /*--------------- tree-like bisection for 2-manifolds */  
    
    typedef M                           mesh_type ;
    typedef P                           pred_type ;
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


    typedef mesh::mesh_params       <
                real_type, 
                iptr_type           >   tree_opts ;

    typedef mesh::tree_timers       <
                real_type ,
                iptr_type           >   tree_stat ;




    /*
    --------------------------------------------------------
     * TREE-MESH: subdivide a 2-manifold mesh. 
    --------------------------------------------------------
     */

    template <
    typename      jlog_file
             >
    __static_call
    __normal_call void_type tree_mesh (
        geom_type &_geom ,
        init_type &_init ,
        hfun_type &_hfun ,
        mesh_type &_mesh ,
        tree_opts &_args ,
        jlog_file &_dump
        )
    {   
    //!!mode_type _mode = null_mode ;
                
    /*------------------------------ push log-file header */
        if (_args.verb() >= 0 )
        {
            _dump.push(
    "#------------------------------------------------------------\n"
    "#    |ITER.|      |CEL-1|      |CEL-2| \n"
    "#------------------------------------------------------------\n"
                    ) ;
        }

    #   ifdef  __use_timers
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttic ;
        typename std ::chrono::
        high_resolution_clock::
            time_point _ttoc ;
        typename std ::chrono::
        high_resolution_clock _time ;

        __unreferenced(_time) ; // why does MSVC need this??
    #   endif//__use_timers

    /*------------------------------ ensure deterministic */  
        std::srand( +1 ) ;

        tree_stat _tcpu  ;


        __unreferenced(_geom) ;
        __unreferenced(_init) ;
        __unreferenced(_hfun) ;
        __unreferenced(_mesh) ;


    /*------------------------------ calc. hfun. at seeds */

        /*
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
        */


    /*-------------------- MAIN: refine edges/faces/cells */

        // push all cells onto queues

        // keep doing passes until queues empty

            // get list of unique edges/faces/etc from active sets

            // determine which edges to refine:
            // actually -- farm this out to a predicate, like rdel-mesh?
            // an edge is refined if it's too big (sample at midpoint + vertices)
            // or if it's hausdorf error is too big, a'la rdel-mesh

            // any elements


    }

    } ;

 
    }
    
#   endif   //__TREE_MESH_2__



