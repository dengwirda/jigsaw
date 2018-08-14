
    /*
    --------------------------------------------------------
     * INIT-MSH_T: init. msh_t for lib_jigsaw.
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
     * Last updated: 10 June, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __INIT_MSH_T__
#   define __INIT_MSH_T__

    /*
    --------------------------------------------------------
     * ALLOC-MSH_T: alloc. msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw
    
    void jigsaw_alloc_vert2 (               // _new-alloc
        jigsaw_VERT2_array_t *_vert2 ,
        indx_t _size
        )
    {
        _vert2->_size = _size ;       
        _vert2->_data = 
        (jigsaw_VERT2_t *) malloc (
            _size * sizeof(jigsaw_VERT2_t)) ;   
    }
    
    void jigsaw_alloc_vert3 (               // _new-alloc
        jigsaw_VERT3_array_t *_vert3 ,
        indx_t _size
        )
    {
        _vert3->_size = _size ;       
        _vert3->_data = 
        (jigsaw_VERT3_t *) malloc (
            _size * sizeof(jigsaw_VERT3_t)) ;   
    }
    
    void jigsaw_alloc_edge2 (               // _new-alloc
        jigsaw_EDGE2_array_t *_edge2 ,
        indx_t _size
        )
    {
        _edge2->_size = _size ;       
        _edge2->_data = 
        (jigsaw_EDGE2_t *) malloc (
            _size * sizeof(jigsaw_EDGE2_t)) ;   
    }
 
    void jigsaw_alloc_tria3 (               // _new-alloc
        jigsaw_TRIA3_array_t *_tria3 ,
        indx_t _size
        )
    {
        _tria3->_size = _size ;       
        _tria3->_data = 
        (jigsaw_TRIA3_t *) malloc (
            _size * sizeof(jigsaw_TRIA3_t)) ;   
    }
    
    void jigsaw_alloc_tria4 (               // _new-alloc
        jigsaw_TRIA4_array_t *_tria4 ,
        indx_t _size
        )
    {
        _tria4->_size = _size ;       
        _tria4->_data = 
        (jigsaw_TRIA4_t *) malloc (
            _size * sizeof(jigsaw_TRIA4_t)) ;   
    }

    void jigsaw_alloc_reals (               // _new-alloc
        jigsaw_REALS_array_t *_reals ,
        indx_t _size
        )
    {
        _reals->_size = _size ;       
        _reals->_data = (real_t *) 
            malloc (_size * sizeof(real_t)) ;   
    }
    
#   endif//__lib_jigsaw 
    
    /*
    --------------------------------------------------------
     * FREE-MSH_T: _free msh_t for lib_jigsaw.
    --------------------------------------------------------
     */
    
#   ifdef  __lib_jigsaw
    
    void jigsaw_free_msh_t (                // free-alloc
        jigsaw_msh_t *_mesh
        )
    {
    /*-------------------------- dealloc. msh_t data type */
        _mesh->_flags = JIGSAW_NULL_FLAG ;
    
        if (_mesh->_vert2._size > +0)
        {
            _mesh->_vert2._size = +0;
             free(
            _mesh->_vert2._data) ;
        }
        
        if (_mesh->_vert3._size > +0)
        {
            _mesh->_vert3._size = +0;
             free(
            _mesh->_vert3._data) ;
        }
        
        if (_mesh->_power._size > +0)
        {
            _mesh->_power._size = +0;
             free(
            _mesh->_power._data) ;
        }
        
        if (_mesh->_edge2._size > +0)
        {
            _mesh->_edge2._size = +0;
             free(
            _mesh->_edge2._data) ;
        }
        
        if (_mesh->_tria3._size > +0)
        {
            _mesh->_tria3._size = +0;
             free(
            _mesh->_tria3._data) ;
        }
        
        if (_mesh->_tria4._size > +0)
        {
            _mesh->_tria4._size = +0;
             free(
            _mesh->_tria4._data) ;
        }
        
        if (_mesh->_xgrid._size > +0)
        {
            _mesh->_xgrid._size = +0;
             free(
            _mesh->_xgrid._data) ;
        }
        
        if (_mesh->_ygrid._size > +0)
        {
            _mesh->_ygrid._size = +0;
             free(
            _mesh->_ygrid._data) ;
        }
        
        if (_mesh->_zgrid._size > +0)
        {
            _mesh->_zgrid._size = +0;
             free(
            _mesh->_zgrid._data) ;
        }

        if (_mesh->_value._size > +0)
        {
            _mesh->_value._size = +0;
             free(
            _mesh->_value._data) ;
        }
        
    }
    
#   endif//__lib_jigsaw

    /*
    --------------------------------------------------------
     * INIT-MSH_T: init. msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw
    
    void jigsaw_init_msh_t (                // init. msh_t
        jigsaw_msh_t *_mesh
        )
    {
        _mesh->_flags = JIGSAW_NULL_FLAG ;
        
        _mesh->_vert2._size = +0 ;
        _mesh->_vert2._data = nullptr ;
        
        _mesh->_vert3._size = +0 ;
        _mesh->_vert3._data = nullptr ;
        
        _mesh->_power._size = +0 ;
        _mesh->_power._data = nullptr ;
        
        _mesh->_edge2._size = +0 ;
        _mesh->_edge2._data = nullptr ;
        
        _mesh->_tria3._size = +0 ;
        _mesh->_tria3._data = nullptr ;
        
        _mesh->_tria4._size = +0 ;
        _mesh->_tria4._data = nullptr ;
        
        _mesh->_value._size = +0 ;
        _mesh->_value._data = nullptr ;
        
        _mesh->_xgrid._size = +0 ;
        _mesh->_xgrid._data = nullptr ;
        
        _mesh->_ygrid._size = +0 ;
        _mesh->_ygrid._data = nullptr ;
        
        _mesh->_zgrid._size = +0 ;
        _mesh->_zgrid._data = nullptr ;
    }

#   endif//__lib_jigsaw

#   endif



