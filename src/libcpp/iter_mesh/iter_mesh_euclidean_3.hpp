
    /*
    --------------------------------------------------------
     * ITER-MESH-EUCLIDEAN-3: mesh obj. for ITER-MESH-3.
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
     * Last updated: 15 Jul., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_MESH_EUCLIDEAN_3__
#   define __ITER_MESH_EUCLIDEAN_3__

    namespace mesh {

template <
    typename R ,
    typename I
             >
    class iter_mesh_euclidean_3d
    {
    public  :

    /*-------------------- euclidean mesh for ITER-MESH-3 */

    typedef R                       real_type ;
    typedef I                       iptr_type ;

    typedef mesh::mesh_complex_edge_2 <
            iptr_type >             edge_base ;

    typedef mesh::mesh_complex_tria_3 <
            iptr_type >             tri3_base ;
    typedef mesh::mesh_complex_quad_4 <
            iptr_type >             quad_base ;

    typedef mesh::mesh_complex_tria_4 <
            iptr_type >             tri4_base ;
    typedef mesh::mesh_complex_hexa_8 <
            iptr_type >             hexa_base ;
    typedef mesh::mesh_complex_wedg_6 <
            iptr_type >             wedg_base ;
    typedef mesh::mesh_complex_pyra_5 <
            iptr_type >             pyra_base ;

    typedef mesh::mesh_complex_node_4 <
            iptr_type ,
            real_type >             node_base ;

    class node_type : public node_base
        {
    /*------------------------- node type for ITER-MESH-3 */
        public  :

        iptr_type                     _hidx = 0 ;
        iptr_type                     _itag = 0 ;

        char_type                     _fdim = 0 ;
        char_type                     _feat = 0 ;

        public  :

        __inline_call iptr_type      & hidx (
            )
        {   return  this->_hidx ;
        }
        __inline_call iptr_type const& hidx (
            ) const
        {   return  this->_hidx ;
        }
        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        __inline_call char_type      & fdim (
            )
        {   return  this->_fdim ;
        }
        __inline_call char_type const& fdim (
            ) const
        {   return  this->_fdim ;
        }
        __inline_call char_type      & feat (
            )
        {   return  this->_feat ;
        }
        __inline_call char_type const& feat (
            ) const
        {   return  this->_feat ;
        }
        } ;

    class edge_type : public edge_base
        {
    /*------------------------- edge type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    class tri3_type : public tri3_base
        {
    /*------------------------- face type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    class quad_type : public quad_base
        {
    /*------------------------- face type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    class tri4_type : public tri4_base
        {
    /*------------------------- cell type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    class hexa_type : public hexa_base
        {
    /*------------------------- cell type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    class wedg_type : public wedg_base
        {
    /*------------------------- cell type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    class pyra_type : public pyra_base
        {
    /*------------------------- cell type for ITER-MESH-3 */
        public  :

        iptr_type                     _itag = 0 ;

        public  :

        __inline_call iptr_type      & itag (
            )
        {   return  this->_itag ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return  this->_itag ;
        }
        } ;

    typedef mesh::mesh_complex_3 <
            node_type ,
            edge_type ,
            tri3_type , quad_type,
            tri4_type , hexa_type,
            wedg_type ,
            pyra_type            >  mesh_type ;

    public  :

    mesh_type                         _mesh ;

    public  :

    /*
    --------------------------------------------------------
     * INIT-MESH: impose (optional) offsets
    --------------------------------------------------------
     */

    __normal_call void_type init (
        float _xoff = + 0.f ,
        float _yoff = + 0.f ,
        float _zoff = + 0.f
        )
    {
        for (auto  _iter  =
             this->_mesh.node().head () ;
                   _iter !=
             this->_mesh.node().tend () ;
                 ++_iter  )
        {
            if (_iter->mark() >= +0)
            {
                _iter->pval(0) -= _xoff ;
                _iter->pval(1) -= _yoff ;
                _iter->pval(2) -= _zoff ;
            }
        }
    }

    } ;


    }

#   endif   //__ITER_MESH_EUCLIDEAN_3__



