
    /*
    --------------------------------------------------------
     * ITER-MESH-EUCLIDEAN-2: mesh obj. for ITER-MESH-2.
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
     * Last updated: 24 Feb., 2021
     *
     * Copyright 2013-2021
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __ITER_MESH_EUCLIDEAN_2__
#   define __ITER_MESH_EUCLIDEAN_2__

    namespace mesh {

    template <
    typename R ,
    typename I
             >
    class iter_mesh_euclidean_2d
    {
    public  :

    /*-------------------- euclidean mesh for ITER-MESH-2 */

    typedef R                       real_type ;
    typedef I                       iptr_type ;

    typedef mesh::mesh_complex_edge_2 <
            iptr_type >             edge_base ;

    typedef mesh::mesh_complex_tria_3 <
            iptr_type >             tri3_base ;
    typedef mesh::mesh_complex_quad_4 <
            iptr_type >             quad_base ;

    typedef mesh::mesh_complex_node_3 <
            iptr_type ,
            real_type >             node_base ;

    class node_type : public node_base
        {
    /*------------------------- node type for ITER-MESH-2 */
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
    /*------------------------- edge type for ITER-MESH-2 */
        public  :

        iptr_type                     _itag ;

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
    /*------------------------- face type for ITER-MESH-2 */
        public  :

        iptr_type                     _itag ;

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
    /*------------------------- face type for ITER-MESH-2 */
        public  :

        iptr_type                     _itag ;

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

    typedef mesh::mesh_complex_2 <
            node_type ,
            edge_type ,
            tri3_type ,
            quad_type            >  mesh_type ;

    public  :

    mesh_type                         _mesh ;

    } ;


    }

#   endif   //__ITER_MESH_EUCLIDEAN_2__



