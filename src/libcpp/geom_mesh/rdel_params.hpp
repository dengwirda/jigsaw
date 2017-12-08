
    /*
    --------------------------------------------------------
     * RDEL-PARAMS: parameters for RDEL-MESH-K.
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
     * Last updated: 04 September, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __RDEL_PARAMS__
#   define __RDEL_PARAMS__

    namespace mesh {

    /*
    --------------------------------------------------------
     * RDEL-PARAMS: user-param's for RDEL-MESH-K
    --------------------------------------------------------
     */
    
    template <
    typename R , 
    typename I
             >
    class rdel_params
        {
        public  :
        
        typedef  R                      real_type ;
        typedef  I                      iptr_type ;
        
        typedef mesh::rdel_params<R, I> self_type ;
        
        enum node_kind{ null_kind , 
             fail_kind, circ_kind , 
             offh_kind, offc_kind , 
             sink_kind, disk_kind , 
             last_kind} ;

        iptr_type           _verb ; // logfile output verbosity

        iptr_type           _seed ; // no. init. "seed" nodes
        
        real_type           _phi1 ; // 1-"hard" angle tolerance
        real_type           _phi2 ; // 2-"hard" angle tolerance
        
        real_type           _eta1 ; // 1-"soft" angle tolerance
        real_type           _eta2 ; // 2-"soft" angle tolerance

        bool_type           _feat ; // true for feature calc.

        iptr_type           _dims ; // topo. dimensions to mesh
        
        real_type           _siz1 ; // 1-dim. element size mul.
        real_type           _siz2 ; // 2-dim. element size mul.
        real_type           _siz3 ; // 3-dim. element size mul.
        
        real_type           _eps1 ; // 1-dim. 1-hausdorff error
        real_type           _eps2 ; // 2-dim. 2-hausdorff error
        
        real_type           _rad2 ; // 2-dim. radius-edge ratio
        real_type           _rad3 ; // 3-dim. radius-edge ratio
        
        real_type           _off2 ; // 2-off. radius-edge ratio
        real_type           _off3 ; // 3-off. radius-edge ratio
        
        real_type           _snk2 ; // 2-dim. sink "safe" ratio
        real_type           _snk3 ; // 3-dim. sink "safe" ratio

        real_type           _vol3 ; // volume-length ratio

        bool_type           _top1 ; // impose "1-manifold-ness"
        bool_type           _top2 ; // impose "2-manifold-ness"

        iptr_type           _iter ; // max. no. refinement iter.
 
        public  :
            
        __static_call
        __inline_call real_type init_siz1 (
            )
        {   return 1. * 4./3. ;
        }
        __static_call
        __inline_call real_type init_siz2 (
            )
        {   return .5 *(4./3. + 
                2./(1.+std::sqrt(1./3.))) ;
        }
        __static_call
        __inline_call real_type init_siz3 (
            )
        {   return .5 *(4./3. + 
                2./(1.+std::sqrt(3./8.))) ;
        }
         
        __static_call
        __inline_call iptr_type init_iter (
            )
        {   return iptr_type (
        std::numeric_limits<iptr_type>::max()) ;
        }
            
        public  : 
        
    /*-------------------------- construct default param. */
        __inline_call rdel_params  (
            ) : 
            _verb(iptr_type(+   0)) ,
            
            _seed(iptr_type(+   8)) ,
                
            _phi1(real_type(+ 60.)) ,
            _phi2(real_type(+ 60.)) ,
            
            _eta1(real_type(+ 45.)) ,
            _eta2(real_type(+ 45.)) ,
            
            _feat(bool_type(false)) ,
            
            _dims(iptr_type(+   3)) ,
            
            _siz1(init_siz1())  ,
            _siz2(init_siz2())  ,
            _siz3(init_siz3())  ,
            
            _eps1(real_type(+.333)) ,
            _eps2(real_type(+.333)) ,
            
            _rad2(real_type(+1.05)) ,
            _rad3(real_type(+2.05)) ,
            
            _off2(real_type(+0.90)) ,
            _off3(real_type(+1.10)) ,
            
            _snk2(real_type(+.200)) ,
            _snk3(real_type(+.333)) ,
            
            _vol3(real_type(+.000)) ,
            
            _top1(bool_type(false)) ,
            _top2(bool_type(false)) ,
     
            _iter(init_iter())              
        {   // load default values
        }
    
    /*------------------------------------ "write" access */
        __inline_call iptr_type      & verb (
            )
        {   return  this->_verb ;
        }
        
        __inline_call iptr_type      & iter (
            )
        {   return  this->_iter ;
        }
        
        __inline_call iptr_type      & seed (
            )
        {   return  this->_seed ;
        }
        
        __inline_call real_type      & phi1 (
            )
        {   return  this->_phi1 ;
        }     
        __inline_call real_type      & phi2 (
            )
        {   return  this->_phi2 ;
        }
        
        __inline_call real_type      & eta1 (
            )
        {   return  this->_eta1 ;
        }
        __inline_call real_type      & eta2 (
            )
        {   return  this->_eta2 ;
        }
        
        __inline_call bool_type      & feat (
            )
        {   return  this->_feat ;
        }
        
        __inline_call iptr_type      & dims (
            )
        {   return  this->_dims ;
        }

        __inline_call real_type      & siz1 (
            )
        {   return  this->_siz1 ;
        }
        __inline_call real_type      & siz2 (
            )
        {   return  this->_siz2 ;
        }
        __inline_call real_type      & siz3 (
            )
        {   return  this->_siz3 ;
        }
        
        __inline_call real_type      & eps1 (
            )
        {   return  this->_eps1 ;
        }
        __inline_call real_type      & eps2 (
            )
        {   return  this->_eps2 ;
        }
        
        __inline_call real_type      & rad2 (
            )
        {   return  this->_rad2 ;
        }
        __inline_call real_type      & rad3 (
            )
        {   return  this->_rad3 ;
        }
        
        __inline_call real_type      & off2 (
            )
        {   return  this->_off2 ;
        }
        __inline_call real_type      & off3 (
            )
        {   return  this->_off3 ;
        }
        
        __inline_call real_type      & snk2 (
            )
        {   return  this->_snk2 ;
        }
        __inline_call real_type      & snk3 (
            )
        {   return  this->_snk3 ;
        }
        
        __inline_call real_type      & vol3 (
            )
        {   return  this->_vol3 ;
        }

        __inline_call bool_type      & top1 (
            )
        {   return  this->_top1 ;
        }
        __inline_call bool_type      & top2 (
            )
        {   return  this->_top2 ;
        }
        
    /*------------------------------------ "const" access */
        __inline_call iptr_type const& verb (
            ) const
        {   return  this->_verb ;
        }
        
        __inline_call iptr_type const& iter (
            ) const
        {   return  this->_iter ;
        }
        
        __inline_call iptr_type const& seed (
            ) const
        {   return  this->_seed ;
        }
        
        __inline_call real_type const& phi1 (
            ) const
        {   return  this->_phi1 ;
        }     
        __inline_call real_type const& phi2 (
            ) const
        {   return  this->_phi2 ;
        }
        
        __inline_call real_type const& eta1 (
            ) const
        {   return  this->_eta1 ;
        }
        __inline_call real_type const& eta2 (
            ) const
        {   return  this->_eta2 ;
        }
        
        __inline_call bool_type const& feat (
            ) const
        {   return  this->_feat ;
        }
        
        __inline_call iptr_type const& dims (
            ) const
        {   return  this->_dims ;
        }

        __inline_call real_type const& siz1 (
            ) const
        {   return  this->_siz1 ;
        }
        __inline_call real_type const& siz2 (
            ) const
        {   return  this->_siz2 ;
        }
        __inline_call real_type const& siz3 (
            ) const
        {   return  this->_siz3 ;
        }
        
        __inline_call real_type const& eps1 (
            ) const
        {   return  this->_eps1 ;
        }
        __inline_call real_type const& eps2 (
            ) const
        {   return  this->_eps2 ;
        }
        
        __inline_call real_type const& rad2 (
            ) const
        {   return  this->_rad2 ;
        }
        __inline_call real_type const& rad3 (
            ) const
        {   return  this->_rad3 ;
        }
        
        __inline_call real_type const& off2 (
            ) const
        {   return  this->_off2 ;
        }
        __inline_call real_type const& off3 (
            ) const
        {   return  this->_off3 ;
        }
        
        __inline_call real_type const& snk2 (
            ) const
        {   return  this->_snk2 ;
        }
        __inline_call real_type const& snk3 (
            ) const
        {   return  this->_snk3 ;
        }
        
        __inline_call real_type const& vol3 (
            ) const
        {   return  this->_vol3 ;
        }

        __inline_call bool_type const& top1 (
            ) const
        {   return  this->_top1 ;
        }
        __inline_call bool_type const& top2 (
            ) const
        {   return  this->_top2 ;
        }
        
        } ;

    }

#   endif // __RDEL_PARAMS__



