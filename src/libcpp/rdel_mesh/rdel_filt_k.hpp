
    /*
    --------------------------------------------------------
     * RDEL-FILT-K: proximity filters for RDEL-MESH-K.
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

#   pragma once

#   ifndef __RDEL_FILT_K__
#   define __RDEL_FILT_K__

    namespace mesh {

    template <
    typename R,
    typename I
             >
    class keep_base_2d
    {
/*--------------------- "base" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;
     
    class node_data
        {
    /*----------------------- internal "co-ord." type */
        public  :
            real_type          _ppos[2]={(real_type)0.} ;
            
            char_type          _hits ;
            
            char_type          _feat ;
            char_type          _topo ;
            
            iptr_type          _itag ;          
        public  :
        __inline_call real_type      & pval (
            iptr_type     _ipos
            )
        {   return this-> _ppos[_ipos] ;
        }
        __inline_call char_type      & hits (
            )
        {   return this-> _hits ;
        }
        __inline_call char_type      & feat (
            )
        {   return this-> _feat ;
        }
        __inline_call char_type      & topo (
            )
        {   return this-> _topo ;
        }
        __inline_call iptr_type      & itag (
            )
        {   return this-> _itag ;
        }
        __inline_call real_type const& pval (
            iptr_type     _ipos
            ) const
        {   return this-> _ppos[_ipos] ;
        }
        __inline_call char_type const& hits (
            ) const
        {   return this-> _hits ;
        }
        __inline_call char_type const& feat (
            ) const
        {   return this-> _feat ;
        }
        __inline_call char_type const& topo (
            ) const
        {   return this-> _topo ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return this-> _itag ;
        }
        } ;  
    } ;

    template <
    typename R,
    typename I
             >
    class keep_base_3d
    {
/*-------------------------- intersection "base" type */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ; 
     
    class node_data
        {
    /*----------------------- internal "co-ord." type */
        public  :
            real_type          _ppos[3]={(real_type)0.} ;
            
            char_type          _hits ;
            
            char_type          _feat ;
            char_type          _topo ;
            
            iptr_type          _itag ;          
        public  :
        __inline_call real_type      & pval (
            iptr_type     _ipos
            )
        {   return this-> _ppos[_ipos] ;
        }
        __inline_call char_type      & hits (
            )
        {   return this-> _hits ;
        }
        __inline_call char_type      & feat (
            )
        {   return this-> _feat ;
        }
        __inline_call char_type      & topo (
            )
        {   return this-> _topo ;
        }
        __inline_call iptr_type      & itag (
            )
        {   return this-> _itag ;
        }
        __inline_call real_type const& pval (
            iptr_type     _ipos
            ) const
        {   return this-> _ppos[_ipos] ;
        }
        __inline_call char_type const& hits (
            ) const
        {   return this-> _hits ;
        }
        __inline_call char_type const& feat (
            ) const
        {   return this-> _feat ;
        }
        __inline_call char_type const& topo (
            ) const
        {   return this-> _topo ;
        }
        __inline_call iptr_type const& itag (
            ) const
        {   return this-> _itag ;
        }
        } ;  
    } ;

    template <
    typename R,
    typename I
             >
    class keep_all_2d : public keep_base_2d<R,I>
    {
/*----------------- "keep-all" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;
    
    typedef keep_base_2d <
            real_type , 
            iptr_type    >  pred_base ;

    public  :
    containers::array <
        typename pred_base::node_data>  _list ;
        
    iptr_type                   _inum ;
    bool_type                   _find ;
        
    public  :
/*------------------------------ construct from _src. */
    __inline_call keep_all_2d (
        ) 
    {
        this->_find = false;
        this->_inum =   +0 ;
    }
/*------------------------------ nullify intersection */
    __inline_call void_type clear (
        )
    {
        this->_list.clear();
    
        this->_find = false;
        this->_inum =   +0 ;
    }
/*------------------------------ process intersection */
    __inline_call 
        void_type operator()  (
    __const_ptr ( real_type) _xpos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag
        )
    {
        this->_inum +=  +1 ;
        
        this->_list.push_tail() ;
        this->_list.tail()
            ->hits () = _hits ;
        this->_list.tail()
            ->feat () = _feat ;
        this->_list.tail()
            ->topo () = _topo ;
        this->_list.tail()
            ->itag () = _itag ;
        
        this->_list.tail()
            ->pval(0) = _xpos[0];
        this->_list.tail()
            ->pval(1) = _xpos[1];
        
        this->_find = true ;
    }  
    } ;
    
    template <
    typename R,
    typename I
             >
    class keep_all_3d : public keep_base_3d<R,I>
    {
/*----------------- "keep-all" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;

    typedef keep_base_3d <
            real_type , 
            iptr_type    >  pred_base ;

    public  :
    containers::array <
        typename pred_base::node_data>  _list ;
        
    iptr_type                   _inum ;
    bool_type                   _find ;
    
    public  :
/*------------------------------ construct from _src. */
    __inline_call keep_all_3d (
        ) 
    {
        this->_find = false;
        this->_inum =   +0 ;
    }
/*------------------------------ nullify intersection */
    __inline_call void_type clear (
        )
    {
        this->_list.clear();
    
        this->_find = false;
        this->_inum =   +0 ;
    }
/*------------------------------ process intersection */
    __inline_call 
        void_type operator()  (
    __const_ptr ( real_type) _xpos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag
        )
    {
        this->_inum +=  +1 ;
        
        this->_list.push_tail() ;
        this->_list.tail()
            ->hits () = _hits ;
        this->_list.tail()
            ->feat () = _feat ;
        this->_list.tail()
            ->topo () = _topo ;
        this->_list.tail()
            ->itag () = _itag ;
        
        this->_list.tail()
            ->pval(0) = _xpos[0];
        this->_list.tail()
            ->pval(1) = _xpos[1];
        this->_list.tail()
            ->pval(2) = _xpos[2];
        
        this->_find = true ;
    }  
    } ;

    template <
    typename R,
    typename I
             >
    class keep_minmax_length_2d 
            : public keep_base_2d<R, I>
    {
/*--------- "min/max-distance" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;

    typedef keep_base_3d <
            real_type  , 
            iptr_type    >  pred_base ;

    public  :
    containers::
    fixed_array<real_type, +2>  _ppos ;
    
    typename 
    pred_base::node_data        _pmin ;
    typename 
    pred_base::node_data        _pmax ;

    real_type                   _dmin ;
    real_type                   _dmax ;

    iptr_type                   _inum ;
    
    bool_type                   _find ;

    public  :
/*------------------------------ construct from _src. */
    __inline_call keep_minmax_length_2d (
    __const_ptr ( real_type) _psrc
        ) 
    {
        this->_ppos[ 0] = _psrc[ 0] ;
        this->_ppos[ 1] = _psrc[ 1] ;

        this->_inum    = +0 ;

        this->_dmin = 
    +std::numeric_limits<real_type>::infinity();
        this->_dmax = 
    -std::numeric_limits<real_type>::infinity();
    
        this->_find    = false ;
    }

/*------------------------------ process intersection */
    __inline_call void_type operator() (
    __const_ptr ( real_type) _xpos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag
        )
    {
        this->_inum += +1 ;
    /*-------------------------------- calc. distance */
        real_type _dist = 
            geometry::lensqr_2d (
                _xpos, &this->_ppos[0]);
    /*------------------------------------- keep min. */
        if (this->_dmin  > _dist)
        {
            this->
           _pmin.pval(0) = _xpos[0];
            this->
           _pmin.pval(1) = _xpos[1];
            
            this->
           _pmin.hits () = _hits;
            this->
           _pmin.feat () = _feat;
            this->
           _pmin.topo () = _topo; 
            this->
           _pmin.itag () = _itag;
        
            this->_dmin  = _dist;
            
            this->_find  =  true;
        }
    /*------------------------------------- keep max. */
        if (this->_dmax  < _dist ||
           (this->_dmax == _dist && 
            this->
           _pmax.hits()  < _hits) )
        {
            this->
           _pmax.pval(0) = _xpos[0];
            this->
           _pmax.pval(1) = _xpos[1];
            
            this->
           _pmax.hits () = _hits;
            this->
           _pmax.feat () = _feat;
            this->
           _pmax.topo () = _topo; 
            this->
           _pmax.itag () = _itag;
        
            this->_dmax  = _dist;
            
            this->_find  =  true;
        }
    }
    } ;

    template <
    typename R,
    typename I
             >
    class keep_minmax_length_3d
            : public keep_base_3d<R, I>
    {
/*--------- "min/max-distance" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;

    typedef keep_base_3d <
            real_type  , 
            iptr_type    >  pred_base ;

    public  :
    containers::
    fixed_array<real_type, +3>  _ppos ;
    
    typename 
    pred_base::node_data        _pmin ;
    typename 
    pred_base::node_data        _pmax ;

    real_type                   _dmin ;
    real_type                   _dmax ;

    iptr_type                   _inum ;
    
    bool_type                   _find ;

    public  :
/*------------------------------ construct from _src. */
    __inline_call keep_minmax_length_3d (
    __const_ptr ( real_type) _psrc
        ) 
    {
        this->_ppos[ 0] = _psrc[ 0] ;
        this->_ppos[ 1] = _psrc[ 1] ;
        this->_ppos[ 2] = _psrc[ 2] ;

        this->_inum    = +0 ;

        this->_dmin = 
    +std::numeric_limits<real_type>::infinity();
        this->_dmax = 
    -std::numeric_limits<real_type>::infinity();
    
        this->_pmin.hits() = 
            geometry::null_hits;
        this->_pmax.hits() = 
            geometry::null_hits;
    
        this->_find    = false ;
    }

/*------------------------------ process intersection */
    __inline_call void_type operator() (
    __const_ptr ( real_type) _xpos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag
        )
    {
        this->_inum += +1 ;
    /*-------------------------------- calc. distance */
        real_type _dist = 
            geometry::lensqr_3d (
                _xpos, &this->_ppos[0]);
    /*------------------------------------- keep min. */
        if (this->_dmin  > _dist)
        {
            this->
           _pmin.pval(0) = _xpos[0];
            this->
           _pmin.pval(1) = _xpos[1];
            this->
           _pmin.pval(2) = _xpos[2];
        
            this->
           _pmin.hits () = _hits;
            this->
           _pmin.feat () = _feat;
            this->
           _pmin.topo () = _topo; 
            this->
           _pmin.itag () = _itag;
        
            this->_dmin  = _dist;
            
            this->_find  =  true;
        }
    /*------------------------------------- keep max. */
        if (this->_dmax  < _dist)
        {
            this->
           _pmax.pval(0) = _xpos[0];
            this->
           _pmax.pval(1) = _xpos[1];
            this->
           _pmax.pval(2) = _xpos[2];
            
            this->
           _pmax.hits () = _hits;
            this->
           _pmax.feat () = _feat;
            this->
           _pmax.topo () = _topo; 
            this->
           _pmax.itag () = _itag;
        
            this->_dmax  = _dist;
            
            this->_find  =  true;
        }
    }
    } ;

    template <
    typename R,
    typename I
             >
    class keep_minmax_cosine_2d
            : public keep_base_2d<R, I>
    {
/*--------------- "max-cosine" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;

    typedef keep_base_3d <
            real_type  , 
            iptr_type    >  pred_base ;

    public  :
    containers::
    fixed_array<real_type, +2>  _ppos ;
    
    containers::
    fixed_array<real_type, +2>  _dvec ;

    typename
    pred_base::node_data        _proj ;

    real_type                   _cmin ;
    real_type                   _best ;

    iptr_type                   _inum ;
    
    bool_type                   _find ;

    public  :
/*------------------------------ construct from _src. */
    __inline_call keep_minmax_cosine_2d (
    __const_ptr ( real_type) _psrc,
    __const_ptr ( real_type) _vsrc,
        real_type const&_cmin = real_type(+.50)
        ) 
    {
        this->_find     = false;
        this->_inum     = +0 ;
        
        this->_ppos[ 0] = _psrc[ 0];
        this->_ppos[ 1] = _psrc[ 1];
        
        this->_dvec[ 0] = _vsrc[ 0];
        this->_dvec[ 1] = _vsrc[ 1];
        
        this->_cmin     = _cmin;
        this->_best     = _cmin;
    }

/*------------------------------ process intersection */
    __inline_call void_type operator() (
    __const_ptr ( real_type) _xpos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag
        )
    {
    /*---------------------------------- calc. cosine */
        real_type _xvec[ 2];
        _xvec[ 0] =
        _xpos[ 0] - this->_ppos[ 0];
        _xvec[ 1] =
        _xpos[ 1] - this->_ppos[ 1];

        real_type _vcos = 
            geometry::cosine_2d (
               &this->_dvec[0], _xvec) ;

    /*---------------------------------- keep maximum */
        if (this->_cmin  < _vcos)
        {
            this->_inum += +1 ;
        }
        if (this->_best  < _vcos)
        {
            this->
           _proj.pval(0) = _xpos[0];
            this->
           _proj.pval(1) = _xpos[1];
            
            this->
           _proj.hits () = _hits;
            this->
           _proj.feat () = _feat;
            this->
           _proj.topo () = _topo; 
            this->
           _proj.itag () = _itag;
        
            this->_best  = _vcos;
            
            this->_find  =  true;
        }
    }
    } ;

    template <
    typename R,
    typename I
             >
    class keep_minmax_cosine_3d
            : public keep_base_3d<R, I>
    {
/*--------------- "max-cosine" intersection predicate */
    public  :   
    typedef R               real_type ;
    typedef I               iptr_type ;

    typedef keep_base_3d <
            real_type  , 
            iptr_type    >  pred_base ;

    public  :
    containers::
    fixed_array<real_type, +3>  _ppos ;
    
    containers::
    fixed_array<real_type, +3>  _dvec ;

    typename
    pred_base::node_data        _proj ;

    real_type                   _cmin ;
    real_type                   _best ;

    iptr_type                   _inum ;
    
    bool_type                   _find ;

    public  :
/*------------------------------ construct from _src. */
    __inline_call keep_minmax_cosine_3d (
    __const_ptr ( real_type) _psrc,
    __const_ptr ( real_type) _vsrc,
        real_type const&_cmin = real_type(+.50)
        ) 
    {
        this->_find     = false;
        this->_inum     = +0 ;
        
        this->_ppos[ 0] = _psrc[ 0];
        this->_ppos[ 1] = _psrc[ 1];
        this->_ppos[ 2] = _psrc[ 2];

        this->_dvec[ 0] = _vsrc[ 0];
        this->_dvec[ 1] = _vsrc[ 1];
        this->_dvec[ 2] = _vsrc[ 2];
        
        this->_cmin     = _cmin;
        this->_best     = _cmin;
    }

/*------------------------------ process intersection */
    __inline_call void_type operator() (
    __const_ptr ( real_type) _xpos,
        char_type _hits,
        char_type _feat,
        char_type _topo,
        iptr_type _itag
        )
    {
    /*---------------------------------- calc. cosine */
        real_type _xvec[ 3];
        _xvec[ 0] =
        _xpos[ 0] - this->_ppos[ 0];
        _xvec[ 1] =
        _xpos[ 1] - this->_ppos[ 1];
        _xvec[ 2] =
        _xpos[ 2] - this->_ppos[ 2];

        real_type _vcos = 
            geometry::cosine_3d (
               &this->_dvec[0], _xvec) ;

    /*---------------------------------- keep maximum */
        if (this->_cmin  < _vcos)
        {
            this->_inum += +1 ;
        }
        if (this->_best  < _vcos)
        {
            this->
           _proj.pval(0) = _xpos[0];
            this->
           _proj.pval(1) = _xpos[1];
            this->
           _proj.pval(2) = _xpos[2];
      
            this->
           _proj.hits () = _hits;  
            this->
           _proj.feat () = _feat;
            this->
           _proj.topo () = _topo; 
            this->
           _proj.itag () = _itag;
        
            this->_best  = _vcos;
            
            this->_find  =  true;
        }
    }
    } ;
    
    
    }

#   endif   //__RDEL_FILT_K__
    
    
    
