    
    /*
    --------------------------------------------------------
     * MATCH-GRAPH: algorithms for graph matching.
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
     * Last updated: 20 August, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
#   pragma once

#   ifndef __MATCH_GRAPH__
#   define __MATCH_GRAPH__

    namespace  graph {
    
    /*------------------------ base graph node/edge types */
    
    class graph_node_base
    {
    public  :
    } ;

    template <
    typename C,
    typename I
             >
    class graph_edge_base
    {
    public  :
        typedef C                       cost_type ;
        typedef I                       iptr_type ;
    
    public  :
        cost_type                      _cost ;
        
        iptr_type                      _inod ;
        iptr_type                      _jnod ;
    } ;
             
    /*
    --------------------------------------------------------
     * MATCH-GRAPH: greedy approx. max-weight matching.
    --------------------------------------------------------
     */
             
    template <
    typename E
             >
    class graph_match
    {
    public  :
    typedef  E                          edge_type ;
    
    typedef typename 
            edge_type::cost_type        cost_type ;
    typedef typename 
            edge_type::iptr_type        iptr_type ;
  
    class edge_comp
        {
        public  :
        __inline_call bool_type operator () (
            edge_type const& _idat,
            edge_type const& _jdat
            ) const
        {   return _idat._cost >
                   _jdat._cost ;
        }
        } ;
    
    typedef containers:: array_list <
            iptr_type               >   node_sets ;              
                    
    typedef containers::prioritymap <
            edge_type , 
            edge_comp               >   edge_heap ;
    
    typedef typename 
            edge_heap::kptr_type        hptr_type ;
    
    typedef containers:: array      < 
            hptr_type               >   heap_ptrs ;
    
    typedef containers:: array      < 
            iptr_type               >   iptr_list ;
    
    typedef containers:: array      < 
            edge_type               >   edge_list ;
    
    /*
    --------------------------------------------------------
     * EDGE-COST: calc. current 'cost' for edge.
    --------------------------------------------------------
     */
         
    __static_call
    __normal_call cost_type edge_cost (
        edge_list &_edge,
        iptr_list &_mark,
        node_sets &_nadj,
        iptr_type  _enum
        )
    {
        cost_type _cost = (cost_type)+0 ;
        
        cost_type _ival = (cost_type)+0 ;
        iptr_type _ideg = (iptr_type)+0 ;
            
        for (auto _epos  = _nadj[
            _edge[_enum]._inod].head();
                  _epos != _nadj[
            _edge[_enum]._inod].tend(); 
                ++_epos  )
        {
            if (_mark[*_epos] >= +0)
            {
            if (_enum !=  *_epos)
            {
                _ideg += +1 ;
                _ival -= 
                    _edge[*_epos]._cost ;
            }
            }
        }
        
        cost_type _jval = (cost_type)+0 ;
        iptr_type _jdeg = (iptr_type)+0 ;
        
        for (auto _epos  = _nadj[
            _edge[_enum]._jnod].head();
                  _epos != _nadj[
            _edge[_enum]._jnod].tend(); 
                ++_epos  )
        {
            if (_mark[*_epos] >= +0)
            {
            if (_enum !=  *_epos)
            {
                _jdeg += +1 ;
                _jval -= 
                    _edge[*_epos]._cost ;
            }
            }
        }

        _cost += _ival ;
        _cost += _jval ;

        return _cost+_edge[_enum]._cost ;            
    } 
    
    /*
    --------------------------------------------------------
     * GET-MATCH: find approx. max-weight matching.
    --------------------------------------------------------
     */
    
    __static_call
    __normal_call void_type get_match (
        edge_list &_edge,
        edge_list &_pair
        )
    {
        iptr_type _nmax =  (iptr_type)0 ;
        node_sets _nadj;
        edge_heap _heap; 
        heap_ptrs _hptr;
        iptr_list _nmrk, _emrk, _nset;
        
        _emrk.set_count(_edge.count(), 
            containers::loose_alloc, +0);
            
        _hptr.set_count(_edge.count(), 
            containers::loose_alloc, +0);
    
        {
        iptr_type _enum  = (iptr_type)0 ;
        for (auto _epos  = _edge.head() ;
                  _epos != _edge.tend() ;
                ++_epos, ++_enum )
        {
            _nadj.push (
                _enum, _epos->_inod) ;
            _nadj.push (
                _enum, _epos->_jnod) ;
            
            _nmax = std::max (
                _nmax, _epos->_inod) ;
            _nmax = std::max (
                _nmax, _epos->_jnod) ;
        }
        }
    
        {
        iptr_type _enum  = (iptr_type)0 ;
        for (auto _epos  = _edge.head() ;
                  _epos != _edge.tend() ;
                ++_epos, ++_enum )
        {
            edge_type  _edat;
            _edat._inod = _epos->_inod  ;
            _edat._jnod = _epos->_jnod  ;
            
            _edat._cost = edge_cost(
            _edge, _emrk, _nadj, _enum) ;  
        
            _hptr[*_epos] 
                = _heap.push( _edat) ;    
        }
        }
    
        iptr_type _pass = (iptr_type)+0 ;
        
        _nmrk.set_count(_nmax, 
            containers::loose_alloc, +0);
        
        for ( ; !_heap.empty() ; )
        {
            edge_type _edat;
            _heap._pop_root(_edat) ;
            
            _pair.push_tail(_edat) ;
            
            _pass +=    +1 ;
            
            _mark[_edat._inod] =  _pass ;
            _mark[_edat._jnod] =  _pass ;
    
    
            _nset.set_count ( +0) ;        
            
            for (auto _epos  = 
                _nadj[_edat._inod].head() ;
                      _epos != 
                _nadj[_edat._inod].tend() ; 
                    ++_epos  )
            {
                if (_emrk[*_epos] >= +0)
                {
                iptr_type _inod = 
                    _edge[*_epos]._inod;
                iptr_type _jnod = 
                    _edge[*_epos]._jnod;
                
                _emrk[*_epos] = 
                    -_emrk[*_epos]) ;
                
                if (_mark[_inod] != _pass)
                {
                    _mark[_inod]  = _pass;
                    _nset.push_tail(_inod);
                }
                if (_mark[_jnod] != _pass)
                {
                    _mark[_jnod]  = _pass;
                    _nset.push_tail(_jnod);
                }
                
                _heap._pop(_hptr[*_epos]);
                }
            }
            
            for (auto _epos  = 
                _nadj[_edat._jnod].head() ;
                      _epos != 
                _nadj[_edat._jnod].tend() ; 
                    ++_epos  )
            {
                if (_emrk[*_epos] >= +0)
                {
                iptr_type _inod = 
                    _edge[*_epos]._inod;
                iptr_type _jnod = 
                    _edge[*_epos]._jnod;
                
                _emrk[*_epos] = 
                    -_emrk[*_epos] ;
                
                if (_mark[_inod] != _pass)
                {
                    _mark[_inod]  = _pass;
                    _nset.push_tail(_inod);
                }
                if (_mark[_jnod] != _pass)
                {
                    _mark[_jnod]  = _pass;
                    _nset.push_tail(_jnod);
                }
                
                _heap._pop(_hptr[*_epos]);
                }
            }
    
    
            for (auto _node  =_nset.head();
                      _node !=_nset.tend();
                    ++_node  )
            {
            for (auto _epos  = 
                      _nadj[*_node].head();
                      _epos != 
                      _nadj[*_node].tend(); 
                    ++_epos  )
            {
                if (_emrk[*_epos] != _pass)
                {  
                edge_type  _enew ;
                _enew._inod = 
                    _edge[*_epos]._inod ;
                _enew._jnod = 
                    _edge[*_epos]._jnod ;
            
                _emrk[*_epos] =   _pass ;
            
                cost_type  _cost =
                    edge_cost (
                _edge, _emrk, _nadj, *_epos);
      
                _enew._cost = _cost;
      
                _heap.
                update(_hptr[*_epos], _enew);
                }
            }      
            }        
            
        }
    
    }
    
    } ;
    
    
    }
    
#   endif   //__MATCH_GRAPH__
    
    
    
