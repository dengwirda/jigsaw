
    template <
    typename F
             >
    class hj_mesh_2
    {
    
    __static_call 
    __normal_call void_type limit_edge_2 (
        )
    {
    }
    
    __static_call 
    __normal_call void_type limit_tria_3 (
        )
    {
    }
    
    __static_call
    __normal_call void_type limit_mesh (
        ffun_type &_ffun,
        real_type  _DFDX
        )
    {
    
        containers::priority_map<> _sort;
    
        // push nodes onto pq
        // pop min at each pass
        // limit any cell neighbours
        // if a cell is limited, update pq of its nodes
        
        // that's it!
    
    }
     
    } ;
    
    
    
