

    template <
        typename  tria_type ,
        typename  tria_list ,
        typename  iptr_list ,
        typename  pred_type
             >
    __normal_call void_type cavity_inflate (
        tria_type &_tria ,
        iptr_type  _nnew ,
        iptr_list &_tset ,
        tria_list &_tnew ,
        iptr_list &_tpop ,
        pred_type &_pred
        )
    {

        while (true)
        {
        
        }
        
    }
    
    template <
        typename  tria_type ,
        typename  tria_list ,
        typename  iptr_list ,
        typename  pred_type
             >
    __normal_call void_type cavity_deflate (
        tria_type &_tria ,
        iptr_type  _nnew ,
        iptr_list &_tset ,
        tria_list &_tnew ,
        iptr_list &_tpop ,
        pred_type &_pred
        )
    {
    }

    template <
        typename  tria_type ,
        typename  tria_list ,
        typename  iptr_list ,
        typename  pred_type
             >
    __normal_call void_type cavity_optimal (
        tria_type &_tria ,
        iptr_type  _star ,
        iptr_list &_tset ,
        tria_list &_tnew ,
        iptr_list &_tpop ,
        pred_type &_pred
        )
    {
    
        cavity_deflate(_tria, _star, _tset);
    
    
    }
    
    
    
