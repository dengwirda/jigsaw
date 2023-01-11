
//  gcc -Wall -Wextra test_all.c -Xlinker -rpath=../lib
//  -L ../lib -ljigsaw -o test_all

//  Call all unit tests and bail-out on any failed runs

#   define __SKIP_MAIN__

#   include "test2d_a.c"
#   include "test2d_b.c"
#   include "test2d_c.c"
#   include "test2d_d.c"
#   include "test2d_e.c"
#   include "test2d_f.c"
#   include "test2d_g.c"
#   include "test2d_h.c"
#   include "test2d_i.c"

#   include "test2s_a.c"
#   include "test2s_b.c"
#   include "test2s_c.c"
//  include "test2s_d.c"
#   include "test2s_e.c"
#   include "test2s_f.c"
#   include "test2s_g.c"
#   include "test2s_h.c"
#   include "test2s_i.c"

#   include "test3d_a.c"
#   include "test3d_b.c"
#   include "test3d_c.c"
#   include "test3d_d.c"
#   include "test3d_e.c"
#   include "test3d_f.c"
#   include "test3d_g.c"
//  include "test3d_h.c"
//  include "test3d_i.c"

    int main ()
    {
        int _retv = 0;
        int _verb =-1;  // verbosity

    /*---------------------------------- run all 2d tests */
        if ((_retv = test2d_a(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_b(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_c(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_d(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_e(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_f(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_g(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_h(_verb)) != 0)
            return _retv ;
        if ((_retv = test2d_i(_verb)) != 0)
            return _retv ;

    /*---------------------------------- run all 2s tests */
        if ((_retv = test2s_a(_verb)) != 0)
            return _retv ;
        if ((_retv = test2s_b(_verb)) != 0)
            return _retv ;
        if ((_retv = test2s_c(_verb)) != 0)
            return _retv ;
    //  if ((_retv = test2s_d(_verb)) != 0)
    //      return _retv ;
        if ((_retv = test2s_e(_verb)) != 0)
            return _retv ;
        if ((_retv = test2s_f(_verb)) != 0)
            return _retv ;
        if ((_retv = test2s_g(_verb)) != 0)
            return _retv ;
        if ((_retv = test2s_h(_verb)) != 0)
            return _retv ;
        if ((_retv = test2s_i(_verb)) != 0)
            return _retv ;

    /*---------------------------------- run all 3d tests */
        if ((_retv = test3d_a(_verb)) != 0)
            return _retv ;
        if ((_retv = test3d_b(_verb)) != 0)
            return _retv ;
        if ((_retv = test3d_c(_verb)) != 0)
            return _retv ;
        if ((_retv = test3d_d(_verb)) != 0)
            return _retv ;
        if ((_retv = test3d_e(_verb)) != 0)
            return _retv ;
        if ((_retv = test3d_f(_verb)) != 0)
            return _retv ;
        if ((_retv = test3d_g(_verb)) != 0)
            return _retv ;
    //  if ((_retv = test3d_h(_verb)) != 0)
    //      return _retv ;
    //  if ((_retv = test3d_i(_verb)) != 0)
    //      return _retv ;

        return _retv ;
    }



