
    /*
    --------------------------------------------------------
     * SAVE-MSH_T: save msh_t for lib_jigsaw.
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
     * Last updated: 06 June, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __SAVE_MSH_T__
#   define __SAVE_MSH_T__

    /*
    --------------------------------------------------------
     * SAVE-MSH_T: save msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

    indx_t jigsaw_save_msh_t (              // save msh_t
        char *_fchr,
        jigsaw_msh_t *_jmsh
        )
    {
        iptr_type _errv = __no_error ;

        __unreferenced(_fchr) ;
        __unreferenced(_jmsh) ;
    
        //!! TODO
    
    
        return  _errv ;
    }

#   endif//__lib_jigsaw

#   endif



