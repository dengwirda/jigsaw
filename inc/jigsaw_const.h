
    /*
    --------------------------------------------------------
     * Constants for the JIGSAW meshing library.
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
     * Last updated: 20 December, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * darren.engwirda@columbia.edu
     * https://github.com/dengwirda
     *
    --------------------------------------------------------
     */
     
    /* 
    --------------------------------------------------------
     * return codes for JIGSAW.
    --------------------------------------------------------
     */ 

#       define JIGSAW_UNKNOWN_ERROR    -1

#       define JIGSAW_NO_ERROR         +0

#       define JIGSAW_FILE_NOT_LOCATED +2
#       define JIGSAW_FILE_NOT_CREATED +3

#       define JIGSAW_INVALID_ARGUMENT +4

    /* 
    --------------------------------------------------------
     * constants for JIGSAW.
    --------------------------------------------------------
     */ 

#       define JIGSAW_NULL_FLAG      -100

#       define JIGSAW_EUCLIDEAN_MESH +100
#       define JIGSAW_EUCLIDEAN_GRID +101
#       define JIGSAW_EUCLIDEAN_DUAL +102
    
#       define JIGSAW_ELLIPSOID_MESH +200
#       define JIGSAW_ELLIPSOID_GRID +201
#       define JIGSAW_ELLIPSOID_DUAL +202

#       define JIGSAW_POINT_TAG      +10
#       define JIGSAW_EDGE2_TAG      +20
#       define JIGSAW_TRIA3_TAG      +30
#       define JIGSAW_QUAD4_TAG      +40
#       define JIGSAW_TRIA4_TAG      +50
#       define JIGSAW_HEXA8_TAG      +60
#       define JIGSAW_WEDG6_TAG      +70
#       define JIGSAW_PYRA5_TAG      +80

#       define JIGSAW_HFUN_RELATIVE  +300
#       define JIGSAW_HFUN_ABSOLUTE  +301
    
#       define JIGSAW_KERN_DELFRONT  +400
#       define JIGSAW_KERN_DELAUNAY  +401

#       define JIGSAW_BNDS_TRIACELL  +402
#       define JIGSAW_BNDS_DUALCELL  +403
     
     
    
