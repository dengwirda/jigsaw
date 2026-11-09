
    /*
    --------------------------------------------------------
     * MSH-READ: parse *.MSH file into mesh data.
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
     * Last updated: 09 Nov., 2026
     *
     * Copyright 2013-2026
     * Darren Engwirda
     * d.engwirda@gmail.com
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MSH_READ__
#   define __MSH_READ__

    class jmsh_reader_base
    {
/*------------------------ inherit from base to read file */
    public  :
    __normal_call void_type initialise() {}
    __normal_call void_type postscript() {}
    __normal_call void_type push_mshid (
        std::int32_t  /*_ftag*/,
        jmsh_kind::enum_data  /*_kind*/
        ) { }
    __normal_call void_type push_ndims (
        std:: size_t  /*_ndim*/
        ) { }
    __normal_call void_type push_radii (
        double      * /*_erad*/
        ) { }
    __normal_call void_type open_point (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_point (
        std:: size_t  /*_irow*/,
        double     *  /*_pval*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_seeds (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_seeds (
        std:: size_t  /*_irow*/,
        double     *  /*_pval*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_power (
        std:: size_t  /*_nrow*/,
        std:: size_t  /*_npwr*/
        ) { }
    __normal_call void_type push_power (
        std:: size_t  /*_irow*/,
        double     *  /*_xpwr*/
        ) { }
    __normal_call void_type open_edge2 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_edge2 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_tria3 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_tria3 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_quad4 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_quad4 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_tria4 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_tria4 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_hexa8 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_hexa8 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_wedg6 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_wedg6 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_pyra5 (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_pyra5 (
        std:: size_t  /*_irow*/,
        std::int32_t* /*_node*/,
        std::int32_t  /*_itag*/
        ) { }
    __normal_call void_type open_value (
        std:: size_t  /*_nrow*/,
        std:: size_t  /*_nval*/
        ) { }
    __normal_call void_type push_value (
        std:: size_t  /*_irow*/,
        float       * /*_vdat*/
        ) { }
    __normal_call void_type open_slope (
        std:: size_t  /*_nrow*/,
        std:: size_t  /*_nval*/
        ) { }
    __normal_call void_type push_slope (
        std:: size_t  /*_irow*/,
        float       * /*_vdat*/
        ) { }
    __normal_call void_type open_bound (
        std:: size_t  /*_nrow*/
        ) { }
    __normal_call void_type push_bound (
        std:: size_t  /*_irow*/,
        std::int32_t  /*_itag*/,
        std::int32_t  /*_inum*/,
        std::int32_t  /*_kind*/
        ) { }
    __normal_call void_type open_coord (
        std:: size_t  /*_idim*/,
        std:: size_t  /*_ndim*/
        ) { }
    __normal_call void_type push_coord (
        std:: size_t  /*_idim*/,
        std:: size_t  /*_irow*/,
        double        /*_ppos*/
        ) { }
    } ;

#   include "msh_read_txt.hpp"
#   include "msh_read_api.hpp"
#   include "msh_read_nc_.hpp"

#   endif   //__MSH_READ__



