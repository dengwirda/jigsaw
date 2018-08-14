
    /* 
    --------------------------------------------------------
     * a bit of template meta-magic on iterator types.
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
     * Last updated: 03 May, 2017
     *
     * Copyright 2013-2017
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */
     
#   pragma once

#   ifndef __ITER_BASE__
#   define __ITER_BASE__

#   include <type_traits>

#   include <cstddef>

    namespace containers {

    /*
    --------------------------------------------------------
     * iterator concepts via typedefs and inheritence.
    --------------------------------------------------------
     */
     
    class null_iterator_kind 
        { /*------ non-iterator types */
        } ;
    class base_iterator_kind
        { /*------ all iterator types */
        } ;
    class single_iterator_kind
        : public   base_iterator_kind
        { /* unidirectional iterators */
        } ;
    class double_iterator_kind 
        : public single_iterator_kind
        { /* bi-directional iterators */
        } ;
    class random_iterator_kind 
        : public double_iterator_kind
        { /* random--access iterators */
        } ;
    class pointer_iterator_kind
        : public random_iterator_kind
        { /* iterators as raw pointer */
        } ;

    class iterator_base
        { /*--------- "non-iterator" base type */
        public  :
        typedef   null_iterator_kind    iter_kind ;
        } ;
    class single_iterator_base
        : public iterator_base
        { /* unidirectional iterator base type */
        public  :
        typedef single_iterator_kind    iter_kind ;
        } ;
    class double_iterator_base
        : public iterator_base
        { /* bi-directional iterator base type */
        public  :
        typedef double_iterator_kind    iter_kind ;
        } ;
    class random_iterator_base
        : public iterator_base
        { /* random--access iterator base type */
        public  :
        typedef random_iterator_kind    iter_kind ;
        } ;

    /*
    --------------------------------------------------------
     * traits class to distinguish iterator types.
    --------------------------------------------------------
     */
     
    template <
    typename iter_type , 
    bool     flag_kind
             >
    class iterator_traits_base
        { 
    /*-------------- real iterator types store their kind */
        public  :
        typedef typename 
                iter_type::iter_kind    iter_kind ;
        typedef typename 
                iter_type::data_type    data_type ;
        typedef typename 
                iter_type::diff_type    diff_type ;
        typedef typename 
                iter_type::size_type    size_type ;
        } ;

    template <
    typename iter_type
             >
    class iterator_traits_base<iter_type, false >
        { 
    /*-------------- noniterator types assigned null kind */
        public  :
        typedef null_iterator_kind      iter_kind ;
        typedef typename 
                std::remove_pointer<
                    iter_type>::type    data_type ;
        typedef std::size_t             size_type ;
        typedef std::ptrdiff_t          diff_type ;
        } ;

    template <
    typename iter_type
             >
    class iterator_traits: public 
            iterator_traits_base <iter_type ,
              std::is_base_of
           <iterator_base,iter_type>::value >
        { 
    /*--------------------------------- just wraps base!! */ 
        } ;

    /*
    --------------------------------------------------------
     * ITER-KIND: traits fun. to return dummy kind.
    --------------------------------------------------------
     */
     
    template <
    typename     iter_type
             >
    __inline_call typename 
        iterator_traits<iter_type>::iter_kind 
        iter_kind (iter_type)
    {   return (  typename 
        iterator_traits<iter_type>::iter_kind()
             ) ;
    }


    }    //} containers

#   endif//__ITER_BASE__



