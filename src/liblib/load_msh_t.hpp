
    /*
    --------------------------------------------------------
     * LOAD-MSH_T: load msh_t for lib_jigsaw.
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
     * Last updated: 31 July, 2018
     *
     * Copyright 2013-2018
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __LOAD_MSH_T__
#   define __LOAD_MSH_T__

    /*
    --------------------------------------------------------
     * LOAD-MSH_T: load msh_t for lib_jigsaw.
    --------------------------------------------------------
     */

#   ifdef  __lib_jigsaw

#   include  "../msh_read.hpp"

    indx_t jigsaw_load_msh_t (              // load msh_t
        char *_fchr,
        jigsaw_msh_t *_jmsh
        )
    {
        class msht_reader: 
            public jmsh_reader_base
        {
        public  :
            jigsaw_msh_t        *_jmsh ;
            
            std::int32_t         _errv ;
            
            std::int32_t         _ftag ;
            jmsh_kind::
            enum_data            _kind ;
            std::int32_t         _ndim ;
            std::int32_t         _nval ;
            
        public  :
        __normal_call msht_reader (
            jigsaw_msh_t *_jsrc =  nullptr
            ) : _jmsh(_jsrc) , 
                _errv(  +0 ) {}            
    /*-------------------------------- read MSHID section */
        __normal_call void_type push_mshid (
            std::int32_t  _ftag ,
            jmsh_kind::enum_data _kind
            )
        {   
            this->_ftag = _ftag ;
            this->_kind = _kind ;
            
            this->_jmsh->_flags  = JIGSAW_NULL_FLAG;
            
            if (_kind == 
                 jmsh_kind::euclidean_mesh)
            this->_jmsh->_flags  = 
                JIGSAW_EUCLIDEAN_MESH ;
            else
            if (_kind == 
                 jmsh_kind::euclidean_grid)
            this->_jmsh->_flags  = 
                JIGSAW_EUCLIDEAN_GRID ;
            else
            if (_kind == 
                 jmsh_kind::ellipsoid_mesh)
            this->_jmsh->_flags  = 
                JIGSAW_ELLIPSOID_MESH ;
            else
            if (_kind == 
                 jmsh_kind::ellipsoid_grid)
            this->_jmsh->_flags  = 
                JIGSAW_ELLIPSOID_GRID ;
        }
    /*-------------------------------- read NDIMS section */
        __normal_call void_type push_ndims (
            std::int32_t  _ndim
            )
        {   
            this->_ndim = _ndim ;
        }
    /*-------------------------------- read RADII section */
        __normal_call void_type push_radii (
            double      *_erad
            ) 
        { 
            jigsaw_alloc_reals (
               &this->_jmsh->_radii, +3  ) ;
        
            this->_jmsh->
                _radii._data[0] = _erad[0] ;
            this->_jmsh->
                _radii._data[1] = _erad[1] ;
            this->_jmsh->
                _radii._data[2] = _erad[2] ;
        } 
    /*-------------------------------- open POINT section */
        __normal_call void_type open_point (
            std::int32_t  _nrow
            ) 
        { 
            if (this->_kind == 
                jmsh_kind:: euclidean_mesh)
            {
            
            if (this->_ndim == +2)
            jigsaw_alloc_vert2 (
               &this->_jmsh->_vert2, _nrow) ;
            else
            if (this->_ndim == +3)
            jigsaw_alloc_vert3 (
               &this->_jmsh->_vert3, _nrow) ;
               
            jigsaw_alloc_reals (
               &this->_jmsh->_power, _nrow) ;
            
            }
            else
            if (this->_kind == 
                jmsh_kind:: ellipsoid_mesh)
            {
            
            if (this->_ndim == +2)
            jigsaw_alloc_vert2 (
               &this->_jmsh->_vert2, _nrow) ;
            else
            if (this->_ndim == +3)
            jigsaw_alloc_vert3 (
               &this->_jmsh->_vert3, _nrow) ;
            
            jigsaw_alloc_reals (
               &this->_jmsh->_power, _nrow) ;
            
            }
        }
    /*-------------------------------- read POINT section */
        __normal_call void_type push_point (
            std::int32_t  _ipos ,
            double       *_pval ,
            std::int32_t  _itag
            )
        {
            if (this->_ndim == +2 )
            {
            
            if (_ipos < this->_jmsh->_vert2._size)
            {
            this->_jmsh->_vert2._data[_ipos].
                _ppos[0] = _pval[0] ;
            this->_jmsh->_vert2._data[_ipos].
                _ppos[1] = _pval[1] ;
                
            this->_jmsh->_vert2.
                _data[_ipos]._itag  = _itag ;
                 
            this->_jmsh->
                _power._data[_ipos] = +0.00 ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            } 
                
            }           
            else
            if (this->_ndim == +3 )
            {
            
            if (_ipos < this->_jmsh->_vert2._size)
            {
            this->_jmsh->_vert3._data[_ipos].
                _ppos[0] = _pval[0] ;
            this->_jmsh->_vert3._data[_ipos].
                _ppos[1] = _pval[1] ;
            this->_jmsh->_vert3._data[_ipos].
                _ppos[2] = _pval[2] ;
                
            this->_jmsh->_vert3.
                _data[_ipos]._itag  = _itag ;
                 
            this->_jmsh->
                _power._data[_ipos] = +0.00 ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            }
                 
            }
        }       
    /*-------------------------------- read POWER section */
        __normal_call void_type push_power (
            std::int32_t _ipos,
            double     * _xpwr
            ) 
        { 
            if (this->_ndim == +2 )
            {
            
            if (_ipos < this->_jmsh->_vert2._size)
            {
            this->_jmsh->_power.
                _data[_ipos] = _xpwr[ +0] ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            } 
                
            }           
            else
            if (this->_ndim == +3 )
            {
            
            if (_ipos < this->_jmsh->_vert3._size)
            {
            this->_jmsh->_power.
                _data[_ipos] = _xpwr[ +0] ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            }
                 
            }
        }    
    /*-------------------------------- open COORD section */
        __normal_call void_type open_coord (
            std::int32_t _idim,
            std::int32_t _nrow
            )
        {   
            if (_idim == +1)
            jigsaw_alloc_reals (
                &this->_jmsh->_xgrid, _nrow) ;
            else
            if (_idim == +2)
            jigsaw_alloc_reals (
                &this->_jmsh->_ygrid, _nrow) ;
            else
            if (_idim == +3)
            jigsaw_alloc_reals (
                &this->_jmsh->_zgrid, _nrow) ;       
        }
    /*-------------------------------- read COORD section */
        __normal_call void_type push_coord (
            std::int32_t _idim,
            std::int32_t _ipos,
            double       _ppos
            )
        {
            if (_idim == +1)
            {
            
            if (_ipos < this->_jmsh->_xgrid._size)
            {
                this->_jmsh->_xgrid.
                    _data[_ipos] = _ppos ;
            }
            else
            {
                this->_errv = __invalid_argument ;
            }
            
            }
            else
            if (_idim == +2)
            {
            
            if (_ipos < this->_jmsh->_ygrid._size)
            {
                this->_jmsh->_ygrid.
                    _data[_ipos] = _ppos ;
            }
            else
            {
                this->_errv = __invalid_argument ;
            }
            
            }
            else
            if (_idim == +3)
            {
            
            if (_ipos < this->_jmsh->_zgrid._size)
            {
                this->_jmsh->_zgrid.
                    _data[_ipos] = _ppos ;
            }
            else
            {
                this->_errv = __invalid_argument ;
            }
            
            }   
        }
    /*-------------------------------- open EDGE2 section */
        __normal_call void_type open_edge2 (
            std::int32_t  _nrow
            ) 
        {
            jigsaw_alloc_edge2 (
               &this->_jmsh->_edge2, _nrow) ;
        }
    /*-------------------------------- push EDGE2 section */
        __normal_call void_type push_edge2 (
            std::int32_t  _ipos,
            std::int32_t* _node,
            std::int32_t  _itag
            ) 
        { 
            if (_ipos < this->_jmsh->_edge2._size)
            {
            this->_jmsh->_edge2._data[_ipos].
                _node[0] = _node[0] ;
            this->_jmsh->_edge2._data[_ipos].
                _node[1] = _node[1] ;
                
            this->_jmsh->_edge2.
                _data[_ipos]._itag =  _itag ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            }
        }
    /*-------------------------------- open TRIA3 section */
        __normal_call void_type open_tria3 (
            std::int32_t  _nrow
            ) 
        {
            jigsaw_alloc_tria3 (
               &this->_jmsh->_tria3, _nrow) ;
        }
    /*-------------------------------- push TRIA3 section */
        __normal_call void_type push_tria3 (
            std::int32_t  _ipos,
            std::int32_t* _node,
            std::int32_t  _itag
            ) 
        { 
            if (_ipos < this->_jmsh->_tria3._size)
            {
            this->_jmsh->_tria3._data[_ipos].
                _node[0] = _node[0] ;
            this->_jmsh->_tria3._data[_ipos].
                _node[1] = _node[1] ;
            this->_jmsh->_tria3._data[_ipos].
                _node[2] = _node[2] ;
                
            this->_jmsh->_tria3.
                _data[_ipos]._itag =  _itag ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            }
        }
    /*-------------------------------- open TRIA4 section */
        __normal_call void_type open_tria4 (
            std::int32_t  _nrow
            ) 
        {
            jigsaw_alloc_tria4 (
               &this->_jmsh->_tria4, _nrow) ;
        }
    /*-------------------------------- push TRIA4 section */
        __normal_call void_type push_tria4 (
            std::int32_t  _ipos,
            std::int32_t* _node,
            std::int32_t  _itag
            ) 
        { 
            if (_ipos < this->_jmsh->_tria4._size)
            {
            this->_jmsh->_tria4._data[_ipos].
                _node[0] = _node[0] ;
            this->_jmsh->_tria4._data[_ipos].
                _node[1] = _node[1] ;
            this->_jmsh->_tria4._data[_ipos].
                _node[2] = _node[2] ;
            this->_jmsh->_tria4._data[_ipos].
                _node[3] = _node[3] ;
                
            this->_jmsh->_tria4.
                _data[_ipos]._itag =  _itag ;
            }
            else
            {
            this->_errv =__invalid_argument ;
            }
        }
    /*-------------------------------- open VALUE section */
        __normal_call void_type open_value (
            std::int32_t  _nrow,
            std::int32_t  _nval
            ) 
        { 
            this->_nval = _nval;

            jigsaw_alloc_reals (
                &this->
            _jmsh->_value , _nrow * _nval ) ;
        }
    /*-------------------------------- push VALUE section */
        __normal_call void_type push_value (
            std::int32_t  _ipos,
            double      * _vdat
            ) 
        { 
            if (_ipos < this->_jmsh->_value._size)
            {
            
            for (auto _ival = +0; 
                _ival< this->_nval; ++_ival)
            {
                this->_jmsh->_value.
            _data[_ipos*(_ival+1)] = _vdat[_ival];
            }
            
            }
            else
            {
                this->_errv =  __invalid_argument;
            }
        }
        
        } ;
    
    /*---------------------------------- parse MSH_T file */
        iptr_type _errv  = __no_error ;
    
        try
        {
            jmsh_reader   _read ;
            std::ifstream _file ;
            
            std::string _fstr(_fchr) ;
            std::string _path ;
            std::string _name ;
            std::string _fext ;
            file_part(_fstr, 
                _path, _name, _fext) ;
             
            _file.open(
                _fstr, std::ifstream:: in) ;

            if (_file.is_open() )
            {
                _read.read_file (
                _file, msht_reader(_jmsh)) ;
            }
            else
            {           
                _errv = __file_not_located ;
            }
            
            _file.close ();

            if(!_read._errs.empty())
            {
                _errv = __invalid_argument ;
            }
        }
        catch (...)
        {
            _errv = __unknown_error ;
        }

    /*---------------------------------- return read flag */
        return  _errv ;
    }

#   endif//__lib_jigsaw

#   endif



