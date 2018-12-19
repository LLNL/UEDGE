#ifndef IDFUSR_H
#define IDFUSR_H

/*

============================================================
       Import of Data from File (IDF) package
     for time-saving programming of data input
              to scientific codes.

                 IDF Version 1.1   
              Date: February 1, 2000


 Copyright 1999 by A. Yu. Pigarov and I.V. Saltanova


                   Developers:
         A. Yu. Pigarov and I.V. Saltanova


         E-mails: apigarov@psfc.mit.edu
                  apigarov@pppl.gov
                  iras@rex.pfc.mit.edu

         IDF Documentation is available at:
  http://www2.psfc.mit.edu/library/preprints.html

          A.Yu. Pigarov, I.V. Saltanova, 
      "Import Data from File (IDF) utilities
  for programming data input to scientific codes",
         Plasma Science and Fusion Center, 
       Massachusetts Institute of Technology,
    Cambridge, MA, PSFC/JA-00-01, January 2000.



    The IDF package is distributed in the hope 
that it will be useful to many computational scientists.

    The IDF package is FREE SOFTWARE. 
You can use, copy, and modify this software for any purpose
and without fee provided that the above copyright
notice appear in all copies.

    The IDF package is distributed "AS IS", 
i.e without any warranty including all implied warranties
of merchantability and fitness.

============================================================

*/

/*
=========================================================
In creating the IDF library,  
macros IDF_CPP_STYLE , IDF_CPP , and IDF_FORTRAN_TYPE
should be defined in "idflib.h" header file.

If you are including "idfusr.h" separately
in your C or C++ code, you must first define
IDF_CPP_STYLE , IDF_CPP , and IDF_FORTRAN_TYPE macros 
in your source code in a string preceeding
#include "idfusr.h".
=========================================================
*/


/*---------------------------------------------------- 
   macro IDF_CPP_STYLE controls the declaration 
of C function arguments as well as function prototype

If IDF_CPP_STYLE is defined as 1
then function prototypes will be declared in 
1989 Standard C style:

type fun(type Arg1, type Arg2)
{ function body }

If IDF_CPP_STYLE is defined as 0 then:
then function prototypes will be declared in 
Pre-Standard C style:

type fun(Arg1,Arg2)
     type Arg1; type Arg2;
{ function body } 
--------------------------------------------------------*/
#ifndef IDF_CPP_
#ifndef IDF_CPP_STYLE
#define IDF_CPP_ 1
#else
#define IDF_CPP_ IDF_CPP_STYLE
#endif
#endif


/*
------------------------------------------
 Specification of naming conventions for
   C functions to be called from FORTRAN
------------------------------------------
 Define the macros IDF_FORTRAN_TYPE:
    for   ANSI_F77  as   1
    for       CRAY  as   1
    for        AIX  as   2
    for       HPUX  as   2
    for        SVS  as   2
    for      other  as   3
*/
#ifndef IDF_FORTRAN_TYPE
#define IDF_FORTRAN_TYPE 3
#endif

#if IDF_FORTRAN_TYPE == 1
#define IDF_FORTRAN(fid_, fid, FID) FID 
#endif
#if IDF_FORTRAN_TYPE == 2
#define IDF_FORTRAN(fid_, fid, FID) fid 
#endif
#if IDF_FORTRAN_TYPE == 3
#define IDF_FORTRAN(fid_, fid, FID) fid_ 
#endif

#ifndef IDF_FORTRAN
#define IDF_FORTRAN(fid_, fid, FID) fid_ 
#endif




/*---------------------------------
  define IDF_CPP as 1 to protect
  C function prototypes under C++ 
-----------------------------------*/

#ifndef IDF_CPP
#define IDF_CPP 0
#endif


#if IDF_CPP == 1
extern "C" {
#endif


#if IDF_CPP_ == 1

/*idf_ini.c*/
extern  int   idf_ini  (char*);
extern  void  idf_end  ();

/*idf_init.c*/
extern  int   idf_init   ();
extern  int   idf_open   (char*);
extern  void  idf_close  ();
extern  void  idf_finish ();

/*idf_do.c*/
extern int   idf_c     (char*, char*);
extern int   idf_uc    (char*, unsigned char*);
extern int   idf_u     (char*, unsigned int*);
extern int   idf_s     (char*, char*, int);
extern int   idf_t     (char*, char*, int);
extern int   idf_i     (char*, int*);
extern int   idf_l     (char*, long*);
extern int   idf_f     (char*, float*);
extern int   idf_d     (char*, double*);
extern int   idf_z     (char*, float*);
extern int   idf_w     (char*, double*);
extern int   idf_iarr  (char*, int*, int);
extern int   idf_larr  (char*, long*, int);
extern int   idf_farr  (char*, float*, int);
extern int   idf_darr  (char*, double*, int);
extern int   idf_zarr  (char*, float*, int);
extern int   idf_warr  (char*, double*, int);

/*idf_dof.c*/
extern int   IDF_FORTRAN(idfinit ,idfinit_ ,IDFINIT ) (int*);
extern void  IDF_FORTRAN(idfinish,idfinish_,IDFINISH) ();
extern int   IDF_FORTRAN(idfini  ,idfini_  ,IDFINI  ) (char*,int*);
extern void  IDF_FORTRAN(idfend  ,idfend_  ,IDFEND  ) ();
extern int   IDF_FORTRAN(idfopen ,idfopen_ ,IDFOPEN ) (char*,int*);
extern void  IDF_FORTRAN(idfclose,idfclose_,IDFCLOSE) ();
extern void  IDF_FORTRAN(idferprn,idferprn_,IDFERPRN) ();
extern int   IDF_FORTRAN(idfc    ,idfc_    ,IDFC    )
                        (char*, int*, void*);
extern int   IDF_FORTRAN(idfs    ,idfs_    ,IDFS    )
                        (char*, int*, void*,int*);
extern int   IDF_FORTRAN(idft    ,idft_    ,IDFT    )
                        (char*, int*, void*,int*);
extern int   IDF_FORTRAN(idfi    ,idfi_    ,IDFI    )
                        (char*, int*, int*);
extern int   IDF_FORTRAN(idff    ,idff_    ,IDFF    )
                        (char*, int*, float*);
extern int   IDF_FORTRAN(idfd    ,idfd_    ,IDFD    )
                        (char*, int*, double*);
extern int   IDF_FORTRAN(idfz    ,idfz_    ,IDFZ    )
                        (char*, int*, float*);
extern int   IDF_FORTRAN(idfw    ,idfw_    ,IDFW    )
                        (char*, int*, double*);
extern int   IDF_FORTRAN(idfiarr ,idfiarr_ ,IDFIARR )
                        (char*, int*, int*,    int*);
extern int   IDF_FORTRAN(idffarr ,idffarr_ ,IDFFARR )
                        (char*, int*, float*,  int*);
extern int   IDF_FORTRAN(idfdarr ,idfdarr_ ,IDFDARR ) 
                        (char*, int*, double*, int*);
extern int   IDF_FORTRAN(idfzarr ,idfzarr_ ,IDFZARR )
                        (char*, int*, float*,  int*);
extern int   IDF_FORTRAN(idfwarr ,idfwarr_ ,IDFWARR ) 
                        (char*, int*, double*, int*);
extern int   IDF_FORTRAN(idfget  ,idfget_  ,IDFGET  ) 
                        (char*, int*, char*, int*, int*, void*);
extern int   IDF_FORTRAN(idfarray,idfarray_,IDFARRAY) 
                        (char*, int*, char*, int*,       void*);

/*idf_getone.c*/
extern int   idf_get_one (char*, char, int, void*);
extern int   idf_get_one_(char*,int, char, int, void*);

/*idf_getarr.c*/
extern int   idf_get_1Darray (char*, char, int, void*);
extern int   idf_get_1Darray_(char*,int, char, int, void*);
extern int   idf_getarray    (char*,int, char, int,int*,int, void*);

/*idf_array.c*/
extern int   idf_get_array (char*, char*, void*);
extern int   idf_get_array_(char*,int, char*,int, void*);

/*idf_get.c*/
extern int   idf_get (char*, char*, int, void*);
extern int   idf_get_(char*,int, char*,int, int, void*);


#else /* IDF_CPP_ */


/*idf_ini.c*/
extern int   idf_ini  ();
extern void  idf_end  ();

/*idf_init.c*/
extern int   idf_init   ();
extern int   idf_open   ();
extern void  idf_close  ();
extern void  idf_finish ();

/*idf_do.c*/
extern int   idf_c     ();
extern int   idf_uc    ();
extern int   idf_u     ();
extern int   idf_s     ();
extern int   idf_t     ();
extern int   idf_i     ();
extern int   idf_l     ();
extern int   idf_f     ();
extern int   idf_d     ();
extern int   idf_z     ();
extern int   idf_w     ();
extern int   idf_iarr  ();
extern int   idf_larr  ();
extern int   idf_farr  ();
extern int   idf_darr  ();
extern int   idf_zarr  ();
extern int   idf_warr  ();

/*idf_dof.c*/
extern int   IDF_FORTRAN(idfinit ,idfinit_ ,IDFINIT ) ();
extern void  IDF_FORTRAN(idfinish,idfinish_,IDFINISH) ();
extern int   IDF_FORTRAN(idfini  ,idfini_  ,IDFINI  ) ();
extern void  IDF_FORTRAN(idfend  ,idfend_  ,IDFEND  ) ();
extern int   IDF_FORTRAN(idfopen ,idfopen_ ,IDFOPEN ) ();
extern void  IDF_FORTRAN(idfclose,idfclose_,IDFCLOSE) ();
extern void  IDF_FORTRAN(idferprn,idferprn_,IDFERPRN) ();
extern int   IDF_FORTRAN(idfc    ,idfc_    ,IDFC    ) ();
extern int   IDF_FORTRAN(idfs    ,idfs_    ,IDFS    ) ();
extern int   IDF_FORTRAN(idft    ,idft_    ,IDFT    ) ();
extern int   IDF_FORTRAN(idfi    ,idfi_    ,IDFI    ) ();
extern int   IDF_FORTRAN(idff    ,idff_    ,IDFF    ) ();
extern int   IDF_FORTRAN(idfd    ,idfd_    ,IDFD    ) ();
extern int   IDF_FORTRAN(idfz    ,idfz_    ,IDFZ    ) ();
extern int   IDF_FORTRAN(idfw    ,idfw_    ,IDFW    ) ();
extern int   IDF_FORTRAN(idfiarr ,idfiarr_ ,IDFIARR ) ();
extern int   IDF_FORTRAN(idffarr ,idffarr_ ,IDFFARR ) ();
extern int   IDF_FORTRAN(idfdarr ,idfdarr_ ,IDFDARR ) ();
extern int   IDF_FORTRAN(idfzarr ,idfzarr_ ,IDFZARR ) ();
extern int   IDF_FORTRAN(idfwarr ,idfwarr_ ,IDFWARR ) ();
extern int   IDF_FORTRAN(idfget  ,idfget_  ,IDFGET  ) ();
extern int   IDF_FORTRAN(idfarray,idfarray_,IDFARRAY) ();

/*idf_getone.c*/
extern int   idf_get_one ();
extern int   idf_get_one_();

/*idf_getarr.c*/
extern int   idf_get_1Darray ();
extern int   idf_get_1Darray_();
extern int   idf_getarray    ();

/*idf_array.c*/
extern int   idf_get_array ();
extern int   idf_get_array_();

/*idf_get.c*/
extern int   idf_get ();
extern int   idf_get_();


#endif /* IDF_CPP_ */


/* protect C function prototypes under C++ */
#if IDF_CPP == 1
           };
#endif


#endif /*IDFUSR_H*/
