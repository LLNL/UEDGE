/*

=======================================================
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
   
=======================================================

*/


#ifndef IDF_H
#define IDF_H


/*--------------*/
#include "idf.i"
/*--------------*/


#ifndef IDF_CPP_

#ifndef IDF_CPP_STYLE
#define IDF_CPP_ 1
#else
#define IDF_CPP_ IDF_CPP_STYLE
#endif

#endif


#if IDF_CPP_ == 1

/*idf_std.c*/
extern int idf_string_over  (unsigned int);
extern int idf_number_over  (unsigned int);
extern int idf_name_seps    (unsigned int);
extern int idf_data_seps    (unsigned int);
extern int idf_data_symbol  (unsigned int);
extern int idf_name_symbol  (unsigned int);
extern int idf_hex_symbol   (unsigned int);
extern int idf_octal_symbol (unsigned int);
extern int idf_exp_symbol   (unsigned int);
extern int idf_form_symbol  (unsigned int);
extern int idf_file_symbol  (unsigned int);

/*idf_keyw.c*/
extern int idf_keyw       (char*);
extern int idf_keyw_void  (int);
extern int idf_keyw_local (int);
extern int idf_keyw_digit (int);
extern int idf_keyw_cnv   (int,int);
extern int idf_keyw_frmt  (int,int);

/*idf_keyF.c*/
extern char* idf_keywF_name (int);
extern int   idf_keywF_lname(int);
extern int   idf_keywF      (char*,int);
extern int   idf_keywF_narg (int);
extern int   idf_keywF_targ (int);
extern int   idf_keywF_tfun (int);

/*idf_float.c*/
#if IDF_FLOAT_DEFAULT == 0
extern int  idf_float_ini();
extern void idf_float_end();
#endif
extern int  idf_float_compose  (double,int,double*);
extern void idf_float_decompose(double,int*,double*,int*);

/*idf_mathd.c*/
extern int idf_mathd_add  (double,double,double*);
extern int idf_mathd_sub  (double,double,double*);
extern int idf_mathd_mul  (double,double,double*);
extern int idf_mathd_div  (double,double,double*);
extern int idf_mathd_pow  (double,double,double*);
extern int idf_mathd      (char, double,double,double*);

/*idf_mathdf.c*/
extern int idf_mathd_abs    (double,double*);
extern int idf_mathd_int    (double,double*);
extern int idf_mathd_sin    (double,double*);
extern int idf_mathd_cos    (double,double*);
extern int idf_mathd_tan    (double,double*);
extern int idf_mathd_sinh   (double,double*);
extern int idf_mathd_cosh   (double,double*);
extern int idf_mathd_tanh   (double,double*);
extern int idf_mathd_log    (double,double*);
extern int idf_mathd_log10  (double,double*);
extern int idf_mathd_sqrt   (double,double*);
extern int idf_mathd_exp    (double,double*);
extern int idf_mathd_asin   (double,double*);
extern int idf_mathd_acos   (double,double*);
extern int idf_mathd_atan   (double,double*);
extern int idf_mathd_asinh  (double,double*);
extern int idf_mathd_acosh  (double,double*);
extern int idf_mathd_atanh  (double,double*);
extern int idf_mathd_hypot  (double,double,double*);
extern int idf_mathd_min    (double,double,double*);
extern int idf_mathd_max    (double,double,double*);
extern int idf_mathd_sq1mz2 (double,double*);
extern int idf_mathd_sq1pz2 (double,double*);
extern int idf_mathd_pow10  (double,double*);

/*idf_mathc.c*/
extern int idf_mathc_add  (double,double,double,double,double*,double*);
extern int idf_mathc_sub  (double,double,double,double,double*,double*);
extern int idf_mathc_mul  (double,double,double,double,double*,double*);
extern int idf_mathc_div  (double,double,double,double,double*,double*);
extern int idf_mathc_pow  (double,double,double,double,double*,double*);

/*idf_mathcf.c*/
extern int idf_mathc_cmplx  (double,double,double*,double*);
extern int idf_mathc_real   (double,double,double*,double*);
extern int idf_mathc_imag   (double,double,double*,double*);
extern int idf_mathc_polar  (double,double,double*,double*);
extern int idf_mathc_abs    (double,double,double*);
extern int idf_mathc_arg    (double,double,double*);
extern int idf_mathc_sin    (double,double,double*,double*);
extern int idf_mathc_cos    (double,double,double*,double*);
extern int idf_mathc_tan    (double,double,double*,double*);
extern int idf_mathc_sinh   (double,double,double*,double*);
extern int idf_mathc_cosh   (double,double,double*,double*);
extern int idf_mathc_tanh   (double,double,double*,double*);
extern int idf_mathc_log    (double,double,double*,double*);
extern int idf_mathc_log10  (double,double,double*,double*);
extern int idf_mathc_sqrt   (double,double,double*,double*);
extern int idf_mathc_exp    (double,double,double*,double*);
extern int idf_mathc_sq1mz2 (double,double,double*,double*);
extern int idf_mathc_sq1pz2 (double,double,double*,double*);
extern int idf_mathc_asin   (double,double,double*,double*);
extern int idf_mathc_acos   (double,double,double*,double*);
extern int idf_mathc_atan   (double,double,double*,double*);
extern int idf_mathc_asinh  (double,double,double*,double*);
extern int idf_mathc_acosh  (double,double,double*,double*);
extern int idf_mathc_atanh  (double,double,double*,double*);
extern int idf_mathc_min    (double,double,double,double,double*,double*);
extern int idf_mathc_max    (double,double,double,double,double*,double*);
extern int idf_mathc_pow10  (double,double,double*,double*);

/*idf_math.c*/
extern void idf_matherr_put(int,int);
extern void idf_matherr_nul();
extern void idf_matherr_prn();
extern int  idf_mathd      (char,double,double,double*);
extern int  idf_mathc      (char,double,double,double,double,double*,double*);
extern int  idf_math       (char, int,double,double, int,double,double,
                            int*,double*,double*);
extern int  idf_mathdf     (int,double,double,double*,double*);
extern int  idf_mathcf     (int,double,double,double,double,double*,double*);
extern int  idf_math_func  (int,int, int,double,double, int,double,double,
                            int*,double*,double*);

/*idf_dtos.c*/
extern char* idf_dtos(double,char*);

/*idf_strtol.c*/
extern int idf_str_to_l(char*,int, long*);

/*idf_strtod.c*/
extern int idf_str_to_d(char*,int, double*);

/*idf_strtoc.c*/
extern int idf_str_to_c(char*,int, char*);

/*idf_strtos.c*/
extern int idf_str_to_s(char*,int, char*,int);

/*idf_strtoz.c*/
extern int idf_str_to_z(char*,int, double*);

/*idf_strtov.c*/
extern int idf_str_to_v(char*,int, double*, int*);

/*idf_jus.c*/
extern void idf_jus_fort();
extern void idf_jus_c   ();
extern int  idf_jus_get ();

/*idf_cnv.c*/
extern int idf_cnv_c_to_val   (int, char,     void*,int);
extern int idf_cnv_s_to_val   (int, char*,int,void*,int);
extern int idf_cnv_s_to_val_c (int, char*,int,void*,int);
extern int idf_cnv_s_to_val_uc(int, char*,int,void*,int);
extern int idf_cnv_i_to_val   (int, int,      void*);
extern int idf_cnv_l_to_val   (int, long,     void*);
extern int idf_cnv_f_to_val   (int, float,    void*);
extern int idf_cnv_d_to_val   (int, double,   void*);
extern int idf_cnv_z_to_val   (int, float*,   void*);
extern int idf_cnv_w_to_val   (int, double*,  void*);

/*idf_sname.c*/
extern void  idf_Sname_put (char*,int);
extern char* idf_Sname_get ();

/*idf_stxt.c*/
extern int idf_stxt_skip_cmnt (int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_form (int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_char (int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_str  (int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_cmplx(int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_crpt (int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_trig (int*,int*,long*, unsigned int*);
extern int idf_stxt_skip_txt  (int*,int*,long*, unsigned int*);

/*idf_sdata.c*/
extern int idf_skip_data       (int*,int*,long*, unsigned int*);
extern int idf_skip_data_period(int*,int*,long*, unsigned int*);

/*idf_err.c*/
extern void idf_err_clear    ();
extern void idf_err_put      (int, int);
extern void idf_err_prn      ();

/*idf_geterr.c*/
extern void idf_get_err_prn  (char*,int, int,int);
extern void idf_list_err_prn ();

/*idf_txterr.c*/
extern void idf_txt_err_nul  ();
extern void idf_txt_err_put  (int, int,int,long, unsigned int);
extern void idf_txt_err_prn  ();

/*idf_dim.c*/
extern int idf_dim        (IDF_NAME_*);
extern int idf_dim_nmax   (int,int*);
extern int idf_dim_offset (int,int*, int, int,int*);

/*idf_align.c*/
extern int   idf_align_bnd  (int,int);
extern void  idf_align      (int,int,int, int,int, int*,int*);

/*idf_frmtype.c*/
extern int idf_frmt_type     (unsigned int);
extern int idf_frmt_ctype    (char);
extern int idf_frmt_table_prp();
extern int idf_frmt_size     (int);
extern int idf_frmt_psize    (int);
extern int idf_frmt_word     (int);
extern int idf_frmt_cns      (int);
extern int idf_frmt_amode    (int);

/*idf_frmt.c*/
extern IDF_FRMT* idf_frmt_address();
extern int       idf_frmt_ini    ();
extern void      idf_frmt_end    ();
extern void      idf_frmt_clean  ();
extern int       idf_frmt_prp    (char*,int);
extern void      idf_frmt_prn    ();
extern int       idf_frmt_noalloc();
extern void      idf_frmt_err_prn();

/*idf_frmtget.c*/
extern IDF_FRMT_TREAT* idf_frmtreat_address();
extern void            idf_frmtreat_nul    ();
extern void            idf_frmtreat_prn    ();
extern int             idf_frmt_begin      (int,int,char*);
extern int             idf_frmt_set_pntr   (void**,int);
extern void            idf_frmt_put_pntr   (void**,void*);
extern char*           idf_frmt_get_pntr   (void**);
extern int             idf_frmtreat_array  (void*);
extern int             idf_frmt_get        ();
extern void            idf_frmtreat_err_prn();

/*idf_name.c*/
extern IDF_NAME* idf_name_address();
extern int       idf_name_ini   ();
extern void      idf_name_end   ();
extern void      idf_name_prn   ();
extern void      idf_name_errprn();
extern void      idf_name_begin (int*,int*,long*);
extern int       idf_name_buf   (int*,int*,long*);
extern int       idf_name_raw   ();
extern int       idf_name       ();
extern int       idf_name_dpos  (int);
extern int       idf_name_get   (int*,int*,long*,int*);

/*idf_extname.c*/
extern int       idf_extname_check (char*,int);

/* idf_fil.c*/
extern IDF_FILE* idf_file_address ();
extern int       idf_file_ini     ();
extern char*     idf_file_Cname   ();
extern void      idf_file_end     ();
extern void      idf_file_prn     ();
extern int       idf_file_open    (char* );
extern void      idf_file_close   ();
extern int       idf_file_name    (char*);
extern void      idf_file_Cname_free();
extern int       idf_file_name_   (char*,int);
extern int       idf_file_opened  ();
extern int       idf_file_ferr    ();
extern int       idf_file_getc    (unsigned int* );
extern int       idf_file_gets    (char* ,int, int*);
extern int       idf_file_ftell   (long*);
extern int       idf_file_fseek   (long);
extern int       idf_file_rewind  ();
extern void      idf_file_errprn  ();
extern void      idf_file_clearerr();

/* idf_list.c*/
extern IDF_NAME_LIST* idf_list_address ();
extern void           idf_list_nul     ();
extern int            idf_list_ini     ();
extern void           idf_list_end     ();
extern void           idf_list_clean   ();
extern int            idf_list_num     ();
extern int            idf_list_numL    ();
extern IDF_NAME_*     idf_list_get     (char*,int,int);
extern IDF_NAME_L*    idf_list_getL    (char*,int,int);
extern int            idf_list_getLV   (char*,int, int*,double*);
extern int            idf_list_put     (IDF_NAME*);
extern int            idf_list_putL    (IDF_NAME*);
extern void           idf_list_prn_    (IDF_NAME_*);
extern void           idf_list_prnL_   (IDF_NAME_L*);
extern void           idf_list_prn     (int);
extern void           idf_list_prnL    (int);
extern int            idf_list_do      ();
extern int            idf_list_check   ();
extern void           idf_list_catalog ();
extern void           idf_list_catalogL();
extern void           idf_list_clearerr();

/*idf_mlist.c*/
extern int idf_mlist  (char*,int, double*);

/*idf_plist.c*/
extern int idf_plist  (char*,int, double*);

/*idf_clist.c*/
extern int idf_clist  (char*,int, double*);

/*idf_dbuf.c*/
extern IDF_DBUF* idf_dbuf_address();
extern int       idf_dbuf_ini    ();
extern void      idf_dbuf_end    ();
extern void      idf_dbuf_nul    ();
extern void      idf_dbuf_prn    ();
extern void      idf_dbuf_begin  (int*,int*,long*, int);
extern int       idf_dbuf_number (int*,int*,long*, unsigned int*);
extern int       idf_dbuf_char   (int*,int*,long*, unsigned int*);
extern int       idf_dbuf_str    (int*,int*,long*, unsigned int*);
extern int       idf_dbuf_txt    (int*,int*,long*, unsigned int*);
extern int       idf_dbuf_cmplx  (int*,int*,long*, unsigned int*);
extern int       idf_dbuf_add    (unsigned int);
extern int       idf_dbuf_to_d_  (double*);
extern int       idf_dbuf_to_f_  (float*);
extern int       idf_dbuf_to_l_  (long*);
extern int       idf_dbuf_to_i_  (int*);
extern int       idf_dbuf_to_d   (double*);
extern int       idf_dbuf_to_f   (float*);
extern int       idf_dbuf_to_l   (long*);
extern int       idf_dbuf_to_i   (int*);
extern int       idf_dbuf_to_s   (char*,int*);
extern int       idf_dbuf_to_t   (char*,int*);
extern int       idf_dbuf_to_c   (char*);
extern int       idf_dbuf_to_z   (float*);
extern int       idf_dbuf_to_w   (double*);
extern int       idf_dbuf_to_v   (double*,int*);

/*idf_fbuf.c*/
extern IDF_FBUF* idf_fbuf_address();
extern int       idf_fbuf_ini    ();
extern void      idf_fbuf_end    ();
extern void      idf_fbuf_nul    ();
extern void      idf_fbuf_prn    ();
extern void      idf_fbuf_begin  (int*,int*,long*, int);
extern int       idf_fbuf_number (int*,int*,long*, unsigned int*);
extern int       idf_fbuf_name   (int*,int*,long*, unsigned int*);
extern int       idf_fbuf_conv   ();

/*idf_form.c*/
extern IDF_FORM* idf_form_address();
extern int       idf_form_ini    ();
extern void      idf_form_end    ();
extern void      idf_form_nul    ();
extern void      idf_form_prn    ();
extern void      idf_form_view   ();
extern int       idf_form_do     (int*,int*,long*, unsigned int*);
extern int       idf_form_up     (int,void*);

/*idf_data.c*/
extern IDF_DATA* idf_data_address();
extern int       idf_data_ini    ();
extern void      idf_data_end    ();
extern void      idf_data_prn0   ();
extern void      idf_data_prn1   ();
extern int       idf_data_begin  (IDF_NAME_*);
extern int       idf_data_beginL (IDF_NAME_L*);
extern int       idf_data_unit   ();
extern int       idf_data_value  ();
extern int       idf_data_up     (int,void*,int);

/*idf_val.c*/
extern IDF_VAL*  idf_val_address();
extern int       idf_val_ini    ();
extern void      idf_val_end    ();
extern void      idf_val_prn    ();
extern int       idf_val_put    (int);
extern int       idf_val_buf    ();
extern int       idf_val_form   ();
extern int       idf_val_get    (int, void*,int);
extern int       idf_val_lcns   (IDF_NAME_L*);

/*idf_fflist.c*/
extern int idf_list_ffl ();

/*idf_catal.c*/
extern int idf_catalog ();

/*idf_order.c*/
extern void idf_order_put (int);
extern int  idf_order_cmp (int);



#else



/*idf_std.c*/
extern int idf_string_over  ();
extern int idf_number_over  ();
extern int idf_name_seps    ();
extern int idf_data_seps    ();
extern int idf_data_symbol  ();
extern int idf_name_symbol  ();
extern int idf_hex_symbol   ();
extern int idf_octal_symbol ();
extern int idf_exp_symbol   ();
extern int idf_form_symbol  ();
extern int idf_file_symbol  ();

/*idf_keyw.c*/
extern int idf_keyw       ();
extern int idf_keyw_void  ();
extern int idf_keyw_local ();
extern int idf_keyw_digit ();
extern int idf_keyw_cnv   ();
extern int idf_keyw_frmt  ();

/*idf_keyF.c*/
extern char* idf_keywF_name ();
extern int   idf_keywF_lname();
extern int   idf_keywF      ();
extern int   idf_keywF_narg ();
extern int   idf_keywF_targ ();
extern int   idf_keywF_tfun ();

/*idf_float.c*/
#if IDF_FLOAT_DEFAULT == 0
extern int  idf_float_ini();
extern void idf_float_end();
#endif
extern int  idf_float_compose  ();
extern void idf_float_decompose();

/*idf_mathd.c*/
extern int idf_mathd_add  ();
extern int idf_mathd_sub  ();
extern int idf_mathd_mul  ();
extern int idf_mathd_div  ();
extern int idf_mathd_pow  ();

/*idf_mathdf.c*/
extern int idf_mathd_abs    ();
extern int idf_mathd_int    ();
extern int idf_mathd_sin    ();
extern int idf_mathd_cos    ();
extern int idf_mathd_tan    ();
extern int idf_mathd_sinh   ();
extern int idf_mathd_cosh   ();
extern int idf_mathd_tanh   ();
extern int idf_mathd_log    ();
extern int idf_mathd_log10  ();
extern int idf_mathd_sqrt   ();
extern int idf_mathd_exp    ();
extern int idf_mathd_asin   ();
extern int idf_mathd_acos   ();
extern int idf_mathd_atan   ();
extern int idf_mathd_asinh  ();
extern int idf_mathd_acosh  ();
extern int idf_mathd_atanh  ();
extern int idf_mathd_hypot  ();
extern int idf_mathd_min    ();
extern int idf_mathd_max    ();
extern int idf_mathd_sq1mz2 ();
extern int idf_mathd_sq1pz2 ();
extern int idf_mathd_pow10  ();

/*idf_mathc.c*/
extern int idf_mathc_add  ();
extern int idf_mathc_sub  ();
extern int idf_mathc_mul  ();
extern int idf_mathc_div  ();
extern int idf_mathc_pow  ();

/*idf_mathcf.c*/
extern int idf_mathc_cmplx  ();
extern int idf_mathc_real   ();
extern int idf_mathc_imag   ();
extern int idf_mathc_polar  ();
extern int idf_mathc_abs    ();
extern int idf_mathc_arg    ();
extern int idf_mathc_sin    ();
extern int idf_mathc_cos    ();
extern int idf_mathc_tan    ();
extern int idf_mathc_sinh   ();
extern int idf_mathc_cosh   ();
extern int idf_mathc_tanh   ();
extern int idf_mathc_log    ();
extern int idf_mathc_log10  ();
extern int idf_mathc_sqrt   ();
extern int idf_mathc_exp    ();
extern int idf_mathc_sq1mz2 ();
extern int idf_mathc_sq1pz2 ();
extern int idf_mathc_asin   ();
extern int idf_mathc_acos   ();
extern int idf_mathc_atan   ();
extern int idf_mathc_asinh  ();
extern int idf_mathc_acosh  ();
extern int idf_mathc_atanh  ();
extern int idf_mathc_min    ();
extern int idf_mathc_max    ();
extern int idf_mathc_pow10  ();

/*idf_math.c*/
extern void idf_matherr_put();
extern void idf_matherr_nul();
extern void idf_matherr_prn();
extern int  idf_mathd      ();
extern int  idf_mathc      ();
extern int  idf_math       ();
extern int  idf_mathdf     ();
extern int  idf_mathcf     ();
extern int  idf_math_func  ();

/*idf_dtos.c*/
extern char* idf_dtos();

/*idf_strtol.c*/
extern int idf_str_to_l();

/*idf_strtod.c*/
extern int idf_str_to_d();

/*idf_strtoc.c*/
extern int idf_str_to_c();

/*idf_strtos.c*/
extern int idf_str_to_s();

/*idf_strtoz.c*/
extern int idf_str_to_z();

/*idf_strtov.c*/
extern int idf_str_to_v();

/*idf_jus.c*/
extern void idf_jus_fort();
extern void idf_jus_c   ();
extern int  idf_jus_get ();

/*idf_cnv.c*/
extern int idf_cnv_c_to_val   ();
extern int idf_cnv_s_to_val   ();
extern int idf_cnv_s_to_val_c ();
extern int idf_cnv_s_to_val_uc();
extern int idf_cnv_i_to_val   ();
extern int idf_cnv_l_to_val   ();
extern int idf_cnv_f_to_val   ();
extern int idf_cnv_d_to_val   ();
extern int idf_cnv_z_to_val   ();
extern int idf_cnv_w_to_val   ();

/*idf_sname.c*/
extern void  idf_Sname_put ();
extern char* idf_Sname_get ();

/*idf_stxt.c*/
extern int idf_stxt_skip_cmnt ();
extern int idf_stxt_skip_form ();
extern int idf_stxt_skip_char ();
extern int idf_stxt_skip_str  ();
extern int idf_stxt_skip_crpt ();
extern int idf_stxt_skip_cmplx();
extern int idf_stxt_skip_trig ();
extern int idf_stxt_skip_txt  ();

/*idf_dim.c*/
extern int idf_dim        ();
extern int idf_dim_nmax   ();
extern int idf_dim_offset ();

/*idf_align.c*/
extern int   idf_align_bnd  ();
extern void  idf_align      ();

/*idf_frmtype.c*/
extern int idf_frmt_type     ();
extern int idf_frmt_ctype    ();
extern int idf_frmt_table_prp();
extern int idf_frmt_size     ();
extern int idf_frmt_word     ();
extern int idf_frmt_psize    ();
extern int idf_frmt_cns      ();
extern int idf_frmt_amode    ();

/*idf_frmt.c*/
extern IDF_FRMT* idf_frmt_address();
extern int       idf_frmt_ini    ();
extern void      idf_frmt_end    ();
extern void      idf_frmt_clean  ();
extern int       idf_frmt_prp    ();
extern void      idf_frmt_prn    ();
extern void      idf_form_view   ();
extern int       idf_frmt_noalloc();
extern void      idf_frmt_err_prn();

/*idf_frmtget.c*/
extern IDF_FRMT_TREAT* idf_frmtreat_address();
extern void            idf_frmtreat_nul    ();
extern void            idf_frmtreat_prn    ();
extern int             idf_frmt_begin      ();
extern int             idf_frmt_set_pntr   ();
extern void            idf_frmt_put_pntr   ();
extern char*           idf_frmt_get_pntr   ();
extern int             idf_frmtreat_array  ();
extern int             idf_frmt_get        ();
extern void            idf_frmtreat_err_prn();

/*idf_sdata.c*/
extern int idf_skip_data       ();
extern int idf_skip_data_period();

/*idf_err.c*/
extern void idf_err_clear    ();
extern void idf_err_put      ();
extern void idf_err_prn      ();

/*idf_geterr.c*/
extern void idf_get_err_prn  ();
extern void idf_list_err_prn ();

/*idf_txterr.c*/
extern void idf_txt_err_nul  ();
extern void idf_txt_err_put  ();
extern void idf_txt_err_prn  ();

/*idf_name.c*/
extern IDF_NAME* idf_name_address();
extern int       idf_name_ini    ();
extern void      idf_name_end    ();
extern void      idf_name_prn    ();
extern void      idf_name_errprn ();
extern void      idf_name_begin  ();
extern int       idf_name_buf    ();
extern int       idf_name_raw    ();
extern int       idf_name        ();
extern int       idf_name_dpos   ();
extern int       idf_name_get    ();

/*idf_extname.c*/
extern int       idf_extname_check ();

/* idf_fil.c*/
extern IDF_FILE* idf_file_address ();
extern char*     idf_file_Cname   ();
extern int       idf_file_ini     ();
extern void      idf_file_end     ();
extern void      idf_file_prn     ();
extern int       idf_file_open    ();
extern int       idf_file_name    ();
extern void      idf_file_Cname_free();
extern int       idf_file_name_   ();
extern void      idf_file_close   ();
extern int       idf_file_opened  ();
extern int       idf_file_ferr    ();
extern int       idf_file_getc    ();
extern int       idf_file_gets    ();
extern int       idf_file_ftell   ();
extern int       idf_file_fseek   ();
extern int       idf_file_rewind  ();
extern void      idf_file_errprn  ();
extern void      idf_file_clearerr();

/* idf_list.c*/
extern IDF_NAME_LIST* idf_list_address();
extern int            idf_list_ini    ();
extern void           idf_list_end    ();
extern void           idf_list_clean  ();
extern int            idf_list_num    ();
extern int            idf_list_numL   ();
extern IDF_NAME_*     idf_list_get    ();
extern IDF_NAME_L*    idf_list_getL   ();
extern int            idf_list_getLV  ();
extern int            idf_list_put    ();
extern int            idf_list_putL   ();
extern void           idf_list_prn_   ();
extern void           idf_list_prn    ();
extern void           idf_list_prnL_  ();
extern void           idf_list_prnL   ();
extern int            idf_list_do     ();
extern int            idf_list_check  ();
extern void           idf_list_catalog ();
extern void           idf_list_catalogL();
extern void           idf_list_clearerr();

/*idf_mlist.c*/
extern int idf_mlist  ();

/*idf_plist.c*/
extern int idf_plist  ();

/*idf_clist.c*/
extern int idf_clist  ();

/*idf_dbuf.c*/
extern IDF_DBUF* idf_dbuf_address();
extern int       idf_dbuf_ini    ();
extern void      idf_dbuf_end    ();
extern void      idf_dbuf_prn    ();
extern void      idf_dbuf_nul    ();
extern void      idf_dbuf_begin  ();
extern int       idf_dbuf_number ();
extern int       idf_dbuf_char   ();
extern int       idf_dbuf_txt    ();
extern int       idf_dbuf_str    ();
extern int       idf_dbuf_cmplx  ();
extern int       idf_dbuf_add    ();
extern int       idf_dbuf_to_d_  ();
extern int       idf_dbuf_to_f_  ();
extern int       idf_dbuf_to_l_  ();
extern int       idf_dbuf_to_i_  ();
extern int       idf_dbuf_to_d   ();
extern int       idf_dbuf_to_f   ();
extern int       idf_dbuf_to_l   ();
extern int       idf_dbuf_to_i   ();
extern int       idf_dbuf_to_s   ();
extern int       idf_dbuf_to_t   ();
extern int       idf_dbuf_to_c   ();
extern int       idf_dbuf_to_z   ();
extern int       idf_dbuf_to_w   ();
extern int       idf_dbuf_to_v   ();

/*idf_fbuf.c*/
extern IDF_FBUF* idf_fbuf_address();
extern int       idf_fbuf_ini    ();
extern void      idf_fbuf_end    ();
extern void      idf_fbuf_nul    ();
extern void      idf_fbuf_prn    ();
extern void      idf_fbuf_begin  ();
extern int       idf_fbuf_number ();
extern int       idf_fbuf_name   ();
extern int       idf_fbuf_conv   ();

/*idf_form.c*/
extern IDF_FORM* idf_form_address();
extern int       idf_form_ini    ();
extern void      idf_form_end    ();
extern void      idf_form_nul    ();
extern void      idf_form_prn    ();
extern int       idf_form_do     ();
extern int       idf_form_up     ();

/*idf_data.c*/
extern IDF_DATA* idf_data_address();
extern int       idf_data_ini    ();
extern void      idf_data_end    ();
extern void      idf_data_prn0   ();
extern void      idf_data_prn1   ();
extern int       idf_data_begin  ();
extern int       idf_data_beginL ();
extern int       idf_data_unit   ();
extern int       idf_data_value  ();
extern int       idf_data_up     ();

/*idf_val.c*/
extern IDF_VAL*  idf_val_address();
extern int       idf_val_ini    ();
extern void      idf_val_end    ();
extern void      idf_val_prn    ();
extern int       idf_val_put    ();
extern int       idf_val_buf    ();
extern int       idf_val_form   ();
extern int       idf_val_get    ();
extern int       idf_val_lcns   ();

/*idf_fflist.c*/
extern int idf_list_ffl ();

/*idf_catal.c*/
extern int idf_catalog ();

/*idf_order.c*/
extern void idf_order_put ();
extern int  idf_order_cmp ();


#endif


#endif /*IDF_H*/
