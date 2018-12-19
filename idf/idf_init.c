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


#include "idflib.h"


#include "idf.h"
#include "idfusr.h"



int idf_init()
{
int ierr=0;
int jerr;

idf_err_clear();

/*-------------------
  init internal data
 --------------------*/

   jerr = idf_file_ini();
if(jerr) {ierr=2; jerr=5+jerr; goto err;} /*err:6,7*/

   jerr = idf_frmt_ini();
if(jerr) {ierr=3; jerr=7+jerr; goto err;} /*err:8,9,10*/

   jerr = idf_name_ini();
if(jerr) {ierr=4; jerr=10+jerr; goto err;}/*err:11,12*/

          idf_list_nul();
   jerr = idf_list_ini();
if(jerr) {ierr=5; jerr=13; goto err;} /*err:13*/

   jerr = idf_data_ini();
if(jerr) {ierr=6; jerr=13+jerr; goto err;}/*err:14,15*/

   jerr = idf_dbuf_ini();
if(jerr) {ierr=7; jerr=15+jerr; goto err;}/*err:16,17*/

   jerr = idf_val_ini();
if(jerr) {ierr=8; jerr=17+jerr; goto err;}/*err:18,19*/

   jerr = idf_fbuf_ini();
if(jerr) {ierr=7; jerr=25+jerr; goto err;}/*err:26,27*/

   jerr = idf_form_ini();
if(jerr) {ierr=8; jerr=27+jerr; goto err;}/*err:28,29*/

#if IDF_FLOAT_DEFAULT == 0
   jerr = idf_form_ini();
if(jerr) {ierr=7; jerr=30; goto err;}/*err:30*/
#endif

idf_order_put( IDF_ORDER_INIT );

err:
if(ierr)
 {
  idf_err_put(ierr,jerr);

#if IDF_MISTAKE == 1
  idf_err_prn();
#endif
 }

return ierr;
}



#if IDF_CPP_ == 1
int idf_open(char* IDFfile_name)
#else
int idf_open(IDFfile_name)
char *IDFfile_name;
#endif
{
/*----------------------------------
  open IDFfile and create name list
 -----------------------------------*/

int ierr=0;
int jerr;

 /* undef error flags */
 idf_txt_err_nul  ();
 idf_file_clearerr();
 idf_list_clearerr();

   jerr = idf_order_cmp( IDF_ORDER_INIT );
if(jerr)
 {
  ierr=10;
  jerr=21;
  goto err;
 }

   jerr = idf_file_name(IDFfile_name);
if(jerr) {ierr=1;goto err;} /*err:1-5*/

   jerr = idf_file_open(IDFfile_name);
if(jerr) {ierr=9; jerr=20; goto err;}/*err:20*/


   jerr = idf_catalog();
if(jerr) {ierr=11; jerr=21+jerr;} /*22-25*/

idf_order_put( IDF_ORDER_OPEN );

idf_jus_c();

err:
if(ierr)
 {
  idf_err_put(ierr,jerr);

#if IDF_MISTAKE == 1
  idf_err_prn();
#endif
 }

return ierr;
}



void idf_close()
{
 idf_file_close();
 idf_list_clean();

#if IDF_MISTAKE == 1
 idf_err_prn();
#endif

 idf_order_put( IDF_ORDER_INIT );
}



void idf_finish()
{
#if IDF_MISTAKE == 1
 idf_err_prn();
#endif

#if IDF_FLOAT_DEFAULT == 0
 idf_float_end();
#endif

 idf_form_end();
 idf_fbuf_end();
 idf_val_end ();
 idf_dbuf_end();
 idf_data_end();
 idf_list_end();
 idf_name_end();
 idf_frmt_end();
 idf_file_end();

 idf_err_clear();

 idf_order_put( IDF_ORDER_START );
}
