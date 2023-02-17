      SUBROUTINE WRITE_C(nout,n_c,c,n_l)
!***********************************************************************
!WRITE_C writes out n_c character variables in fields of length n_l
!Input:
!  nout-output file unit number (-)
!  n_c-number of character variables
!  c()-character array
!  n_l-length of field
!  W.A. Houlberg 2/99
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  c(*)
      INTEGER        n_c,                     n_l,
     #               nout
!Declaration of local variables
      CHARACTER*30   char
      CHARACTER*2    c_c,                     c_l
      INTEGER        i
      WRITE(c_c,'(i2)') n_c
      WRITE(c_l,'(i2)') n_l
      char='('//c_c//'a'//c_l//')'
      WRITE(nout,char) (c(i),i=1,n_c)
      RETURN
      END
      SUBROUTINE WRITE_IR(nout,n_i,i,n_r,r,k_format)
!***********************************************************************
!WRITE_IR writes a writes n_i integer variables followed by n_r real
!  variables to the unit nout
!Input:
!  nout-output file unit number (-)
!  n_i-number of integer variables (-)
!  i()-integer array
!  n_r-number of real variables (-)
!  r()-real array
!  k_format-real format option
!          =1 use f
!          =2 use 1pe
!          =else use e
!  W.A. Houlberg 2/99
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        k_format,                n_i,
     #               n_r,                     nout
      INTEGER        i(*)
      REAL           r(*)
!Declaration of local variables
      CHARACTER*30   char
      CHARACTER*2    c_i,                     c_r
      INTEGER        j
      WRITE(c_i,'(i2)') n_i
      WRITE(c_r,'(i2)') n_r
      IF(n_i.eq.0) THEN
!Real data only
        IF(k_format.eq.1) THEN
          char='('//c_r//'(f12.6))'
        ELSEIF(k_format.eq.2) THEN
          char='('//c_r//'(1pe12.4))'
        ELSE
          char='('//c_r//'(e12.4))'
        ENDIF
        WRITE(nout,char) (r(j),j=1,n_r)
      ELSEIF(n_r.eq.0) THEN
!Integer data only
        char='('//c_i//'(i12))'
        WRITE(nout,char) (i(j),j=1,n_i)
      ELSE
!Both integer and real data
        IF(k_format.eq.1) THEN
          char='('//c_i//'(i12),'//c_r//'(f12.6))'
        ELSEIF(k_format.eq.2) THEN
          char='('//c_i//'(i12),'//c_r//'(1pe12.4))'
        ELSE
          char='('//c_i//'(i12),'//c_r//'(e12.4))'
        ENDIF
        WRITE(nout,char) (i(j),j=1,n_i),(r(j),j=1,n_r)
      ENDIF
      RETURN
      END
      SUBROUTINE WRITE_LINE(nout,label,nabove,nbelow)
!***********************************************************************
!WRITE_LINE writes a label preceeded by nabove blank lines and
!  followed by nbelow blank lines to the unit nout
!Input:
!  nout-output file unit number (-)
!  label-label to be printed
!  nabove-number of blank lines above label
!  nbelow-number of blanklines below label
!  W.A. Houlberg 2/99
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  label
      INTEGER        nabove,                  nbelow,
     #               nout
!Declaration of local variables
      INTEGER        j
      IF(nabove.gt.0) THEN
        DO j=1,nabove
          WRITE(nout,'( )')
        ENDDO
      ENDIF
      WRITE(nout,'(a)') label
      IF(nbelow.gt.0) THEN
        DO j=1,nbelow
          WRITE(nout,'( )')
        ENDDO
      ENDIF
      RETURN
      END
      SUBROUTINE WRITE_LINE_IR(nout,label,n_i,i,n_r,r,k_format)
!***********************************************************************
!WRITE_LINE_IR writes a writes a label followed by n_i integer numbers
!  and n_r real numbers to the unit nout
!Input:
!  nout-output file unit number
!  label-label to be printed
!  n_i-number of integer variables
!  i()-integer array
!  n_r-number of real variables
!  r()-real array
!  k_format-format option for real variables
!          =1 use f
!          =2 use 1pe
!          =else use e
!  W.A. Houlberg 12/98
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      CHARACTER*(*)  label
      INTEGER        k_format,                n_i,
     #               n_r,                     nout
      INTEGER        i(*)
      REAL           r(*)
!Declaration of local variables
      CHARACTER*30   char
      CHARACTER*2    c_i,                     c_r
      INTEGER        j
      WRITE(c_i,'(i2)') n_i
      WRITE(c_r,'(i2)') n_r
      IF(n_i.eq.0) THEN
!Real data only
        IF(k_format.eq.1) THEN
          char='(a48,'//c_r//'(f12.6))'
        ELSEIF(k_format.eq.2) THEN
          char='(a48,'//c_r//'(1pe12.4))'
        ELSE
          char='(a48,'//c_r//'(e12.4))'
        ENDIF
        WRITE(nout,char) label(1:48),(r(j),j=1,n_r)
      ELSEIF(n_r.eq.0) THEN
!Integer data only
        char='(a48,'//c_i//'(i12))'
        WRITE(nout,char) label(1:48),(i(j),j=1,n_i)
      ELSE
!Both integer and real data
        IF(k_format.eq.1) THEN
          char='(a48,'//c_i//'(i12),'//c_r//'(f12.6))'
        ELSEIF(k_format.eq.2) THEN
          char='(a48,'//c_i//'(i12),'//c_r//'(1pe12.4))'
        ELSE
          char='(a48,'//c_i//'(i12),'//c_r//'(e12.4))'
        ENDIF
        WRITE(nout,char) label(1:48),(i(j),j=1,n_i),(r(j),j=1,n_r)
      ENDIF
      RETURN
      END
