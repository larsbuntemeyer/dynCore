      BLOCK DATA VARN
      IMPLICIT NONE
      INCLUDE "corg.h"
C
      INTEGER :: NTAB
C
      DATA    YDVARN      , YTVARN     , YMVARN     , YNVARN
     1      / NZMXVD*'  ' , NZMXVT*'  ', NZMXVM*'  ', NZMXVN*'  ' /
C
      DATA        (YAVARN(NTAB), NTAB=1,NZMXVA) /
     1           'U       ','V       ','PS      ','T       ','QD      ',
     2           'QW      ','TSLECH  ','TSWECH  ','TSIECH  ','TSN     ',
     3           'TD3     ','TD4     ','TD5     ','TD      ','TDCL    ',
     4           'WSECH   ','SN      ','WL      ','FIB     ','BLA     ',
     5           'GLAC    ','AZ0     ','VGRAT   ','FOREST  ','ALBECH  ',
     6           'WSMX    ','VLT     ','FAO     ','QDBL    ','SEAICE  ',
     7           'VAROR   ','BETA    ','WMINLOK ','WMAXLOK ','VPOR    ',
     8           'WS1     ','WS2     ','WS3     ','WS4     ','WS5     ',
     9           'DZR     ','DZS     ','FKSAT   ','FMPOT   ','BCLAPP  ',
     A           'FTKVM   ','FTKVH   ','FI      ','VERVEL  ','APRL    ',
     1           'APRC    ','APRS    ','ACLC    ','ACLCV   ','ALBEDO  ',
     2           'EVAP    ','BFLHS   ','T2MAX   ','T2MIN   ','TSMAX   ',
     3           'TSMIN   ','WIMAX   ','SICED   ','TOPMAX  ','TKE     ',
     4           'RGCGN   ','EVAPM   ','DSNAC   ','PHI     ','TMCH    ',
     5           'RLA     ','ACLCAC  ','DEW2    ','TEMP2   ','USTAR3  ',
     6           'USTR    ','VSTR    ','U10     ','V10     ','VDIS    ',
     7           'WIND10  ','ACLCOV  ','ALWCVI  ','QVI     ','SNMEL   ',
     8           'TSLIN   ','SRADS   ','SRADSU  ','SRAD0   ','SRAD0U  ',
     9           'TRADS   ','TRADSU  ','TRAD0   ','SCLF0   ','SCLFS   ',
     B           'SRAF0   ','SRAFS   ','TCLF0   ','TCLFS   ','TRAF0   ',
     1           'TRAFS   ','EMTER   ','TRSOL   ','EMTEF   ','TRSOF   ',
     2           'TSURF   ','USTRGW  ','VSTRGW  ','VDISGW  ','AHFL    ',
     3           'AHFS    ','TEFF    ','DRAIN   ','USTRL   ','USTRW   ',
     4           'USTRI   ','VSTRL   ','VSTRW   ','VSTRI   ','EVAPL   ',
     5           'EVAPW   ','EVAPI   ','AHFSL   ','AHFSW   ','AHFSI   ',
     6           'AZ0L    ','AZ0W    ','AZ0I    ','ALSOL   ','ALSOW   ',
     7           'ALSOI   ','AHFICE  ','QRES    ','TMCHL   ','TMCHW   ',
     8           'TMCHI   ','QDB     ','QDBW    ','QDBI    ','TSECH   ',
     9           'BFLHSL  ','BFLHSW  ','BFLHSI  ','SRFL    ','BFLQDSL ',
     C           'BFLQDSW ','BFLQDSI ','QDBOXS  ','QWBOXS  ','EKBOXS  ',
     1           'FHBOXS  ','FIBOXS  ','TLAMBDA ','DLAMBDA ','PORVOL  ',
     2           'FCAP    ','WI3     ','WI4     ','WI5     ','WI      ',
     3           'WICL    ','BFLQDS  ','RUNOFF  ','TMCM    ','CAPE    ',
     4           'VBM10M  ','ETRANS  ','EBSOIL  ','ESNOW   ','ESKIN   ',
     5           'ERES    ','QI      ','QIVI    ','QIBOXS  ','W       ',
     6           'DWDT    ','PINT    ','RPRAC   ' /
C
      DATA        (YRVARN(NTAB), NTAB=1,9) /
     1           'U       ','V       ','PS      ','T       ','QD      ',
     2           'TSWECH  ','TSIECH  ','SEAICE  ','QDBL    '/
C
      DATA        (YFVARN(NTAB), NTAB=1,NZMXVF) /
     1           'U       ','V       ','PS      ','T       ','QD      ',
     2           'QW      ','TSLECH  ','TSWECH  ','TSIECH  ','TSN     ',
     3           'TD3     ','TD4     ','TD5     ','TD      ','TDCL    ',
     4           'WSECH   ','SN      ','WL      ','FIB     ','BLA     ',
     5           'GLAC    ','AZ0     ','VGRAT   ','FOREST  ','ALBECH  ',
     6           'WSMX    ','VLT     ','FAO     ','QDBL    ','SEAICE  ',
     7           'VAROR   ','BETA    ','WMINLOK ','WMAXLOK ','VBM10M  ',
     8           'FTKVM   ','FTKVH   ','FI      ','VERVEL  ','APRL    ',
     9           'APRC    ','APRS    ','ACLC    ','ACLCV   ','ALBEDO  ',
     A           'EVAP    ','BFLHS   ','T2MAX   ','T2MIN   ','TSMAX   ',
     1           'TSMIN   ','WIMAX   ','SICED   ','TOPMAX  ','TKE     ',
     2           'RGCGN   ','EVAPM   ','DSNAC   ','PHI     ','TMCH    ',
     3           'RLA     ','ACLCAC  ','DEW2    ','TEMP2   ','USTAR3  ',
     4           'USTR    ','VSTR    ','U10     ','V10     ','VDIS    ',
     5           'WIND10  ','ACLCOV  ','ALWCVI  ','QVI     ','SNMEL   ',
     6           'TSLIN   ','SRADS   ','SRADSU  ','SRAD0   ','SRAD0U  ',
     7           'TRADS   ','TRADSU  ','TRAD0   ','SCLF0   ','SCLFS   ',
     8           'SRAF0   ','SRAFS   ','TCLF0   ','TCLFS   ','TRAF0   ',
     9           'TRAFS   ','EMTER   ','TRSOL   ','EMTEF   ','TRSOF   ',
     B           'TSURF   ','USTRGW  ','VSTRGW  ','VDISGW  ','AHFL    ',
     1           'AHFS    ','TEFF    ','DRAIN   ','USTRL   ','USTRW   ',
     2           'USTRI   ','VSTRL   ','VSTRW   ','VSTRI   ','EVAPL   ',
     3           'EVAPW   ','EVAPI   ','AHFSL   ','AHFSW   ','AHFSI   ',
     4           'AZ0L    ','AZ0W    ','AZ0I    ','ALSOL   ','ALSOW   ',
     5           'ALSOI   ','AHFICE  ','QRES    ','TMCHL   ','TMCHW   ',
     6           'TMCHI   ','QDB     ','QDBW    ','QDBI    ','TSECH   ',
     7           'BFLHSL  ','BFLHSW  ','BFLHSI  ','SRFL    ','BFLQDSL ',
     8           'BFLQDSW ','BFLQDSI ','QDBOXS  ','QWBOXS  ','EKBOXS  ',
     9           'FHBOXS  ','FIBOXS  ','TLAMBDA ','DLAMBDA ','PORVOL  ',
     C           'FCAP    ','WI3     ','WI4     ','WI5     ','WI      ',
     1           'WICL    ','BFLQDS  ','RUNOFF  ','TMCM    ','CAPE    ',
     2           'WS1     ','WS2     ','WS3     ','WS4     ','WS5     ',
     3           'DZR     ','DZS     ','FKSAT   ','FMPOT   ','BCLAPP  ',
     4           'VPOR    ','ETRANS  ','EBSOIL  ','ESNOW   ','ESKIN   ',
     5           'ERES    ','QI      ','QIVI    ','QIBOXS  ','W       ',
     6           'DWDT    ','PINT    ','RPRAC   ' /
C
      END BLOCK DATA VARN
