C
C     COMMON-BLOCKS FROM DWD-PHYSICS
C
C     ORG: NZT, NHORPH
C     COMDYN: DT
C     COMPHY: LAKKU
C
      SUBROUTINE ECACCU(NAKTION,
     &    KE    , IEJE   ,
     &    ACLCAC, TEFF  , BFLHS , BFLQDS , AHFL   , AHFS   , APRC   ,
     &    APRL  , APRS  , EVAP  , DEW2   , TEMP2  , USTAR3 , USTR   ,
     &    U10   , VDIS  , VSTR  , V10    , WIND10 , ACLCOV , ALWCVI ,
     &    QVI   , SCLF0 , SCLFS , SRAF0  , SRAFS  , TCLF0  , TCLFS  ,
     &    TRAF0 , TRAFS , SRADS , SRADSU , SRAD0  , SRAD0U , TRADS  ,
     &    TRADSU, TRAD0 , DSNAC , RUNOFF , SNMEL  , TSLIN  , TSURF  ,
     &    USTRGW, VDISGW, VSTRGW, T2MAX  , T2MIN  , TSMAX  , TSMIN  ,
     &    WIMAX , TOPMAX, DRAIN , USTRL  , USTRW  , USTRI  , VSTRL  ,
     &    VSTRW , VSTRI , EVAPL , EVAPW  , EVAPI  ,
     &    AHFSL , AHFSW , AHFSI , BFLQDSL, BFLQDSW, BFLQDSI,
     &    BFLHSL, BFLHSW, BFLHSI, AHFICE , QRES   , QDBOXS , QWBOXS ,
     &    EKBOXS, FHBOXS, FIBOXS, VBM10M , CAPE   , ETRANS , EBSOIL ,
     &    ESNOW , ESKIN , ERES  , QIVI   , QIBOXS , RPRAC)
C 
      IMPLICIT NONE
C
      INCLUDE "org.h"
      INCLUDE "comdyn.h"
      INCLUDE "comphy.h"
C
      INTEGER, INTENT(IN) :: IEJE, NAKTION, KE
      REAL, DIMENSION(IEJE,KE), INTENT(INOUT) :: ACLCAC
      REAL, DIMENSION(IEJE),    INTENT(INOUT) :: 
     &  ACLCOV,AHFICE,AHFL,AHFS,
     &  AHFSI, AHFSL, AHFSW, ALWCVI, APRC, APRL, APRS, BFLHS,
     &  BFLHSI, BFLHSL, BFLHSW, BFLQDS, BFLQDSI, BFLQDSL, BFLQDSW
      REAL, DIMENSION(IEJE),    INTENT(INOUT) :: 
     &  CAPE,DEW2,DRAIN,DSNAC, 
     &  EBSOIL, EKBOXS, ERES, ESKIN, ESNOW, ETRANS, EVAP, EVAPI, EVAPL, 
     &  EVAPW, FHBOXS, FIBOXS
      REAL, DIMENSION(IEJE),    INTENT(INOUT) :: 
     &  QDBOXS, QIBOXS, QIVI, 
     &  QRES, QVI, QWBOXS, RUNOFF, SCLF0, SCLFS, SNMEL, SRAD0, SRAD0U, 
     &  SRADS, SRADSU, SRAF0, SRAFS, T2MAX, T2MIN, TCLF0, TCLFS
      REAL, DIMENSION(IEJE),    INTENT(INOUT) :: 
     &  TEFF, TEMP2, TOPMAX, 
     &  TRAD0, TRADS, TRADSU, TRAF0, TRAFS, TSLIN, TSMAX, TSMIN, 
     &  TSURF, U10, USTAR3, USTR, USTRGW, USTRI, USTRL, USTRW, V10
      REAL, DIMENSION(IEJE),    INTENT(INOUT) :: 
     &  VBM10M, VDIS, VDISGW, 
     &  VSTR, VSTRGW, VSTRI, VSTRL, VSTRW, WIMAX, WIND10 
      REAL, INTENT(INOUT)    :: RPRAC   (IEJE,KE)
C
      INTEGER :: I, K
      REAL    :: FAK, FAKDIA, FAKDIM, FAKFL, ZDT, ZTIMEA
C
!      DIMENSION ACLCAC  (IEJE,KE)
!      DIMENSION TEFF    (IEJE)
!      DIMENSION BFLHS   (IEJE)
!      DIMENSION BFLQDS  (IEJE)
!      DIMENSION AHFL    (IEJE)
!      DIMENSION AHFS    (IEJE)
!      DIMENSION APRC    (IEJE)
!      DIMENSION APRL    (IEJE)
!      DIMENSION APRS    (IEJE)
!      DIMENSION EVAP    (IEJE)
!      DIMENSION DEW2    (IEJE)
!      DIMENSION TEMP2   (IEJE)
!      DIMENSION USTAR3  (IEJE)
!      DIMENSION USTR    (IEJE)
!      DIMENSION U10     (IEJE)
!      DIMENSION VDIS    (IEJE)
!      DIMENSION VSTR    (IEJE)
!      DIMENSION V10     (IEJE)
!      DIMENSION WIND10  (IEJE)
!      DIMENSION ACLCOV  (IEJE)
!      DIMENSION ALWCVI  (IEJE)
!      DIMENSION QVI     (IEJE)
!      DIMENSION SCLF0   (IEJE)
!      DIMENSION SCLFS   (IEJE)
!      DIMENSION SRAF0   (IEJE)
!      DIMENSION SRAFS   (IEJE)
!      DIMENSION TCLF0   (IEJE)
!      DIMENSION TCLFS   (IEJE)
!      DIMENSION TRAF0   (IEJE)
!      DIMENSION TRAFS   (IEJE)
!      DIMENSION SRADS   (IEJE)
!      DIMENSION SRADSU  (IEJE)
!      DIMENSION SRAD0   (IEJE)
!      DIMENSION SRAD0U  (IEJE)
!      DIMENSION TRADS   (IEJE)
!      DIMENSION TRADSU  (IEJE)
!      DIMENSION TRAD0   (IEJE)
!      DIMENSION DSNAC   (IEJE)
!      DIMENSION RUNOFF  (IEJE)
!      DIMENSION DRAIN   (IEJE)
!      DIMENSION SNMEL   (IEJE)
!      DIMENSION TSLIN   (IEJE)
!      DIMENSION TSURF   (IEJE)
!      DIMENSION USTRGW  (IEJE)
!      DIMENSION VDISGW  (IEJE)
!      DIMENSION VSTRGW  (IEJE)
!C
!      DIMENSION T2MAX   (IEJE)
!      DIMENSION T2MIN   (IEJE)
!      DIMENSION TSMAX   (IEJE)
!      DIMENSION TSMIN   (IEJE)
!      DIMENSION WIMAX   (IEJE)
!      DIMENSION TOPMAX  (IEJE)
!      DIMENSION USTRL   (IEJE)
!      DIMENSION USTRW   (IEJE)
!      DIMENSION USTRI   (IEJE)
!      DIMENSION VSTRL   (IEJE)
!      DIMENSION VSTRW   (IEJE)
!      DIMENSION VSTRI   (IEJE)
!      DIMENSION EVAPL   (IEJE)
!      DIMENSION EVAPW   (IEJE)
!      DIMENSION EVAPI   (IEJE)
!      DIMENSION AHFSL   (IEJE)
!      DIMENSION AHFSW   (IEJE)
!      DIMENSION AHFSI   (IEJE)
!      DIMENSION BFLQDSL (IEJE)
!      DIMENSION BFLQDSW (IEJE)
!      DIMENSION BFLQDSI (IEJE)
!      DIMENSION BFLHSL  (IEJE)
!      DIMENSION BFLHSW  (IEJE)
!      DIMENSION BFLHSI  (IEJE)
!      DIMENSION AHFICE  (IEJE)
!      DIMENSION QRES    (IEJE)
!      DIMENSION QDBOXS  (IEJE)
!      DIMENSION QWBOXS  (IEJE)
!      DIMENSION EKBOXS  (IEJE)
!      DIMENSION FHBOXS  (IEJE)
!      DIMENSION FIBOXS  (IEJE)
!      DIMENSION VBM10M  (IEJE)
!      DIMENSION CAPE    (IEJE)
!      DIMENSION ETRANS  (IEJE)
!      DIMENSION EBSOIL  (IEJE)
!      DIMENSION ESNOW   (IEJE)
!      DIMENSION ESKIN   (IEJE)
!      DIMENSION ERES    (IEJE)
!      DIMENSION QIVI    (IEJE)
!      DIMENSION QIBOXS  (IEJE)
C
      IF (NZT.EQ.0) THEN
         ZDT=2.*DT
      ELSE
         ZDT=DT
      ENDIF
      ZTIMEA=FLOAT(NDMXN)
C
      IF (NZT.EQ.0 .AND. NEAA.LE.1) THEN
         ZTIMEA=1.
      ENDIF
C
      FAKFL=1.
      IF (NAKTION.EQ.1) THEN
         FAK=1.0/(ZTIMEA*ZDT)
         FAKDIA=FAK
         IF (LNEAR) FAKDIA=1.
         FAKDIM=1000.
      ENDIF
      IF (NAKTION.EQ.2) THEN
         FAK=(ZTIMEA*ZDT)
         FAKDIA=FAK
         IF (LNEAR) FAKDIA=1.
         FAKDIM=0.001
      ENDIF
      IF (NAKTION.EQ.3) THEN
         FAKDIA=0.
         FAK=0.
         IF (LAKKU) THEN
            FAKDIM=1.
         ELSE
            FAKDIM=0.
         ENDIF
         FAKFL=0.
      ENDIF
C
      DO I=1,IEJE
         TEFF   (I) = TEFF    (I) * FAK
         BFLHS  (I) = BFLHS   (I) * FAKFL
         BFLHSL (I) = BFLHSL  (I) * FAKFL
         BFLHSW (I) = BFLHSW  (I) * FAKFL
         BFLHSI (I) = BFLHSI  (I) * FAKFL
         BFLQDS (I) = BFLQDS  (I) * FAKFL
         BFLQDSL(I) = BFLQDSL (I) * FAKFL
         BFLQDSW(I) = BFLQDSW (I) * FAKFL
         BFLQDSI(I) = BFLQDSI (I) * FAKFL
         AHFS   (I) = AHFS    (I) * FAKFL
         AHFSL  (I) = AHFSL   (I) * FAKFL
         AHFSW  (I) = AHFSW   (I) * FAKFL
         AHFSI  (I) = AHFSI   (I) * FAKFL
         AHFL   (I) = AHFL    (I) * FAKFL
         APRC   (I) = APRC    (I) * FAKDIM
         APRL   (I) = APRL    (I) * FAKDIM
         APRS   (I) = APRS    (I) * FAKDIM
         EVAP   (I) = EVAP    (I) * FAKDIM
         EVAPL  (I) = EVAPL   (I) * FAKDIM
         EVAPW  (I) = EVAPW   (I) * FAKDIM
         EVAPI  (I) = EVAPI   (I) * FAKDIM
         DEW2   (I) = DEW2    (I) * FAKDIA
         TEMP2  (I) = TEMP2   (I) * FAKDIA
         USTAR3 (I) = USTAR3  (I) * FAK
         USTR   (I) = USTR    (I) * FAK
         USTRL  (I) = USTRL   (I) * FAK
         USTRW  (I) = USTRW   (I) * FAK
         USTRI  (I) = USTRI   (I) * FAK
         U10    (I) = U10     (I) * FAKDIA
         VDIS   (I) = VDIS    (I) * FAK
         VSTR   (I) = VSTR    (I) * FAK
         VSTRL  (I) = VSTRL   (I) * FAK
         VSTRW  (I) = VSTRW   (I) * FAK
         VSTRI  (I) = VSTRI   (I) * FAK
         V10    (I) = V10     (I) * FAKDIA
         WIND10 (I) = WIND10  (I) * FAKDIA
         ACLCOV (I) = ACLCOV  (I) * FAK
         ALWCVI (I) = ALWCVI  (I) * FAK
         QVI    (I) = QVI     (I) * FAK
         SCLF0  (I) = SCLF0   (I) * FAK
         SCLFS  (I) = SCLFS   (I) * FAK
         SRAF0  (I) = SRAF0   (I) * FAK
         SRAFS  (I) = SRAFS   (I) * FAK
         TCLF0  (I) = TCLF0   (I) * FAK
         TCLFS  (I) = TCLFS   (I) * FAK
         TRAF0  (I) = TRAF0   (I) * FAK
         TRAFS  (I) = TRAFS   (I) * FAK
         SRADS  (I) = SRADS   (I) * FAK
         SRADSU (I) = SRADSU  (I) * FAK
         SRAD0  (I) = SRAD0   (I) * FAK
         SRAD0U (I) = SRAD0U  (I) * FAK
         TRADS  (I) = TRADS   (I) * FAK
         TRADSU (I) = TRADSU  (I) * FAK
         TRAD0  (I) = TRAD0   (I) * FAK
         DSNAC  (I) = DSNAC   (I) * FAKDIM
         RUNOFF (I) = RUNOFF  (I) * FAKDIM
         DRAIN  (I) = DRAIN   (I) * FAKDIM
         SNMEL  (I) = SNMEL   (I) * FAKDIM
         TSLIN  (I) = TSLIN   (I) * FAK
         AHFICE (I) = AHFICE  (I) * FAK
         QRES   (I) = QRES    (I) * FAK
         TSURF  (I) = TSURF   (I) * FAK
         USTRGW (I) = USTRGW  (I) * FAK
         VDISGW (I) = VDISGW  (I) * FAK
         VSTRGW (I) = VSTRGW  (I) * FAK
         QDBOXS (I) = QDBOXS  (I) * FAKFL
         QWBOXS (I) = QWBOXS  (I) * FAKFL
         EKBOXS (I) = EKBOXS  (I) * FAKFL
         FHBOXS (I) = FHBOXS  (I) * FAKFL
         FIBOXS (I) = FIBOXS  (I) * FAKFL
         CAPE   (I) = CAPE    (I) * FAKDIA
         ETRANS (I) = ETRANS  (I) * FAKDIM
         EBSOIL (I) = EBSOIL  (I) * FAKDIM
         ESNOW  (I) = ESNOW   (I) * FAKDIM
         ESKIN  (I) = ESKIN   (I) * FAKDIM
         ERES   (I) = ERES    (I) * FAKDIM
         QIVI   (I) = QIVI    (I) * FAK
         QIBOXS (I) = QIBOXS  (I) * FAKFL
      ENDDO
C
      DO K=1,KE
         DO I=1,IEJE
            ACLCAC  (I,K) = ACLCAC  (I,K) * FAK
            RPRAC   (I,K) = RPRAC   (I,K) * FAK
         ENDDO
      ENDDO
C
      IF (NAKTION.EQ.3) THEN
         DO I=1,IEJE
            T2MAX(I)=-99999.
            T2MIN(I)=99999.
            TSMAX(I)=-99999.
            TSMIN(I)=99999.
            WIMAX(I)=0.
            VBM10M(I)=0.
            TOPMAX(I)=99999.
         ENDDO
      ENDIF
C
      RETURN
C
      END SUBROUTINE ECACCU
