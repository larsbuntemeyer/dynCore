!

! Code converted using TO_F90 by Alan Miller
! Date: 2015-09-07  Time: 14:27:24

!     SUBROUTINE GAUSS

!@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
!@(#) MODUL GAUSS.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94

!**** GAUSS    -   UP:LOESUNG TRIDIAGONALER GLEICHUNGSSYSTEME
!**   AUFRUF   :   CALL GAUSS (IANF,IEND,M,IAGA)
!**   ENTRIES  :      ---
!**   ZWECK    :   PARALLELE LOESUNG VON IANF-IEND+1 SYSTEMEN VON KE
!**                SIMULTANEN DREIPUNKT-DIFFERENZENGLEICHUNGEN (EIN-
!**                DIMENSIONALE RANDWERTAUFGABE)
!**   VERSIONS-
!**   DATUM    :   20. 2. 89
!**
!**   EXTERNALS:   KEINE
!**
!**   EINGABE-
!**   PARAMETER:   IANF: INDEX DES ERSTEN TRIDIAGONALEN SYSTEMS
!**                IEND: INDEX DES LETZTEN     ,,         ,,
!**                M   : ( = 1,3) BEZEICHNET EINEN BESTIMMTEN SATZ VON
!**                      GLEICHUNGSSYSTEMEN
!**
!**   PARAMETER:      ---
!**
!**   COMMON-
!**   BLOECKE  :   PARAM
!**
!**   METHODE  :   ZEILENWEISE INVERSION ( ELIMINATIONSVERFAHREN VON
!**                GAUSS )
!**   FEHLERBE-
!**   HANDLUNG :      ---
!**   VERFASSER:   G. DOMS

SUBROUTINE gauss (ianf,iend,m,aga,agb,agc,agd,age)

USE Hydro
USE Grid

IMPLICIT NONE

INTEGER, INTENT(IN)        :: ianf
INTEGER, INTENT(IN)        :: iend
INTEGER, INTENT(IN OUT)    :: m
REAL, INTENT(IN)           :: aga(ie,ke,4)
REAL, INTENT(IN)           :: agb(ie,ke,4)
REAL, INTENT(OUT)       :: agc(ie,ke,4)
REAL, INTENT(OUT)       :: agd(ie,ke,4)
REAL, INTENT(OUT)       :: age(ie,ke,4)


!     ERWARTET WERDEN DIE KOEFFIZIENTENMATRIZEN AGA(I,K,M), AGB(I,K,M)
!     UND AGC(I,K,M) SOWIE DIE INHOMOGENEN VEKTOREN AGD(I,K,M).
!     DIE LOESUNGSVEKTOREN ERSCHEINEN AUF AGE(I,K,M).



!     DIMENSIONIERUNG DER ERFORDERLICHEN FELDER





!     Local Variables

INTEGER :: i,k
REAL :: zzz

!     TRANSFORMATION DER KOEFFIZINENTENMATRIZEN IN EINE OBERE DREIECKS-
!     FORM MIT 1 IN DER DIAGONALEN

DO i  = ianf , iend
  agc(i,1,m) = agc(i,1,m)/agb(i,1,m)
  agd(i,1,m) = agd(i,1,m)/agb(i,1,m)
END DO
DO k  =   2 , ke-1
  DO i  = ianf, iend
    zzz        = 1. / ( agb(i,k,m) - aga(i,k,m)*agc(i,k-1,m) )
    agc(i,k,m) = agc(i,k,m)*zzz
    agd(i,k,m) = ( agd(i,k,m) - aga(i,k,m)*agd(i,k-1,m) )*zzz
  END DO
END DO
DO i  = ianf, iend
  age(i,ke,m) = ( agd(i,ke,m) - aga(i,ke,m)*agd(i,ke-1,m) ) /  &
      ( agb(i,ke,m) - aga(i,ke,m)*agc(i,ke-1,m) )
END DO

!     LOESUNG DURCH RUECKWAERTSSUBSTITUTION

DO k  = ke-1 , 1 , -1
  DO i  = ianf , iend
    age(i,k,m) = agd(i,k,m) - agc(i,k,m)*age(i,k+1,m)
  END DO
END DO

RETURN
END SUBROUTINE gauss
