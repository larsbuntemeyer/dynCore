!*==remorg.f    processed by SPAG 6.71Rc at 15:45 on 12 Dec 2014
      !
      PROGRAM REMORG
      !
      IMPLICIT NONE
!*--REMORG7
      !
      ! REMORG   -   HP:ORGANISATION DER REMO-PROGNOSE
      ! AUFRUF   :    --
      ! ENTRIES  :   KEINE
      ! ZWECK    :   ORGANISATION DER GESAMTEN REMO-PROGNOSE VOM HOLEN DER
      !              ANFANGS- UND RANDDATEN BIS ZUM SCHREIBEN DER ERGEB-
      !              NISSE
      ! VERSIONS-
      ! DATUM    :   29.01.04
      !              2007
      !
      ! EXTERNALS:   PLATTE, MEDEA , KONST, DATUTC, EC4ORG
      !              INIMIT, INIZR, CLOSEZR
      ! EINGABE-
      ! PARAMETER:   INPUT FUER REMORG UEBER NAMELIST
      ! AUSGABE-
      ! PARAMETER:   KEINE
      !
      ! COMMON-
      ! BLOECKE  :   PARAM, ORG, CORG, COMDYN, COMDIA
      !
      ! METHODE  :   SETZEN VON STEUERPARAMETERN ABHAENGIG VOM ZEITSCHRITT
      !              SEQUENTIELLER AUFRUF VON UP'S ZUR PROGNOSE-DURCH-
      !              FUEHRUNG
      ! FEHLERBE-
      ! HANDLUNG :   STOP IM FEHLERFALL
      ! VERFASSER:   R.PODZUN
      !
      !----------------------------------------------------------------
      ! Local Varibales
      !----------------------------------------------------------------
      !
      ! This seems to define the MPI-ID that should
      ! read in boundary data
      !
      INTEGER :: MYIDTMP
      !
      ! Number of Processors for the planar grid decomposition
      ! in x- and y-direction
      !
      INTEGER  :: NPROCXM, NPROCYM
      !
      ! Some 'date'
      !
      INTEGER  :: NZTANF
      !
      ! Simple timer for total runtime
      !
      REAL     :: T0,T1,TTOTAL,TMIN,TMAX,TMEAN
      !
      ! External 'block data'. These will initialize
      ! some common block data.
      !
      EXTERNAL :: ECMGBDT, UNITDT, VARN, VARN2
      EXTERNAL :: SUAER, SUCFC, SULW, SURAD, SUSW
      !
      ! Include header files
      !
      INCLUDE "parorg.h" ! CONTAINS ORGANIZATIONAL DATA OF THE PARALLEL (NODE-) PROGRAM
      INCLUDE "param.h"  ! DIMENSIONEN DES REMO_MPI
      INCLUDE "org.h"    ! GROESSEN FUER DIE ORGANISATION DES MODELLS
      INCLUDE "corg.h"   ! CHARACTER-GROESSEN FUER DIE DATEIENVERWALTUNG
      INCLUDE "comdyn.h" ! STEUERPARAMETER DER ADIABATISCHEN PROGNOSE
      !not used INCLUDE "comdia.h" ! STEUERPARAMETER DER MODELL-DIAGNOSE
      !not used INCLUDE "comphy.h" ! STEUERPARAMETER DER DIABATISCHEN PROGNOSE
      INCLUDE "comphy.h"
      !
      CALL CPU_TIME(T0)
      !
      ! MPI Initialization
      !
      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERROR)
      !
      ! DATA-STATEMENTS ERSETZT DURCH ZUWEISUNGEN
      !
      NA = 3
      NJ = 1
      NE = 2
      NA2 = 2
      NJ2 = 1
      NRD1 = 1
      NRD2 = 2
      !
      ! Some logical stuff
      !
      LWRITEE = .FALSE.
      LWRITED = .FALSE.
      LWRITEF = .FALSE.
      LWRITET = .FALSE.
      !
      ! KONSTANTEN UMSPEICHERN
      !
      NHORPH = NHORPHD
      !
      ! Get some global MPI information
      !
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERROR)
      !
      REQUESTNR = 1
      REQUESTNRR = 1
      !
      ! IN REMO VERWENDETEN DATEIEN DEFINIEREN UND EROEFFNEN
      !
      MYIDTMP = 0
      !
      ! This loops over all MPI task, so each MPI task
      ! will read in its data one after the other.
      ! This obviously should prevent all MPI tasks to read
      ! from the harddisk at the same time.
      !
      DO
         IF (MYID==MYIDTMP) THEN
            !
            CALL PLATTE(MYID)
            !
            ! EINLESEN DER NAMELIST-AUFTRAGSKARTEN
            ! This routine reads in all data from the
            ! INPUT file
            !
            CALL MEDEA(MOIE,MOJE,MOKE,MOIEJE,MOIEJEKE,MOIEKE,MOKE1,MYID,
     &                 NPROCXM,NPROCYM)
            !
            ! SICHEHEITS-CHECK FUER MITTELWERTBILDUNG
            !
            IF (MOIEJE/=IDIMSUM) THEN
               PRINT *, 'PARAMETER IDIMSUM IS NOT EQUAL TO MOIE*MOJE !'
               PRINT *, 'IDIMSUM=', IDIMSUM, ' MOIEJE=', MOIEJE
               PRINT *, 'ERROR: CHECK COMMON-BLOCK PARAM.H'
               STOP
            ENDIF
            !
         ENDIF
         !
         ! All MPI tasks will wait for the current task
         ! to be finished with reading in data
         !
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         !
         ! This is a weird way to exit a loop.
         ! Should be "do MYIDTMP=0,NPROC-1" from the beginning
         ! (See above)
         !
         IF (MYIDTMP==(NPROC-1)) EXIT
         MYIDTMP = MYIDTMP + 1
         !
      ENDDO
      !
      !
      ! COMMON-BLOECKE DER MITTELUNG ANLEGEN
      !
      !
      IF (MYID==0) CALL INIMIT
      !
      ! KONSTANTEN BESETZEN
      !
      NPROCX = NPROCXM
      NPROCY = NPROCYM
      !
      ! Grid Decomposition
      !
      CALL DECOMP
      !
      ! KONSTANTEN BESETZEN
      !
      CALL KONST
      !
      ! BERECHNUNG DES AKTUELLEN DATUMS
      !
      NZT = NANF
      NZTANF = NANF
      IF (NANF>0) NZTANF = NANF + 1
      !
      CALL DATUTC(NZTANF,YADAT,DT,YAKDAT1,YAKDAT2,NAKJATA,AKHH)
      YANDAT = YAKDAT1
      !
      ! ZEITREIHEN-DATEIEN ANLEGEN
      !
      IF (MYID==0) CALL INIZR(NZTANF)
      !
      !================================================================
      ! SCHALTER FUER DIE JEWEILIGE PHYSIK
      ! Here comes the main part: Evolve Remo
      !
      CALL EC4ORG(NOZ)
      !
      !================================================================
      !
      ! Close files?
      !
      IF (MYID==0) CALL CLOSEZR
      !
      ! Calls MPI_Wait so all tasks can synchronize
      !
      CALL PWAIT
      !
      ! Runtime
      !
      CALL CPU_TIME(T1)
      TTOTAL = T1-T0
      TMAX   = TTOTAL
      TMIN   = TTOTAL
      TMEAN  = TTOTAL
      !
      COUNT=1
      CALL PALLREDUCER(TMAX,MPI_MAX)
      CALL PALLREDUCER(TMIN,MPI_MIN)
      CALL PALLREDUCER(TMEAN,MPI_SUM)
      TMEAN = TMEAN/NPROC
      !
      IF (MYID==0) THEN
        WRITE(*,*) ''
        WRITE(*,*) '-------------------------------------------------'
        WRITE(*,*) 'Master  Runtime: ', TTOTAL
        WRITE(*,*) 'Maximum Runtime: ', TMAX
        WRITE(*,*) 'Minimum Runtime: ', TMIN
        WRITE(*,*) 'Mean    Runtime: ', TMEAN
        WRITE(*,*) '-------------------------------------------------'
        WRITE(*,*) ''
      ENDIF
      !
      ! Finalize MPI
      !
      CALL MPI_FINALIZE(IERROR)
      !
      !
      END PROGRAM REMORG
      !
