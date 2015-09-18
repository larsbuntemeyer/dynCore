      SUBROUTINE GETBUF(VAR01, VAR02,
     &     VAR03, VAR04, VAR05, VAR06, VAR07, VAR08, VAR09, VAR10,
     &     REBUF, KSIZE, ILO, IUP, JLO, JUP, IMESLEN, IBUF, LPRES,
     &     IBUFLEN)
!
!     GET ALL VALUES FROM SENDBUF
!
      IMPLICIT NONE
!
      INCLUDE "parorg.h"
!
      INTEGER, INTENT(IN) ::
     &     KSIZE(10)

      REAL, INTENT(OUT) ::
     &     VAR01(IE, JE, KSIZE( 1))

      REAL, OPTIONAL, INTENT(OUT) ::
     &     VAR02(IE, JE, KSIZE( 2)),
     &     VAR03(IE, JE, KSIZE( 3)),
     &     VAR04(IE, JE, KSIZE( 4)),
     &     VAR05(IE, JE, KSIZE( 5)),
     &     VAR06(IE, JE, KSIZE( 6)),
     &     VAR07(IE, JE, KSIZE( 7)),
     &     VAR08(IE, JE, KSIZE( 8)),
     &     VAR09(IE, JE, KSIZE( 9)),
     &     VAR10(IE, JE, KSIZE(10))

      INTEGER, INTENT(IN) ::
     &     IBUFLEN

      REAL, INTENT(INOUT) ::
     &     REBUF(IBUFLEN,6)

      INTEGER, INTENT(IN) ::
     &     ILO, IUP, JLO, JUP

      INTEGER, INTENT(OUT) ::
     &     IMESLEN

      INTEGER, INTENT(IN) ::
     &     IBUF

      LOGICAL, INTENT(IN) ::
     &     LPRES(9)
!
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K

      IMESLEN = 0

      DO K = 1, KSIZE( 1)
        DO J = JLO, JUP
          DO I = ILO, IUP
            IMESLEN = IMESLEN + 1
            VAR01(I,J,K) = REBUF(IMESLEN,IBUF)
          END DO
        END DO
      END DO

      IF (LPRES(1) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 2)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR02(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      END IF

      IF (LPRES(2) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 3)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR03(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(3) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 4)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR04(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(4) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 5)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR05(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(5) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 6)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR06(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(6) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 7)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR07(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(7) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 8)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR08(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(8) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE( 9)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR09(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      IF (LPRES(9) .EQV. .TRUE.) THEN
        DO K = 1, KSIZE(10)
          DO J = JLO, JUP
            DO I = ILO, IUP
              IMESLEN = IMESLEN + 1
              VAR10(I,J,K) = REBUF(IMESLEN,IBUF)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE GETBUF
