C     SRIF Data Processing using Householder Tranfsormations
C
C     Observations are merged into the SRIF state by rewriting
C     the problem in an interesting form, then using Householder
C     transformations on a matrix, and drawing the results out
C     of the transformed matrix.
C
C     We arrange our equations in this form:
C
C     .      ||  |  Im |       |  Iv | ||2
C     J(x) = ||  |     | X0 -  |     | ||
C     .      ||  |  H  |       |  z  | ||
C
C
C     where the state of the sytem "x" minimizes J. We want
C     of course to transform this to give us a new Im,Iv to
C     replace the old one, which includes the H,z information.
C
C     Applying an orthogonal transformation T to the vector being
C     minimized will not change its length, and J will be minimized by
C     the same X vector:
C
C     .      ||  |  Im |       |  Iv | ||2
C     J(x) = || T|     | X0 - T|     | ||
C     .      ||  |  H  |       |  z  | ||
C
C     We carefully select the T matrix to triangularize the left hand
C     block, leaving us with an equation in this form:
C
C     .      ||  |  Im |       |  Iv | ||2
C     J(x) = ||  |     | X0 -  |     | ||
C     .      ||  |  0  |       |  e  | ||
C
C     Apologies for re-using Im and Iv here, expanding the names to
C     differentiate between prior estimates and updated estimates makes
C     it pretty unreadable. Note that the new Im is upper triangular;
C     the new Im and Iv are appropriate for use going forward to
C     represent a refined estimate of the state X, and the E value left
C     in the lower part of the right box represents a cumulative amount
C     of known error that can not be reduced out of J by any X.
C
C     Note that we never actually construct T. Our only interest in it
C     is the assertion that it is an orthogonal transform!
C
C     This function operates "in place" on all data items, updating Im
C     and Iv in place, and destroying H and Z.
C
C     @param    N       (in)     dimension of state vector
C     @param    M       (in)     dimension of sense vector
C     @param    Im      (in/out) estimate information matrix
C     @param    Iv      (in/out) estimate information vector
C     @param    H       (in/out) sensor coefficients
C     @param    Z       (in/out) sensor readings
C
C     Runtime is O(N * N * (N+M))

      SUBROUTINE SRIF_DATA(N,M,IM,IV,H,Z)
      IMPLICIT NONE
      INTEGER N,M
      DOUBLE PRECISION IM(N,N),IV(N),H(M,N),Z(M)
C
      INTEGER DIA,ROW,COL
      DOUBLE PRECISION S,BETA,GAMMA,IMD,IMD2,UTOP
C
      DO DIA=1,N
C
         IMD  = IM(DIA,DIA)
         IMD2 = IMD**2
         S = IMD2
         DO ROW=DIA+1,N
            S = S + IM(ROW,DIA)**2
         END DO
         DO ROW=1,M
            S = S + H(ROW,DIA)**2
         END DO
C
         IF (S .GT. IMD2) THEN
C
            S =   SQRT(S)
            IF (IMD .GT. 0.) S = -S
C
            UTOP = IMD - S
            BETA = 1. / (S * UTOP)
C
            IM(DIA,DIA) = S
C
            DO COL=DIA+1,N
C
               GAMMA =         UTOP        *IM(DIA,COL)
               DO ROW=DIA+1,N
                  GAMMA = GAMMA + IM(ROW,DIA)*IM(ROW,COL)
               END DO
               DO ROW=1,M
                  GAMMA = GAMMA + H(ROW,DIA)*H(ROW,COL)
               END DO
               GAMMA = BETA * GAMMA
C
               IM(DIA,COL) = IM(DIA,COL) + GAMMA * UTOP
               DO ROW=DIA+1,N
                  IM(ROW,COL) = IM(ROW,COL) + GAMMA * IM(ROW,DIA)
               END DO
               DO ROW=1,M
                  H(ROW,COL) = H(ROW,COL) + GAMMA * H(ROW,DIA)
               END DO
            END DO
C
            GAMMA =         UTOP        *IV(DIA)
            DO ROW=DIA+1,N
               GAMMA = GAMMA + IM(ROW,DIA)*IV(ROW)
            END DO
            DO ROW=1,M
               GAMMA = GAMMA + H(ROW,DIA)*Z(ROW)
            END DO
            GAMMA = BETA * GAMMA
C
            IV(DIA) = IV(DIA) + GAMMA * UTOP
            DO ROW=DIA+1,N
               IV(ROW) = IV(ROW) + GAMMA * IM(ROW,DIA)
            END DO
            DO ROW=1,M
               Z(ROW) = Z(ROW) + GAMMA * H(ROW,DIA)
            END DO
C
         END IF
C
         DO ROW=DIA+1,N
            IM(ROW,DIA) = 0
         END DO
         DO ROW=1,M
            H(ROW,DIA) = 0
         END DO
C
      END DO
C
      END
