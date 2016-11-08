C     Upper Triangular udut Factorizaton
C
C     P(N,N)    matrix to be factored
C     U(N,N)    where to store U factor
C     D(N)      where to store D factor
C
C     Given a positive definite matrix P,
C     construct the diagonal matrix D
C     and the unit upper triangular matrix U
C     such that P = U D U**T
C
C     The content of P is destroyed.
C
      SUBROUTINE UDUT_FACT(N,P,U,D)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION P(N,N),U(N,N),D(N)
C
      DOUBLE PRECISION ALPHA,BETA
      INTEGER I,J,K
C
      DO J = N, 2, -1
         D(J) = P(J,J)
         ALPHA = 1. / D(J)
         DO K = 1, J-1
            BETA = P(K,J)
            U(K,J) = ALPHA * BETA
            DO I = 1, K
               P(I,K) = P(I,K) - BETA * U(I,J)
            END DO
         END DO
C
         U(J,J) = 1.
         DO K = J+1, N
            U(K,J) = 0.
         END DO
C
      END DO
C
      D(1) = P(1,1)
      U(1,1) = 1.
C
      END
