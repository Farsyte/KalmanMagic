C     Square Root Information Filter Estimate Production
C
C     This function examines the state of the SRIF and
C     provides the current state estimate and covariance;
C     as a byproduct, also provides S (where S*trans(S) = P).
C
C     The system is described by
C     .     Iv = Im x + v
C
C     NOTE: Im is upper triangular,
C     which simplifies this tremendously.
C
      SUBROUTINE SRIF_READ(N,X,P,S,IM,IV)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(N),P(N,N),S(N,N),IM(N,N),IV(N)
C
      INTEGER I,J,K
      DOUBLE PRECISION SUM
C
C     X = Im \ IV
C
      DO J=N,1,-1
         SUM = 0.
         DO K=J+1,N
            SUM = SUM + IM(J,K)*X(K)
         END DO
         X(J) = (IV(J) - SUM) / IM(J,J)
      END DO
C
C     S = inv(Im)
C
      S(1,1) = 1. / IM(1,1)
      DO J = 2, N
         S(J,J) = 1. / IM(J,J)
         DO K = 1, J - 1
            SUM = 0.
            DO I = K, J - 1
               SUM = SUM - S(K,I) * IM(I,J)
            END DO
            S(K,J) = SUM * S(J,J)
         END DO
      END DO
C
C     P = S * S';
C
      DO J=1,N
         DO I=1,J-1
            SUM = 0.
            DO K=J,N
               SUM = SUM + S(I,K) * S(J,K)
            END DO
            P(I,J) = SUM
            P(J,I) = SUM
         END DO
         SUM = 0.
         DO K=J,N
            SUM = SUM + S(J,K)**2
         END DO
         P(J,J) = SUM
      END DO
C
      END
