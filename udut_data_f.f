C     U-D Update
C
C     Inputs:
C     X(n,1)    Prior Estimate
C     U(n,n)    Prior Covariance, U factor
C     D(n,1)    Prior Covariance, D factor
C     Z(m,1)    Observation
C     H(m,n)    Observation Coefficients
C     R(m)      Observation Error Variances
C
C       Note: this procedure requires
C       that the R matrix is diagonal.
C
C     Outputs:
C     X(n,1)    Updated Estimate
C     U(n,n)    Upated Covariance, U factor
C     D(n,1)    Upated Covariance, D factor
C
C     The covariance is P = U D U'
C     where D is diagonal, and U is unit upper triangular.
C
C     Elements on and below the diagonal of U are neither
C     referenced nor updated.
C
      SUBROUTINE UDUT_DATA(N,M,X,U,D,H,R,Z)
      IMPLICIT NONE
      INTEGER N,M
      DOUBLE PRECISION X(N),U(N,N),D(N),H(M,N),R(M),Z(M)
C
      INTEGER I,J,K,L
      DOUBLE PRECISION ALPHA,BETA,GAMMA,LAMDA
      DOUBLE PRECISION B(N)
C
      DO J = 1, M
C
         DO K = 1, N
            Z(J) = Z(J) - H(J,K) * X(K)
         END DO
C
         DO K = N, 2, -1
            DO L = 1, K-1
               H(J,K) = H(J,K) + U(L,K)*H(J,L)
            END DO
            B(K) = D(K) * H(J,K)
         END DO
         B(1) = D(1) * H(J,1)
C
C     Z(J) := Z(J) - H'X (the computed residual),
C     B = D U' H, and H := U'H have all been computed.
C
         ALPHA = R(J) + B(1)*H(J,1)
         GAMMA = 1. / ALPHA
         D(1) = R(J) * GAMMA * D(1)
C
         DO K = 2, N
            BETA = ALPHA
            ALPHA = ALPHA + B(K)*H(J,K)
            LAMDA = -H(J,K) * GAMMA
            GAMMA = 1. / ALPHA
            D(K) = BETA * GAMMA * D(K)
            DO I = 1, K-1
               BETA = U(I,K)
               U(I,K) = BETA + B(I) * LAMDA
               B(I) = B(I) + B(K) * BETA
            END DO
         END DO
C
         Z(J) = Z(J) * GAMMA
         DO K = 1, N
            X(K) = X(K) + B(K) * Z(J)
         END DO
C
      END DO
C
      END
