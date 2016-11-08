C     Data Processing usiug the Potter Mechanization
C
C     This function updates the estimate of the state x and the state
C     error covariance square root estimate S based on an observation
C     vector z derived from the state via the observation coefficient
C     matrix H plus a noise vector.
C
C     The observation coefficient vector and observation value must be
C     normalized so that the variance of the observation noise is one.
C
C     PARAMETERS:
C     N  - state vector dimension
C     M  - sensor vector dimension
C     X  - (N   input)  state estimate (in-place update)
C     S  - (NxN input)  error covariance (in-place update)
C     H  - (MxN input)  normalized observation coefficients
C     Z  - (M   input)  normalized noisy observation
C
      SUBROUTINE POTTER_DATA(N,M,X,S,H,Z)
      IMPLICIT NONE
      INTEGER N,M
      DOUBLE PRECISION X(N),S(N,N),H(M,N),Z(M)
C
      INTEGER I,J,K
      DOUBLE PRECISION SIGMA,DELTA,GAMMA,ALPHA
      DOUBLE PRECISION V(N)
C
      DO J = 1, M
C
         SIGMA = 1.
         DELTA = Z(J)
         DO I = 1, N
            V(I) = 0.
            DO K = 1, N
               V(I) = V(I) + S(K,I) * H(J,K)
            END DO
            DELTA = DELTA - H(J,I) * X(I)
            SIGMA = SIGMA + V(I) ** 2
         END DO
C
         SIGMA = 1. / SIGMA
         DELTA = DELTA * SIGMA
         GAMMA = SIGMA / (1. + SQRT(SIGMA))
         DO I = 1, N
            ALPHA = 0.
            DO K = 1, N
               ALPHA = ALPHA + S(I,K) * V(K)
            END DO
            X(I) = X(I) + ALPHA * DELTA
            ALPHA = ALPHA * GAMMA
            DO K = 1, N
               S(I,K) = S(I,K) - ALPHA * V(K)
            END DO
         END DO
C
      END DO
C
      END
