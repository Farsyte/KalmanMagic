C     Time Updating of U-D Factors
C     using the MWG-S algorithm.
C
C     See baseline p132 mechanization
C     for additional comments.
C
C     NOTE: this subroutine MUST NOT modify "F"
C
C     NOTE: this subroutine MAY modify "G"
C
      SUBROUTINE UDUT_TIME(N,M,X,U,D,F,G)
      IMPLICIT NONE
      INTEGER N
      INTEGER M
      DOUBLE PRECISION X(N)
      DOUBLE PRECISION U(N,N)
      DOUBLE PRECISION D(N)
      DOUBLE PRECISION F(N,N)
      DOUBLE PRECISION G(N,M)
C
      INTEGER I
      INTEGER J
      INTEGER K
      DOUBLE PRECISION SUM
      DOUBLE PRECISION DINV
      DOUBLE PRECISION SIGMA
C
C     The "f2c" utility, which is used here only to construct the
C     correct function prototype to allow C to call Fortran, does not
C     like it when we declare local arrays whose dimensions are
C     subroutine parameters. Nonetheless, it still generates a useful
C     function prototype.
C
C     Fortunately, G77 accepts this extension, so I do not have to pass
C     these areas in from outside, and can use this form during my
C     prototyping and experimentation state.
C
      DOUBLE PRECISION TX(N)                    ! generic temp N-vector
C
      DOUBLE PRECISION VJ(N)                    ! Vj
      DOUBLE PRECISION GJ(M)                    ! Gj
      DOUBLE PRECISION DVJ(N)                   ! D .* Vj
C
C     Calculate U := F*U
C     using TX to hold U(:,J)
C     while it is being rewritten.
C
      DO J = N, 1, -1
         DO I = 1, J-1
            TX(I) = U(I,J)
         END DO
         DO I = 1, N
            SUM = F(I,J)
            DO K = 1, J - 1
               SUM = SUM + F(I,K) * TX(K)
            END DO
            U(I,J) = SUM
         END DO
      END DO
C
C     Calculate X := F*X
C     using TX to hold the old X.
C
      DO J = 1, N
         TX(J) = X(J)
      END DO
      DO J = 1, N
         SUM = 0.
         DO K = 1, N
            SUM = SUM + F(J,K) * TX(K)
         END DO
         X(J) = SUM
      END DO
C
C
C     Modified Weighted Gram-Schmidt Orthogonalization
C
      DO J = N, 2, -1
C
C     Ujj := Vj D Vj + Gj Q Gj
C
         SUM = 0.
         DO K = 1, N
            VJ(K) = U(J,K)
            DVJ(K) = D(K) * VJ(K)
            SUM = SUM + VJ(K)*DVJ(K)
         END DO
         DO K = 1, M
            GJ(K) = G(J,K)
            SUM = SUM + GJ(K)**2
         END DO
         U(J,J) = SUM
C
         DINV = 1. / SUM                        ! DINV is 1/Dj
         DO K = 1, J - 1
C
C     Ujk := ( Uk D Vj + Gk Q Gj ) / Ujj
C
            SUM = 0.
            DO I = 1, N
               SUM = SUM + U(K,I)*DVJ(I)
            END DO
            DO I = 1, M
               SUM = SUM + G(K,I)*GJ(I)
            END DO
C
            SIGMA = SUM * DINV
            U(J,K) = SIGMA
C
C     Uk := Uk - Vj Ujk
C
            DO I = 1, N
               U(K,I) = U(K,I) - SIGMA*VJ(I)
            END DO
C
C     Gk := Gk - Gj Ujk
C
C     This destroys the G parameter, but skipping
C     this step introduces errors.
C
            DO I = 1, M
               G(K,I) = G(K,I) - SIGMA*GJ(I)
            END DO
         END DO
      END DO
C
C     Ujj := Vj D Vj + Gj Q Gj   (specific to j=1)
C
      SUM = 0.
      DO K = 1, N
         SUM = SUM + U(1,K)**2 * D(K)
      END DO
      DO K = 1, M
         SUM = SUM + G(1,K)**2
      END DO
C
C
C     WAIT! Instead of stuffing SUM back into U(1,1) ...
C
C     Move all results into their final locations:
C     D data is found along the diagonal of the U storage, and
C     U data is found in the lower triangular part.
C     Need to move D data out, copy U data into place,
C     then set U diagonal elements to 1 and elements below
C     the diagonal to 0.
C
C     NOTE: be careful to copy the resuts out of the lower triangular
C     area BEFORE clearing that area. Do not fall into the trap of
C     reordeing the loops to optimally pass over the U memory.
C
      D(1) = SUM
      U(1,1) = 1.
C
      DO J = 2, N
         D(J) = U(J,J)
         DO I = 1, J - 1
            U(I,J) = U(J,I)
            U(J,I) = 0.
         END DO
         U(J,J) = 1.
      END DO
C
      END
