C     Triangularize a four-part matrix using Householder Transforms
C
C     We are given the four quadrants of a matrix A, presented as:
C
C     |        |
C     | Aa  Ab |
C     |        |
C     | Ac  Ad |
C     |        |
C
C     Aa is a DxD matrix.
C     Ab is a DxN matrix.
C     Ac is a MxD matrix.
C     Ad is a MxN matrix.
C
C     This function partially triangularizes this matrix by walking down
C     the diagonal of Aa, selecting a householder transformation that
C     clears elements below that diagonal element, and applying that
C     transformation to the remainder of the matrix.
C
C     The Householder transform can be viewed geometrically as producing
C     the mirror image of a vector Y using a plane that is defined by
C     the transforming vector U.
C
C     At each iteration, we use the transform that takes the column
C     vector from the selected element to the bottom of the matrix, and
C     transforms it so that all elements but the first are zero.
C
C     The sign of "s" is selected to avoid subtracting two large numbers
C     with a small difference, and more importantly to avoid inverting a
C     very tiny number dominated by numerical errors.
C
C     Note that this function never explicitly stores the T matrix and
C     the U vector is overwritten by result data. Frequently the actual
C     T matrix is not needed. If the cumulative T matrix is required,
C     append an identity matrix to the right side of the input; the T
C     matrix will appear in that area.
C
C     Runtime is O( D*(D+M)*(D+N) )

      SUBROUTINE hh_tri_QUAD(D,M,N,AA,AB,AC,AD)
      IMPLICIT NONE
      INTEGER D,M,N
      DOUBLE PRECISION AA(D,D),AB(D,N),AC(M,D),AD(M,N)
C
      INTEGER L,I,K
      DOUBLE PRECISION S,BETA,GAMMA
C
      DO L=1,D
C
         S = 0.
         DO I=L,D
            S = S + AA(I,L)**2
         END DO                                 ! I=L,D
         DO I=1,M
            S = S + AC(I,L)**2
         END DO                                 ! I=1,M
         S = SQRT(S)
C
         IF (AA(L,L) .GT. 0) S = -S
C
         AA(L,L) = AA(L,L) - S
         BETA = 1. / (S * AA(L,L))

         DO K=L+1,D
C
            GAMMA = 0.
            DO I=L,D
               GAMMA = GAMMA + AA(I,L)*AA(I,K)
            END DO                              ! I=L,D
            DO I=1,M
               GAMMA = GAMMA + AC(I,L)*AC(I,K)
            END DO                              ! I=1,M
            GAMMA = BETA * GAMMA
C
            DO I=L,D
               AA(I,K) = AA(I,K) + GAMMA * AA(I,L)
            END DO                              ! I=L,D
            DO I=1,M
               AC(I,K) = AC(I,K) + GAMMA * AC(I,L)
            END DO                              ! I=1,M
C
            CONTINUE
         END DO                                 ! K=L+1,D
C
         DO K=1,N
C
            GAMMA = 0.
            DO I=L,D
               GAMMA = GAMMA + AA(I,L)*AB(I,K)
            END DO                              ! I=L,D
            DO I=1,M
               GAMMA = GAMMA + AC(I,L)*AD(I,K)
            END DO                              ! I=1,M
            GAMMA = BETA * GAMMA
C
            DO I=L,D
               AB(I,K) = AB(I,K) + GAMMA * AA(I,L)
            END DO                              ! I=L,D
            DO I=1,M
               AD(I,K) = AD(I,K) + GAMMA * AC(I,L)
            END DO                              ! I=1,M
C
            CONTINUE
         END DO                                 ! K=1,N
C
         AA(L,L) = S
         DO I=L+1,D
            AA(I,L) = 0
         END DO                                 ! I=L+1,D
         DO I=1,M
            AC(I,L) = 0
         END DO                                 ! I=1,M
C
         CONTINUE
      END DO                                    ! L=1,D
C
      END

