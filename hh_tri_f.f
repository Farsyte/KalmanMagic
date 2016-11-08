C     Triangularize (part of) a matrix using Householder Transforms
C
C     For each of the first "d" diagonal elements of A, a Householder
C     transform that produces zeros in all entries below that element is
C     applied to the entire submatrix from that element to the lower
C     right corner of the matrix, inclusive.
C
C     The Householder transform can be viewed geometrically as producing
C     the mirror image of a vector V using a plane that is defined by
C     the transforming vector U.
C
C     At each iteration, we use the transform that examines another
C     column of the matrix, and transforms it so that all elements below
C     the diagonal are zero.
C
C     If all elements below the diagonal are already zero, then no
C     transform is required, of couse.
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
C     Runtime is O(dmn), with O(d) square roots and O(d) inversions.

      SUBROUTINE HH_TRI(M,N,A,D)
      IMPLICIT NONE
      INTEGER M,N,D
      DOUBLE PRECISION A(*)
C
      INTEGER I,K,L,MK,ML
      DOUBLE PRECISION S,BETA,GAMMA
C
      ML = 0
      DO L=1, D
         S = 0.
         DO I = L, M
            S = S + A(I+ML) ** 2
         END DO
         S = SQRT(S)
C
         IF (A(L+ML) .GT. 0) S = -S
C
         A(L+ML) = A(L+ML) - S
         BETA = 1. / (S * A(L+ML))
C
         MK = ML + M
         DO K = L+1, N
            GAMMA = 0.
            DO I = L, M
               GAMMA = GAMMA + A(I+ML) * A(I+MK)
            END DO
            GAMMA = GAMMA * BETA
            DO I = L, M
               A(I+MK) = A(I+MK) + GAMMA * A(I+ML)
            END DO
            MK = MK + M
         END DO
C
         A(L+ML) = S
         DO I = L + 1, M
            A(I+ML) = 0.
         END DO
         ML = ML + M
      END DO
C
      END
