C     Triangularize a four-part matrix using Householder Transforms
C
C     We are given the six sections of a matrix A, presented as:
C
C     |            |
C     | Aa  Ab  Ac |
C     |            |
C     | Ad  Ae  Af |
C     |            |
C
C     where Aa and Ae are square.
C
C     This function partially triangularizes this matrix by walking down
C     the diagonal of A, selecting a householder transformation that
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
C     The sign of "s" is selected to avoid subtradting two large numbers
C     with a small difference, and more importantly to avoid inverting a
C     very tiny number dominated by numerical errors.
C
C     Runtime is O(D1*(D1+D2+D3)*(D1+D2))

      SUBROUTINE hh_tri_HEX(D1,D2,D3,AA,AB,AC,AD,AE,AF)
      IMPLICIT NONE
      INTEGER D1,D2,D3
      DOUBLE PRECISION AA(D1,D1),AB(D1,D2),AC(D1,D3)
      DOUBLE PRECISION AD(D2,D1),AE(D2,D2),AF(D2,D3)
C
      INTEGER I,J,K
      DOUBLE PRECISION S,BETA,GAMMA
C
      DO I=1,D1
C
         S = 0.
         DO K=I+1,D1
            S = S + AA(K,I)**2
         END DO                                 ! K=I+1,D1
         DO K=1,D2
            S = S + AD(K,I)**2
         END DO                                 ! K=1,D2
C
         IF (S .GT. 0) THEN
C
            S = -DSIGN(SQRT(AA(I,I)**2+S), AA(I,I))
C
            AA(I,I) = AA(I,I) - S
            BETA = 1. / (S * AA(I,I))

            DO J=I+1,D1
C
               GAMMA = 0.
               DO K=I,D1
                  GAMMA = GAMMA + AA(K,I)*AA(K,J)
               END DO                           ! K=I,D1
               DO K=1,D2
                  GAMMA = GAMMA + AD(K,I)*AD(K,J)
               END DO                           ! K=1,D2
C
               GAMMA = BETA * GAMMA
C
               DO K=I,D1
                  AA(K,J) = AA(K,J) + GAMMA * AA(K,I)
               END DO                           ! K=I,D1
               DO K=1,D2
                  AD(K,J) = AD(K,J) + GAMMA * AD(K,I)
               END DO                           ! K=1,D2
            END DO                              ! J=I+1,D1
C
            DO J=1,D2
C
               GAMMA = 0.
               DO K=I,D1
                  GAMMA = GAMMA + AA(K,I)*AB(K,J)
               END DO                           ! K=I,D1
               DO K=1,D2
                  GAMMA = GAMMA + AD(K,I)*AE(K,J)
               END DO                           ! K=1,D2
C
               GAMMA = BETA * GAMMA
C
               DO K=I,D1
                  AB(K,J) = AB(K,J) + GAMMA * AA(K,I)
               END DO                           ! K=I,D1
               DO K=1,D2
                  AE(K,J) = AE(K,J) + GAMMA * AD(K,I)
               END DO                           ! K=1,D2
            END DO                              ! J=1,D2
C
            DO J=1,D3
C
               GAMMA = 0.
               DO K=I,D1
                  GAMMA = GAMMA + AA(K,I)*AC(K,J)
               END DO                           ! K=I,D1
               DO K=1,D2
                  GAMMA = GAMMA + AD(K,I)*AF(K,J)
               END DO                           ! K=1,D2
C
               GAMMA = BETA * GAMMA
C
               DO K=I,D1
                  AC(K,J) = AC(K,J) + GAMMA * AA(K,I)
               END DO                           ! K=I,D1
               DO K=1,D2
                  AF(K,J) = AF(K,J) + GAMMA * AD(K,I)
               END DO                           ! K=1,D2
            END DO                              ! J=1,D3
C
            AA(I,I) = S
            DO K=I+1,D1
               AA(K,I) = 0
            END DO                              ! K=I+1,D1
            DO K=1,D2
               AD(K,I) = 0
            END DO                              ! K=1,D2
         END IF
      END DO                                    ! I=1,D1
C
      DO I=1,D2-1
C
         S = 0.
         DO K=I+1,D2
            S = S + AE(K,I)**2
         END DO                                 ! K=I+1,D2
C
         IF (S .GT. 0) THEN
C
            S = -DSIGN(SQRT(AE(I,I)**2+S), AE(I,I))
C
            AE(I,I) = AE(I,I) - S
            BETA = 1. / (S * AE(I,I))

            DO J=I+1,D2
C
               GAMMA = 0.
               DO K=I,D2
                  GAMMA = GAMMA + AE(K,I)*AE(K,J)
               END DO                           ! K=I,D2
C
               GAMMA = BETA * GAMMA
C
               DO K=I,D2
                  AE(K,J) = AE(K,J) + GAMMA * AE(K,I)
               END DO                           ! K=I,D2
            END DO                              ! J=I+1,D2
C
            DO J=1,D3
C
               GAMMA = 0.
               DO K=I,D2
                  GAMMA = GAMMA + AE(K,I)*AF(K,J)
               END DO                           ! K=I,D2
C
               GAMMA = BETA * GAMMA
C
               DO K=I,D2
                  AF(K,J) = AF(K,J) + GAMMA * AE(K,I)
               END DO                           ! K=I,D2
            END DO                              ! J=1,D3
C
            AE(I,I) = S
            DO K=I+1,D2
               AE(K,I) = 0
            END DO                              ! K=I+1,D2
         END IF
      END DO                                    ! I=1,D2-1
C
      END
