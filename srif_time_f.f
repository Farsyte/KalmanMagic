C     SQUARE ROOT INFORMATION FILTERING: TIME PROCESSING
C
C     We are going to use a trick here, similar to what we just did in
C     the SRIF Data processing above: transform the whole problem into a
C     single matrix, and triangularize the matrix to get the results
C     that we are seeking.
C
C     The equation of interest is as follows:
C
C     .    new X = F X + G W
C
C     where
C     .    F   is the state update transform,
C     .    G   is the process noise mapping matrix,
C     .    W   is the process noise vector from N(0,I)
C
C     We are going to use G to represent our process noise sources: it
C     is merely the mapping of a collection of independent gaussian
C     noise sources with zero mean and unit variance (represented by the
C     W vector in the equation) into the state space. The covariance
C     matrix for the process noise Q is G*G'.
C
C     We could factor Q into GG' but in practice, it is as easy to
C     initially construct G for the system as it is to construct Q (in
C     fact, when I look back at my constructions, I tended to construct
C     Q by first building a matrix that could have been calleed G).
C
C     We have the following equations coming into the procedure:
C
C     .    Iv0 = Im0 X0 +   Vx,   Vx from N(0,I)
C     .     X1 =  F  X0 + G W0
C     and describing "W0" we have:
C     .    Wv0 = Wm0 W0 +   Vw,   Vw from N(0,I)
C
C     ( Note that this actualy gives us a lot more freedom in
C     . what is going on with W0 than we need; in practice,
C     . we can select G so that W0 is from N(0,I) but that
C     . does not actually simplify the math that is coming. )
C
C     Similar to before, we want an (Iv1,Im1) that will determine the
C     unknown X1 state, as always, minimizing the error in a
C     least-squares sense.
C
C     First, express it all in terms of X1. This will of course require
C     us to have on hand Fi, the inverse of F ...
C
C     Solving the second equation for X0 in terms of X1,
C     .     X0 = Fi (X1 - G w)
C     then substituting that back into the first,
C     .    Iv0 = Im0 Fi (X1 - G W0) + v
C
C     Combining the equation describing W0 with
C     the new equation describing X1 and W0,
C
C     .    Wv0 =                    Wm0   W0 + Vw
C     .    Iv0 = (Im0 Fi) X1 - (Im0 Fi G) W0 + Vx
C
C     rewriting this back into our matrix forms and turning it
C     into a performance functional,
C
C     .          ||  |  Wm0         0    |   | W0 |    | Wv0 | ||2
C     J(W0,X1) = ||  |                   | * |    | -  |     | ||
C     .          ||  | -Im0*Fi*G  Im0*Fi |   | X1 |    | Iv0 | ||
C
C     Now hit the stuff inside the squared-two-norm with a nice
C     Householder diagonalization transform ... and we get:
C
C     .          ||  |  Wm0         0    |   | W0 |    | Wv0 | ||2
C     J(W0,X1) = || T|                   | * |    | - T|     | ||
C     .          ||  | -Im0*Fi*G  Im0*Fi |   | X1 |    | Iv0 | ||
C
C     .          ||  |  ^T1        ^T2   |   | W0 |    | ^T3 | ||2
C     .        = ||  |                   | * |    | -  |     | ||
C     .          ||  |   0         Im1   |   | X1 |    | Iv1 | ||
C
C     which separates into minimizing
C     .    ^T1 W0 + ^T2 X1 - ^T3
C     and
C     .             Im1 X1 - Iv1
C
C     Note that this last line is precisely the form we want. All we
C     lack is the proof that J(W0,X1) is minimized by X1 when the last
C     line is minimized; Bierman goes through the magic that proves this
C     result. Essentially, ^T1 is invertable diagonal, so the top part
C     can always minimize to zero when W0 = (^T2 X1 - ^T3) / ^T1 which
C     means that the value of X1 is determined only by the last line.
C
C     So once again, we marshal our data and apply Householder to the
C     result, giving us a nice upper triangular Im1 matrix and an Iv1
C     vector that we can use as the filter state going forward.
C
      SUBROUTINE SRIF_TIME(N,M,IM,IV,FI,G,WM,WV)
      IMPLICIT NONE
      INTEGER N,M
      DOUBLE PRECISION IM(N,N),IV(N),FI(N,N),G(N,M),WM(N,N),WV(M)
C
      INTEGER I,J,K,L
      DOUBLE PRECISION WX(M,N),IW(N,M)
      DOUBLE PRECISION SUM,TN(N)
      DOUBLE PRECISION OLDD,OLDD2,S2,S,TOPU
      DOUBLE PRECISION BETA,GAMMA
C
C     Im = Im*Fi
C
      DO I = 1, N
         DO J = 1, N
            SUM = 0.
            DO K = 1, N
               SUM = SUM + IM(I,K) * FI(K,J)
            END DO
            TN(J) = SUM
         END DO
         DO J = 1, N
            IM(I,J) = TN(J)
         END DO
      END DO
C
C     iw = -Im*G
C
      DO J = 1, M
         DO I = 1, N
            SUM = 0.
            DO K = 1, N
               SUM = SUM + IM(I,K) * G(K,J)
            END DO
            IW(I,J) = SUM
         END DO
      END DO
C
C     wx = zeros(m,n)
C
      DO J = 1, N
         DO I = 1, M
            WX(I,J) = 0.
         END DO
      END DO
C
C     Matrix to be triangularized is
C
C     .     | Wm  wx  Wv |  (m)
C     .     |            |
C     .     | iw  Im  Iv |  (n)
C     .       (m) (n) (1)
C
C     "wx" is zero when we start.
C     "iw" is zero when we finish.
C
C     This function partially triangularizes this matrix by walking down
C     the diagonal of Wm then Im, selecting a householder transformation
C     that clears elements below that diagonal element, and applying
C     that transformation to the remainder of the matrix. This will
C     result in Wm triangular, iw zero, and Im triangular.
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
C     Runtime is O(m*(m+n+1)*(m+n))
C
      DO L = 1, M
C
         OLDD = WM(L,L)
         OLDD2 = OLDD ** 2
C
         S2 = OLDD2
         DO I = L + 1, M
            S2 = S2 + WM(I,L) ** 2
         END DO
         DO I = 1, N
            S2 = S2 + IW(I,L) ** 2
         END DO
C
         IF (S2 > OLDD2) THEN
            S = SQRT(S2)
            IF (OLDD > 0) S = -S
C
            TOPU = OLDD - S
            BETA = 1 / (S * TOPU)
C
            DO J = L + 1, M
C
               GAMMA = TOPU * WM(L,J)
               DO I = L + 1, M
                  GAMMA = GAMMA + WM(I,L) * WM(I,J)
               END DO
               DO I = 1, N
                  GAMMA = GAMMA + IW(I,L) * IW(I,J)
               END DO
               GAMMA = GAMMA * BETA
C
               WM(L,J) = WM(L,J) + GAMMA * TOPU
               DO I = L + 1, M
                  WM(I,J) = WM(I,J) + GAMMA * WM(I,L)
               END DO
               DO I = 1, N
                  IW(I,J) = IW(I,J) + GAMMA * IW(I,L)
               END DO
            END DO
C
            DO J = 1, N
               GAMMA = TOPU * WX(L,J)
               DO I = L + 1, M
                  GAMMA = GAMMA + WM(I,L) * WX(I,J)
               END DO
               DO I = 1, N
                  GAMMA = GAMMA + IW(I,L) * IM(I,J)
               END DO
               GAMMA = GAMMA * BETA
C
               WX(L,J) = WX(L,J) + GAMMA * TOPU
               DO I = L + 1, M
                  WX(I,J) = WX(I,J) + GAMMA * WM(I,L)
               END DO
               DO I = 1, N
                  IM(I,J) = IM(I,J) + GAMMA * IW(I,L)
               END DO
            END DO
C
            GAMMA = TOPU * WV(L)
            DO I = L + 1, M
               GAMMA = GAMMA + WM(I,L) * WV(I)
            END DO
            DO I = 1, N
               GAMMA = GAMMA + IW(I,L) * IV(I)
            END DO
            GAMMA = GAMMA * BETA
C
            WV(L) = WV(L) + GAMMA * TOPU
            DO I = L + 1, M
               WV(I) = WV(I) + GAMMA * WM(I,L)
            END DO
            DO I = 1, N
               IV(I) = IV(I) + GAMMA * IW(I,L)
            END DO
C
            WM(L,L) = S
         END IF
         DO I = L + 1, M
            WM(I,L) = 0
         END DO
         DO I = 1, N
            IW(I,L) = 0
         END DO
      END DO
C
      DO L = 1, N - 1
C
         OLDD = IM(L,L)
         OLDD2 = OLDD ** 2
         S2 = OLDD2
         DO I = L + 1, M
            S2 = S2 + IM(I,L) ** 2
         END DO
         IF (S2 > OLDD2) then
            S = SQRT(S2)
            IF (OLDD > 0) S = -S
C
            TOPU = OLDD - S
            BETA = 1 / (S * TOPU)
C
            DO J = L + 1, N
               GAMMA = TOPU * IM(L,J)
               DO I = L + 1, M
                  GAMMA = GAMMA + IM(I,L) * IM(I,J)
               END DO
               GAMMA = GAMMA * BETA
C
               IM(L,J) = IM(L,J) + GAMMA * TOPU
               DO I = L + 1, M
                  IM(I,J) = IM(I,J) + GAMMA * IM(I,L)
               END DO
            END DO
C
            GAMMA = TOPU * IV(L)
            DO I = L + 1, M
               GAMMA = GAMMA + IM(I,L) * IV(I)
            END DO
            GAMMA = GAMMA * BETA
C
            IV(L) = IV(L) + GAMMA * TOPU
            DO I = L + 1, M
               IV(I) = IV(I) + GAMMA * IM(I,L)
            END DO
C
            IM(L,L) = S
         END IF
         DO I = L + 1, M
            IM(I,L) = 0
         END DO
      END DO
C
      END
