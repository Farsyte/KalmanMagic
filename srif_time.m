%% SRIF Time Processing [REALLY SLOW IN MATLAB]
%
% We are going to use a trick here, similar to what we just did
% in the SRIF Data processing above: transform the whole problem
% into a single matrix, and triangularize the matrix to get the
% results that we are seeking.
%
% The equation of interest is as follows:
%
%    new X = F X + G W
%
% where
%     F   is the state update transform,
%     G   is the process noise mapping matrix,
%     W   is the process noise vector from N(0,I)
%
% We are going to use G to represent our process noise sources: it is
% merely the mapping of a collection of independent gaussian noise
% sources with zero mean and unit variance (represented by the W vector
% in the equation) into the state space. The covariance matrix for the
% process noise Q is G*G'.
%
% We could factor Q into GG' but in practice, it is as easy to initially
% construct G for the system as it is to construct Q (in fact, when I
% look back at my constructions, I tended to construct Q by first
% building a matrix that could have been calleed G).
%
% We have the following equations coming into the procedure:
%
%   Iv0 = Im0 X0 +   Vx,   Vx from N(0,I)
%    X1 =  F  X0 + G W0
% and describing "W0" we have:
%   Wv0 = Wm0 W0 +   Vw,   Vw from N(0,I)
%
% ( Note that this actualy gives us a lot more freedom in
%   what is going on with W0 than we need; in practice,
%   we can select G so that W0 is from N(0,I) but that
%   does not actually simplify the math that is coming. )
%
% Similar to before, we want an (Iv1,Im1) that will
% determine the unknown X1 state, as always, minimizing
% the error in a least-squares sense.
%
% First, express it all in terms of X1. This will of course
% require us to have on hand Fi, the inverse of F.
%
% Solving the second equation for X0 in terms of X1,
%     X0 = Fi (X1 - G w)
% then substituting that back into the first,
%   Iv0 = Im0 Fi (X1 - G W0) + v
%
% Combining the equation describing W0 with
% the new equation describing X1 and W0,
%
%   Wv0 =                    Wm0   W0 + Vw
%   Iv0 = (Im0 Fi) X1 - (Im0 Fi G) W0 + Vx
%
% rewriting this back into our matrix forms and turning it
% into a performance functional,
%
%            ||  |  Wm0         0    |   | W0 |    | Wv0 | ||2
% J(W0,X1) = ||  |                   | * |    | -  |     | ||
%            ||  | -Im0*Fi*G  Im0*Fi |   | X1 |    | Iv0 | ||
%
% Now hit the stuff inside the squared-two-norm with a nice
% Householder diagonalization transform ... and we get:
%
%            ||  |  Wm0         0    |   | W0 |    | Wv0 | ||2
% J(W0,X1) = || T|                   | * |    | - T|     | ||
%            ||  | -Im0*Fi*G  Im0*Fi |   | X1 |    | Iv0 | ||
%
%            ||  |  ^T1        ^T2   |   | W0 |    | ^T3 | ||2
%          = ||  |                   | * |    | -  |     | ||
%            ||  |   0         Im1   |   | X1 |    | Iv1 | ||
%
% which separates into minimizing
%   ^T1 W0 + ^T2 X1 - ^T3
% and
%            Im1 X1 - Iv1
%
% Note that this last line is precisely the form we want. All we lack is
% the proof that J(W0,X1) is minimized by X1 when the last line is
% minimized; Bierman goes through the magic that proves this
% result. Essentially, ^T1 is invertable diagonal, so the top part can
% always minimize to zero when W0 = (^T2 X1 - ^T3) / ^T1 which means
% that the value of X1 is determined only by the last line.
%
% So once again, we marshal our data and apply Householder to the
% result, giving us a nice upper triangular Im1 matrix and an Iv1
% vector that we can use as the filter state going forward.

function [Im Iv] = srif_time(Im,Iv,Fi,G,Wm,Wv)

    [n m] = size(G);

    Im =  Im*Fi;
    iw = -Im*G;
    wx = zeros(m,n);

    % Matrix to be triangularized is
    %
    %       | Wm  wx  Wv |  (m)
    %       |            |
    %       | iw  Im  Iv |  (n)
    %         (m) (n) (1)
    %
    % "wx" is zero when we start.
    % "iw" is zero when we finish.

    % This function partially triangularizes this matrix by walking down
    % the diagonal of Wm then Im, selecting a householder
    % transformation that clears elements below that diagonal element,
    % and applying that transformation to the remainder of the
    % matrix. This will result in Wm triangular, iw zero, and Im
    % triangular.
    %
    % The Householder transform can be viewed geometrically as
    % producing the mirror image of a vector Y using a plane that is
    % defined by the transforming vector U.
    %
    % At each iteration, we use the transform that takes the column
    % vector from the selected element to the bottom of the matrix,
    % and transforms it so that all elements but the first are zero.
    %
    % The sign of "s" is selected to avoid subtracting two large
    % numbers with a small difference, and more importantly to avoid
    % inverting a very tiny number dominated by numerical errors.
    %
    % Runtime is O(m*(m+n+1)*(m+n))

    for l=1:m

        oldD = Wm(l,l);
        oldD2 = oldD ^ 2;
        s2 = oldD2 + n2(Wm(l+1:m,l)) + n2(iw(1:n,l));
        if (s2 > oldD2)
            s = sqrt(s2);
            if (oldD > 0) s = -s; end;

            topU = oldD - s;
            beta = 1 / (s * topU);

            for c1=l+1:m
                gamma =              topU       * Wm(l    ,c1)   ;
                gamma = gamma + dot( Wm(l+1:m,l), Wm(l+1:m,c1) ) ;
                gamma = gamma + dot( iw(  1:n,l), iw(  1:n,c1) ) ;
                gamma = gamma * beta;

                Wm(l    ,c1) = Wm(l    ,c1) + gamma * topU        ;
                Wm(l+1:m,c1) = Wm(l+1:m,c1) + gamma * Wm(l+1:m,l) ;
                iw(  1:n,c1) = iw(  1:n,c1) + gamma * iw(  1:n,l) ;
            end

            for c2=1:n
                gamma =              topU       * wx(l    ,c2)   ;
                gamma = gamma + dot( Wm(l+1:m,l), wx(l+1:m,c2) ) ;
                gamma = gamma + dot( iw(  1:n,l), Im(  1:n,c2) ) ;
                gamma = gamma * beta;

                wx(l    ,c2) = wx(l    ,c2) + gamma * topU        ;
                wx(l+1:m,c2) = wx(l+1:m,c2) + gamma * Wm(l+1:m,l) ;
                Im(  1:n,c2) = Im(  1:n,c2) + gamma * iw(  1:n,l) ;
            end

            gamma =              topU       * Wv(l    , 1)   ;
            gamma = gamma + dot( Wm(l+1:m,l), Wv(l+1:m, 1) ) ;
            gamma = gamma + dot( iw(  1:n,l), Iv(  1:n, 1) ) ;
            gamma = gamma * beta;

            Wv(l    , 1) = Wv(l    , 1) + gamma * topU        ;
            Wv(l+1:m, 1) = Wv(l+1:m, 1) + gamma * Wm(l+1:m,l) ;
            Iv(  1:n, 1) = Iv(  1:n, 1) + gamma * iw(  1:n,l) ;

            Wm(l    ,l) = s;
        end
        Wm(l+1:m,l) = 0;
        iw(  1:n,l) = 0;
    end

    for l=1:n-1

        oldD = Im(l,l);
        oldD2 = oldD ^ 2;
        s2 = oldD2 + n2(Im(l+1:m,l)) ;
        if (s2 > oldD2)
            s = sqrt(s2);
            if (oldD > 0) s = -s; end;

            topU = oldD - s;
            beta = 1 / (s * topU);

            for c1=l+1:n
                gamma =              topU       * Im(l    ,c1)   ;
                gamma = gamma + dot( Im(l+1:m,l), Im(l+1:m,c1) ) ;
                gamma = gamma * beta;

                Im(l    ,c1) = Im(l    ,c1) + gamma * topU        ;
                Im(l+1:m,c1) = Im(l+1:m,c1) + gamma * Im(l+1:m,l) ;
            end

            gamma =              topU       * Iv(l    , 1)   ;
            gamma = gamma + dot( Im(l+1:m,l), Iv(l+1:m, 1) ) ;
            gamma = gamma * beta;

            Iv(l    , 1) = Iv(l    , 1) + gamma * topU        ;
            Iv(l+1:m, 1) = Iv(l+1:m, 1) + gamma * Im(l+1:m,l) ;

            Im(l    ,l) = s;
        end
        Im(l+1:m,l) = 0;
    end
