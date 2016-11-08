%% Triangularize a four-part matrix using Householder Transforms
%
% We are given the four quadrants of a matrix A, presented as:
%
%       |  Aa    Ab | (Mt)
%       |           |
%       |  Ac    Ad | (Mb)
%         (Nl)  (Nr)
%
% This function partially triangularizes this matrix by walking
% down the diagonal of Aa, selecting a householder transformation
% that clears elements below that diagonal element, and applying
% that transformation to the remainder of the matrix.
%
% The Householder transform can be viewed geometrically as producing
% the mirror image of a vector Y using a plane that is defined by the
% transforming vector U.
%
% At each iteration, we use the transform that takes the column vector
% from the selected element to the bottom of the matrix, and
% transforms it so that all elements but the first are zero.
%
% The sign of "s" is selected to avoid subtracting two large numbers
% with a small difference, and more importantly to avoid inverting a
% very tiny number dominated by numerical errors.
%
% Runtime is O( D*(D+M)*(D+N) )

function [Aa Ab Ac Ad] = hh_tri_quad(Aa, Ab, Ac, Ad)

    Mt = max([size(Aa,1) size(Ab,1)]);
    Mb = max([size(Ac,1) size(Ad,1)]);

    Nl = max(size(Aa,2), size(Ac,2));
    Nr = max(size(Ab,2), size(Ad,2));

    d = min(Mt, Nl);

    % m = max(size(Ac,1), size(Ad,1));
    % n = max(size(Ab,2), size(Ad,2));

    for l=1:d

        % grab the sum of the squares of the elements
        % below the diagonal. If this is zero, then we
        % can skip this column entirely.

        s = n2(Aa(l+1:Mt,l)) + n2(Ac(1:Mb,l));

        if (s > 0)
            s = sqrt(Aa(l,l)^2 + s);
            if (Aa(l,l) > 0) s = -s; end;

            Aa(l,l) = Aa(l,l) - s;
            beta = 1 / (s * Aa(l,l));

            for k=l+1:Nl
                gamma = beta * ( dot( Aa(l:Mt,l), Aa(l:Mt,k) ) +     ...
                                 dot( Ac(1:Mb,l), Ac(1:Mb,k) ) );
                Aa(l:Mt,k) = Aa(l:Mt,k) + gamma * Aa(l:Mt,l);
                Ac(1:Mb,k) = Ac(1:Mb,k) + gamma * Ac(1:Mb,l);
            end

            for k=1:Nr
                gamma = beta * ( dot( Aa(l:Mt,l),Ab(l:Mt,k) ) +      ...
                                 dot( Ac(1:Mb,l),Ad(1:Mb,k) ) );
                Ab(l:Mt,k) = Ab(l:Mt,k) + gamma * Aa(l:Mt,l);
                Ad(1:Mb,k) = Ad(1:Mb,k) + gamma * Ac(1:Mb,l);
            end

            Aa(l,l) = s;
            Aa(l+1:d,l) = 0;
            Ac(1:Mb,l) = 0;
        end
    end
