%% Triangularize a six-part matrix using Householder Transforms
%
% We are given the six quadrants of a matrix A, presented as:
%
%       |             |
%       | Aa  Ab  Ac  |
%       |             |
%       | Ad  Ae  Af  |
%       |             |
%
% This function partially triangularizes this matrix by walking
% down the diagonal of A, selecting a householder transformation
% that clears elements below that diagonal element, and applying
% that transformation to the remainder of the matrix. Assuming that
% Aa and Ae are square, this will give results in Aa and Ae that
% are upper triangular, and will leave Ad entirely zero.
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
% Runtime is O(d1*(d1+d2+d3)*(d1+d2))

function [Aa Ab Ac Ad Ae Af] = hh_tri_hex(Aa,Ab,Ac,Ad,Ae,Af)

    d1 = max([size(Aa,1) size(Ab,1) size(Ac,1) size(Aa,2) size(Ad,2) ]);
    d2 = max([size(Ad,1) size(Ae,1) size(Af,1) size(Ab,2) size(Ae,2) ]);
    d3 = max(size(Ac,2), size(Af,2));

    for l=1:d1

        % grab the sum of the squares of the elements
        % below the diagonal. If this is zero, then we
        % can skip this column entirely.

        s = n2(Aa(l+1:d1,l)) + n2(Ad(1:d2,l));

        if (s > 0)
            s = sqrt( Aa(l,l)^2 + s );
            if (Aa(l,l) > 0) s = -s; end;

            Aa(l,l) = Aa(l,l) - s;
            beta = 1 / (s * Aa(l,l));

            for c1=l+1:d1
                gamma = beta * ( dot( Aa(l:d1,l), Aa(l:d1,c1) ) +    ...
                                 dot( Ad(1:d2,l), Ad(1:d2,c1) ) );
                Aa(l:d1,c1) = Aa(l:d1,c1) + gamma * Aa(l:d1,l);
                Ad(1:d2,c1) = Ad(1:d2,c1) + gamma * Ad(1:d2,l);
            end

            for c2=1:d2
                gamma = beta * ( dot( Aa(l:d1,l),Ab(l:d1,c2) ) +     ...
                                 dot( Ad(1:d2,l),Ae(1:d2,c2) ) );
                Ab(l:d1,c2) = Ab(l:d1,c2) + gamma * Aa(l:d1,l);
                Ae(1:d2,c2) = Ae(1:d2,c2) + gamma * Ad(1:d2,l);
            end

            for c3=1:d3
                gamma = beta * ( dot( Aa(l:d1,l),Ac(l:d1,c3) ) +     ...
                                 dot( Ad(1:d2,l),Af(1:d2,c3) ) );
                Ac(l:d1,c3) = Ac(l:d1,c3) + gamma * Aa(l:d1,l);
                Af(1:d2,c3) = Af(1:d2,c3) + gamma * Ad(1:d2,l);
            end

            Aa(l,l) = s;
            Aa(l+1:d1,l) = 0;
            Ad(1:d2,l) = 0;
        end
    end

    for l=1:d2-1

        % grab the sum of the squares of the elements
        % below the diagonal. If this is zero, then we
        % can skip this column entirely.
        %
        % XXX: comparing S with Ae(l,l)^2*eps might be
        % better than comparing it to zero, so we do
        % not have to perform a transform using a U that
        % is almost exactly along the Lth coordinate.

        s = n2(Ae(l+1:d2,l));
        if (s > 0)
            s = sqrt( Ae(l,l)^2 + s );
            if (Ae(l,l) > 0) s = -s; end;

            Ae(l,l) = Ae(l,l) - s;
            beta = 1 / (s * Ae(l,l));

            for c2=l+1:d2
                gamma = beta * ( dot( Ae(l:d2,l), Ae(l:d2,c2) ) );
                Ae(l:d2,c2) = Ae(l:d2,c2) + gamma * Ae(l:d2,l);
            end

            for c3=1:d3
                gamma = beta * ( dot( Ae(l:d2,l), Af(l:d2,c3) ) );
                Af(l:d2,c3) = Af(l:d2,c3) + gamma * Ae(l:d2,l);
            end

            Ae(l,l) = s;
            Ae(l+1:d2,l) = 0;

            % If the comparison was not against zero,
            % then we will need an "else" clause here
            % to take appropriate action to manage
            % the approximation.
        end
    end
