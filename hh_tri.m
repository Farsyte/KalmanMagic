%% Triangularize (part of) a matrix using Householder Transforms
%
% For each of the first "d" diagonal elements of A, a Householder
% transform that produces zeros in all entries below that element is
% applied to the entire submatrix from that element to the lower right
% corner of the matrix, inclusive.
%
% The Householder transform can be viewed geometrically as producing
% the mirror image of a vector V using a plane that is defined by the
% transforming vector U.
%
% At each iteration, we use the transform that examines another column
% of the matrix, and transforms it so that all elements below the
% diagonal are zero.
%
% If all elements below the diagonal are already zero, then no
% transform is required, of couse.
%
% The sign of "s" is selected to avoid subtracting two large numbers
% with a small difference, and more importantly to avoid inverting a
% very tiny number dominated by numerical errors.
%
% Note that this function never explicitly stores the T matrix and the
% U vector is overwritten by result data. Frequently the actual T
% matrix is not needed. If the cumulative T matrix is required, append
% an identity matrix to the right side of the input; the T matrix will
% appear in that area.
%
% Runtime is O(dmn), with O(d) square roots and O(d) scalar
% inversions.

function A = hh_tri(A,d)

    [ m n ] = size(A);

    % Walk down the first "d" diagonal elements.
    for l=1:d

        s = sqrt(n2( A(l:m,l) ));

        if (A(l,l) > 0) s = -s; end;

        A(l,l) = A(l,l) - s;
        beta = 1 / (s * A(l,l));

        for k=l+1:n
            gamma = beta * ( dot(A(l:m,l), A(l:m,k)) );
            A(l:m,k) = A(l:m,k) + gamma * A(l:m,l);
        end

        A(l,l) = s;
        A(l+1:end,l) = 0;
    end

