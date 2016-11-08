%% Upper Triangular Factorization: UU' = P
%
% If P>0 it has an upper triangular factorization U, UU'=P.
% The following procedure constructs the factorization
% with positive diagonal elements.
%
% NOTE: this is not the Cholesky factorization,
% which would be
%       U = chol(P) giving U'U = P

function U = uut_fact(P)

    n = size(P,1);

    U = zeros(n,n);

    for j=n:-1:2
        if (P(j,j) > 0)
            U(j,j) = sqrt(P(j,j));
            for k=1:j-1
                U(k,j) = P(k,j) / U(j,j);
            end
            for k=1:j-1
                for i=1:k
                    P(i,k) = P(i,k) - U(i,j) * U(k,j);
                end
            end
        end
    end
    if (P(1,1) > 0)
        U(1,1) = sqrt(P(1,1));
    end
