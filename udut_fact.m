%% Upper Triangular Square Root Free Cholesky Decomposition
%
% If P>0 then there exist U, D such that
%     P = U * diag(D) * U'
% where U is Unit Upper Triangular.

function [U D] = udut_fact(P)

    n = size(P,1);
    U = eye(n);
    D = zeros(n,1);

    for j=n:-1:2
        D(j) = P(j,j);
        for k=1:j-1
            U(k,j) = P(k,j) / D(j);
        end
        for k=1:j-1
            for i=1:k
                P(i,k) = P(i,k) - U(i,j) * P(k,j);
            end
        end
    end
    D(1) = P(1,1);
