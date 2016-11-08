%% Agee-Turner PD Factorization Update
%
% Based on the algorithm as presented by Bierman (p45)
%
% Given the U,D factors of a PD matrix P,
%     P = U*D*U'
% calculate the updated U,D factors:
%     U*D*U' = P + c*A*A'
% where "c" is a scalar and "A" is an N-vector.
%
% This algorithm may have some issues when c is negative,
% where rounding errors may drive output D(j) negative,
% but should be just fine for positive c.
%
% Note that delaying the storage of the new D(j) value
% until the end of the loop allows in-place operation.

function [U D] = udut_caat(U,D,c,A)

    n = size(A,1);

    for j=n:-1:2

        nDj = D(j) + c*A(j)^2;

        c = c / nDj;

        k=1:j-1;

        A(k) = A(k) - A(j)*U(k,j);

        U(k,j) = U(k,j) + c*A(j)*A(k);

        c = c * D(j);

        D(j) = nDj;
    end

    D(1) = D(1) + c*A(1)^2;
