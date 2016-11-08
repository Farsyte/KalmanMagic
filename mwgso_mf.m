%% MWG-S Orthogonalization and Matrix Factorization
%
% This code performs a Modified Weighted Gram-Schmidt
% Orthogonalization of the vectors Wj presented in
% the rows of the W matrix, relative to the weights
% presented in the D vector, generating Vj vectors
% into the V matrix result.
%
% The resulting V matrix has the property that
% the matrix product Dbar = V D V' is diagonal, and
% the transformation used to generate V from W
% was unit upper triangular.
%
% This gives us the following relationships:
%
%       W D W' = (Ubar V) D (Ubar V)'
%              =  Ubar V  D  V' Ubar'
%              =  Ubar (V D V') Ubar'
%              =  Ubar   Dbar   Ubar'
%
% This completes the realization of our initial value
% as a (UDU') factored form.

function [Ubar Dbar V] = mwgso_mf(W,D)

    [n N] = size(W);
    V = W;

    % Note: W(n,N) and D(N), where n < N,
    % and we want Ubar(n,n) and Dbar(n).

    Ubar = eye(n);
    Dbar = zeros(n,1);

    for j=n:-1:1

        Vj = V(j,:);
        DVj = D' .* Vj;

        Dbar(j) = sum(Vj .* DVj);

        k = 1:j-1;

        Ubar(k,j) = V(k,:) * DVj' / Dbar(j);

        V(k,:) = V(k,:) - Ubar(k,j) * Vj;
    end
