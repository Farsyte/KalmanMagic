%% Modified Weighted Gram-Schmidt Orthogonalization
%
% The input to the algorithm is a collection
% of linearly independent vectors Wj presented
% in the rows of the W matrix, and a collection
% of positive scale factors Dj presented in
% the vector D, which is treated in the math as
% the diagonal of a matrix D.
%
% The result of the algorithm is a collection
% of vectors Vj presented in the rows of the V
% matrix, that are (when weighted by D) orthogonal:
% that is,
%     Vi' D Vj = 0  whenever i != j.
%
% Generate the V vectors by:
%     Vn = Wn
%               k=n   Wj' D Vk
%     Vj = Wj - Sum   -------- Vk
%              k=j+1  Vk' D Vk
%
% This formulation of the algorithm is presented in
% a backward recursive form because we will be wanting
% to use it to generate an upper triangular transform.
%
% Actual implementations may, of course, want to keep
% track of the actual transformation being applied.

function V = mwgso(W,D)

    [n N] = size(W);

    for j=n:-1:1
        Wj = W(j,:);
        Vj = Wj;
        for k=j+1:n
            Vj = Vj - (sum(Wj .* DV(k,:)) / VDV(k)) * V(k,:);
        end
        V(j,:) = Vj;
        DV(j,:) = D' .* Vj;
        VDV(j) = dot(Vj, DV(j,:));
    end
