%% Generate a random P-D matrix of specified dimension.
%
% Note that this routine MIGHT not produce
% a positive definite matrix.

function P = randpd(d)
    S = randn(d);
    P = S * S';
