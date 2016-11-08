%% Square Root Information Filter Estimate Production
%
% This function examines the state of the SRIF and
% provides the current state estimate and covariance.
%
% The system is described by
%       Im x = Iv - v
%
% NOTE: Im is upper triangular.
%
% Note the use here of "matrix predivision"
% which can give more useful results than
% attempting to work with inv(S).
%
% The "uut_fact" routine turns out an upper
% triangular matrix.

function [Im Iv] = srif_init(x,P)

    Nx = size(x,1);
    S = uut_fact(P);
    Im = S \ eye(Nx);
    Iv = S \ x;

    Im(isnan(Im)) = 0;
    Im(isinf(Im)) = 0;

    Iv(isnan(Iv)) = 0;
    Iv(isinf(Iv)) = 0;

