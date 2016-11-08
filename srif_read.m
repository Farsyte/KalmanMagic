%% Square Root Information Filter Estimate Production
%
% This function examines the state of the SRIF and
% provides the current state estimate and covariance;
% as a byproduct, also provides S (where S S' = P).
%
% The system is described by
%       Iv = Im x + v
%
% NOTE: Im is upper triangular.

function [Xest Pest Sest] = srif_read(Im,Iv)

    Nx = size(Iv,1);

    Xest = Im \ Iv;
    Sest = Im \ eye(Nx);
    Pest = Sest * Sest';
