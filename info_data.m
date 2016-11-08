%% Data Processing using the Information Matrix
%
% This function merges a measurement and its noise into the
% information matrix estimating the state of the observed system
% using the Im matrix and Iv vector.
%
% The H matrix describes the relationship between the system being
% measured and the observation in the z vector via the equation
%
%       z = H*x + v
%
% where
%     z(m,1)     observation vector
%     H(m,n)     observation coefficient matrix
%     x(n,1)     the true state of the system
%     v(n,1)     observation noise,
%                gaussian with zero mean, and covariance R.
%
% The algorithm would normally do some expensive matrix operations
% which involve data that is usually invariant. Given the sensor
% observation coefficient matrix H and the sensor noise covariance
% matrix R, these are:
%
%    HtRi    = H' / R;
%    HtRiH   = HtRi * H;
%
% Parameters and results are:
%   Im(n,n)     information matrix
%   Iv(n,1)     information vector
%   HtRi(n,m)   Precalculated matrix, (H'/R)
%   HtRiH(n,n)  Precalcualted matrix, (H'/R)*H
%   z(m,1)      observation vector
%
% Runtime for this code is O( n * (n+m) )

function [Im, Iv] = info_data(Im, Iv, HtRi, HtRiH, z)

    Im = Im + HtRiH;
    Iv = Iv + HtRi * z;
