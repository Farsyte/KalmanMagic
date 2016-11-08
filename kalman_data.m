%% Generic Kalman Filter: data processing
%
% Baseline implementation of the Data Processing phase of a Kalman
% Filter estimating a linear discrete system.
%
%   x(n,1)   state estimate
%   P(n,n)   state covariance
%   H(m,n)   observation coefficients
%   R(m,m)   observation covariance
%   z(m,1)   observation

function [x P] = kalman_data(x,P,H,R,z)

    y = z - H*x;                % innovation residual
    S = H*P*H' + R;             % innovation covariance
    K =   P*H' / S;             % optimal kalman gain
    x = x + K*y;                % updated state estimate
    P = P - K*H*P;              % updated state covariance
