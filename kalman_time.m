%% Generic Kalman Filter: time processing
%
% Based on the article in Wikipedia.
%
% NOTE: no "u" term here. Real implementations will
% need to expand the code to include control terms.
%
% Parameters and return values:
%   x(n,1)      state estimate
%   P(n,n)      error covaraince
%   Q(n,n)      process noise covariance
%   F(n,n)      state update transform matrix

function [x,P] = kalman_time(x,P,Q,F)

    x = F*x;            % update state estimate
    P = F*P*F' + Q;     % update error covariance
