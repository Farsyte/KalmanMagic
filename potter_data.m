%% Data Processing using the Potter Algorithm
%
%   x(n,1)      state estimate
%   S(n,n)      error covariance square root
%   RdH(m,n)    normalized observation transfer matrix
%   Rdz(m)      normalized observation vector
%
% The potter mechanization requires that the sensor
% data and its transfer matrix have been normalized
% so that the sensor noise is independent, zero mean
% and unit variance.

function [x S] = potter_data(x,S,RdH,Rdz)

    m = size(Rdz,1);

    for j=1:m
        Hj = RdH(j,:);
        Zj = Rdz(j);

        v = Hj * S;
        sigma = 1 / (v*v' + 1);         % residual variance inverse
        K = S*v';                       % unweighted kalman gain
        d = Zj - Hj*x;                  % predicted residuals
        x = x + K*(d*sigma);            % updated state estimate
        gamma = sigma / (1 + sqrt(sigma));
        S = S - gamma*K*v;              % square root covariance update

    end
