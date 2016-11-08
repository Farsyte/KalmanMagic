%% Time Processing using the Information Matrix
%
% This function updates the information matrix based on a state
% transition of the system through the state update function,
%       X = F X + W
% where W is gaussian noise with zero mean, and covariance Q.
%
% Parametes and return data:
%
%   Im(n,n)     information matrix
%   Iv(n,1)     information vector
%   F(n,n)      state upate transformation matrix
%   Fi(n,n)     precalculated inverse of the F matrix
%   Q(n,n)      process noise covariance
%   Qi(n,n)     precalculated inverse of the Q matrix

function [Im Iv] = info_time(Im,Iv,F,Fi,Q,Qi)

    Nx = size(Iv,1);
    M = Fi' * Im * Fi;
    C = M / (M + Qi);
    L = eye(Nx) - C;
    Im = L*M*L' + C*Qi*C';
    Iv = L*Fi'*Iv;
