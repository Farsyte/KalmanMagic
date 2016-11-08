%% Information Filter: Read Estimate
%
% Given the state of the information filter, calculate the
% corresponding estimate of the state of the system and
% the covariance of the error in that estimate.
%
% For this to work, the Im matrix must be invertable.

function [x P] = info_read(Im, Iv)

    P = inv(Im);
    x = P * Iv;
