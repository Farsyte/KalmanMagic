%% Information Filter: Initialize
%
% Given an estimate of the state of the system and the covariance of
% the error, create the corresponding state variables for the
% Information Filter.
%
% The parameters describe the actual state X via the equation:
%        X =  x + N(0,P)
%
% The filter state describes the actual state X via the equation:
%     Im*X = Iv + N(0,I)
%
% It is better to construct (Im,Iv) by hand, carefully, than to
% attempt to mess with generating them from P, especially if
% your initial P matrix is not invertable.
%
% The best way to initialize Im and Iv is to fill them with zeros
% representing a complete lack of information, then update the
% filter with the initial observations.

function [Im Iv] = info_init(x, P)

    Im = inv(P);
    Iv = Im * x;
