%% Set up to model System One.
%
% This system has three state variables,
% and four sensors. It does not evolve.
%
% z = H x + w


Nx = 3;                         % number of states
Nz = 4;                         % number of sensors

% ... for now, force Nw == Nx and Nv == Nz,
% the reduced rank noise matrices seem to cause problems.

Nw = Nx;
Nv = Nz;

% Noise Properties

Tw = eye(Nx) * 1e-3;            % process noise transfer matrix
Tv = diag([0.1 0.2 0.3 0.4]);   % sensor noise transfer matrix

% System Model: never changes.

F = eye(Nx);

% Sensor Model

H  = [ 1 2 3 ; ...
       4 5 6 ; ...
       1 3 5 ; ...
       2 5 9 ];

% Process Noise

Q = outer(Tw);
w = zeros(Nw,Nt);

% System True State

x0 = [ 10 20 30 ]';             % initial true state
x_ = x0 * ones(1,Nt);           % constant true state

%% Initial State Estimates
%
% NOTE: covariance methods do not work well
% if the initial covariance is too large,
% and there seems to be a problem with Potter
% when we start with any S other than eye(Nx).

S = eye(Nx) * 1e3;
P = S * S';
Pu = eye(Nx);
Pd = diag(P);

x = zeros(Nx,1);
