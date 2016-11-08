%% System Model Two
%
% This matlab script generates the initial data
% for System Two, of the form
%
%       x := F x + w
%       z := H x + v

%% System Two: oscillations in the plane.
%
% This system models the position and velocity
% of a particle on a plane. The particle starts
% at the origin, with initial velocity selected
% so that the range of the position along each
% axis is equal to the range of the velocity along
% the same axis.
%
% States:
%   x(1) = 10 sin (w1 t) = X position
%   x(2) = w1 cos (w1 t) = X velocity
%   x(3) = 10 sin (w2 t) = Y position
%   x(4) = w2 cos (w2 t) = Y velocity
%
%   x0(1) =  0
%   x0(2) = w1
%   x0(3) =  0
%   x0(4) = w2
%
%   dot x(1) =  w1   cos(w1 t) =     x(2)
%   dot x(2) = -w1^2 sin(w1 t) = -w1 x(1)
%   dot x(3) =  w2   cos(w2 t) =     x(4)
%   dot x(4) = -w2^2 sin(w2 t) = -w2 x(3)
%
% Writing the "dot x(*)" terms as a matrix, we get:
%
%       +   0     1   0     0  +
%       |  -w1^2  0   0     0  |
%       |   0     0   0     1  |
%       +   0     0  -w2^2  0  +
%
% So: "w1" is the square root of the spring constant in X,
% and "w2" is the square root of the spring constant in Y.
%
% Oscillation frequency in X is 2*pi*w1
% Oscillation frequency in Y is 2*pi*w2
%
% Call this "logF" for the moment ...
%
% For very tiny dt,
%       F = I + logF*dt
% but only when dt is very very small.
%
% As it turns out, for any dt,
%       F(dt) = expm(logF*dt)
% which is where the name logF comes from.
%
% Damping: in addition to the above, augmented
% the matrix to include a friction term which
% gives an acceleration proportional to
% the velocity.

%% Specify Fixed Properties
%
% These values can not be changed without
% making corresponding source code changes.
%

Nx = 4;         % dimension of state vector
Nw = 2;         % number of independent system noise sources
Nz = 5;         % number of sensors
Nv = 2;         % number of independent sensor noise sources

% ... for now, force Nw == Nx and Nv == Nz,
% the reduced rank noise matrices seem to cause problems.

Nw = Nx;
Nv = Nz;

%% Specify Variable Properties
%
% These values could be set externally to any
% desried value (within reason).

dt = 0.01;      % elapsed simulation time per update
k1 = 4.0;       % spring constant in X direction
k2 = 9.0;       % spring constant in Y direction
k3 = -0.00;     % velocity damping

Tw = randn(Nx,Nw) * 1e-2;        % process noise transfer matrix

% Tv = randn(Nz,Nv) * 1e-1;        % sensor noise transfer matrix
% R = Tv*Tv' needs to be diagonal for our UD data update to work.
Tv = eye(Nz) * 1e-1;

% Decrease the noise on the Position channels.
Tw(1,:) = 1e-3 * Tw(1,:);
Tw(3,:) = 1e-3 * Tw(3,:);


%% Specify System Model
%
% We use Matlab's Matrix Exponential here
% to generate an exact transfer matrix for
% the current delta-time, corresponding to
% the differential equations:
%
%     dot x(1) =     x(2)
%     dot x(2) = -k1 x(1)
%     dot x(3) =     x(4)
%     dot x(4) = -k2 x(3)
%
% Oscillation frequency in X is 2pi*sqrt(k1)
% Oscillation frequency in Y is 2pi*sqrt(k2)

logF = [ 0   1   0   0 ;
         -k1 k3   0   0 ;
         0   0   0   1 ;
         0   0  -k2 k3 ];

F = expm(logF*dt);

%% Specify Sensor Model
%
% z = H*x + v

H = randn(Nz,Nx);

%% Derive Process Noise
%
% w sampled from N(0,Q)
%
% Note: if Nw < Nx, we may end up
% with complex data -- or NaN -- in
% the otherwise unused parts of lQ.

Q = outer(Tw);                  % noise covariance matrix
w = Tw * randn(Nw,Nt);          % complete process noise

% To verify noise statistics,
% eWmean = mean(w,2);
% eWcov = cov(w') - Q;


%% Derive System True State

x0 = [ 1 0 0 1 ]';              % initial state
x_ = zeros(Nx,Nt);              % preallocate storage
x = x0;
for t=1:Nt
    x = F*x + w(:,t);
    x_(:,t) = x;
end

%% Initial State Estimates

x = zeros(Nx,1);

S = eye(Nx) * 1e3;
P = S * S';
Pu = eye(Nx);
Pd = diag(P);

%% Phase Subspace Plots
%
% Disabled for now.
% Besides, moved "z" production
% out to the common code.
%
% fig = figure('Position',[1920 600 600 400], ...
%              'Name','System Two State Truth');
% subplot(2,2,1); plot(x_(1,:),x_(2,:));
% subplot(2,2,2); plot(x_(1,:),x_(3,:));
% subplot(2,2,3); plot(x_(2,:),x_(4,:));
% subplot(2,2,4); plot(x_(3,:),x_(4,:));
%
% fig = figure('Position',[2560  24 600 400], ...
%              'Name','System Two State Sense');
% subplot(2,2,1); plot(z(1,:),z(2,:));
% subplot(2,2,2); plot(z(1,:),z(3,:));
% subplot(2,2,3); plot(z(2,:),z(4,:));
% subplot(2,2,4); plot(z(3,:),z(4,:));
