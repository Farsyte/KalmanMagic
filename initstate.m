%% Generate the initial state.

%% Observation Noise.
%
% We are going to need Tv below,
% so we can bring pure uses of Tv
% down into system-independent code.
%
% v will be sampled from N(0,R), where
%     R = Tv*Tv'

R = outer(Tv);                  % noise covariance matrix
v = Tv * randn(Nv,Nt);          % complete sensor noise

% To verify noise statistics,
% eVmean = mean(v,2);
% eVcov = cov(v') - R;

%% Derive Observations.
%
% NOTE: x_(:,i) is the true state at the end of the Ith interval,
% so the phasing is to do the time step to the time of the sensor
% observation, process the sensor observation, and generate the
% estimate of the state at the time the sensor was active.

z_ = H*x_     ;
z  =   z_ + v ;

%% Whiten Observation Noise.

% Some data processing routines require
% the equations to be normalized so that
% the sensor noise variance is 1, so we
% provide pretransformed H and z matrices
% that have that property.

% RdH = diag(Rsd)\H;
% Rdz = diag(Rsd)\z;

RdH = Tv\H;
Rdz = Tv\z;

% H'H is a frequent term. Precalculate it,
% and the normalized version of the same.

HtH   = H' * H;
RdHtH = RdH' * RdH;

% establish an error reporting high water mark.

eco = 1e-16;

% NOTE: each estimation method will have
% to provide its own "initial estimate"
% corresponding either to "no information"
% or seeded from the initial observation.


% XXX: maybe preload these for information filters?
% [Wm Wv] = p070_srif_init(zeros(Nx,1),Q);
% [Im Iv] = p070_srif_init(x,P);

HtRi    = H' / R;
HtRiH   = HtRi * H;

