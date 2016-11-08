%% Construct a clean workspace with data needed for Estimation.

clear;

Nt = 1000;              % number of iterations
rand('twister',0);      % repeatable random numbers
randn('state',0);       % repeatable random numbers

systemone;

FASTreps = 10;          % how many times to run FAST implementations
SLOWreps =  1;          % how many times to run SLOW implementations

initstate;

disp(sprintf('\nestablishing data procesing runtime\n'));

fmt = '%-24s  %16.6f s/rep, %16.3f us/iter';

[tdts dtus] = time_kalman_data(x,P,H,R,z,FASTreps);
disp(sprintf(fmt, 'kalman_data',                tdts, dtus));

[tdts dtus] = time_info_data(x,H,P,R,z,FASTreps);
disp(sprintf(fmt, 'info_data',                  tdts, dtus));

% was commented out ... too slow?
[tdts dtus] = time_potter_data(x,S,RdH,Rdz,SLOWreps);
disp(sprintf(fmt, 'potter_data',                tdts, dtus));

[tdts dtus] = time_potter_data_mex(x,S,RdH,Rdz,FASTreps);
disp(sprintf(fmt, 'potter_data_mex',            tdts, dtus));

[tdts dtus] = time_udut_data_mex(x,Pu,Pd,H,R,z,FASTreps);
disp(sprintf(fmt, 'udut_data_mex',              tdts, dtus));

[tdts dtus] = time_srif_data(x,P,RdH,Rdz,SLOWreps);
disp(sprintf(fmt, 'srif_data',                  tdts, dtus));

[tdts dtus] = time_srif_data_mex(x,P,RdH,Rdz,FASTreps);
disp(sprintf(fmt, 'srif_data_mex',              tdts, dtus));

[tdts dtus] = time_srif_read(x,P,RdH,Rdz,FASTreps);
disp(sprintf(fmt, 'srif_read',                  tdts, dtus));

[tdts dtus] = time_srif_read_mex(x,P,RdH,Rdz,SLOWreps);
disp(sprintf(fmt, 'srif_read_mex',              tdts, dtus));

