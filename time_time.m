%% Construct a clean workspace with data needed for Estimation.

clear;

Nt = 1000;              % number of iterations
rand('twister',0);      % repeatable random numbers
randn('state',0);       % repeatable random numbers

systemtwo;

FASTreps = 10;          % how many times to run FAST implementations
SLOWreps =  1;          % how many times to run SLOW implementations

initstate;

fmt = '%-24s  %16.6f s/rep, %16.3f us/iter';

disp(sprintf('\nestablishing time processing runtime\n'));

[tdts dtus] = time_kalman_time(x,z,F,H,P,Q,R,FASTreps);
disp(sprintf(fmt, 'kalman_time',                tdts, dtus));

[tdts dtus] = time_info_time(x,Rdz,F,RdH,P,Q,R,FASTreps);
disp(sprintf(fmt, 'info_time',                  tdts, dtus));

[tdts dtus] = time_udut_time_mex(x,z,F,Tw,H,Pu,Pd,R,FASTreps);
disp(sprintf(fmt, 'udut_time_mex',              tdts, dtus));

[tdts dtus] = time_srif_time(x,Rdz,F,RdH,P,Q,R,SLOWreps);
disp(sprintf(fmt, 'srif_time',                  tdts, dtus));

[tdts dtus] = time_srif_time_mex(x,Rdz,F,RdH,P,Q,R,FASTreps);
disp(sprintf(fmt, 'srif_time_mex',              tdts, dtus));

