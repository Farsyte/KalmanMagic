%% Construct a clean workspace with data needed for Estimation.

clear;

Nt = 1000;              % number of iterations
rand('twister',0);      % repeatable random numbers
randn('state',0);       % repeatable random numbers

systemtwo;

initstate;

disp(' ');
disp('comparing system update functions');
disp(' ');

% baseline system that provides reference data

[xref vref] = test_kalman_time(x,z,F,H,P,Q,R,Nt,Nx);
plot_zxe('System Two',z,H,xref,vref,z_,x_);

[xchk vchk] = test_info_time(x,z,F,H,P,Q,R);
test_disp_xe('info_time',                xchk,vchk,xref,vref,z,H,z_,x_);

[xchk vchk] = test_srif_time(x,Rdz,F,RdH,P,Q);
test_disp_xe('srif_time',                xchk,vchk,xref,vref,z,H,z_,x_);

[xchk vchk] = test_srif_time_mex(x,Rdz,F,RdH,P,Q);
test_disp_xe('srif_time_mex',            xchk,vchk,xref,vref,z,H,z_,x_);

[xchk vchk] = test_udut_time_mex(x,z,F,Tw,H,Pu,Pd,R);
test_disp_xe('udut_time_mex',            xchk,vchk,xref,vref,z,H,z_,x_);

