%% Construct a clean workspace with data needed for Estimation.

clear;

Nt = 1000;              % number of iterations
rand('twister',0);      % repeatable random numbers
randn('state',0);       % repeatable random numbers

systemone;

initstate;

disp(' ');
disp('comparing data processing functions');
disp(' ');

% baseline system that provides reference data

[xref vref] = test_kalman_data(x,P,H,R,z);
plot_zxe('System One',z,H,xref,vref,z_,x_);

% crosscheck INFO vs KF ...

[xchk vchk] = test_info_read(x,P,HtRi,HtRiH,z);
test_disp_xe('info_data',                xchk,vchk,xref,vref,z,H,z_,x_);

% crosscheck POTTER vs KF ...

[xchk vchk] = test_potter_data(x,S,RdH,Rdz);
test_disp_xe('potter_data',              xchk,vchk,xref,vref,z,H,z_,x_);

[xchk vchk] = test_potter_data_mex(x,S,RdH,Rdz);
test_disp_xe('potter_data_mex',          xchk,vchk,xref,vref,z,H,z_,x_);

% crosscheck UDUT vs KF ...

[xchk vchk] = test_udut_data_mex(x,Pu,Pd,R,H,z);
test_disp_xe('udut_data_mex',            xchk,vchk,xref,vref,z,H,z_,x_);

% crosscheck SRIF vs KF ...

[xchk vchk] = test_srif_read(x,P,RdH,Rdz);
test_disp_xe('srif_read',                xchk,vchk,xref,vref,z,H,z_,x_);

[xchk vchk] = test_srif_data_mex(x,P,RdH,Rdz);
test_disp_xe('srif_data_mex',            xchk,vchk,xref,vref,z,H,z_,x_);

[xchk vchk] = test_srif_read_mex(x,P,RdH,Rdz);
test_disp_xe('srif_read_mex',            xchk,vchk,xref,vref,z,H,z_,x_);

