function [tdts dtus] = time_udut_data_mex(x,Pu,Pd,H,R,z,reps)

    Nx = size(x,1);
    Nz = size(z,1);
    Nt = size(z,2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    % udut_data assumes that there is
    % no correlation between sensors.

    Rd = diag(R);

    for i=1:Nt
        [x Pu Pd] = udut_data_mex(x,Pu,Pd,H,Rd,z(:,i));
    end

    tic;
        for r=1:reps
            for i=1:Nt
                [x Pu Pd] = udut_data_mex(x,Pu,Pd,H,Rd,z(:,i));
            end
        end
    tdts = toc / reps;

    iters = Nt;
    dtus = tdts*1e6/iters;
