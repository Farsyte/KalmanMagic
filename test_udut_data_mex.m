function [ xo vo ] = test_udut_data_mex(x,Pu,Pd,R,H,z)

    Nx = size(x,1);
    Nz = size(z,1);
    Nt = size(z,2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    % p100_udut_data assumes that there is
    % no correlation between sensors.

    Rd = diag(R);

    for i=1:Nt

        [x Pu Pd] = udut_data_mex(x,Pu,Pd,H,Rd,z(:,i));

        xo(:,i) = x;

        for k=1:Nx
            vo(k,i) = Pd(k) +                                        ...
                      sum( Pu(k,k+1:Nx) .*                           ...
                           Pd(k+1:Nx)' .*                            ...
                           Pu(k,k+1:Nx) );
        end
    end
