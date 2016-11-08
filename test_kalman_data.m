function [ xo, vo ] = test_kalman_data(x,P,H,R,z)

    Nx = size(x, 1);
    Nt = size(z, 2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    for i=1:Nt
        [x P] = kalman_data(x,P,H,R,z(:,i));
        xo(:,i) = x;
        vo(:,i) = diag(P);
    end
