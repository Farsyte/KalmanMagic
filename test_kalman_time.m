function [ xo, vo ] = test_kalman_time(x,z,F,H,P,Q,R,Nt,Nx)

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    for j=1:Nt

        [x P] = kalman_time(x,P,Q,F);
        [x P] = kalman_data(x,P,H,R,z(:,j));

        xo(:,j) = x;
        vo(:,j) = diag(P);

    end
