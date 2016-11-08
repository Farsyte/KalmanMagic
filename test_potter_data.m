function [ xo, vo ] = test_potter_data(x,S,RdH,Rdz)

    Nx = size(x,1);
    Nz = size(Rdz,1);
    Nt = size(Rdz,2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    for i=1:Nt
        [x S] = potter_data(x,S,RdH,Rdz(:,i));
        P = S * S';
        xo(:,i) = x;
        vo(:,i) = diag(P);
    end
