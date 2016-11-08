function [ xo vo ] = test_srif_read(x,P,RdH,Rdz)

    Nx = size(x,1);
    Nt = size(Rdz,2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    [Im Iv] = srif_init(x,P);

    for i=1:Nt
        [Im Iv] = srif_data(Im,Iv,RdH,Rdz(:,i));
        [x P S] = srif_read(Im,Iv);
        xo(:,i) = x;
        vo(:,i) = diag(P);
    end
