function [ xo, vo ] = test_info_read(x,P,HtRi,HtRiH,z)

    Nx = size(x,1);
    Nt = size(z,2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    [Im Iv] = info_init(x,P);

    for i=1:Nt
        [Im Iv] = info_data(Im,Iv,HtRi,HtRiH,z(:,i));
        [x P]   = info_read(Im,Iv);
        xo(:,i) = x;
        vo(:,i) = diag(P);
    end
