function [ xo, vo ] = test_info_time(x,z,F,H,P,Q,R)

    Nx = size(x,1);
    Nt = size(z,2);

    Fi = inv(F);
    Qi = inv(Q);


    HtRi  = H' / R;
    HtRiH = HtRi * H;

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    [Im Iv] = info_init(x,P);

    for i=1:Nt

        [Im Iv] = info_time(Im,Iv,F,Fi,Q,Qi);
        [Im Iv] = info_data(Im,Iv,HtRi,HtRiH,z(:,i));
        [x P]   = info_read(Im,Iv);

        xo(:,i) = x;
        vo(:,i) = diag(P);

    end
