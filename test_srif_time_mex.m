function [ xo, vo ] = test_srif_time_mex(x,Rdz,F,RdH,P,Q)

    Nx = size(x,1);
    Nt = size(Rdz,2);

    xo = zeros(Nx,Nt);
    vo = zeros(Nx,Nt);

    Fi = inv(F);
    G  = eye(Nx);

    [Wm Wv] = srif_init(zeros(Nx,1),Q);
    [Im Iv] = srif_init(x,P);

    for i=1:Nt

        [Im Iv] = srif_time_mex(Im,Iv,Fi,G,Wm,Wv);
        [Im Iv] = srif_data_mex(Im,Iv,RdH,Rdz(:,i));
        [x P S] = srif_read_mex(Im,Iv);

        xo(:,i) = x;
        vo(:,i) = diag(P);

    end
