function [tdts dtus] = time_srif_time_mex(x,z,F,H,P,Q,R,reps)

    Nx = size(x,1);
    Nt = size(z,2);

    Fi = inv(F);
    G  = eye(Nx);

    [Wm Wv] = srif_init(zeros(Nx,1),Q);
    [Im Iv] = srif_init(x,P);

    for i=1:Nt
        [Im Iv] = srif_time_mex(Im,Iv,F,G,Wm,Wv);
        [Im Iv] = srif_data_mex(Im,Iv,H,z(:,i));
        [x P S] = srif_read_mex(Im,Iv);
    end

    tic;
        for r=1:reps
            for i=1:Nt
                [Im Iv] = srif_time_mex(Im,Iv,F,G,Wm,Wv);
                [Im Iv] = srif_data_mex(Im,Iv,H,z(:,i));
                [x P S] = srif_read_mex(Im,Iv);
            end
        end
        tdts = toc / reps;

        iters = Nt;
        dtus = tdts*1e6/iters;
