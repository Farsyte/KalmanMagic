function [tdts dtus] = time_srif_data_mex(x,P,RdH,Rdz,reps)

    Nt = size(Rdz,2);

    [Im Iv] = srif_init(x, P);
    for i=1:Nt
        [Im Iv] = srif_data_mex(Im,Iv,RdH,Rdz(:,i));
    end

    tic;
        for r=1:reps
            for i=1:Nt
                [Im Iv] = srif_data_mex(Im,Iv,RdH,Rdz(:,i));
            end
        end
        tdts = toc / reps;

        iters = Nt;
        dtus = tdts*1e6/iters;
