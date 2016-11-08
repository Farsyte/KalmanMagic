function [tdts dtus] = time_srif_data(x,P,RdH,Rdz,reps)

    Nt = size(Rdz,2);

    [Im Iv] = srif_init(x, P);
    for i=1:Nt
        [Im Iv] = srif_data(Im,Iv,RdH,Rdz(:,i));
    end

    tic;
        for r=1:reps
            for i=1:Nt
                [Im Iv] = srif_data(Im,Iv,RdH,Rdz(:,i));
            end
        end
    tdts = toc / reps;

    iters = Nt;
    dtus = tdts*1e6/iters;
