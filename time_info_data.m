function [tdts dtus] = time_info_data(x,H,P,R,z,reps)

    Nt = size(z,2);

    HtRi  = H' / R;
    HtRiH = HtRi * H;

    [Im Iv] = info_init(x,P);

    for i=1:Nt
        [Im Iv] = info_data(Im,Iv,HtRi,HtRiH,z(:,i));
        [x P]   = info_read(Im,Iv);
    end

    tic;
    for r=1:reps
        for i=1:Nt
            [Im Iv] = info_data(Im,Iv,HtRi,HtRiH,z(:,i));
            [x P]   = info_read(Im,Iv);
        end
    end
    tdts = toc / reps;

    iters = Nt;
    dtus = tdts*1e6/iters;
