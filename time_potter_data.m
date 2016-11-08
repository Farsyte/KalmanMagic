function [tdts dtus] = time_potter_data(x,S,RdH,Rdz,reps)

    Nt = size(Rdz,2);

    for i=1:Nt
        [x S] = potter_data(x,S,RdH,Rdz(:,i));
    end

    tic;
        for r=1:reps
            for i=1:Nt
                [x S] = potter_data(x,S,RdH,Rdz(:,i));
            end
        end
    tdts = toc / reps;

    iters = Nt;
    dtus = tdts*1e6/iters;
