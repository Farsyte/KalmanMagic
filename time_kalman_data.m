function [tdts dtus] = time_kalman_data(x,P,H,R,z,reps)

    Nt = size(z,2);

    for i=1:Nt
        [x P] = kalman_data(x,P,H,R,z(:,i));
    end

    tic;
        for r=1:reps
            for i=1:Nt
                [x P] = kalman_data(x,P,H,R,z(:,i));
            end
        end
    tdts = toc / reps;

    iters = Nt;
    dtus = tdts*1e6/iters;
