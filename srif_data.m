%% SRIF Data Processing Update [REALLY SLOW IN MATLAB]
%
% This code implements SRIF data processing by forming the
% filter state and observation into a single matrix, then
% performaing Householder transforms to triangularize it.
%
% Why do we do this? It all comes back to minimizing a
% performance function. The incoming filter contains
% matrix Im and vector Iv such that our estimated state
% minimizes the performance function:
%
% ~J(x) = || ~Im x - ~Iv ||^2
%
% The incoming observation also carries information about
% the system, in that the system x is expected to minimize
% a sensor-related performance function:
%
% zJ(x) = ||  H  x -  z  ||^2
%
% Knowing both our incoming filter state (Im,Iv) and the
% recent observation (H,z) we refine our estimate of the
% true system state x; we think that it is the state that
% minimizes the sum of the above two functions:
%
% ^J(x) = || ~Im x - ~Iv ||^2 + ||  H  x -  z  ||^2
%
% This procedure will transform our incoming data to
% provide updated values ^Im and ^Iv that represent
% this revised understanding of the true state x.
%
% First, rewrite ^J(x) as:
%
%         ||  | ~Im |      | ~Iv | ||2
% ^J(x) = ||  |     | x -  |     | ||
%         ||  |  H  |      |  z  | ||
%
% ( Sorry about the Bad ASCII Graphics: concatinate ~Im and H
%   vertically into a matrix, and multiply that by X; then
%   concatinate vector ~Iv with observation z, and subtract
%   from the result of prior multiply. Take the sum of the
%   squares of all resulting components. )
%
% Then we note that we can apply any orthogonal transform T to the
% matrices inside, without changing the value of X that minimizes
% the performance function.
%
%         ||  | ~Im |      | ~Iv | ||2
% ^J(x) = || T|     | x - T|     | ||
%         ||  |  H  |      |  z  | ||
%
% So carefully select the T matrix to change H to 0, and we
% get this result:
%
%         ||  | ^Im |      | ^Iv | ||2
% ^J(x) = ||  |     | x -  |     | ||
%         ||  |  0  |      |  e  | ||
%
% We do not need to obtain T, we just need to obtain the results
% of applying it to our data.
%
% One stable and efficient way to do this is to use Householder
% transformations to triangularize the matrix.
%
% The basic unit of work is to use a Householder transform to
% change the first column of a matrix so that it is all zeros
% except for the first element. Geometrically viewed, we select
% a vector "U" to be used as the mirror direction, which will
% map the first column (as a vector) onto the desired result.
%
% The Householder transform of a vector V using the mirror
% direction U is done by taking the component of V that is
% parallel to U, and reversing it:
%
%   V = V - 2 * ( V'U / U'U )
%
% We select the U vector so that the first column of the matrix
% has all but its first element cleared, then apply the same
% transform to all other columns of the matrix.
%
% The above procedure is repeated for submatrices of decreasing
% size, ranging from the next diagonal element to be processed
% down to the lower right corner of the matrix. Note that we do
% not need to apply the transform above or to the left of the
% diagonal element being considered: to the left is all zeros,
% and the U vector is implicitly zero in all elements above the
% area we are manipulating.
%
%
% The implementation here is slightly complicated by the fact
% that the matrix on which we are operating is kept in four
% separate matrices, causing chunks of code to be repeated
% for the different parts of the big matrix.
%
% At each iteration, we use the transform that takes the column vector
% from the selected element to the bottom of the matrix, and
% transforms it so that all elements but the first are zero.
%
% The sign of "s" is selected to avoid subtracting two large numbers
% with a small difference, and more importantly to avoid inverting a
% very tiny number dominated by numerical errors.
%
% Runtime is O( D*(D+M)*(D+N) )

function [Im Iv] = srif_data(Im,Iv,H,z)

%  Im(n,n)      information matrix
%  Iv(n,1)      information vector
%  H (m,n)      observation transform   (normalized)
%  z (m,1)      observation vector      (normalized)
%
% Note that H and Z must be normalized so that the sensor
% noise is a collection of independent gaussian variables
% with zero mean and unit variance.
%
% We will perform N householder transforms on the matrix
%
%       |  ~Im  ~Iv |  (n)
%       |   H    z  |  (m)
%          (n)  (1)
%
% giving us a result matrix
%
%       |  ^Im  ^Iv |  (n)
%       |   0    e  |  (m)
%          (n)  (1)
%
% The actula transform T is never constructed, just the
% results. We discard the final value of "e".
%
% As it turns out, we do not actually need to write
% the "zero" values back into the H area, saving us a
% little time, but H is still destroyed by the procedure.

    [m n] = size(H);

    for l=1:n

        % Use a Householder transform on the submatrix (from this
        % diagonal element through the lower right corner) to set
        % all elements below the diagonal in this column to zero.
        %
        % If they are already "close enough" to zero, just pushing
        % zero values into them gives the same result as doing the
        % full transform.
        %
        % Detect "close enough" by comparing the square of the
        % diagonal element value with the sum of the squares of all
        % values; if adding in the squares of the data below the
        % diagonal does not increase the total, they are small
        % enough.

        oldD = Im(l,l);
        oldD2 = oldD ^ 2;
        s2 = oldD2 + n2(Im(l+1:n,l)) + n2(H(1:m,l));

        if (s2 > oldD2)

            % "s" will be the final value we write back to the
            % diagonal element. The magnitude is dictated by the
            % requirement to maintain the two-norm of the column
            % vector, and we select the sign to avoid numerical
            % problems when calculating the value of the first
            % entry in the transformation's U vector.

            s = sqrt(s2);
            if (oldD > 0) s = -s; end;

            topU = oldD - s;
            beta = 1 / (s * topU);

            % Now we perform the transform, where the first
            % element of the U vector is in topU, and the rest
            % of it is in Im(l+1:n,l) and H(1:m,l).

            for k=l+1:n
                gamma = topU * Im(l,k);
                gamma = gamma +  dot( Im(l+1:n,l), Im(l+1:n,k) );
                gamma = gamma +  dot( H (  1:m,l), H (  1:m,k) );
                gamma = beta * gamma;

                Im(l    ,k) = Im(l    ,k) + gamma * topU       ;
                Im(l+1:n,k) = Im(l+1:n,k) + gamma * Im(l+1:n,l);
                H (  1:m,k) = H (  1:m,k) + gamma * H (  1:m,l);
            end

                gamma = topU * Iv(l,1);
                gamma = gamma +  dot( Im(l+1:n,l),Iv(l+1:n,1) );
                gamma = gamma +  dot( H (  1:m,l),z (  1:m,1) );
                gamma = beta * gamma;

                Iv(l    ,1) = Iv(l    ,1) + gamma * topU;
                Iv(l+1:n,1) = Iv(l+1:n,1) + gamma * Im(l+1:n,l);
                z(   1:m,1) = z(   1:m,1) + gamma * H (  1:m,l);

            Im(l,l) = s;
        end

        % Explicitly set elements below the diagonal to zero,
        % even if they were so small that we did not need
        % to perform the transformation.
        %
        % Writing zeros to H(1:m,l) turns out to be unnecessary.

        Im(l+1:n,l) = 0;

    end
