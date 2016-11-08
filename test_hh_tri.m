function Essq = test_hh_tri(m,n,d)

    A = [ randn(m,n) eye(m) ];

    B = hh_tri(A,d);

    % The augmented part of the matrix, which was the identity
    % on input, should contain the transform matrix on output.

    T = B(:,n+1:end);

    % Test that T is orthogonal.

    Eorth = T*T' - eye(m);

    % Test that output is input, transformed by T.

    Etab  = T*A  - B;

    % Test that output is zero below the diagonal.

    Ezbd = B(1:m,1:d);
    for i=1:d
        Ezbd(1:i,i) = 0;
    end;

    Essq = n2(Eorth) + n2(Etab) + n2(Ezbd);
