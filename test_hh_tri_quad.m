function Essq = test_hh_tri_quad(m,n,d)

    Ar = eye(m);

    Aa = randpd(d);
    Ab = [ randn(d,n-d)    Ar(1:d,:) ];
    Ac = randn(m-d,d);
    Ad = [ randn(m-d,n-d)  Ar(d+1:m,:) ];

    A = [ Aa Ab ; Ac Ad ];

    [Ba Bb Bc Bd] = hh_tri_quad(Aa,Ab,Ac,Ad);

    B = [Ba Bb ; Bc Bd];

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
