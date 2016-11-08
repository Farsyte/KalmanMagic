function Essq = test_hh_tri_hex_mex(d1,d2,d3)

    d = d1 + d2 - 1;
    m = d1 + d2;
    n = d1 + d2 + d3;

    Ar = eye(m);

    Aa = randpd(d1);
    Ab = randn(d1,d2);
    Ac = [ randn(d1,d3)    Ar(1:d1,:) ];

    Ad = randn(d2,d1);
    Ae = randpd(d2);
    Af = [ randn(d2,d3)   Ar(d1+1:d1+d2,:) ];

    A = [ Aa Ab Ac ; Ad Ae Af ];

    [Ba Bb Bc Bd Be Bf] = hh_tri_hex_mex(Aa,Ab,Ac,Ad,Ae,Af);

    B = [Ba Bb Bc ; Bd Be Bf];


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
