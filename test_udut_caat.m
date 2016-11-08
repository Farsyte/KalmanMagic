function ssqe = test_udut_caat(d)

    P = randpd(d);
    c = exp(randn);
    a = randn(d,1);

    Pchk = P + c*a*a';

    [U    D   ] = udut_fact_mex(P   );
    [Uchk Dchk] = udut_fact_mex(Pchk);
    [Ubar Dbar] = udut_caat(U,D,c,a);

    Uerr = Ubar - Uchk;
    Derr = Dbar - Dchk;
    ssqe = n2(Uerr) + n2(Derr);
