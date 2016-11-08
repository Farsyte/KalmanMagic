function e = test_udut_fact_mex(d)

    I = eye(d);
    S = randn(d,d);
    P = S * S';

    Ptmp = P;

    [U D] = udut_fact_mex(Ptmp);

    Pchg = Ptmp - P;

    Perr = U * diag(D) * U' - P;

    e = n2(Pchg) + n2(Perr);

