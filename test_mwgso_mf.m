function e = test_mwgso_mf(n,N)

    % usually called with n<N

    W = randn(n,N);
    D = exp(randn(N,1));

    % NOTE: assumes the rows of W are linearly independent.

    WDWt = W * diag(D) * W';

    [Ubar Dbar V] = mwgso_mf(W,D);

    VDVt = V * diag(D) * V';
    e1 = VDVt - diag(Dbar);

    udut = Ubar * diag(Dbar) * Ubar';
    e2 = udut - WDWt;

    e = n2(e1) + n2(e2);
