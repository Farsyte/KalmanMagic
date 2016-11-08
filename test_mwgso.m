function e = test_mwgso(n,N)

    W = randn(n,N);
    D = exp(randn(N,1));

    % REQUIRED:
    %     W is a collection of n N-vectors (N>n).
    %     D is a collection of N positive scalars.
    %     W vectors are linearly independent.

    V = mwgso(W,D);

    % RESULT:
    %     Vi*D*Vj is zero when (i != j).

    VDVt = V * diag(D) * V';

    dVDVt = diag(VDVt);

    odVDVt = VDVt - diag(dVDVt);

    e = n2(odVDVt);
