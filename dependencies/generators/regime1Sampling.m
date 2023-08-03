function S = regime1Sampling(n, m, T, p)
    if m > n
        error('m cannot be larger than n');
    end
    S = zeros(T*m, T*n);
    randmSamp = datasample(1:n, m, 'Weights', p);
    %randmSamp = datasample(1:n, m);
    sampM = zeros(m,n);
    for j = 1:m
        sampM(j, randmSamp(j)) = 1;
    end
    for i = 1:T
        S(((i-1)*m + 1) : i*m, (i-1)*n+1:i*n) = sampM;
    end
    S = sparse(S);
end