function S = regime2Sampling(n, M, T,p)
    Ms = sum(M);
    if Ms > T*n
        error('M cannot be larger than T*n');
    end
    S = zeros(Ms, T*n);
    index = 1;
    for i = 1:T
        probDist = p((i-1)*n + 1 : i*n);
        %for j = 1:n
        %    probDist(j) = unifrnd(0,1);
        %end
        %s = sum(probDist);
        %for j = 1:n
        %    probDist(j) = probDist(j)/s;
        %end
        m = M(i);
        randmSamp = datasample(1:n, m, 'Weights', probDist);
        %randmSamp = datasample(1:n, m);
        sampM = zeros(m,n);
        for j = 1:m
            sampM(j, randmSamp(j)) = 1;
        end
        S(index : index + m - 1, (i-1)*n+1 : i*n) = sampM;
        index = index + m;
    end
end