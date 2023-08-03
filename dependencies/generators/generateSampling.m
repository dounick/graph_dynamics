function S = generateSampling(regime, n, m, T, p)
    switch(regime)
        case 0
            S = fullSampling(n,T);
        case 1
            S = regime1Sampling(n,m,T, p);
        case 2
            S = regime2Sampling(n,m,T, p);
        case 3
            S = regime3Sampling(n,m,T, p);
    end
end