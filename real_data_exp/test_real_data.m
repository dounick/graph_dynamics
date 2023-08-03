function errs = test_real_data(ptype, sparsity, T)
    %ptype = 'sparsity', 'bandlimited', 'unif'

    % include dependencies
    addpath(genpath('..')); 

    %load graph and signals
    load('full_graph_pressure.mat', 'G');
    load('Sea_level_Pressure.mat', 'Data')
    
    n = G.N;
    bwidth = 200;
    
    L = G.L;
    A = expm(-4*L);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [V,D] = eigs(A,n, bwidth);
    V = real(V);
    D = real(D);    
    
    Vk = V(:,1:bwidth);
    
    fT = zeros(n, 1);
    for i = 1:n
        for t = 0:T-1
            fT(i) = fT(i) + D(i,i)^(2*t);
        end
        fT(i) = sqrt(fT(i));
    end
    
    Utilde = zeros(T*n, bwidth);
    for i = 0:T-1
        for j = 1:bwidth
            coeff = D(j,j)^i/fT(j);
            Utilde(i*n + 1:(i+1)*n, j) = coeff*Vk(:,j);
        end
    end
    
    f = Data(:, 1);
    
    ratio = ones(T,1)/T;
    switch ptype
        case 'unif'
            p = ones(T*n, 1)/n;
        case 'sparsity'
            [ratio, p] = compute_sparsity_dist(Utilde, T, sparsity, 2);
        case 'bandlimited'
            p = compute_bandlimited_dist(Utilde, T, bwidth, 2);
    end
    pdiag = diag(p);
    normF = vecnorm(f);
    embF = embed(f,A,T);
    
    stepSize = 20;
    minSamples = 20;
    maxSamples = 200;
    trials = 100;
    
    errs = zeros(trials, (maxSamples - minSamples)/stepSize + 1);
    
    for t = 1:trials
        disp(t)
        for m = minSamples:stepSize:maxSamples
            S = sparse(regime2Sampling(n,round(m*ratio),T,p));
            POmega = diag(diag(S*pdiag*transpose(S)));
            Phi = POmega^(1/2)*S;
            u = Phi*embF;
            freqapprox = cosamp(Phi*Utilde, u, sparsity, bwidth, 10);
            approx = Utilde*freqapprox;
            err = vecnorm(approx(1:n) - f)/normF;
            errs(t,m/stepSize) = err;
        end
    end