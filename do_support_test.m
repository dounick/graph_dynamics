function errors = do_support_test(graph, tau, k, s, ptype, algo)

    %basic recovery error experiment for a fixed graph ("m", "r"), 
    %   ptype ("sparsity", "bandlimited"), and algo ("cosamp", "lasso")
    % example usage: do_support_test("m", 5, 100, 2, "sparsity", "cosamp")
    addpath(genpath("."))
    
    switch graph
        case "m"
            G = gsp_minnesota();
        case "r"
            G = gsp_ring(2642);
    end    
    N = G.N;
    W = G.A;
    D = diag(sum(W,2));
    L = D - W;
    
    Deltat = 4;
    A = expm(-Deltat*L);
    [U, E] = eigs(A, N);
    
    fT = zeros(N, 1);
    for i = 1:N
        fT(i) = sqrt((E(i,i)^(2*tau) - 1)/(E(i,i)^2 - 1));
    end
    
    Ukt = zeros(tau*N, k);
    for i = 0:tau-1
        for j = 1:k
            coeff = E(j,j)^i/fT(j);
            Ukt(i*N + 1:(i+1)*N, j) = coeff*U(:,j);
        end
    end

    trial = 100;
    switch ptype
        case "sparsity"
            [~, p] = compute_sparsity_dist(Ukt, tau, s, 2);
        case "bandlimited"
            [~, p] = compute_bandlimited_dist(Ukt, tau, k, 2);
    end
    diag_p = diag(p);
    ratio = ones(tau, 1)/tau;
    errors = zeros(tau, trial);
    for i = 1:100
        i
        x = generateSparse(tau*N, s, k, Ukt);
        truefreq = Ukt(:, 1:k)'*x;
        [~, trueSupp] = maxk(truefreq, s, 'ComparisonMethod','abs');
        trueSupp = sort(trueSupp);
        for m = 10:10:100
            S = sparse(regime2Sampling(N,ceil(m*ratio),tau,p));
            POmega = diag(diag(S*diag_p*transpose(S)));
            Phi = POmega^(1/2)*diag(1/sqrt(m/10)*ones(m,1))*S;
            u = Phi*x;
            switch algo
                case "cosamp"
                    freqapprox = cosamp(Phi*Ukt, u, s, k, 10);
                    [~, supp] = maxk(freqapprox, s, 'ComparisonMethod','abs');
                    supp = sort(supp);
                    if supp == trueSupp
                        errors(m/10, i) = 1;
                    end
                case "lasso"
                    B = lasso(Phi*Ukt, u);
                    sz = size(B,2);
                    err = 1e5;
                    bestIndex = 0;
                    for j = 1:sz
                        approx = Ukt*B(:,j);
                        newerr = vecnorm(approx(1:N) - x(1:N));
                        if newerr < err
                            err = newerr;
                            bestIndex = j;
                        end
                    end
                    freqapprox = B(:, bestIndex);
                    [~, supp] = maxk(freqapprox, s, 'ComparisonMethod','abs');
                    supp = sort(supp);
                    if supp == trueSupp
                        errors(m/10, i) = 1;
                    end
            end
        end
    end

