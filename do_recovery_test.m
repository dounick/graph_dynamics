function errors = do_recovery_test(graph, tau, k, s, ptype, algo)

    %basic recovery error experiment for a fixed graph ("m", "r"), 
    %   ptype ("sparsity", "bandlimited"), and algo ("cosamp", "lasso")
    %sample usage: do_recovery_test("m", 5, 100, 2, "sparsity", "cosamp")
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
    for i = 1:trial
        i
        x = generateSparse(tau*N, s, k, Ukt);
        for m = tau:tau:tau^2
            S = sparse(regime2Sampling(N,ceil(m*ratio),tau,p));
            POmega = diag(diag(S*diag_p*transpose(S)));
            Phi = POmega^(1/2)*diag(1/sqrt(m/10)*ones(m,1))*S;
            u = Phi*x;
            switch algo
                case "cosamp"
                    freqapprox = cosamp(Phi*Ukt, u, s, k, 10);
                    approx = Ukt(1:N, :)*freqapprox;
                    err = vecnorm(approx-x(1:N))/vecnorm(x(1:N));
                case "lasso"
                    B = lasso(Phi*Ukt, u);
                    sz = size(B,2);
                    err = 1e5;
                    for j = 1:sz
                        approx = Ukt(1:N, :)*B(:,j);
                        newerr = vecnorm(approx - x(1:N))/vecnorm(x(1:N));
                        if newerr < err
                            err = newerr;
                        end
                    end
            end
            errors(m/tau, i) = err;
        end
    end

