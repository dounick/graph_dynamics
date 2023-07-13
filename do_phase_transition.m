function success_rate = do_phase_transition(graph, regime, ptype, tau)
    addpath(genpath('..'))
    switch graph
        case "m"
            G = gsp_minnesota();
            N = G.N;
        case "r"
            N = 2642;
            G = gsp_ring(N);
    end
    bwidth = 400;
    A = exp(-4*G.L);
    [U, Eigs] = eigs(A, bwidth);
    disp("Finished initializing graphs and operators")
    s_min = 2;
    s_max = 2;
    s_step = 1;
    switch regime
        case 2
            p = ones(tau*N,1)/N;
            ratio = ones(tau,1)/tau;
        case 3
            p = ones(tau*N,1)/(tau*N);
            ratio = 1;
    end
    large_diag = diag(p);
    

    fT = zeros(bwidth, 1);
    for i = 1:bwidth
        fT(i) = sqrt((Eigs(i,i)^(2*tau) - 1)/(Eigs(i,i)^2 - 1));
    end
    
    Ukt = zeros(tau*N, bwidth);
    for i = 0:tau-1
        for j = 1:bwidth
            coeff = Eigs(j,j)^i/fT(j);
            Ukt(i*N + 1:(i+1)*N, j) = coeff*U(:,j);
        end
    end
    switch regime
        case 1
            p = compute_opt_dist(Utilde, k);
            ratio = 1;
    end
    success_rate = zeros(100, (s_max- s_min)/s_step + 1);
    write_path = strcat(graph, "_", ptype, "_", int2str(regime), "_", int2str(tau), ".csv");
    trials = 1;
    error_tol = 0.01;

    disp("Finished setting up for sampling, beginning recovery tests")
    for s = s_min:s_step:s_max
        s
        if ptype == "opt"
                [ratio, p] = compute_new_dist(Ukt, tau, s, regime);
                large_diag = diag(p);
        end
                
        for m = 10:10:1000
            successes = 0;
            for i = 1:trials
                x = generateSparse(tau*N, s, bwidth, Ukt);
                switch regime
                    case 1
                        S = sparse(regime2Sampling(N,ceil(m*ratio),tau,p));
                    case 2
                        S = sparse(regime2Sampling(N,ceil(m*ratio),tau,p));
                    case 3
                        S = sparse(regime3Sampling(N,m,tau,p));
                end
                POmega = diag(diag(S*large_diag*transpose(S)));
                Phi = POmega^(1/2)*S;
                u = Phi*x;
                freqapprox = cosamp(Phi*Ukt, u, s, bwidth, 10);
                approx = Ukt(1:N, :)*freqapprox;
                err = vecnorm(approx-x(1:N));
                if err < error_tol
                    successes = successes + 1;
                end
            end
            success_rate(m/10, s/s_step) = successes;
	    writematrix(success_rate, write_path);
        end
    end

