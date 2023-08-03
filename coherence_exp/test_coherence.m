function output = test_coherence(graph, regime, ptype, tau)
    %graph = 'm', 'r'
    %regime = 2, 3
    %ptype = 'unif', 'opt'
    %tau >= 1

    % include dependencies
    addpath(genpath('..')); 
    
    %initialize graph
    N = 2642;
    switch graph
        case "m"
            G = gsp_minnesota(); 
        case "r"
            G = gsp_ring(N);
    end

    %iterate over s and k
    bwidth_max = 2600;
    [U, Eigs] = eigs(expm(-4*G.L), bwidth_max);
    bwidth_min = 50;
    bwidth_step = 50;
    s_min = 10;
    s_max = 200;
    s_step = 10;
    
    fT = zeros(N, 1);
    for i = 1:bwidth_max
        fT(i) = sqrt((Eigs(i,i)^(2*tau) - 1)/(Eigs(i,i)^2 - 1));
    end
    
    Ukt = zeros(tau*N, bwidth_max);
    for i = 0:tau-1
        for j = 1:bwidth_max
            coeff = Eigs(j,j)^i/fT(j);
            Ukt(i*N + 1:(i+1)*N, j) = coeff*U(:,j);
        end
    end
    
    output = zeros((s_max-s_min)/s_step + 1,(bwidth_max-bwidth_min)/bwidth_step+1);    
    for k = bwidth_min:bwidth_step:bwidth_max
        disp(k)
        for s = s_min:s_step:s_max
            output(s/s_step,k/bwidth_step) = compute_sparsity_coherence(Ukt(:, 1:k), s, tau, regime, ptype);
        end
    end




