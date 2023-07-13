clear all, close all, clc;
%%%
addpath(genpath('.')); 

graph_name = "minnesota";

switch graph_name
    case "minnesota"
        G = gsp_minnesota();
        N = G.N;
        bwidth = 400;
        [U, Eigs] = eigs(G.L, bwidth, 'smallestabs');
        popt = compute_opt_dist(U, bwidth);
        p2 = compute_2_dist(U, bwidth);
        s_min = 2;
        s_max = 100;
        s_step = 2;
    case "ring"
        N = 2642;
        G = gsp_ring(N);
        bwidth = 400;
        [U, Eigs] = eigs(G.L, bwidth, 'smallestabs');
        popt = compute_opt_dist(U, bwidth);
        p2 = compute_2_dist(U, bwidth);
        s_min = 2;
        s_max = 100;
        s_step = 2;
end

trials = 100;
error_tol = 0.01;

diag_opt = diag(popt);
diag_2 = diag(p2);

success_rate = zeros(1000/10, (s_max- s_min)/s_step + 1);

%required_samples = zeros((s_max - s_min)/s_step + 1);

for s = s_min:s_step:s_max
    display(s)
    for m = 10:10:1000
        successes = 0;
        for i = 1:trials
            x = generateSparse(N, s, bwidth, U);
            S = sparse(regime3Sampling(N,m,1,p2));
            POmega = sparse(diag(diag(S*diag_2*transpose(S))));
            %S = sparse(regime3Sampling(N,m,1,popt));
            %POmega = sparse(diag(diag(S*diag_opt*transpose(S))));
            Phi = sparse(POmega^(-1/2)*S);
            u = Phi*x;
  
            approx = U*cosamp(Phi*U, u, s, bwidth, 10);
            err = vecnorm(approx-x);
            if err < error_tol
                successes = successes + 1;
            end
        end
        success_rate(m/10, s/s_step) = successes;
    end
end
writematrix(success_rate, "m_2.csv")






