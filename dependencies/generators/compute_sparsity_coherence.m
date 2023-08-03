function c = compute_sparsity_coherence(U, s, T, regime, ptype)
    N = size(U, 1)/T;
    %compute coherence
    c = 0;
    switch regime
        case 1
            for t = 1:T
                twoNorms = zeros(N,1);
                p = ptype((t-1)*N + 1:t*N);
                for i = 1:N
                    bestNums = maxk(U((t-1)*N + i, :), s, 'ComparisonMethod', 'abs');
                    twoNorms(i) = max(twoNorms(i), vecnorm(bestNums)^2/(p(i)));
                end
            end
            c = c + sum(twoNorms);
        case 2
            for t = 1:T
                twoNorms = zeros(N, 1);
                for i = 1:N
                    bestNums = maxk(U((t-1)*N + i, :), s, 'ComparisonMethod', 'abs');
                    twoNorms(i) = vecnorm(bestNums)^2;
                end
                if ptype == "opt"
                    c = c + sum(twoNorms);
                elseif ptype == "unif"
                    c = c + N*max(twoNorms);
                end
            end
        case 3
            twoNorms = zeros(T*N, 1);
            for i = 1:T*N
                bestNums = maxk(U(i, :), s, 'ComparisonMethod', 'abs');
                twoNorms(i) = vecnorm(bestNums)^2;
            end
            if ptype == "opt"
                c = sum(twoNorms);
            elseif ptype == "unif"
                c = T*N*max(twoNorms);
            end
    end