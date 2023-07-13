function [ratio, p] = compute_sparsity_dist(Ukt, T, s, regime)
    N = size(Ukt, 1)/T;
    p = zeros(T*N, 1);
    switch(regime)
        case 1
            ratio = zeros(T, 1);
            for t = 1:T
                twoNorms = zeros(N, 1);
                for i = 1:N
                    twoNorms(i) = vecnorm(Ukt((t-1)*N + i, :))^2;
                end
                sumBest = sum(twoNorms);
                for i = 1:N
                    p((t-1)*N+i) = twoNorms(i)/sumBest;
                end
                ratio(t) = sumBest;
            end
            ratio = ratio/sum(ratio);
        case 2
            ratio = zeros(T, 1);
            for t = 1:T
                twoNorms = zeros(N, 1);
                for i = 1:N
                    bestNums = maxk(Ukt((t-1)*N + i, :), s, 'ComparisonMethod', 'abs');
                    twoNorms(i) = vecnorm(bestNums)^2;
                end
                sumBest = sum(twoNorms);
                for i = 1:N
                    p((t-1)*N+i) = twoNorms(i)/sumBest;
                end
                ratio(t) = sumBest;
            end
            ratio = ratio/sum(ratio);
        case 3
            ratio = 1;
            twoNorms = zeros(T*N, 1);
            for i = 1:T*N
                bestNums = maxk(Ukt(i, :), s, 'ComparisonMethod', 'abs');
                twoNorms(i) = vecnorm(bestNums)^2;
            end
            sumBest = sum(twoNorms);
            for i = 1:T*N
                p(i) = twoNorms(i)/sumBest;
            end
    end