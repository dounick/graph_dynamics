function [ratio, p] = compute_bandlimited_dist(Utilde, T, k, regime)
    N = size(Ukt, 1)/T;
    p = zeros(T*N, 1);
    switch regime
        case 2
            ratio = zeros(T, 1);
            for t = 1:T
                twoNorms = zeros(N, 1);
                for i = 1:N
                    twoNorms(i) = vecnorm(Utilde((t-1)*N + i, :))^2;
                end
                sumBest = sum(twoNorms);
                for i = 1:N
                    p((t-1)*N+i) = twoNorms(i)/sumBest;
                end
                ratio(t) = sumBest;
            end
            ratio = ratio/sum(ratio);
        case 3
            p = zeros(T*N, 1);
            ratio = 1;
            for i = 1:T*N
                curr = vecnorm(Utilde(i, 1:k))^2;
                p(i) = curr/k;
            end
    end
end