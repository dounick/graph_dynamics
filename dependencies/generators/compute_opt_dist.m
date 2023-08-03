function p = compute_opt_dist(Utilde, T, k, regime)
    N = size(Utilde, 1)/T;
    p = zeros(T*N, 1);
    switch regime
        case 2                 
            infNorms = zeros(T*N, 1);
            for t = 1:T
                sumInfs = 0;
                for i = (t-1)*N+1:t*N
                    max = 0;
                    for j = 1:k
                        if abs(Utilde(i,j)) > max
                            max = abs(Utilde(i,j));
                        end
                    end
                    infNorms(i) = max;
                end
                for it = (t-1)*N+1:t*N
                    sumInfs = sumInfs + infNorms(it)^2;
                end
                for it = (t-1)*N+1:t*N
                    p(it) = (infNorms(it)^2)/sumInfs;
                end   
            end
        case 3
        infNorms = zeros(T*N, 1);
        sumInfs = 0;
    
        for i = 1:T*N
            max = 0;
            for j = 1:k
                if abs(Utilde(i,j)) > max
                    max = abs(Utilde(i,j));
                end
            end
            infNorms(i) = max;
        end
    
        for it = 1:T*N
            sumInfs = sumInfs + infNorms(it)^2;
        end
        for it = 1:T*N
            p(it) = (infNorms(it)^2)/sumInfs;
        end   
    end
end