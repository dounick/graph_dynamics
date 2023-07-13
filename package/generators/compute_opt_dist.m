function p = compute_opt_dist(Utilde, k)
    TN = size(Utilde, 1);
    p = zeros(TN, 1);
    infNorms = zeros(TN, 1);

    sumInfs = 0;

    for i = 1:TN
        max = 0;
        for j = 1:k
            if abs(Utilde(i,j)) > max
                max = abs(Utilde(i,j));
            end
        end
        infNorms(i) = max;
    end

    for it = 1:TN
        sumInfs = sumInfs + infNorms(it)^2;
    end
    for it = 1:TN
        p(it) = (infNorms(it)^2)/sumInfs;
    end   
end