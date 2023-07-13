function approximation = cosamp(Phi, u, s, size, iteration)
    error = 0.0000001;
    approx = zeros(size,1);
    hermPhi = transpose(Phi);
    for k = 1:iteration
        v = u - Phi*approx;
        y = hermPhi * v;
        T = union(topk(y, 2*s, error), topk(approx, size, error));
        b = zeros(size,1);
        b(T) = pinv(Phi(:, T)) * u;
        sHighest = topk(b, s, error);
        approx = zeros(size,1);
        for i = 1:length(sHighest)
            approx(sHighest(i)) = b(sHighest(i));
        end
        approximation = approx;
    end
end
