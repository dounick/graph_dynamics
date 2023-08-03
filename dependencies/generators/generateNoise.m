function e = generateNoise(n, amp)
    e = zeros(n, 1);
    for i = 1:n
        e(i) = unifrnd(-1, 1);
    end
    normE = vecnorm(e);
    for i = 1:n
        e(i) = e(i)*amp/normE;
    end
end
