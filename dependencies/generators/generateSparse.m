function signal = generateSparse(n, s, k, U)
    if s > k
        error('s must be smaller than k');
    end
    sampleSet = datasample(1:k, s);
    signal = zeros(n,1);
    for i = 1:length(sampleSet)
        signal = signal + normrnd(0,1)*U(:, sampleSet(i));
    end
    normX = vecnorm(signal);
    for i = 1:length(signal)
        signal(i) = signal(i)/normX;
    end
end
