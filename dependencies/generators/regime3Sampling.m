function S = regime3Sampling(n, M, T, p)
    if M > T*n
        error('M cannot be larger than T*n');
    end
    randmSamp = datasample(1:T*n, M, 'Weights', p);
    %randmSamp = datasample(1:T*n, T*m);
    S = zeros(M,T*n);
    for j = 1:M
        S(j, randmSamp(j)) = 1;
    end
    %
    %for i = 1:M
    %    S(i, :) = S(i, :) * evolA(:,:,floor((randmSamp(i)-1)/n)+1);
    %end
    S = sparse(S);
end