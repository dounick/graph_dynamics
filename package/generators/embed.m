function vec = embed(x, T, tau)
    n = size(x,1);
    vec = zeros(tau*n, 1);
    powA = eye(n,n);
    for i = 0:tau-1
        vec(i*n + 1 : (i+1)*n, 1) = powA*x;
        powA = powA*T;
    end
end