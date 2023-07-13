function indices = topk(a, k, error)
    [~, indices] = maxk(a, k,'ComparisonMethod','abs');
    i = k;
    while i > 0
        if abs(a(indices(i))) > error
            break
        else
            i = i - 1;
        end
    end
    indices = indices(1:i);
end