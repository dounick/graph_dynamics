function c2 = compute_2_coherence(U, s, p)
    n = size(U, 1);

    %compute c2
    c2 = 0;
    for i = 1:n
        sorted_row = sort(U(i,:), 'descend', 'ComparisonMethod', 'abs');
        row_sum = 0;
        for j = 1:s
            row_sum = row_sum + sorted_row(j)^2;
        end
        row_sum = sqrt(row_sum);
        row_coherence = row_sum/sqrt(p(i));
        if row_coherence > c2
            c2 = row_coherence;
        end
    end