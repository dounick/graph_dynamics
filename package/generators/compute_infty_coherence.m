function cinfty = compute_infty_coherence(U, p)
    n = size(U, 1);
    k = size(U, 2);
    %compute cinfty:
    cinfty = 0;
    for i = 1:n
        max_row = 0;
        for j = 1:k
            if abs(U(i,j)) > max_row
                max_row = abs(U(i,j));
            end
        end
        row_coherence = max_row/sqrt(p(i));
        if row_coherence > cinfty
            cinfty = row_coherence;
        end
    end
