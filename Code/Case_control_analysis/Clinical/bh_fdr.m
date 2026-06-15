function adj_p = bh_fdr(p)
    valid_idx = ~isnan(p);
    p_valid = p(valid_idx);
    
    [p_sorted, sort_idx] = sort(p_valid);
    m = length(p_valid);
    adj_p_sorted = zeros(m, 1);
    
    for i = m:-1:1
        if i == m
            adj_p_sorted(i) = p_sorted(i);
        else 
            adj_p_sorted(i) = min(adj_p_sorted(i+1), p_sorted(i) * m / i);
        end
    end
    
    adj_p_valid = zeros(m, 1);
    adj_p_valid(sort_idx) = adj_p_sorted;
    
    adj_p = nan(length(p), 1);
    adj_p(valid_idx) = min(adj_p_valid, 1);
end

