function [est, se, tval, pval] = linear_contrast(B, CovB, c_idx, df)
    if isscalar(c_idx)
        c = zeros(length(B),1);
        c(c_idx) = 1;
    else
        c = c_idx(:);
    end
    est = c' * B;
    se = sqrt(c' * CovB * c);
    tval = est / se;
    pval = 2 * (1 - tcdf(abs(tval), df));
end