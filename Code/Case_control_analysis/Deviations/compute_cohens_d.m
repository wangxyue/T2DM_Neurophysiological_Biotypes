function d = compute_cohens_d(x1, x2)
    n1 = length(x1);
    n2 = length(x2);
    m1 = mean(x1);
    m2 = mean(x2);
    v1 = var(x1);
    v2 = var(x2);
    
    sd_pooled = sqrt(((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2));
    d = (m1 - m2) / sd_pooled;
end