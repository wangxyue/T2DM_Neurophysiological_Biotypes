function perm_results = permutation_test_plsc(X, Y, nPerm)
    n = size(X,1);
    nLV = min(size(X,2), size(Y,2));
    perm_singvals = zeros(nPerm, nLV);

    for i = 1:nPerm
        perm_idx = randperm(n);
        Y_perm = Y(perm_idx, :); 
        s = svd(X' * Y_perm, 'econ');
        perm_singvals(i, 1:length(s)) = s(:)';
    end
    perm_results.perm_singvals = perm_singvals;
end