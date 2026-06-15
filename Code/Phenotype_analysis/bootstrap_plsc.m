function boot_results = bootstrap_plsc(X, Y, U_ref, V_ref, nBoot)
    n = size(X,1);
    nLV = size(V_ref,2);
    boot_U = zeros(size(X,2), nLV, nBoot);
    boot_V = zeros(size(Y,2), nLV, nBoot);

    for b = 1:nBoot
        idx = randsample(n, n, true);
        
        Xb = zscore(X(idx, :));
        Yb = zscore(Y(idx, :));

        [Ub, ~, Vb] = svd(Xb' * Yb, 'econ');

        for lv = 1:nLV
            if corr(Ub(:,lv), U_ref(:,lv), 'rows', 'complete') < 0
                Ub(:,lv) = -Ub(:,lv);
                Vb(:,lv) = -Vb(:,lv);
            end
        end

        boot_U(:,:,b) = Ub;
        boot_V(:,:,b) = Vb;
    end

    X_se = std(boot_U, 0, 3);
    Y_se = std(boot_V, 0, 3);
    
    X_se(X_se == 0) = eps;
    Y_se(Y_se == 0) = eps;
    
    boot_results.boot_U = boot_U;
    boot_results.boot_V = boot_V;
    
    boot_results.X_se = X_se;
    boot_results.Y_se = Y_se;
end