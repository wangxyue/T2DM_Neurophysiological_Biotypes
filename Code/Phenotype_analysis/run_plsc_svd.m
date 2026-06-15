function plsc = run_plsc_svd(X, Y)
    R = X' * Y;  
    [U, S, V] = svd(R, 'econ');
    singvals = diag(S);
    explained = (singvals.^2) / sum(singvals.^2);
    
    plsc.U = U;
    plsc.V = V;
    plsc.S = S;
    plsc.singvals = singvals;
    plsc.explained = explained;
    plsc.X_scores = X * U; 
    plsc.Y_scores = Y * V; 
    plsc.crosscov = R;
end