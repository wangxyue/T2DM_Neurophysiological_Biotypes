function bestAcc = cluster_agreement(label1, label2, k)
    % Evaluates the maximum agreement between two sets of cluster labels 
    % by solving the linear assignment problem (Hungarian algorithm).
    
    W = zeros(k, k);
    for i = 1:k
        for j = 1:k
            W(i,j) = sum(label1 == i & label2 == j);
        end
    end
    
    Cost = -W;
   
    assignment = matchpairs(Cost, 10000);
    
    correct_matches = 0;
    for i = 1:size(assignment, 1)
        correct_matches = correct_matches + W(assignment(i,1), assignment(i,2));
    end
    bestAcc = correct_matches / length(label1);
end