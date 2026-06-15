function results = run_system_level_case_control_analysis(Z_group, Z_HC, parcel_to_economo, ve_names, group_name, save_path)

    fprintf('\n--- Running statistical analysis: %s vs HC (system-level) ---\n', group_name);

    numClasses = 7;
    alpha = 0.05;

    % 1. Calculate average deviation across 7 cytoarchitectonic classes
    fprintf('Calculating system-level scores for 7 cytoarchitectonic classes...\n');

    ClassScore_HC   = zeros(size(Z_HC, 1), numClasses);
    ClassScore_Group = zeros(size(Z_group, 1), numClasses);

    for k = 1:numClasses
        region_idx = (parcel_to_economo == k);
        ClassScore_HC(:, k)    = mean(Z_HC(:, region_idx), 2);
        ClassScore_Group(:, k) = mean(Z_group(:, region_idx), 2);
    end

    % 2. Statistical testing (Rank-sum) and Effect size (Cohen's d)
    p_Class = zeros(numClasses, 1);
    d_Class = zeros(numClasses, 1);

    for k = 1:numClasses
        p_Class(k) = ranksum(ClassScore_Group(:, k), ClassScore_HC(:, k));
        d_Class(k) = compute_cohens_d(ClassScore_Group(:, k), ClassScore_HC(:, k));
    end

    % 3. FDR Correction
    fprintf('Applying Benjamini-Hochberg FDR correction...\n');
    adj_p_Class = bh_fdr(p_Class);
    adj_p_Class = adj_p_Class(:);

    % 4. Generate map of significant effect sizes
    sig_mask = adj_p_Class < alpha;
    final_system_map = d_Class;
    final_system_map(~sig_mask) = 0;

    fprintf('--- %s vs HC: System-level Cohen''s d ---\n', group_name);
    disp(table((1:numClasses)', ve_names(:), d_Class, p_Class, adj_p_Class, sig_mask, ...
        'VariableNames', {'System_ID', 'System_Name', 'Cohens_d', 'P_raw', 'P_FDR', 'Significant'}));

    % 5. save results
    save_file = fullfile(save_path, ['System_Level_EffectSize_' group_name '_vs_HC.mat']);
    save(save_file, ...
        'ClassScore_HC', 'ClassScore_Group', ...
        'p_Class', 'adj_p_Class', ...
        'd_Class', 'sig_mask', 'final_system_map', ...
        've_names', 'group_name');


    results = struct();
    results.group_name = group_name;
    results.ClassScore_HC = ClassScore_HC;                  
    results.ClassScore_Group = ClassScore_Group;            
    results.p_Class = p_Class;                             
    results.adj_p_Class = adj_p_Class;                     
    results.d_Class = d_Class;                             
    results.sig_mask = sig_mask;                           
    results.final_system_map = final_system_map;           
    results.ve_names = ve_names;                           
    results.save_file = save_file;

    fprintf('Results successfully saved to: %s\n', save_file);
end