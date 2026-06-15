function results = run_regional_level_case_control_analysis(Z_group, Z_HC, group_name, save_path)
    fprintf('\n--- Running statistical analysis: %s vs HC (regional-level) ---\n', group_name);
    
    numRegions = size(Z_HC, 2);
    alpha = 0.05;

    if size(Z_group, 2) ~= numRegions
        error('Dimension mismatch');
    end

    % 1. Initialize result arrays
    p_region = zeros(numRegions, 1);
    d_region = zeros(numRegions, 1);

    % 2. Perform non-parametric tests and calculate Cohen's d for each region
    fprintf('Calculating statistical differences across %d regions...\n', numRegions);
    for i = 1:numRegions
        data_hc = Z_HC(:, i);
        data_group = Z_group(:, i);

        % Mann-Whitney U test
        p_region(i) = ranksum(data_group, data_hc);
        d_region(i) = compute_cohens_d(data_group, data_hc);
    end

    % 3. FDR Correction
    fprintf('Applying Benjamini-Hochberg FDR correction...\n');
    adj_p_region = bh_fdr(p_region);
    adj_p_region = adj_p_region(:);

    % 4. Generate final Cohen's d map with significance mask
    sig_mask = adj_p_region < alpha;
    final_map = d_region;
    final_map(~sig_mask) = 0;

    % 5. Save results
    save_file = fullfile(save_path, ['Regional_Level_EffectSize_' group_name '_vs_HC.mat']);
    save(save_file, 'p_region', 'adj_p_region', 'd_region', 'sig_mask', 'final_map', 'group_name');

    fprintf('Analysis complete! Results successfully saved to:\n%s\n', save_file);

    results = struct();
    results.group_name = group_name;
    results.p_region = p_region;
    results.adj_p_region = adj_p_region;
    results.d_region = d_region;
    results.sig_mask = sig_mask;
    results.final_map = final_map;
    results.save_file = save_file;
end