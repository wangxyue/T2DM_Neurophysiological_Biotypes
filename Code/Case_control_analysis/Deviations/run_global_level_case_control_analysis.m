function results = run_global_level_case_control_analysis(Z_group, Z_HC, group_name, save_path)

    group_all_mean = mean(Z_group, 2);
    hc_all_mean    = mean(Z_HC, 2);
    fprintf('\n--- Running statistical analysis: %s vs HC (global) ---\n', group_name);

    % 1. Mann-Whitney U test (Non-parametric)
    [p_raw, ~, stats] = ranksum(group_all_mean, hc_all_mean);
    fprintf('1. Mann-Whitney U test: p = %.6f, z-val = %.3f\n', p_raw, stats.zval);

    % 2. Cohen's d
    d_global = compute_cohens_d(group_all_mean, hc_all_mean);
    fprintf('2. Cohen''s d = %.4f\n', d_global);

    % 3. Permutation Test
    n_perm = 10000;
    rng(123);

    all_data = [group_all_mean; hc_all_mean];
    n_group = length(group_all_mean);

    fprintf('3. Running permutation test (%d permutations)...\n', n_perm);
    obs_diff = abs(mean(group_all_mean) - mean(hc_all_mean));

    count = 0;
    for i = 1:n_perm
        shuffled = all_data(randperm(length(all_data)));
        perm_diff = abs(mean(shuffled(1:n_group)) - mean(shuffled(n_group+1:end)));
        if perm_diff >= obs_diff
            count = count + 1;
        end
    end

    global_p_perm = (count + 1) / (n_perm + 1);
    fprintf('Permutation test P_perm = %.6f\n', global_p_perm);

    save_file = fullfile(save_path, ['Global_Level_EffectSize_' group_name '_vs_HC.mat']);
    save(save_file, ...
        'group_all_mean', 'hc_all_mean', ...
        'p_raw', 'stats', ...
        'd_global', 'global_p_perm');

   
    results = struct();
    results.group_name = group_name;
    results.group_all_mean = group_all_mean;
    results.hc_all_mean = hc_all_mean;
    results.p_raw = p_raw;
    results.zval = stats.zval;
    results.d_global = d_global;
    results.global_p_perm = global_p_perm;
    results.save_file = save_file;

    fprintf('Results successfully saved to: %s\n', save_file);
end
