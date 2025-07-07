function summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, save_path)
% summarize_ripple_models - Extract summary stats for selected models
%
% INPUTS:
%   output         - 1xN struct array from bootstrap results (each with fields: t_stat, pval, R2, model)
%   model_indices  - vector of model indices to summarize (e.g., [13 14 29])
%   save_path      - optional full path to save CSV (e.g., 'Ripple_Model_Summary.csv')
%
% OUTPUT:
%   summary_table  - MATLAB table with summary statistics per model

    if nargin < 3
        save_path = ''; % No save by default
    end

    model_indices = unique(model_indices(:))'; % ensure row vector
    nBoot = numel(output);

    % Preallocate
    summary.model_index = model_indices(:);
    summary.model_name = strings(numel(model_indices),1);
    summary.mean_tstat = nan(numel(model_indices),1);
    summary.ci_low = nan(numel(model_indices),1);
    summary.ci_high = nan(numel(model_indices),1);
    summary.pval = nan(numel(model_indices),1);
    summary.R2 = nan(numel(model_indices),1);

    % Loop through selected models
    for i = 1:numel(model_indices)
        idx = model_indices(i);

        t_vals = nan(nBoot,1);
        p_vals = nan(nBoot,1);
        R2_vals = nan(nBoot,1);

        for b = 1:nBoot
            t_vals(b) = output(b).t_stat{idx};
            p_vals(b) = output(b).pval{idx};
            R2_vals(b) = output(b).R2(idx);
        end

        summary.model_name(i) = output(1).model{idx};
        summary.mean_tstat(i) = mean(t_vals, 'omitnan');
        summary.ci_low(i) = prctile(t_vals, 2.5);
        summary.ci_high(i) = prctile(t_vals, 97.5);
        summary.pval(i) = prctile(p_vals, 50)';
        summary.R2(i) = prctile(R2_vals, 50)';
    end

    % Convert to table
    summary_table = struct2table(summary);

    % Save to CSV if path given
    if ~isempty(save_path)
        writetable(summary_table, save_path);
    end
end
