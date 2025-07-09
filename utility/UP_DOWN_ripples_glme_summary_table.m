function summary_table = UP_DOWN_ripples_glme_summary_table(output, model_indices, save_path)
% UP_DOWN_ripples_glme_summary_table - Summarize bootstrap stats for GLME models
% Handles models with multiple predictors, and includes beta coefficients.
%
% INPUTS:
%   output         - 1xN struct array from bootstrap results (each with fields: t_stat, pval, R2, b, model)
%   model_indices  - vector of model indices to summarize (e.g., [13 14 29])
%   save_path      - optional full path to save CSV (e.g., 'Ripple_Model_Summary.csv')
%
% OUTPUT:
%   summary_table  - MATLAB table with summary statistics per predictor

    if nargin < 3
        save_path = ''; % No save by default
    end

    model_indices = unique(model_indices(:))';
    nBoot = numel(output);

    % Temporary containers
    model_index_all = [];
    model_name_all = {};
    variable_all = {};
    mean_tstat_all = [];
    ci_low_all = [];
    ci_high_all = [];
    pval_all = [];
    R2_all = [];
    mean_beta_all = [];
    beta_ci_low_all = [];
    beta_ci_high_all = [];

    for i = 1:numel(model_indices)
        idx = model_indices(i);

        model_str = output(1).model{idx};
        nVars = numel(output(1).t_stat{idx});  % number of predictors

        % Preallocate bootstrap matrices: [nBoot x nVars]
        t_vals = nan(nBoot, nVars);
        p_vals = nan(nBoot, nVars);
        b_vals = nan(nBoot, nVars);
        R2_vals = nan(nBoot, 1);

        for b = 1:nBoot
            t_vals(b, :) = output(b).t_stat{idx};
            p_vals(b, :) = output(b).pval{idx};
            b_vals(b, :) = output(b).b{idx};
            R2_vals(b) = output(b).R2(idx);  % same for all variables
        end

        % Parse variable names from formula string
        var_names = get_predictor_names(model_str);

        for v = 1:nVars
            model_index_all(end+1,1) = idx;
            model_name_all{end+1,1} = model_str;
            variable_all{end+1,1} = var_names{v};

            mean_tstat_all(end+1,1) = mean(t_vals(:,v), 'omitnan');
            ci = prctile(t_vals(:,v), [2.5 97.5]);
            ci_low_all(end+1,1) = ci(1);
            ci_high_all(end+1,1) = ci(2);
            pval_all(end+1,1) = prctile(p_vals(:,v), 50);
            R2_all(end+1,1) = prctile(R2_vals, 50);

            % Beta stats
            mean_beta_all(end+1,1) = mean(b_vals(:,v), 'omitnan');
            ci_beta = prctile(b_vals(:,v), [2.5 97.5]);
            beta_ci_low_all(end+1,1) = ci_beta(1);
            beta_ci_high_all(end+1,1) = ci_beta(2);
        end
    end

    % Create output table
    summary_table = table( ...
        model_index_all, model_name_all, variable_all, ...
        mean_tstat_all, ci_low_all, ci_high_all, ...
        pval_all, R2_all, ...
        mean_beta_all, beta_ci_low_all, beta_ci_high_all, ...
        'VariableNames', {'model_index', 'model_name', 'variable', ...
                          'mean_tstat', 'ci_low', 'ci_high', 'pval', 'R2', ...
                          'mean_beta', 'beta_ci_low', 'beta_ci_high'});

    % Save to CSV
    if ~isempty(save_path)
        writetable(summary_table, save_path);
    end
end

function predictors = get_predictor_names(model_str)
% Extract predictor names from model string: e.g., "A ~ B + C + D + (1|subjectID)"
    parts = split(model_str, '~');
    if numel(parts) < 2
        predictors = {''};
        return;
    end
    rhs = parts{2};
    rhs = regexprep(rhs, '\(.*?\)', '');         % remove random effects (e.g., (1|subjectID))
    tokens = strtrim(split(rhs, '+'));           % split by '+'
    tokens = tokens(~cellfun(@isempty, tokens)); % remove empty
    predictors = tokens;
end
