function outputStruct = predict_HPC_MUA_DOWN_UP(ipsi_V1_MUA, contra_V1_MUA, ipsi_HPC_MUA, contra_HPC_MUA, varargin)
% predict_HPC_MUA_DOWN_UP
% Predict hippocampal MUA during DOWN-UP transitions using ipsi and contra V1 MUA
%
% OUTPUT:
%   outputStruct : structure containing all betas, permutation p-values, and bootstrap CIs

% SETTINGS
p = inputParser;
addParameter(p, 'nPerm', 1000);
addParameter(p, 'nBoot', 1000);
addParameter(p, 'outputStruct', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
addParameter(p, 'DOWN_UP_lag', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
outputStruct = p.Results.outputStruct;
subject_id = p.Results.subject_id;
DOWN_UP_lag = p.Results.DOWN_UP_lag;

if isempty(outputStruct)
    % Time settings
    time_windows = [-1,1];
    time_bin = 0.01;
    timebin_edge = time_windows(1):time_bin:time_windows(end);
    timeVec = timebin_edge(1:end-1) + time_bin/2; % centers

    nTimes = length(timeVec);

    %% 1. Non-time-resolved GLM [-0.1 0.2s]
%     key_window_idx = find(timeVec >= 0 & timeVec <= 0.05);
    key_window_idx = (timeVec >= -0.1) & (timeVec <= 0.1);
% 
%     ipsiV1_key = ipsi_V1_MUA(:, key_window_idx);
%     contraV1_key = contra_V1_MUA(:, key_window_idx);
%     ipsiHPC_key = ipsi_HPC_MUA(:, key_window_idx);
%     contraHPC_key = contra_HPC_MUA(:, key_window_idx);


%     ipsi_V1_MUA_norm = (ipsi_V1_MUA - min(ipsi_V1_MUA(:))) ./ (prctile(ipsi_V1_MUA(:),99)-min(ipsi_V1_MUA(:)));
%     contra_V1_MUA_norm = (contra_V1_MUA - min(contra_V1_MUA(:))) ./ (prctile(contra_V1_MUA(:),99)-min(contra_V1_MUA(:)));
%     ipsi_HPC_MUA_norm = (ipsi_HPC_MUA - min(ipsi_HPC_MUA(:))) ./ (prctile(ipsi_HPC_MUA(:),99)-min(ipsi_HPC_MUA(:)));
%     contra_HPC_MUA_norm = (contra_HPC_MUA - min(contra_HPC_MUA(:))) ./ (prctile(contra_HPC_MUA(:),99)-min(contra_HPC_MUA(:)));

    ipsiV1_key = ipsi_V1_MUA(:, key_window_idx(end)) - ipsi_V1_MUA(:, key_window_idx(1));
    contraV1_key = contra_V1_MUA(:, key_window_idx(end)) - contra_V1_MUA(:, key_window_idx(1));
%     combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

    ipsiHPC_key = ipsi_HPC_MUA(:, key_window_idx(end)) - ipsi_HPC_MUA(:, key_window_idx(1));
    contraHPC_key = contra_HPC_MUA(:, key_window_idx(end)) - contra_HPC_MUA(:, key_window_idx(1));

    ipsiV1_cat = ipsiV1_key(:);
    contraV1_cat = contraV1_key(:);
    ipsiHPC_cat = ipsiHPC_key(:);
    contraHPC_cat = contraHPC_key(:);
%     combinedV1_cat = combinedV1_key(:);

    validIdx = ~(isnan(ipsiV1_cat) | isnan(contraV1_cat) | isnan(ipsiHPC_cat) | isnan(contraHPC_cat));
    ipsiV1_cat = ipsiV1_cat(validIdx);
    contraV1_cat = contraV1_cat(validIdx);
    ipsiHPC_cat = ipsiHPC_cat(validIdx);
    contraHPC_cat = contraHPC_cat(validIdx);

    % Predictors
    X_raw = DOWN_UP_lag;
%     X = X_raw;

    X_raw = [ipsiV1_cat, contraV1_cat];
    X = (X_raw - mean(X_raw,1)) ./ std(X_raw,[],1); % z-score
    ipsiHPC_cat = (ipsiHPC_cat - mean(ipsiHPC_cat,1)) ./ std(ipsiHPC_cat,[],1); % z-score
    contraHPC_cat = (contraHPC_cat - mean(contraHPC_cat,1)) ./ std(contraHPC_cat,[],1); % z-score
    
    %% 2. Non-time-resolved GLM Fit
    [b_full_ipsi, ~, stats_full_ipsi] = glmfit(X, ipsiHPC_cat, 'normal');
    [b_full_contra, ~, stats_full_contra] = glmfit(X, contraHPC_cat, 'normal');

    disp('Non-time-resolved GLM Results: ipsiHPC');
    disp(table(stats_full_ipsi.beta(2:end), stats_full_ipsi.p(2:end), 'VariableNames', {'Beta','pVal'}, 'RowNames', {'ipsiV1','contraV1','ipsi+contraV1'}));

    disp('Non-time-resolved GLM Results: contraHPC');
    disp(table(stats_full_contra.beta(2:end), stats_full_contra.p(2:end), 'VariableNames', {'Beta','pVal'}, 'RowNames', {'ipsiV1','contraV1','ipsi+contraV1'}));

    %% 3. Non-time-resolved Permutation Test + Bootstrap
    beta_null_ipsi = nan(nPerm, 3);
    beta_null_contra = nan(nPerm, 3);
    beta_boot_ipsi = nan(nBoot, 3);
    beta_boot_contra = nan(nBoot, 3);

    n = length(ipsiHPC_cat);

    for i = 1:nPerm
        idx_perm = randperm(n);
        b_perm = glmfit(X, ipsiHPC_cat(idx_perm), 'normal');
        beta_null_ipsi(i,:) = b_perm(2:end)';

        b_perm = glmfit(X, contraHPC_cat(idx_perm), 'normal');
        beta_null_contra(i,:) = b_perm(2:end)';
    end

    for i = 1:nBoot
        idx_boot = randsample(n, n, true);
        b_boot = glmfit(X(idx_boot,:), ipsiHPC_cat(idx_boot), 'normal');
        beta_boot_ipsi(i,:) = b_boot(2:end)';

        b_boot = glmfit(X(idx_boot,:), contraHPC_cat(idx_boot), 'normal');
        beta_boot_contra(i,:) = b_boot(2:end)';
    end

    % Non-time-resolved permutation p-values
    pval_perm_ipsi = mean(abs(beta_null_ipsi) >= abs(repmat(b_full_ipsi(2:end)', nPerm,1)));
    pval_perm_contra = mean(abs(beta_null_contra) >= abs(repmat(b_full_contra(2:end)', nPerm,1)));

    CI_boot_ipsi = prctile(beta_boot_ipsi, [2.5 97.5]);
    CI_boot_contra = prctile(beta_boot_contra, [2.5 97.5]);

    %% 4. Time-resolved Regression, Permutation and Bootstrap
    beta_time_ipsi = nan(nTimes,3);
    beta_time_contra = nan(nTimes,3);
    beta_null_time_ipsi = nan(nTimes,3,nPerm);
    beta_null_time_contra = nan(nTimes,3,nPerm);
    beta_boot_time_ipsi = nan(nTimes,3,nBoot);
    beta_boot_time_contra = nan(nTimes,3,nBoot);

    for t = 1:nTimes
        X1 = ipsi_V1_MUA(:,t);
        X2 = contra_V1_MUA(:,t);
        Y_ipsi = ipsi_HPC_MUA(:,t);
        Y_contra = contra_HPC_MUA(:,t);

        valid = ~(isnan(X1) | isnan(X2) | isnan(Y_ipsi) | isnan(Y_contra));

        if sum(valid) > 10
            X_stack = [X1(valid), X2(valid), X1(valid)+X2(valid)];
            X_z = (X_stack - mean(X_stack,1)) ./ std(X_stack,[],1);

            B = regress(Y_ipsi(valid), [ones(sum(valid),1), X_z]);
            beta_time_ipsi(t,:) = B(2:end);

            B = regress(Y_contra(valid), [ones(sum(valid),1), X_z]);
            beta_time_contra(t,:) = B(2:end);

            for i = 1:nPerm
                idx_perm = randperm(sum(valid));
                B = regress(Y_ipsi(valid(idx_perm)), [ones(sum(valid),1), X_z]);
                beta_null_time_ipsi(t,:,i) = B(2:end);

                B = regress(Y_contra(valid(idx_perm)), [ones(sum(valid),1), X_z]);
                beta_null_time_contra(t,:,i) = B(2:end);
            end

            for i = 1:nBoot
                idx_boot = randsample(sum(valid), sum(valid), true);
                B = regress(Y_ipsi(valid(idx_boot)), [ones(sum(valid),1), X_z(idx_boot,:)]);
                beta_boot_time_ipsi(t,:,i) = B(2:end);

                B = regress(Y_contra(valid(idx_boot)), [ones(sum(valid),1), X_z(idx_boot,:)]);
                beta_boot_time_contra(t,:,i) = B(2:end);
            end
        end
    end

    %% 5. Organize outputs into structure
    outputStruct.ipsiHPC.betas = stats_full_ipsi.beta(2:end);
    outputStruct.ipsiHPC.pvals_perm = pval_perm_ipsi;
    outputStruct.ipsiHPC.CI_boot = CI_boot_ipsi';

    outputStruct.contraHPC.betas = stats_full_contra.beta(2:end);
    outputStruct.contraHPC.pvals_perm = pval_perm_contra;
    outputStruct.contraHPC.CI_boot = CI_boot_contra';

    outputStruct.time_resolved.ipsiHPC.betas = beta_time_ipsi;
    outputStruct.time_resolved.ipsiHPC.CI_boot = beta_boot_time_ipsi;
    outputStruct.time_resolved.ipsiHPC.perm_null = beta_null_time_ipsi;

    outputStruct.time_resolved.contraHPC.betas = beta_time_contra;
    outputStruct.time_resolved.contraHPC.CI_boot = beta_boot_time_contra;
    outputStruct.time_resolved.contraHPC.perm_null = beta_null_time_contra;

    outputStruct.timeVec = timeVec;
    outputStruct.nPerm = nPerm;
    outputStruct.nBoot = nBoot;
end

if plot_option ~=0
    %% 6. PLOTTING
    customColors = [74,20,134; 228,42,168; 0,90,50]/256; % dark purple, magenta, dark green
    col_order = [3,1,2]; % ipsiV1 (green), contraV1 (purple), sum (magenta)
    timeVec = outputStruct.timeVec;
    time_bin = timeVec(2) - timeVec(1);

    titles = {'ipsiHPC prediction', 'contraHPC prediction'};
    data_beta = {outputStruct.time_resolved.ipsiHPC.betas, outputStruct.time_resolved.contraHPC.betas};
    data_boot = {outputStruct.time_resolved.ipsiHPC.CI_boot, outputStruct.time_resolved.contraHPC.CI_boot};
    data_null = {outputStruct.time_resolved.ipsiHPC.perm_null, outputStruct.time_resolved.contraHPC.perm_null};

    figure;
    for hpc = 1:2
        subplot(3,2,2*hpc-1); % Time-resolved beta curves
        hold on;

        max_beta = max(abs(data_beta{hpc}(:)));
        base_y = 1.2 * max_beta;
        offset = 0.05;
        y_sig = [base_y, base_y + offset, base_y + 2*offset];

        for pred = 1:3
            y = data_beta{hpc}(:,pred);
            boot = squeeze(data_boot{hpc}(:,pred,:));

            % Mean ± CI
            LCI = prctile(boot,2.5,2);
            UCI = prctile(boot,97.5,2);

            % Plot main curve
            x = timeVec;
            plot(x, y, 'Color', customColors(col_order(pred),:), 'LineWidth',2);

            % Shaded error
            patch([x fliplr(x)], [UCI' fliplr(LCI')], customColors(col_order(pred),:),...
                'FaceAlpha',0.3, 'EdgeColor','none');

            % Significance
            null_dist = squeeze(data_null{hpc}(:,pred,:));
            real_beta = y;
            pvals = mean(abs(null_dist) >= abs(real_beta),2);
            sig_idx = find(pvals < 0.05);

            if ~isempty(sig_idx)
                for si = 1:length(sig_idx)
                    s_time = x(sig_idx(si));
                    plot([s_time-time_bin/2 s_time+time_bin/2], [y_sig(pred) y_sig(pred)],...
                        '-', 'Color', customColors(col_order(pred),:), 'LineWidth',4);
                end
            end
        end

        xlabel('Time from DOWN-UP transition (s)');
        ylabel('Beta (standardized)');
        title(titles{hpc});
        xline(0,'k--');
        ylim([-1.2*max_beta, y_sig(3)+0.1]);
    end

    %% New subplot for Non-time-resolved bar plot
    for hpc = 1:2
        subplot(3,2,2*hpc); % Bar plots

        hold on;

        if hpc == 1
            real_beta = outputStruct.ipsiHPC.betas;
            real_CI = outputStruct.ipsiHPC.CI_boot;
            shuffled_dist = outputStruct.non_time_resolved_permutation.ipsiHPC; % assume this field exists
        else
            real_beta = outputStruct.contraHPC.betas;
            real_CI = outputStruct.contraHPC.CI_boot;
            shuffled_dist = outputStruct.non_time_resolved_permutation.contraHPC;
        end

        % Calculate shuffled mean and 95% CI
        shuffled_mean = mean(shuffled_dist,1);
        shuffled_CI = prctile(shuffled_dist, [2.5 97.5]);

        x_real = [1,2,3];
        x_shuff = [1.3,2.3,3.3];

        for i = 1:3
            % Real beta bar
            bar(x_real(i), real_beta(i), 'FaceColor', customColors(col_order(i),:), 'FaceAlpha', 0.8);

            % Error bar (bootstrap CI)
            errorbar(x_real(i), real_beta(i), ...
                real_beta(i)-real_CI(i,1), real_CI(i,2)-real_beta(i), 'k', 'LineWidth',1.5);

            % Shuffled beta bar
            bar(x_shuff(i), shuffled_mean(i), 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);

            % Error bar (shuffled CI)
            errorbar(x_shuff(i), shuffled_mean(i), ...
                shuffled_mean(i)-shuffled_CI(1,i), shuffled_CI(2,i)-shuffled_mean(i), 'k', 'LineWidth',1);
        end

        xlim([0.5 3.8]);
        xticks([1 2 3]);
        xticklabels({'ipsiV1','contraV1','ipsi+contraV1'});
        ylabel('Beta');
        title(['Non-time-resolved Betas: ', titles{hpc}(1:end-11)]); % remove "prediction"
    end

    sgtitle('Prediction of HPC MUA from V1 MUA');
end

end

