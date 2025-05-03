function outputStruct = predict_ripples_by_DOWN_UP_V1_MUA(ipsi_V1_MUA, contra_V1_MUA, ipsi_ripples, contra_ripples, varargin)
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
addParameter(p, 'ripples_lag', []);
addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
outputStruct = p.Results.outputStruct;
subject_id = p.Results.subject_id;
DOWN_UP_lag = p.Results.DOWN_UP_lag;
ripples_lag = p.Results.ripples_lag;
ripples_power = p.Results.ripples_power;


if isempty(outputStruct)
    % Time settings
    time_windows = [-1,1];
    time_bin = 0.01;
    timebin_edge = time_windows(1):time_bin:time_windows(end);
    timeVec = timebin_edge(1:end-1) + time_bin/2; % centers

    nTimes = length(timeVec);

    %% 1. Linear regression
    key_window_idx1 = find(timeVec < 0 & timeVec >= -0.05);
    key_window_idx2 = find(timeVec >= 0 & timeVec <= 0.1);

    ipsiV1_key = max(ipsi_V1_MUA(:, key_window_idx2)')'-min(ipsi_V1_MUA(:, key_window_idx1)')' ;
    contraV1_key = max(contra_V1_MUA(:, key_window_idx2)')' - min(contra_V1_MUA(:, key_window_idx1)')';
    %     combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

    ipsiHPC_key = ipsi_ripples(:,key_window_idx2) ;
    contraHPC_key = contra_ripples(:,key_window_idx2);

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

    %     validIdx = find(ipsiHPC_cat>0 | contraHPC_cat>0);
    %
    %     ipsiV1_cat = ipsiV1_cat(validIdx);
    %     contraV1_cat = contraV1_cat(validIdx);
    %     ipsiHPC_cat = ipsiHPC_cat(validIdx);
    %     contraHPC_cat = contraHPC_cat(validIdx);
    %     DOWN_UP_lag = DOWN_UP_lag(validIdx);
    %     subject_id =  subject_id(validIdx);
    % Predictors
    X_raw = DOWN_UP_lag;
    %     X = X_raw;

    %     X_raw = [ipsiV1_cat, mean([ipsiV1_cat contraV1_cat],2)];
    %     X = (X_raw - mean(X_raw,1)) ./ std(X_raw,[],1); % z-score
%     ipsiHPC_cat = (ipsiHPC_cat - mean(ipsiHPC_cat,1)) ./ std(ipsiHPC_cat,[],1); % z-score
%     contraHPC_cat = (contraHPC_cat - mean(contraHPC_cat,1)) ./ std(contraHPC_cat,[],1); % z-score
%     ipsiV1_cat = (ipsiV1_cat - mean(ipsiV1_cat,1)) ./ std(ipsiV1_cat,[],1); % z-score
%     contraV1_cat = (contraV1_cat - mean(contraV1_cat,1)) ./ std(contraV1_cat,[],1); % z-score

    %     HPC_diff = zscore([ipsiHPC_cat-contraHPC_cat]);


    % Create table
    tbl = table(ipsiHPC_cat, contraHPC_cat,ipsiV1_cat,contraV1_cat,zscore(DOWN_UP_lag), categorical(subject_id), 'VariableNames', {'ipsiHPC','contraHPC','ipsiV1','contraV1','DOWN_UP_lag','subjectID'});

    ipsi_R2 = nan(1000,1);
    ipsi_t_stat = nan(1000,1);
    ipsi_pval = nan(1000,1);

    contra_R2 = nan(1000,1);
    contra_t_stat = nan(1000,1);
    contra_pval = nan(1000,1);
    
    tic
    parfor iBoot = 1:1000
        s = RandStream('philox4x32_10','Seed',iBoot);
        index = randsample(s, size(DOWN_UP_lag,1), size(DOWN_UP_lag,1), true);
        tbl = table(zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)),zscore(ipsiV1_cat(index)),zscore(contraV1_cat(index)),zscore(DOWN_UP_lag(index)), categorical(subject_id(index)), 'VariableNames', {'ipsiHPC','contraHPC','ipsiV1','contraV1','DOWN_UP_lag','subjectID'});

        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1*DOWN_UP_lag + contraV1*DOWN_UP_lag + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1 + contraV1 + DOWN_UP_lag + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1  + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ contraV1  + (1|subjectID)');
        ipsi_glme = fitglme(tbl, 'ipsiHPC ~ DOWN_UP_lag  + (1|subjectID)');

        %         ipsi_glme.Rsquared.Ordinary
        ipsi_R2(iBoot) = ipsi_glme.Rsquared.Ordinary;
        ipsi_t_stat(iBoot) = ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
        ipsi_pval(iBoot) = ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
        %         ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'ipsiV1'));
        %         ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'contraV1'));

        %         contra_glme = fitglme(tbl, 'contraHPC ~ ipsiV1*DOWN_UP_lag + contraV1*DOWN_UP_lag + (1|subjectID)');
        %         contra_glme = fitglme(tbl, 'contraHPC ~ ipsiV1 + contraV1 + DOWN_UP_lag + (1|subjectID)');
        contra_glme = fitglme(tbl, 'contraHPC ~ DOWN_UP_lag  + (1|subjectID)');
        contra_R2(iBoot) = contra_glme.Rsquared.Ordinary;
        contra_t_stat(iBoot) = contra_glme.Coefficients.tStat(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));
        contra_pval(iBoot) = contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));
        %         contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'ipsiV1'));
        %         contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'contraV1'));
    end
    toc
%     disp('')
    output.ipsi_R2 =     ipsi_R2;
    output.ipsi_t_stat =    ipsi_t_stat;
    output.ipsi_pval =    ipsi_pval;

    output.contra_R2 =     contra_R2;
    output.contra_t_stat =    contra_t_stat;
    output.contra_pval =     contra_pval;

    %%%% Permutation test
    ipsi_R2 = nan(1000,1);
    ipsi_t_stat = nan(1000,1);
    ipsi_pval = nan(1000,1);

    contra_R2 = nan(1000,1);
    contra_t_stat = nan(1000,1);
    contra_pval = nan(1000,1);
    tic
    parfor iBoot = 1:1000
        s = RandStream('philox4x32_10','Seed',iBoot);
        index = randsample(s, size(DOWN_UP_lag,1), size(DOWN_UP_lag,1), false);
        tbl = table(zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)),zscore(ipsiV1_cat(:)),zscore(contraV1_cat(:)),zscore(DOWN_UP_lag(:)), categorical(subject_id(:)), 'VariableNames', {'ipsiHPC','contraHPC','ipsiV1','contraV1','DOWN_UP_lag','subjectID'});

        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1*DOWN_UP_lag + contraV1*DOWN_UP_lag + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1 + contraV1 + DOWN_UP_lag + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1  + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ contraV1  + (1|subjectID)');
        ipsi_glme = fitglme(tbl, 'ipsiHPC ~ DOWN_UP_lag  + (1|subjectID)');

        %         ipsi_glme.Rsquared.Ordinary
        ipsi_R2(iBoot) = ipsi_glme.Rsquared.Ordinary;
        ipsi_t_stat(iBoot) = ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
        ipsi_pval(iBoot) = ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
        %         ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'ipsiV1'));
        %         ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'contraV1'));

        %         contra_glme = fitglme(tbl, 'contraHPC ~ ipsiV1*DOWN_UP_lag + contraV1*DOWN_UP_lag + (1|subjectID)');
        %         contra_glme = fitglme(tbl, 'contraHPC ~ ipsiV1 + contraV1 + DOWN_UP_lag + (1|subjectID)');
        contra_glme = fitglme(tbl, 'contraHPC ~ DOWN_UP_lag  + (1|subjectID)');
        contra_R2(iBoot) = contra_glme.Rsquared.Ordinary;
        contra_t_stat(iBoot) = contra_glme.Coefficients.tStat(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));
        contra_pval(iBoot) = contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));
        %         contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'ipsiV1'));
        %         contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'contraV1'));
    end
    toc
    output.ipsi_R2_shuffled =   ipsi_R2;
    output.ipsi_t_stat_shuffled =  ipsi_t_stat;
    output.ipsi_pval_shuffled =  ipsi_pval;

    output.contra_R2_shuffled =   contra_R2;
    output.contra_t_stat_shuffled =  contra_t_stat;
    output.contra_pval_shuffled =   contra_pval;




    %%% lag resolved mixed effect regression
    lag_thresholds = prctile(DOWN_UP_lag,[0:10:100])

    ipsi_R2 = nan(length(lag_thresholds)-1,1000);
    ipsi_t_stat = nan(length(lag_thresholds)-1,1000);
    ipsi_pval = nan(length(lag_thresholds)-1,1000);

    contra_R2 = nan(length(lag_thresholds)-1,1000);
    contra_t_stat = nan(length(lag_thresholds)-1,1000);
    contra_pval = nan(length(lag_thresholds)-1,1000);

    for ngroup = 1:length(lag_thresholds)-1
        tic
        % Find indices for this bin
        index = find(DOWN_UP_lag < lag_thresholds(ngroup+1) & DOWN_UP_lag > lag_thresholds(ngroup));

        % Bootstrap arrays
        ipsi_t_boot = nan(1000,1);
        ipsi_R2_boot = nan(1000,1);
        ipsi_pval_boot = nan(1000,1);

        contra_t_boot = nan(1000,1);
        contra_R2_boot = nan(1000,1);
        contra_pval_boot = nan(1000,1);

        % nBins = [];
        % Bootstrap loop
        N = length(index);
        parfor iBoot = 1:1000
            % Use a specific random seed based on iboot
            s = RandStream('philox4x32_10','Seed',iBoot);
            boot_idx = randsample(s, index, N, true);

            % Build sub-table for current bin
            boot_tbl = table(zscore(ipsiHPC_cat(boot_idx)), zscore(contraHPC_cat(boot_idx)), ...
                zscore(ipsiV1_cat(boot_idx)), zscore(contraV1_cat(boot_idx)), ...
                zscore(DOWN_UP_lag(boot_idx)), categorical(subject_id(boot_idx)), ...
                'VariableNames', {'ipsiHPC', 'contraHPC', 'ipsiV1', 'contraV1', 'DOWN_UP_lag', 'subjectID'});

            ipsi_glme = fitglme(boot_tbl, 'ipsiHPC ~ DOWN_UP_lag + (1|subjectID)');
            contra_glme = fitglme(boot_tbl, 'contraHPC ~ DOWN_UP_lag + (1|subjectID)');

            % Store t-statistics
            ipsi_t_boot(iBoot) = ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
            contra_t_boot(iBoot) = contra_glme.Coefficients.tStat(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));

            ipsi_R2_boot(iBoot) = ipsi_glme.Rsquared.Ordinary;
            contra_R2_boot(iBoot) = contra_glme.Rsquared.Ordinary;

            ipsi_pval_boot(iBoot) = ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
            contra_pval_boot(iBoot) = contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));

        end

        % Compute mean and 95% CI for t-statistics
        ipsi_t_stat(ngroup,:) = ipsi_t_boot;
        contra_t_stat(ngroup,:) = contra_t_boot;

        ipsi_R2(ngroup,:) = ipsi_R2_boot;
        contra_R2(ngroup,:) = contra_R2_boot;

        ipsi_pval(ngroup,:) = ipsi_pval_boot;
        contra_pval(ngroup,:) = contra_pval_boot;

        % Save bin center
        lag_bins(ngroup) = mean([lag_thresholds(ngroup) lag_thresholds(ngroup+1)]);
        toc
    end

    output.lag_bins =     lag_bins;
    output.lag_resolved_ipsi_R2 =     ipsi_R2;
    output.lag_resolved_ipsi_t_stat =    ipsi_t_stat;
    output.lag_resolved_ipsi_pval =    ipsi_pval;

    output.lag_resolved_contra_R2 =     contra_R2;
    output.lag_resolved_contra_t_stat =    contra_t_stat;
    output.lag_resolved_contra_pval =     contra_pval;



    %%% Permutation Test for lag resolved mixed effect regression
    lag_thresholds = prctile(DOWN_UP_lag,[0:10:100])

    ipsi_R2 = nan(length(lag_thresholds)-1,1000);
    ipsi_t_stat = nan(length(lag_thresholds)-1,1000);
    ipsi_pval = nan(length(lag_thresholds)-1,1000);

    contra_R2 = nan(length(lag_thresholds)-1,1000);
    contra_t_stat = nan(length(lag_thresholds)-1,1000);
    contra_pval = nan(length(lag_thresholds)-1,1000);

    for ngroup = 1:length(lag_thresholds)-1
        tic
        % Find indices for this bin
        index = find(DOWN_UP_lag < lag_thresholds(ngroup+1) & DOWN_UP_lag > lag_thresholds(ngroup));

        % Bootstrap arrays
        ipsi_t_boot = nan(1000,1);
        ipsi_R2_boot = nan(1000,1);
        ipsi_pval_boot = nan(1000,1);

        contra_t_boot = nan(1000,1);
        contra_R2_boot = nan(1000,1);
        contra_pval_boot = nan(1000,1);

        % nBins = [];
        % Bootstrap loop
        N = length(index);
        parfor iBoot = 1:1000
            % Use a specific random seed based on iboot
            s = RandStream('philox4x32_10','Seed',iBoot);
            boot_idx = randsample(s, index, N, false);

            % Build sub-table for current bin
            boot_tbl = table(zscore(ipsiHPC_cat(boot_idx)), zscore(contraHPC_cat(boot_idx)), ...
                zscore(ipsiV1_cat(index)), zscore(contraV1_cat(index)), ...
                zscore(DOWN_UP_lag(index)), categorical(subject_id(index)), ...
                'VariableNames', {'ipsiHPC', 'contraHPC', 'ipsiV1', 'contraV1', 'DOWN_UP_lag', 'subjectID'});

            ipsi_glme = fitglme(boot_tbl, 'ipsiHPC ~ DOWN_UP_lag + (1|subjectID)');
            contra_glme = fitglme(boot_tbl, 'contraHPC ~ DOWN_UP_lag + (1|subjectID)');

            % Store t-statistics
            ipsi_t_boot(iBoot) = ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
            contra_t_boot(iBoot) = contra_glme.Coefficients.tStat(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));

            ipsi_R2_boot(iBoot) = ipsi_glme.Rsquared.Ordinary;
            contra_R2_boot(iBoot) = contra_glme.Rsquared.Ordinary;

            ipsi_pval_boot(iBoot) = ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
            contra_pval_boot(iBoot) = contra_glme.Coefficients.pValue(strcmp(contra_glme.Coefficients.Name,'DOWN_UP_lag'));

        end

        % Compute mean and 95% CI for t-statistics
        ipsi_t_stat(ngroup,:) = ipsi_t_boot;
        contra_t_stat(ngroup,:) = contra_t_boot;

        ipsi_R2(ngroup,:) = ipsi_R2_boot;
        contra_R2(ngroup,:) = contra_R2_boot;

        ipsi_pval(ngroup,:) = ipsi_pval_boot;
        contra_pval(ngroup,:) = contra_pval_boot;

        % Save bin center
        lag_bins(ngroup) = mean([lag_thresholds(ngroup) lag_thresholds(ngroup+1)]);
        toc
    end

    output.lag_resolved_ipsi_R2_shuffled =     ipsi_R2;
    output.lag_resolved_ipsi_t_stat_shuffled =    ipsi_t_stat;
    output.lag_resolved_ipsi_pval_shuffled =    ipsi_pval;

    output.lag_resolved_contra_R2_shuffled =     contra_R2;
    output.lag_resolved_contra_t_stat_shuffled =    contra_t_stat;
    output.lag_resolved_contra_pval_shuffled =     contra_pval;



end

fig = figure('Color','w');
fig.Position = [350 59 1650 930];
fig.Name ='HPC MUA predicted by DOWN UP synchrony';

customColors = [0,90,50;74,20,134; 228,42,168 ]/256; % dark purple, dark green, magenta

nexttile
scatter(DOWN_UP_lag,ipsiHPC_cat,'filled','MarkerFaceColor',customColors(1,:),'MarkerFaceAlpha',0.1)
ylim([-2 12])
xlabel('DOWN UP transition ipsi-contra lags')
ylabel('ipsi HPC rebound excitation')

set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

nexttile
scatter(DOWN_UP_lag,contraHPC_cat,'filled','MarkerFaceColor',customColors(2,:),'MarkerFaceAlpha',0.1)
ylim([-2 12])
xlabel('DOWN UP transition ipsi-contra lags')
ylabel('contra HPC rebound excitation')
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

nexttile
scatter(ipsiHPC_cat,contraHPC_cat,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.1)
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
xlabel('ipsi HPC rebound excitation')
ylabel('contra HPC rebound excitation')
xlim([-2 12])
ylim([-2 12])


mean_ipsi = mean(output.ipsi_t_stat);
ci_ipsi = prctile(output.ipsi_t_stat, [2.5 97.5]);

mean_ipsi_shuffled = mean(output.ipsi_t_stat_shuffled);
ci_ipsi_shuffled = prctile(output.ipsi_t_stat_shuffled, [2.5 97.5]);

mean_contra = mean(output.contra_t_stat);
ci_contra = prctile(output.contra_t_stat, [2.5 97.5]);

mean_contra_shuffled = mean(output.contra_t_stat_shuffled);
ci_contra_shuffled = prctile(output.contra_t_stat_shuffled, [2.5 97.5]);

% Organize data for bar plot
barData = [mean_ipsi, mean_ipsi_shuffled, mean_contra, mean_contra_shuffled];
lowerCI = [mean_ipsi(1)-ci_ipsi(1), mean_ipsi_shuffled(1)-ci_ipsi_shuffled(1), mean_contra(1)-ci_contra(1), mean_contra_shuffled(1)-ci_contra_shuffled(1)];
upperCI = [mean_ipsi(1)-ci_ipsi(2), mean_ipsi_shuffled(1)-ci_ipsi_shuffled(2), mean_contra(1)-ci_contra(2), mean_contra_shuffled(1)-ci_contra_shuffled(2)];

nexttile
barColors = [0,90,50;65,171,93;74,20,134;128,125,186]/256;
% Define custom x-positions
x_pos = [1, 2, 4, 5]; % spacing

%     b = bar(1:4, barData, 'FaceColor', 'flat');
for i = 1:4
    hold on
    bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', 'CData', barColors(i,:), 'EdgeColor', 'none');
    errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

% Customize plot
xlim([0.5 5.5])
xticks(x_pos)
xticklabels({'ipsi HPC','ipsi HPC shuffled','contra HPC','contra HPC shuffled'})
ylabel('Bootstrapped t-statistic')
title('t-statistic of HPC excitation')
%     legend({'Data','Shuffled'}, 'Box','off')
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
grid off
hold off

end

