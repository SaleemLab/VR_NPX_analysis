function output = predict_UP_DOWN_V1_MUA_by_bilateral_ripples(ipsi_V1_MUA,contra_V1_MUA,ipsi_HPC_MUA,ipsi_ripples,UP_DOWN_info,varargin)


% SETTINGS
p = inputParser;
addParameter(p, 'nPerm', 1000);
addParameter(p, 'nBoot', 1000);
addParameter(p, 'output', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
addParameter(p, 'session_id', 0);
addParameter(p, 'UP_DOWN_lag', []);
addParameter(p, 'UP_DOWN_index', []);
addParameter(p, 'time_window', -0.1);
% addParameter(p, 'V1_time_window', 0.1);
% addParameter(p, 'ripples_info', []);
% addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
output = p.Results.output;
subject_id = p.Results.subject_id;
session_id = p.Results.session_id;
UP_DOWN_lag = p.Results.UP_DOWN_lag;
UP_DOWN_index = p.Results.UP_DOWN_index;
plot_option = p.Results.plot_option;
% ripples_info = p.Results.ripples_info;
% ripples_power = p.Results.ripples_power;
time_window = p.Results.time_window;
% V1_time_window = p.Results.V1_time_window;



% Time settings
time_windows = [-1,1];
time_bin = 0.01;
timebin_edge = time_windows(1):time_bin:time_windows(end);
timeVec = timebin_edge(1:end-1) + time_bin/2; % centers

time_bin = 0.02;
timebin_edge_prob = time_windows(1):time_bin:time_windows(end);
timeVec_prob = timebin_edge_prob(1:end-1) + time_bin/2; % centers


nTimes = length(timeVec);

%% 1. Linear regression
key_window_idx2 = find(timeVec < 0 & timeVec >= time_window);
key_window_idx1 = find(timeVec > 0 & timeVec <= 0.1);

ipsiV1_key = max(ipsi_V1_MUA(:, key_window_idx2)')'-min(ipsi_V1_MUA(:, key_window_idx1)')' ;
contraV1_key = max(contra_V1_MUA(:, key_window_idx2)')' - min(contra_V1_MUA(:, key_window_idx1)')';
% ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')' -min(ipsi_HPC_MUA(:, key_window_idx1)')' ;
% contraHPC_key = max(contra_HPC_MUA(:, key_window_idx2)')'-min(contra_HPC_MUA(:, key_window_idx1)')' ;

ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')' ;
% ipsiV1_key = mean(ipsi_V1_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_V1_MUA(:, key_window_idx1)',"omitnan")' ;
% contraV1_key = mean(contra_V1_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_V1_MUA(:, key_window_idx1)',"omitnan")';
% ipsiHPC_key = mean(ipsi_HPC_MUA(:, key_window_idx2)',"omitnan")' ;
% contraHPC_key = mean(contra_HPC_MUA(:, key_window_idx2)',"omitnan")';
% ipsiHPC_key = mean(ipsi_HPC_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_HPC_MUA(:, key_window_idx1)',"omitnan")' ;
% contraHPC_key = mean(contra_HPC_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_HPC_MUA(:, key_window_idx1)',"omitnan")';
    % combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

%%%% HPC ripples
key_window_idx2= find(timeVec_prob < 0 & timeVec_prob >= time_window);
%     key_window_idx2 = find(timeVec_prob >= 0 & timeVec_prob <= 0.1);

ipsiRipples_key = sum(ipsi_ripples(:,key_window_idx2),2)>0 ;
ipsiV1_cat = ipsiV1_key(:);
contraV1_cat = contraV1_key(:);
ipsiHPC_cat = ipsiHPC_key(:);
ipsiRipples_cat = ipsiRipples_key(:);
%     combinedV1_cat = combinedV1_key(:);

validIdx = ~(isnan(ipsiV1_cat) | isnan(contraV1_cat) | isnan(ipsiHPC_cat));
ipsiV1_cat = ipsiV1_cat(validIdx);
contraV1_cat = contraV1_cat(validIdx);
ipsiHPC_cat = ipsiHPC_cat(validIdx);
ipsiRipples_cat = ipsiRipples_cat(validIdx);
UP_DOWN_index = 1:length(ipsiRipples_cat);

% UP_DOWN_index = UP_DOWN_index(validIdx);

%     ipsi_ripple_power(ipsiRipples_cat) = ripples_info.ipsi_first_ripples_power_UP(UP_DOWN_index(ipsiRipples_cat))';
%     contra_ripple_power(contraRipples_cat) = ripples_info.contra_first_ripples_power_UP(UP_DOWN_index(contraRipples_cat))';
%     ipsi_ripple_lag(ipsiRipples_cat) = abs(ripples_info.ipsi_first_ripples_lag_UP(UP_DOWN_index(ipsiRipples_cat))');
%     contra_ripple_lag(contraRipples_cat) = abs(ripples_info.contra_first_ripples_lag_UP(UP_DOWN_index(contraRipples_cat))');


ipsi_ripple_power = nan(length(UP_DOWN_index),1);
% contra_ripple_power = nan(length(UP_DOWN_index),1);
ipsi_ripple_duration = nan(length(UP_DOWN_index),1);
% contra_ripple_duration = nan(length(UP_DOWN_index),1);
% ripple_lag = nan(length(DOWN_UP_index),1);
ipsi_ripple_lag = nan(length(UP_DOWN_index),1);
% contra_ripple_lag = nan(length(UP_DOWN_index),1);
ipsi_time_from_last_ripple = nan(length(UP_DOWN_index),1);
% contra_time_from_last_ripple = nan(length(UP_DOWN_index),1);
ipsi_Delta_peaks_zscore_UD = nan(length(UP_DOWN_index),1);
contra_Delta_peaks_zscore_UD = nan(length(UP_DOWN_index),1);


ipsi_ripple_power(ipsiRipples_cat) = UP_DOWN_info.last_ripples_power_UP(UP_DOWN_index(ipsiRipples_cat))';

ipsi_ripple_duration(ipsiRipples_cat) = UP_DOWN_info.last_ripples_duration_UP(UP_DOWN_index(ipsiRipples_cat))';

ipsi_ripple_lag(ipsiRipples_cat) = abs(UP_DOWN_info.last_ripples_lag_diff_UP(UP_DOWN_index(ipsiRipples_cat))');

ipsi_time_from_last_ripple(ipsiRipples_cat) = abs(UP_DOWN_info.time_from_last_ripples_UP(UP_DOWN_index(ipsiRipples_cat))');

ipsi_Delta_peaks_zscore_UD = (UP_DOWN_info.SWpeakmag_UD(UP_DOWN_index)');
contra_Delta_peaks_zscore_UD = (UP_DOWN_info.contra_Delta_peaks_zscore_UD(UP_DOWN_index)');
DOWN_duration = UP_DOWN_info.next_DOWN_duration(UP_DOWN_index);
ipsi_spindles = double(sum(UP_DOWN_info.ipsi_spindles_DOWN(UP_DOWN_index,:),2,'omitnan')>0);
contra_spindles = double(sum(UP_DOWN_info.contra_spindles_DOWN(UP_DOWN_index,:),2,'omitnan')>0);


ripples_index = find(ipsiRipples_cat==1);

% [b, logl, H, stats] = coxphfit(UP_DOWN_lag(index), ipsi_time_from_last_ripple(index), 'Strata', subject_id(index));
if isempty(output)
clear output
%%%%%%%% All other models
    parfor iBoot = 1:1000
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        % index = randsample(s, ripples_index, size(ripples_index,1), true);
        % index = randsample(s, ripples_index, size(ripples_index,1), true);
        index = randsample(s, size(UP_DOWN_index,2), size(UP_DOWN_index,2), true);
        tbl = table(normalize(ipsi_ripple_duration(index)), ...
            normalize(ipsi_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), ...
            ipsiRipples_cat(index), ...
            zscore(ipsiHPC_cat(index)), ...
            zscore(ipsiV1_cat(index)), ...
            zscore(contraV1_cat(index)), ...
            normalize(ipsi_Delta_peaks_zscore_UD(index))', ...
            normalize(contra_Delta_peaks_zscore_UD(index))', ...
            normalize(DOWN_duration(index)), ...
            categorical(ipsi_spindles(index)), ... % Categorical for Grand Average coding
            categorical(contra_spindles(index)), ...
            categorical(subject_id(index)), ...
            categorical(session_id(index)), ...    % Added sessionID
            'VariableNames', {'ipsiDuration', 'ipsiLag', 'ipsiPower', 'ipsiOccurance', ...
            'ipsiHPC', 'ipsiV1', 'contraV1', 'ipsiDelta', 'contraDelta', ...
            'downDuration', 'ipsiSpindle', 'contraSpindle', ...
            'subjectID', 'sessionID'});

        modelList = {
            % Predicting ipsiV1
            'ipsiV1 ~ ipsiPower + (1|subjectID) + (1|sessionID)';
            'ipsiV1 ~ ipsiHPC + (1|subjectID) + (1|sessionID)';
            'ipsiV1 ~ ipsiOccurance + (1|subjectID) + (1|sessionID)';
            'ipsiV1 ~ ipsiDuration + (1|subjectID) + (1|sessionID)';

            % Predicting ipsiDelta
            'ipsiDelta ~ ipsiHPC + (1|subjectID) + (1|sessionID)';
            'ipsiDelta ~ ipsiPower + (1|subjectID) + (1|sessionID)';
            'ipsiDelta ~ ipsiOccurance + (1|subjectID) + (1|sessionID)';
            'ipsiDelta ~ ipsiDuration + (1|subjectID) + (1|sessionID)';
            };

        nModels = numel(modelList);

        % Preallocate local outputs
        local_b = cell(nModels,1);
        local_R2 = zeros(nModels,1);
        local_tstat = cell(nModels,1);
        local_pval = cell(nModels,1);
        local_model = modelList; % store the model formulas
        local_variable = cell(nModels,1);

        for m = 1:nModels
            glme = fitglme(tbl, modelList{m});

            % Save results
            local_b{m} = glme.Coefficients.Estimate(2:end);
            local_R2(m) = glme.Rsquared.Ordinary;
            local_tstat{m} = glme.Coefficients.tStat(2:end);
            local_pval{m} = glme.Coefficients.pValue(2:end);

            % Save variable names -- exclude intercept
            local_variable{m} = glme.CoefficientNames(2:end);
        end

        % Save into output struct
        output(iBoot).b = local_b;
        output(iBoot).R2 = local_R2;
        output(iBoot).t_stat = local_tstat;
        output(iBoot).pval = local_pval;
        output(iBoot).model = local_model;
        output(iBoot).variable = local_variable;
        output(iBoot).type = repmat({'Ripples'}, nModels, 1);
        toc
    end


    % 
    % %%%% Ripple occurance
    % parfor iBoot = 1:1000
    %     tic
    %     s = RandStream('philox4x32_10', 'Seed', iBoot);
    %     index = randsample(s, size(UP_DOWN_index,2), size(UP_DOWN_index,2), true);
    %     %         index = intersect(index,ripples_index)
    %     tbl = table(normalize(ipsi_ripple_duration(index)), ...
    %         normalize(ipsi_ripple_lag(index)), ...
    %         normalize(ipsi_ripple_power(index)), ...
    %         ipsiRipples_cat(index), ...
    %         zscore(ipsiHPC_cat(index)), ...
    %         zscore(ipsiV1_cat(index)), ...
    %         zscore(contraV1_cat(index)), ...
    %         normalize(ipsi_Delta_peaks_zscore_UD(index))', ...
    %         normalize(contra_Delta_peaks_zscore_UD(index))', ...
    %         normalize(DOWN_duration(index)), ...
    %         categorical(ipsi_spindles(index)), ... % Categorical for Grand Average coding
    %         categorical(contra_spindles(index)), ...
    %         categorical(subject_id(index)), ...
    %         categorical(session_id(index)), ...    % Added sessionID
    %         'VariableNames', {'ipsiDuration', 'ipsiLag', 'ipsiPower', 'ipsiOccurance', ...
    %         'ipsiHPC', 'ipsiV1', 'contraV1', 'ipsiDelta', 'contraDelta', ...
    %         'downDuration', 'ipsiSpindle', 'contraSpindle', ...
    %         'subjectID', 'sessionID'});
    % 
    %     modelList = {
    %         'ipsiV1 ~ ipsiOccurance + (1|subjectID) + (1|sessionID)';
    %         'ipsiV1 ~ ipsiHPC + (1|subjectID) + (1|sessionID)';
    %         'ipsiDelta ~ ipsiOccurance + (1|subjectID) + (1|sessionID)';
    %         'ipsiDelta ~ ipsiHPC + (1|subjectID) + (1|sessionID)';
    %         };
    % 
    %     nModels = numel(modelList);
    % 
    %     % Preallocate local outputs
    %     local_b = cell(nModels,1);
    %     local_R2 = zeros(nModels,1);
    %     local_tstat = cell(nModels,1);
    %     local_pval = cell(nModels,1);
    %     local_model = modelList; % store model formulas
    %     local_variable = cell(nModels,1);
    % 
    %     for m = 1:nModels
    %         glme = fitglme(tbl, modelList{m});
    % 
    %         % Save results
    %         local_b{m} = glme.Coefficients.Estimate(2:end);
    %         local_R2(m) = glme.Rsquared.Ordinary;
    %         local_tstat{m} = glme.Coefficients.tStat(2:end);
    %         local_pval{m} = glme.Coefficients.pValue(2:end);
    % 
    %         % Save variable names -- exclude intercept
    %         local_variable{m} = glme.CoefficientNames(2:end);
    %     end
    % 
    %     % Save into output struct
    %     output2(iBoot).b = local_b;
    %     output2(iBoot).R2 = local_R2;
    %     output2(iBoot).t_stat = local_tstat;
    %     output2(iBoot).pval = local_pval;
    %     output2(iBoot).model = local_model;
    %     output2(iBoot).variable = local_variable;
    %     output2(iBoot).type = repmat({'All'}, nModels, 1);
    %     toc
    % end

   
    % nBoot = numel(output2); % number of bootstraps
    % 
    % % Get field names automatically
    % fieldNames = fieldnames(output);
    % 
    % % Initialize combined structure
    % combined_output = struct();
    % 
    % for iBoot = 1:nBoot
    %     for f = 1:numel(fieldNames)
    %         thisField = fieldNames{f};
    % 
    %         % Combine the fields
    %         combined_output(iBoot).(thisField) = [output(iBoot).(thisField); output2(iBoot).(thisField)];
    %     end
    % end
    % 
    % output = combined_output; clear combined_output output2
end

% 
% 
if plot_option == 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% V1 depression by HPC activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    index = find(contains(output(1).model,'ipsiV1 ~ ipsiHPC + (1'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(ipsi_b);
    ci_contra = prctile(ipsi_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    % ------------------ Plot Setup ------------------
    fig = figure('Color', 'w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name = 'V1 depression predicted by HPC excitation';
    customColors = [0,90,50; 74,20,134]/256;
    x_pos = [1 2];

    % ------------------ Scatter: ipsi ------------------
    nexttile
    x = ipsiHPC_cat(ripples_index); x= x(:);
    y = ipsiV1_cat(ripples_index);y = y(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.05)
    xlabel('ipsi HPC excitation'); ylabel('ipsi V1 depression')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12); hold on

    % Visual regression line (polyfit)
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'k', 'LineWidth', 2)

    % Annotation using bootstrapped stats
    medR2 = prctile(R2_ipsi,50); medP = prctile(pval_ipsi,50);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.3f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    xlim([-2 15])
    ylim([-2 10])

    % ------------------ Scatter: contra ------------------
    nexttile
    nexttile
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
        if i ==1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    xlim([0.5 2.5]); xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b'); title('HPC excitation')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off






    %%%%%%%%%%%%%%%%
    % ------------------ Extract bootstrapped stats ------------------
    index = find(contains(output(1).model,'ipsiV1 ~ ipsiPower + (1'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(ipsi_b);
    ci_contra = prctile(ipsi_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];


    customColors = [0,90,50; 74,20,134]/256;
    x_pos = [1 2];

    fig = figure('Color','w');
    fig.Position = [350 59 800 930];
    fig.Name = 'V1 depression predicted by ripples';

    % ------------------ Scatter: ripple power vs V1 depression ------------------
    nexttile
    x = ipsi_ripple_power(:); y = ipsiV1_cat(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2)
    xlabel('ipsi ripple power'); ylabel('ipsi V1 depression')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);
    medR2 = prctile(R2_ipsi,50); medP = prctile(pval_ipsi,50);
    % medR2 = median(R2_ipsi_power); medP = median(pval_ipsi_power);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.3f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    ylim([-1 10])


    nexttile
    nexttile
    histogram(ipsiV1_cat(ipsiRipples_cat == 1), -1:0.1:6, 'Normalization', 'probability'); hold on
    histogram(ipsiV1_cat(ipsiRipples_cat == 0), -1:0.1:6, 'Normalization', 'probability');
    [~, pVal_ks_ipsi] = kstest2(ipsiV1_cat(ipsiRipples_cat == 1), ipsiV1_cat(ipsiRipples_cat == 0));
    legend('With ipsi ripple', 'Without ipsi ripple', 'Box','off')
    xlabel('ipsi V1 depression'); ylabel('Proportion');
    text(6, 0.04, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    nexttile
    nexttile
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
        if i ==1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    xticks(x_pos); xticklabels({'ipsi','contra'});
    ylabel('Bootstrapped b'); title('Ripple power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    % ------------------ Bar plot: t-statistics (ripple occurrence) ------------------
    index = find(contains(output(1).model,'ipsiV1 ~ ipsiOccurance + (1'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(ipsi_b);
    ci_contra = prctile(ipsi_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    nexttile
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
        if i == 1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        % else
        %     text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_contra,50)), ...
        %         'FontSize', 10, 'Color', 'k');
        end
    end
    xticks(x_pos); xticklabels({'ipsi','contra'});
    ylabel('Bootstrapped b'); title('Ripple occurrence')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% V1 delta power predicted by HPC excitation
    index = find(contains(output(1).model,'ipsiDelta ~ ipsiHPC + (1'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(ipsi_b);
    ci_contra = prctile(ipsi_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];



    fig = figure('Color', 'w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name = 'V1 delta power predicted by HPC excitation';
    customColors = [0,90,50; 74,20,134]/256;

    % IPSI
    nexttile
    x = ipsiHPC_cat(ripples_index);x = x(:);
    y = ipsi_Delta_peaks_zscore_UD(ripples_index); y = y(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.05)
    xlabel('ipsi HPC excitation'); ylabel('ipsi V1 delta power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    valid = ~isnan(x) & ~isnan(y);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(polyfit(x(valid), y(valid), 1), x_fit);
    % medR2 = median(R2_ipsi); medP = median(pval_ipsi);
    medR2 = prctile(R2_ipsi,50); medP = prctile(pval_ipsi,50);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.3f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)

    % CONTRA
    nexttile
    % T-stat bar
    nexttile
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5)
        if i == 1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b'); title('HPC excitation')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% V1 delta power predicted by ripple features
    index = find(contains(output(1).model,'ipsiDelta ~ ipsiPower + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(ipsi_b);
    ci_contra = prctile(ipsi_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    fig = figure('Color','w');
    fig.Position = [350 59 800 930];
    fig.Name = 'V1 delta power predicted by ripples';

    % Scatter: ripple power
    nexttile
    x = ipsi_ripple_power(:); y = ipsi_Delta_peaks_zscore_UD(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2)
    xlabel('ipsi ripple power'); ylabel('ipsi V1 delta power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    valid = ~isnan(x) & ~isnan(y);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(polyfit(x(valid), y(valid), 1), x_fit);
    R2 = prctile(R2_ipsi, 50); p = prctile(pval_ipsi, 50);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.3f\np = %.2e', R2, p), ...
        'FontSize', 10, 'Color', col)

    nexttile
    % Histogram: ripple occurrence
    nexttile
    histogram(ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 1), 'BinWidth', 0.1, 'Normalization','probability'); hold on
    histogram(ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 0), 'BinWidth', 0.1, 'Normalization','probability');
    [~, pVal_ks_ipsi] = kstest2(ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 1), ...
        ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 0));
    legend('With ipsi ripple','Without','Box','off')
    xlabel('ipsi V1 delta power'); ylabel('Proportion');
    text(min(ipsi_Delta_peaks_zscore_UD), 0.08, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize',10)
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off
    xlim([0 5])

    nexttile
    % Bar: t-stats
    nexttile
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5)
        if i == 1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b'); title('Ripple power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    nexttile
    'ipsiDelta ~ ipsiOccurance + (1|subjectID)+ (1|sessionID)'

    index = find(contains(output(1).model,'ipsiDelta ~ ipsiOccurance + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(ipsi_b);
    ci_contra = prctile(ipsi_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5)
        if i == 1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b'); title('Ripple occurance')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


end
end
