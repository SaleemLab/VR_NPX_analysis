function output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA,contra_V1_MUA,ipsi_HPC_MUA,contra_HPC_MUA,ipsi_ripples,contra_ripples,UP_DOWN_info,varargin)
% SETTINGS
p = inputParser;
addParameter(p, 'nPerm', 1000);
addParameter(p, 'nBoot', 1000);
addParameter(p, 'output', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
addParameter(p, 'session_id', 0);
addParameter(p, 's', 0);
addParameter(p, 'time_window', 0.1);

% addParameter(p, 'DOWN_UP_lag', []);
% addParameter(p, 'DOWN_UP_index', []);
% addParameter(p, 'ripples_info', []);
% addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
output = p.Results.output;
subject_id = p.Results.subject_id;
session_id = p.Results.session_id;
% DOWN_UP_lag = p.Results.DOWN_UP_lag;
% DOWN_UP_index = p.Results.DOWN_UP_index;
plot_option = p.Results.plot_option;
time_window = p.Results.time_window;


% ripples_info = p.Results.ripples_info;
% ripples_power = p.Results.ripples_power;


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
%%%%%% Windows
key_window_idx1 = find(timeVec < 0 & timeVec >= -0.05);

if length(time_window)==2
    key_window_idx2 = find(timeVec > time_window(1) & timeVec <= time_window(end));
else
    key_window_idx2 = find(timeVec > 0 & timeVec <= time_window(end));
end


ipsiV1_key = max(ipsi_V1_MUA(:, key_window_idx2)')'-min(ipsi_V1_MUA(:, key_window_idx1)')' ;
contraV1_key = max(contra_V1_MUA(:, key_window_idx2)')' - min(contra_V1_MUA(:, key_window_idx1)')';
ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')'-min(ipsi_HPC_MUA(:, key_window_idx1)')' ;
contraHPC_key = max(contra_HPC_MUA(:, key_window_idx2)')' - min(contra_HPC_MUA(:, key_window_idx1)')';

% ipsiV1_key = mean(ipsi_V1_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_V1_MUA(:, key_window_idx1)',"omitnan")' ;
% contraV1_key = mean(contra_V1_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_V1_MUA(:, key_window_idx1)',"omitnan")';
% ipsiHPC_key = mean(ipsi_HPC_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_HPC_MUA(:, key_window_idx1)',"omitnan")' ;
% contraHPC_key = mean(contra_HPC_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_HPC_MUA(:, key_window_idx1)',"omitnan")';
%     combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

%%%% HPC ripples
%     key_window_idx2= find(timeVec_prob < 0 & timeVec_prob >= -0.05);
if length(time_window)==2
    key_window_idx2 = find(timeVec_prob > time_window(1) & timeVec_prob <= time_window(end));
else
    key_window_idx2 = find(timeVec_prob > 0 & timeVec_prob <= time_window(end));
end

ipsiRipples_key = sum(ipsi_ripples(:,key_window_idx2),2)>0 ;
contraRipples_key =sum(contra_ripples(:,key_window_idx2),2)>0 ;

ipsiV1_cat = ipsiV1_key(:);
contraV1_cat = contraV1_key(:);
ipsiHPC_cat = ipsiHPC_key(:);
contraHPC_cat = contraHPC_key(:);
ipsiRipples_cat = ipsiRipples_key(:);
contraRipples_cat = contraRipples_key(:);
%     combinedV1_cat = combinedV1_key(:);

% validIdx = ~(isnan(ipsiV1_cat) | isnan(contraV1_cat) | isnan(ipsiHPC_cat) | isnan(contraHPC_cat));
% ipsiV1_cat = ipsiV1_cat(validIdx);
% contraV1_cat = contraV1_cat(validIdx);
% ipsiHPC_cat = ipsiHPC_cat(validIdx);
% contraHPC_cat = contraHPC_cat(validIdx);
% ipsiRipples_cat = ipsiRipples_cat(validIdx);
% contraRipples_cat = contraRipples_cat(validIdx);
DOWN_UP_index = 1:length(ipsi_ripples);

%     ipsi_ripple_power(ipsiRipples_cat) = ripples_info.ipsi_first_ripples_power_UP(UP_DOWN_index(ipsiRipples_cat))';
%     contra_ripple_power(contraRipples_cat) = ripples_info.contra_first_ripples_power_UP(UP_DOWN_index(contraRipples_cat))';
%     ipsi_ripple_lag(ipsiRipples_cat) = abs(ripples_info.ipsi_first_ripples_lag_UP(UP_DOWN_index(ipsiRipples_cat))');
%     contra_ripple_lag(contraRipples_cat) = abs(ripples_info.contra_first_ripples_lag_UP(UP_DOWN_index(contraRipples_cat))');


ipsi_ripple_power = nan(length(DOWN_UP_index),1);
contra_ripple_power = nan(length(DOWN_UP_index),1);
ipsi_ripple_duration = nan(length(DOWN_UP_index),1);
contra_ripple_duration = nan(length(DOWN_UP_index),1);
% ripple_lag = nan(length(DOWN_UP_index),1);
ipsi_ripple_lag = nan(length(DOWN_UP_index),1);
contra_ripple_lag = nan(length(DOWN_UP_index),1);
ipsi_time_from_last_ripple = nan(length(DOWN_UP_index),1);
contra_time_from_last_ripple = nan(length(DOWN_UP_index),1);
ipsi_Delta_peaks_zscore_DU = nan(length(DOWN_UP_index),1);
contra_Delta_peaks_zscore_DU = nan(length(DOWN_UP_index),1);

ipsi_ripple_power(ipsiRipples_cat) = UP_DOWN_info.ipsi_first_ripples_power_UP(DOWN_UP_index(ipsiRipples_cat))';
contra_ripple_power(contraRipples_cat) = UP_DOWN_info.contra_first_ripples_power_UP(DOWN_UP_index(contraRipples_cat))';
ipsi_ripple_duration(ipsiRipples_cat) = UP_DOWN_info.ipsi_first_ripples_duration_UP(DOWN_UP_index(ipsiRipples_cat))';
contra_ripple_duration(contraRipples_cat) = UP_DOWN_info.contra_first_ripples_duration_UP(DOWN_UP_index(contraRipples_cat))';
ipsi_ripple_lag(ipsiRipples_cat) = abs(UP_DOWN_info.ipsi_first_ripples_lag_UP(DOWN_UP_index(ipsiRipples_cat))');
contra_ripple_lag(contraRipples_cat) = abs(UP_DOWN_info.contra_first_ripples_lag_UP(DOWN_UP_index(contraRipples_cat))');
ipsi_time_from_last_ripple(ipsiRipples_cat) = abs(UP_DOWN_info.ipsi_time_from_last_ripples_UP(DOWN_UP_index(ipsiRipples_cat))');
contra_time_from_last_ripple(contraRipples_cat) = abs(UP_DOWN_info.contra_time_from_last_ripples_UP(DOWN_UP_index(contraRipples_cat))');

ipsi_Delta_peaks_zscore_DU = (UP_DOWN_info.SWpeakmag_DU(DOWN_UP_index)');
contra_Delta_peaks_zscore_DU = (UP_DOWN_info.contra_Delta_peaks_zscore_DU(DOWN_UP_index)');
DOWN_duration = UP_DOWN_info.previous_DOWN_duration(DOWN_UP_index);
ipsi_spindles = double(sum(UP_DOWN_info.ipsi_spindles_UP(DOWN_UP_index,:),2,'omitnan')>0);
contra_spindles = double(sum(UP_DOWN_info.contra_spindles_UP(DOWN_UP_index,:),2,'omitnan')>0);

ripples_index = find(ipsiRipples_cat==1 | contraRipples_cat==1);

% [b, logl, H, stats] = coxphfit(UP_DOWN_lag(index), ipsi_time_from_last_ripple(index), 'Strata', subject_id(index));

%     clear output


if isempty(output)
 
    %%%% Ripple occurance
    for iBoot = 1:nBoot
        
        tic
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        index = randsample(s, size(DOWN_UP_index,2), size(DOWN_UP_index,2), true);
        %         index = intersect(index,ripples_index)
        tbl = table(normalize(ipsi_ripple_duration(index)),normalize(contra_ripple_duration(index)),normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), normalize(ipsiHPC_cat(index)),normalize(contraHPC_cat(index)), ...
            normalize(ipsi_Delta_peaks_zscore_DU(index))',  normalize(contra_Delta_peaks_zscore_DU(index))', ...
            normalize(DOWN_duration(index)),ipsi_spindles(index),contra_spindles(index), ...
            categorical(subject_id(index)),categorical(session_id(index)), ...
            'VariableNames', {'ipsiDuration','contraDuration','ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC' ...
            'ipsiDelta','contraDelta','downDuration','ipsiSpindle','contraSpindle','subjectID','sessionID'});

        % Define model list (only one for now)
        modelList = {
            % --- Ipsi Models ---
            'ipsiHPC ~ ipsiDelta + (1|subjectID) + (1|sessionID)';
            'ipsiHPC ~ downDuration + (1|subjectID) + (1|sessionID)';
            'ipsiHPC ~ ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'ipsiPower ~ ipsiDelta + (1|subjectID) + (1|sessionID)';
            'ipsiPower ~ downDuration + (1|subjectID) + (1|sessionID)';
            'ipsiPower ~ ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'ipsiOccurance ~ ipsiDelta + (1|subjectID) + (1|sessionID)';
            'ipsiOccurance ~ downDuration + (1|subjectID) + (1|sessionID)';
            'ipsiOccurance ~ ipsiSpindle + (1|subjectID) + (1|sessionID)';

            % --- Contra Models ---
            'contraHPC ~ ipsiDelta + (1|subjectID) + (1|sessionID)';
            'contraHPC ~ downDuration + (1|subjectID) + (1|sessionID)';
            'contraHPC ~ ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'contraPower ~ ipsiDelta + (1|subjectID) + (1|sessionID)';
            'contraPower ~ downDuration + (1|subjectID) + (1|sessionID)';
            'contraPower ~ ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'contraOccurance ~ ipsiDelta + (1|subjectID) + (1|sessionID)';
            'contraOccurance ~ downDuration + (1|subjectID) + (1|sessionID)';
            'contraOccurance ~ ipsiSpindle + (1|subjectID) + (1|sessionID)';

            % --- Multi-Predictor Models ---
            'ipsiOccurance ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'ipsiHPC ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'ipsiPower ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'contraOccurance ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'contraHPC ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|sessionID)';
            'contraPower ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|sessionID)';
            };

        nModels = numel(modelList);

        % Preallocate local outputs
        local_b = cell(nModels,1);
        local_R2 = zeros(nModels,1);
        % local_R2_McFadden = zeros(nModels,1);
        local_tstat = cell(nModels,1);
        local_pval = cell(nModels,1);
        local_model = modelList; % store model formulas
        local_variable = cell(nModels,1);


        %
        % % Fit null model (only intercept)
        % glme_null = fitglme(tbl, 'ipsiOccurance ~ 1 + (1|subjectID)', ...
        %     'Distribution', 'Binomial', ...
        %     'Link', 'logit');
        % ipsi_LL_null = glme_null.LogLikelihood;
        %
        % glme_null = fitglme(tbl, 'contraOccurance ~ 1 + (1|subjectID)', ...
        %     'Distribution', 'Binomial', ...
        %     'Link', 'logit');
        % contra_LL_null = glme_null.LogLikelihood;

        for m = 1:nModels
            if contains(modelList{m},'Occurance ~')
                glme = fitglme(tbl, modelList{m},'Distribution', 'Binomial', 'Link', 'logit');
            else
                glme = fitglme(tbl, modelList{m});
            end
            % Save results
            local_b{m} = glme.Coefficients.Estimate(2:end);
            local_R2(m) = glme.Rsquared.Ordinary;

            % % McFadden's pseudo-R²
            % if contains(modelList{m},'ipsiOccurance')
            %     local_R2_McFadden(m) = 1 - (glme.LogLikelihood / ipsi_LL_null);
            % elseif contains(modelList{m},'contraOccurance')
            %     local_R2_McFadden(m) = 1 - (glme.LogLikelihood / contra_LL_null);
            % end


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


end

if plot_option == 1


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract bootstrapped stats for delta power vs HPC excitation
    index = find(contains(output(1).model,'ipsiHPC ~ ipsiDelta + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    index = find(contains(output(1).model,'contraHPC ~ ipsiDelta + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        contra_b(nBoot) = output(nBoot).b{index};
        R2_contra(nBoot) = output(nBoot).R2(index);
        pval_contra(nBoot) = output(nBoot).pval{index};
    end


    % Compute means and confidence intervals for t-statistics
    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);

    mean_contra = mean(contra_b);
    ci_contra = prctile(contra_b, [2.5 97.5]);

    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    % Prepare figure
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name = 'HPC excitation predicted by DOWN UP delta power';

    customColors = [0,90,50; 74,20,134; 228,42,168]/256;

    % SCATTER: ipsi
    nexttile
    x = ipsi_Delta_peaks_zscore_DU(:);  % force column
    y = ipsiHPC_cat(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.05)
    ylim([-2 15])
    xlim([0 5])
    xlabel('DOWN UP transition delta power')
    ylabel('ipsi HPC rebound excitation')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    hold on
    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);
    medianR2_ipsi = prctile(R2_ipsi, 50);
    medianP_ipsi = prctile(pval_ipsi, 50);
    col = 'k'; if medianP_ipsi < 0.05, col = 'r'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.5*max(y), ...
        sprintf('R^2 = %.3f\np = %.2e', medianR2_ipsi, medianP_ipsi), ...
        'FontSize', 10, 'Color', col)

    % SCATTER: contra
    nexttile
    x = ipsi_Delta_peaks_zscore_DU(:);  % still same x predictor
    y = contraHPC_cat(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.05)
    ylim([-2 15])
    xlim([0 5])
    xlabel('DOWN UP transition delta power')
    ylabel('contra HPC rebound excitation')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    hold on
    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);
    medianR2_contra = prctile(R2_contra, 50);
    medianP_contra = prctile(pval_contra, 50);
    col = 'k'; if medianP_contra < 0.05, col = 'r'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.5*max(y), ...
        sprintf('R^2 = %.3f\np = %.2e', medianR2_contra, medianP_contra), ...
        'FontSize', 10, 'Color', col)

    % BAR PLOT: t-stats
    nexttile
    barColors = [0,90,50; 74,20,134]/256;
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', barColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
    end
    xlim([0.5 2.5])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b')
    title('b of HPC excitation')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off
    hold off



    % Bootstrapped stats for ripple power (delta-based)
    index = find(contains(output(1).model,'ipsiPower ~ ipsiDelta + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    index = find(contains(output(1).model,'contraPower ~ ipsiDelta + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        contra_b(nBoot) = output(nBoot).b{index};
        R2_contra(nBoot) = output(nBoot).R2(index);
        pval_contra(nBoot) = output(nBoot).pval{index};
    end

    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(contra_b);
    ci_contra = prctile(contra_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];



    %%%%%%%%%%%%%%%%%%
    fig = figure('Color','w');
    fig.Position = [350 59 800 930];
    fig.Name = 'Ripples predicted by DOWN UP delta power';

    customColors = [0,90,50;74,20,134; 228,42,168]/256;

    % --- Scatter: ipsi ripple power
    nexttile
    x = ipsi_Delta_peaks_zscore_DU(:);
    y = ipsi_ripple_power(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2)
    ylim([5 25])
    xlim([0 5])
    xlabel('DOWN UP transition delta power')
    ylabel('ipsi ripple power')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    hold on
    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);
    medianR2_ipsi = prctile(R2_ipsi, 50);
    medianP_ipsi = prctile(pval_ipsi, 50);
    col = 'k'; if medianP_ipsi < 0.05, col = 'r'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.5 * max(y), ...
        sprintf('R^2 = %.2e\np = %.2e', medianR2_ipsi, medianP_ipsi), ...
        'FontSize', 10, 'Color', col)

    % --- Scatter: contra ripple power
    nexttile
    x = ipsi_Delta_peaks_zscore_DU(:);
    y = contra_ripple_power(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.2)
    ylim([5 25])
    xlim([0 5])
    xlabel('DOWN UP transition delta power')
    ylabel('contra ripple power')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    hold on
    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);
    medianR2_contra = prctile(R2_contra, 50);
    medianP_contra = prctile(pval_contra, 50);
    col = 'k'; if medianP_contra < 0.05, col = 'r'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.5 * max(y), ...
        sprintf('R^2 = %.2e\np = %.2e', medianR2_contra, medianP_contra), ...
        'FontSize', 10, 'Color', col)

    % --- Histogram: delta power vs ipsi ripple occurrence
    nexttile
    histogram(ipsi_Delta_peaks_zscore_DU(ipsiRipples_cat == 1), 0:0.1:5, 'Normalization', 'probability'); hold on
    histogram(ipsi_Delta_peaks_zscore_DU(ipsiRipples_cat == 0), 0:0.1:5, 'Normalization', 'probability')
    [~, pVal_ks_ipsi] = kstest2(ipsi_Delta_peaks_zscore_DU(ipsiRipples_cat == 1), ipsi_Delta_peaks_zscore_DU(ipsiRipples_cat == 0));
    text(4.5, 0.05, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    legend('With ipsi ripples', 'Without ipsi ripples', 'box', 'off')
    xlabel('DOWN UP transition delta power')
    ylabel('Proportion of events')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)


    % --- Histogram: delta power vs contra ripple occurrence
    nexttile
    histogram(ipsi_Delta_peaks_zscore_DU(contraRipples_cat == 1), 0:0.1:5, 'Normalization', 'probability'); hold on
    histogram(ipsi_Delta_peaks_zscore_DU(contraRipples_cat == 0), 0:0.1:5, 'Normalization', 'probability')
    [~, pVal_ks_ipsi] = kstest2(ipsi_Delta_peaks_zscore_DU(contraRipples_cat == 1), ipsi_Delta_peaks_zscore_DU(contraRipples_cat == 0));
    text(4.5, 0.05, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    legend('With contra ripples', 'Without contra ripples', 'box', 'off')
    xlabel('DOWN UP transition delta power')
    ylabel('Proportion of events')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)

    % --- Bar plot: b ripple power
    % mean_ipsi = mean(ipsi_t);
    % ci_ipsi = prctile(ipsi_t, [2.5 97.5]);
    % mean_contra = mean(contra_t);
    % ci_contra = prctile(contra_t, [2.5 97.5]);
    % barData = [mean_ipsi, mean_contra];
    % lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    % upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    nexttile
    barColors = [0,90,50;74,20,134]/256;
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', barColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
    end
    xlim([0.5 2.5])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b')
    title('b of Ripple power')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off
    hold off


    % --- Bar plot: b ripple occurrence
    index = find(contains(output(1).model,'ipsiOccurance ~ ipsiDelta + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
    index = find(contains(output(1).model,'contraOccurance ~ ipsiDelta + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        contra_b(nBoot) = output(nBoot).b{index};
        R2_contra(nBoot) = output(nBoot).R2(index);
        pval_contra(nBoot) = output(nBoot).pval{index};
    end

    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(contra_b);
    ci_contra = prctile(contra_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];


    nexttile
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', barColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
        if i == 1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        else
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_contra,50)), ...
                'FontSize', 10, 'Color', 'k');
        end

    end
    xlim([0.5 2.5])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b')
    title('b of Ripple occurrence')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off
    hold off







    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name ='HPC excitation predicted by DOWN UP + spindles';

    customColors = [0,90,50;74,20,134; 228,42,168 ]/256;

    % --- Histogram: ipsi HPC excitation with/without spindle
    nexttile
    histogram(ipsiHPC_cat(ipsi_spindles==1), 0:0.1:5, 'Normalization', 'probability'); hold on
    histogram(ipsiHPC_cat(ipsi_spindles==0), 0:0.1:5, 'Normalization', 'probability')
    [~, pVal_ks_ipsi] = kstest2(ipsiHPC_cat(ipsi_spindles==1), ipsiHPC_cat(ipsi_spindles==0));
    text(4.5, 0.05, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    legend('With ipsi spindles','Without ipsi spindles','box','off')
    xlabel('ipsi HPC MUA excitation (z)')
    ylabel('Proportion of events')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    % --- Histogram: contra HPC excitation with/without spindle
    nexttile
    histogram(contraHPC_cat(ipsi_spindles==1), 0:0.1:5, 'Normalization', 'probability'); hold on
    histogram(contraHPC_cat(ipsi_spindles==0), 0:0.1:5, 'Normalization', 'probability')
    [~, pVal_ks_contra] = kstest2(contraHPC_cat(ipsi_spindles==1), contraHPC_cat(ipsi_spindles==0));
    text(4.5, 0.05, sprintf('p = %.2e', pVal_ks_contra), 'FontSize', 10, 'Color', 'k')
    legend('With ipsi spindles','Without ipsi spindles','box','off')
    xlabel('contra HPC MUA excitation (z)')
    ylabel('Proportion of events')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    % --- Bar plot: b spindle
    index = find(contains(output(1).model,'ipsiHPC ~ ipsiSpindle + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
        index = find(contains(output(1).model,'contraHPC ~ ipsiSpindle + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        contra_b(nBoot) = output(nBoot).b{index};
        R2_contra(nBoot) = output(nBoot).R2(index);
        pval_contra(nBoot) = output(nBoot).pval{index};
    end

    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(contra_b);
    ci_contra = prctile(contra_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];


    nexttile
    barColors = [0,90,50;74,20,134]/256;
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', barColors(i,:), 'EdgeColor', 'none','FaceAlpha',0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
    end
    xlim([0.5 2.5])
    ylim([0 0.3])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b')
    title('b of HPC excitation')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
    grid off
    hold off





    % Initialize figure
    fig = figure('Color','w');
    fig.Position = [350 59 800 930*2/3];
    fig.Name = 'Ripples predicted by spindles';

    customColors = [0,90,50;74,20,134]/256;

    % Histogram & KS test: ripple power (ipsi)
    nexttile
    edges = 5:0.5:25;
    histogram(ipsi_ripple_power(ipsi_spindles==1), edges, 'Normalization', 'probability', ...
        'FaceAlpha', 0.5); hold on
    histogram(ipsi_ripple_power(ipsi_spindles==0), edges, 'Normalization', 'probability', ...
        'FaceAlpha', 0.5);

    [~, pVal_ks_ipsi] = kstest2(ipsi_ripple_power(ipsi_spindles==1), ipsi_ripple_power(ipsi_spindles==0));
    text(20, 0.008, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')

    legend('With ipsi spindles', 'Without ipsi spindles', 'Box', 'off')
    xlabel('Ipsi ripple power')
    ylabel('Proportion of events')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    % Histogram & KS test: ripple power (contra)
    nexttile
    histogram(contra_ripple_power(ipsi_spindles==1), edges, 'Normalization', 'probability', ...
        'FaceAlpha', 0.5); hold on
    histogram(contra_ripple_power(ipsi_spindles==0), edges, 'Normalization', 'probability', ...
        'FaceAlpha', 0.5);

    [~, pVal_ks_contra] = kstest2(contra_ripple_power(ipsi_spindles==1), contra_ripple_power(ipsi_spindles==0));
    text(20, 0.008, sprintf('p = %.2e', pVal_ks_contra), 'FontSize', 10, 'Color', 'k')

    legend('With ipsi spindles', 'Without ipsi spindles', 'Box', 'off')
    xlabel('Contra ripple power')
    ylabel('Proportion of events')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    % Bootstrapped t-stat: ripple power (ipsi/contra)
        % --- Bar plot: b spindle
    index = find(contains(output(1).model,'ipsiPower ~ ipsiSpindle + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
        index = find(contains(output(1).model,'contraPower ~ ipsiSpindle + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        contra_b(nBoot) = output(nBoot).b{index};
        R2_contra(nBoot) = output(nBoot).R2(index);
        pval_contra(nBoot) = output(nBoot).pval{index};
    end

    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(contra_b);
    ci_contra = prctile(contra_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];
    

    nexttile
    x_pos = [1, 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none','FaceAlpha',0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
        if i == 1
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
                'FontSize', 10, 'Color', 'k');
        else
            text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_contra,50)), ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    ylim([-0.1 0.2])
    xlim([0.5 2.5])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b')
    title('Ripple power b (spindle effect)')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
    grid off

    % %% Ripple occurrence proportion & t-stat
    % % Ripple proportion: ipsi
    % p_with = sum(ipsiRipples_cat(ipsi_spindles==1)) / sum(ipsi_spindles==1);
    % p_without = sum(ipsiRipples_cat(ipsi_spindles==0)) / sum(ipsi_spindles==0);

    % KS test-like binary comparison (proportions)
    % Bootstrapped t-stats (occurrence)
    index = find(contains(output(1).model,'ipsiOccurance ~ ipsiSpindle + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        ipsi_b(nBoot) = output(nBoot).b{index};
        R2_ipsi(nBoot) = output(nBoot).R2(index);
        pval_ipsi(nBoot) = output(nBoot).pval{index};
    end
        index = find(contains(output(1).model,'contraOccurance ~ ipsiSpindle + (1|subjectID) + (1|sessionID)'));
    for nBoot = 1:1000
        contra_b(nBoot) = output(nBoot).b{index};
        R2_contra(nBoot) = output(nBoot).R2(index);
        pval_contra(nBoot) = output(nBoot).pval{index};
    end

    mean_ipsi = mean(ipsi_b);
    ci_ipsi = prctile(ipsi_b, [2.5 97.5]);
    mean_contra = mean(contra_b);
    ci_contra = prctile(contra_b, [2.5 97.5]);
    barData = [mean_ipsi, mean_contra];
    lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    nexttile
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none','FaceAlpha',0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
    end
    ylim([0 0.6])
    xlim([0.5 2.5])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped b')
    title('Ripple occurrence b (spindle effect)')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
    grid off

    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fig = figure('Color','w');
    % fig.Position = [350 59 800 930/3*2];
    % fig.Name = 'HPC excitation predicted by previous DOWN duration';
    % 
    % customColors = [0,90,50; 74,20,134; 228,42,168]/256;
    % 
    % % Bootstrapped stats
    % for nBoot = 1:1000
    %     ipsi_t(nBoot) = output(nBoot).t_stat{3};
    %     contra_t(nBoot) = output(nBoot).t_stat{4};
    %     R2_ipsi(nBoot) = output(nBoot).R2(3);
    %     R2_contra(nBoot) = output(nBoot).R2(4);
    %     pval_ipsi(nBoot) = output(nBoot).pval{3};
    %     pval_contra(nBoot) = output(nBoot).pval{4};
    % end
    % 
    % % --- Ipsi HPC excitation
    % nexttile
    % x = DOWN_duration(:); y = ipsiHPC_cat(:);
    % scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.05)
    % xlabel('DOWN duration (s)'); ylabel('ipsi HPC excitation')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    % valid = ~isnan(x) & ~isnan(y) & x < 1;
    % coeffs = polyfit(x(valid), y(valid), 1);
    % x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    % y_fit = polyval(coeffs, x_fit);
    % medR2 = median(R2_ipsi); medP = median(pval_ipsi);
    % col = 'k'; if medP < 0.05, col = 'r'; end
    % plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    % text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
    %     sprintf('R^2 = %.2e\np = %.2e', medR2, medP), 'FontSize', 10, 'Color', col)
    % xlim([0 1])
    % 
    % % --- Contra HPC excitation
    % nexttile
    % x = DOWN_duration(:); y = contraHPC_cat(:);
    % scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.05)
    % xlabel('DOWN duration (s)'); ylabel('contra HPC excitation')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    % valid = ~isnan(x) & ~isnan(y) & x < 1;
    % coeffs = polyfit(x(valid), y(valid), 1);
    % x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    % y_fit = polyval(coeffs, x_fit);
    % medR2 = median(R2_contra); medP = median(pval_contra);
    % col = 'k'; if medP < 0.05, col = 'r'; end
    % plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    % text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
    %     sprintf('R^2 = %.2e\np = %.2e', medR2, medP), 'FontSize', 10, 'Color', col)
    % xlim([0 1])
    % 
    % % --- t-stat bar plot
    % mean_ipsi = mean(ipsi_t); ci_ipsi = prctile(ipsi_t, [2.5 97.5]);
    % mean_contra = mean(contra_t); ci_contra = prctile(contra_t, [2.5 97.5]);
    % barData = [mean_ipsi, mean_contra];
    % lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    % upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];
    % 
    % nexttile
    % x_pos = [1 2];
    % for i = 1:2
    %     hold on
    %     bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
    %         'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    %     errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
    %         'k', 'linestyle', 'none', 'linewidth', 1.5);
    % end
    % xlim([0.5 2.5]); xticks(x_pos); xticklabels({'ipsi','contra'})
    % ylabel('Bootstrapped t-statistic'); title('HPC excitation')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); grid off
    % hold off
    % 
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%
    % fig = figure('Color','w');
    % fig.Position = [350 59 800 930];
    % fig.Name = 'Ripples predicted by previous DOWN duration';
    % 
    % % Bootstrapped stats
    % for nBoot = 1:1000
    %     ipsi_t(nBoot) = output(nBoot).t_stat{9};
    %     contra_t(nBoot) = output(nBoot).t_stat{10};
    %     R2_ipsi(nBoot) = output(nBoot).R2(9);
    %     R2_contra(nBoot) = output(nBoot).R2(10);
    %     pval_ipsi(nBoot) = output(nBoot).pval{9};
    %     pval_contra(nBoot) = output(nBoot).pval{10};
    % end
    % 
    % % --- Ipsi ripple power
    % nexttile
    % x = DOWN_duration(:); y = ipsi_ripple_power(:);
    % scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2)
    % xlabel('DOWN duration (s)'); ylabel('ipsi ripple power')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    % valid = ~isnan(x) & ~isnan(y);
    % coeffs = polyfit(x(valid), y(valid), 1);
    % x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    % y_fit = polyval(coeffs, x_fit);
    % medR2 = median(R2_ipsi); medP = median(pval_ipsi);
    % col = 'k'; if medP < 0.05, col = 'r'; end
    % plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    % text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
    %     sprintf('R^2 = %.2e\np = %.2e', medR2, medP), 'FontSize', 10, 'Color', col)
    % xlim([0 1])
    % ylim([0 25])
    % 
    % % --- Contra ripple power
    % nexttile
    % x = DOWN_duration(:); y = contra_ripple_power(:);
    % scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.2)
    % xlabel('DOWN duration (s)'); ylabel('ipsi ripple power')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    % valid = ~isnan(x) & ~isnan(y);
    % coeffs = polyfit(x(valid), y(valid), 1);
    % x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    % y_fit = polyval(coeffs, x_fit);
    % medR2 = median(R2_contra); medP = median(pval_contra);
    % col = 'k'; if medP < 0.05, col = 'r'; end
    % plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    % text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
    %     sprintf('R^2 = %.2e\np = %.2e', medR2, medP), 'FontSize', 10, 'Color', col)
    % xlim([0 1])
    % ylim([0 25])
    % 
    % % --- Histogram: DOWN duration vs. ipsi ripple presence
    % nexttile
    % histogram(DOWN_duration(ipsiRipples_cat==1), 0:0.01:1, 'Normalization','probability'); hold on
    % histogram(DOWN_duration(ipsiRipples_cat==0), 0:0.01:1, 'Normalization','probability')
    % [~, pVal_ks_ipsi] = kstest2(DOWN_duration(ipsiRipples_cat==1), DOWN_duration(ipsiRipples_cat==0));
    % text(0.8, 0.05, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    % legend('With ipsi ripple', 'Without ipsi ripple', 'box','off')
    % xlabel('DOWN duration (s)'); ylabel('Proportion')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12)
    % xlim([0 1])
    % % ylim([0 25])
    % 
    % % --- Histogram: DOWN duration vs. contra ripple presence
    % nexttile
    % histogram(DOWN_duration(contraRipples_cat==1), 0:0.01:1, 'Normalization','probability'); hold on
    % histogram(DOWN_duration(contraRipples_cat==0), 0:0.01:1, 'Normalization','probability')
    % [~, pVal_ks_ipsi] = kstest2(DOWN_duration(contraRipples_cat==1), DOWN_duration(contraRipples_cat==0));
    % text(0.8, 0.05, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    % legend('With contra ripple', 'Without contra ripple', 'box','off')
    % xlabel('DOWN duration (s)'); ylabel('Proportion')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12)
    % xlim([0 1])
    % 
    % % --- t-stat bar plots for ripple power
    % mean_ipsi = mean(ipsi_t); ci_ipsi = prctile(ipsi_t, [2.5 97.5]);
    % mean_contra = mean(contra_t); ci_contra = prctile(contra_t, [2.5 97.5]);
    % barData = [mean_ipsi, mean_contra];
    % lowerCI = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    % upperCI = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];
    % 
    % nexttile
    % for i = 1:2
    %     hold on
    %     bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
    %         'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    %     errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
    %         'k', 'linestyle', 'none', 'linewidth', 1.5);
    % end
    % xlim([0.5 2.5]); xticks(x_pos); xticklabels({'ipsi','contra'})
    % ylabel('Bootstrapped t-statistic'); title('Ripple power')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); grid off
    % hold off
    % 
    % 
    % for nBoot = 1:1000
    %     ipsi_t(nBoot) = output(nBoot).t_stat{19};
    %     contra_t(nBoot) = output(nBoot).t_stat{20};
    %     pval_ipsi(nBoot) = output(nBoot).pval{19};
    %     pval_contra(nBoot) = output(nBoot).pval{20};
    % end
    % 
    % 
    % nexttile
    % barData = [mean(ipsi_t), mean(contra_t)];
    % ci_ipsi = prctile(ipsi_t, [2.5 97.5]);
    % ci_contra = prctile(contra_t, [2.5 97.5]);
    % lowerCI = [barData(1) - ci_ipsi(1), barData(2) - ci_contra(1)];
    % upperCI = [barData(1) - ci_ipsi(2), barData(2) - ci_contra(2)];
    % for i = 1:2
    %     hold on
    %     bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
    %         'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    %     errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
    %         'k', 'linestyle', 'none', 'linewidth', 1.5);
    %     if i == 1
    %         text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_ipsi,50)), ...
    %             'FontSize', 10, 'Color', 'k');
    %     else
    %         text(x_pos(i)+0.25, barData(i), sprintf('p = %.2e', prctile(pval_contra,50)), ...
    %             'FontSize', 10, 'Color', 'k');
    %     end
    % end
    % xlim([0.5 2.5]); xticks(x_pos); xticklabels({'ipsi','contra'})
    % ylabel('Bootstrapped t-statistic'); title('t-statistic of Ripple occurrence')
    % set(gca, 'TickDir','out','Box','off','Color','none','FontSize',12); grid off
    % 



end
end

