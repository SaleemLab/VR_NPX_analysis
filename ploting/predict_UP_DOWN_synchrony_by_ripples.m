function output = predict_UP_DOWN_synchrony_by_ripples(ipsi_V1_MUA,contra_V1_MUA,ipsi_HPC_MUA,contra_HPC_MUA,ipsi_ripples,contra_ripples,UP_DOWN_info,varargin)
% predict_HPC_MUA_DOWN_UP
% Predict hippocampal ripple occurance and ripple power and ripple
% synchrony during DOWN-UP transitions using ipsi and contra V1 MUA and
% ipsi-contra lag
%
% OUTPUT:
%   outputStruct : structure containing all betas, permutation p-values, and bootstrap CIs

% SETTINGS
p = inputParser;
addParameter(p, 'nPerm', 1000);
addParameter(p, 'nBoot', 1000);
addParameter(p, 'output', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
addParameter(p, 'UP_DOWN_lag', []);
addParameter(p, 'UP_DOWN_index', []);
% addParameter(p, 'ripples_info', []);
% addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
output = p.Results.output;
subject_id = p.Results.subject_id;
UP_DOWN_lag = p.Results.UP_DOWN_lag;
UP_DOWN_index = p.Results.UP_DOWN_index;
plot_option = p.Results.plot_option;
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
key_window_idx2 = find(timeVec < 0 & timeVec >= -0.1);
key_window_idx1 = find(timeVec > 0 & timeVec <= 0.05);

ipsiV1_key = max(ipsi_V1_MUA(:, key_window_idx2)')'-min(ipsi_V1_MUA(:, key_window_idx1)')' ;
contraV1_key = max(contra_V1_MUA(:, key_window_idx2)')' - min(contra_V1_MUA(:, key_window_idx1)')';
% ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')' -min(ipsi_HPC_MUA(:, key_window_idx1)')' ;
% contraHPC_key = max(contra_HPC_MUA(:, key_window_idx2)')'-min(contra_HPC_MUA(:, key_window_idx1)')' ;

ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')' ;
contraHPC_key = max(contra_HPC_MUA(:, key_window_idx2)')';

% ipsiV1_key = mean(ipsi_V1_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_V1_MUA(:, key_window_idx1)',"omitnan")' ;
% contraV1_key = mean(contra_V1_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_V1_MUA(:, key_window_idx1)',"omitnan")';
% ipsiHPC_key = mean(ipsi_HPC_MUA(:, key_window_idx2)',"omitnan")' ;
% contraHPC_key = mean(contra_HPC_MUA(:, key_window_idx2)',"omitnan")';
% ipsiHPC_key = mean(ipsi_HPC_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_HPC_MUA(:, key_window_idx1)',"omitnan")' ;
% contraHPC_key = mean(contra_HPC_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_HPC_MUA(:, key_window_idx1)',"omitnan")';
    % combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

%%%% HPC ripples
key_window_idx2= find(timeVec_prob < 0 & timeVec_prob >= -0.1);
%     key_window_idx2 = find(timeVec_prob >= 0 & timeVec_prob <= 0.1);

ipsiRipples_key = sum(ipsi_ripples(:,key_window_idx2),2)>0 ;
contraRipples_key =sum(contra_ripples(:,key_window_idx2),2)>0 ;

ipsiV1_cat = ipsiV1_key(:);
contraV1_cat = contraV1_key(:);
ipsiHPC_cat = ipsiHPC_key(:);
contraHPC_cat = contraHPC_key(:);
ipsiRipples_cat = ipsiRipples_key(:);
contraRipples_cat = contraRipples_key(:);
%     combinedV1_cat = combinedV1_key(:);

validIdx = ~(isnan(ipsiV1_cat) | isnan(contraV1_cat) | isnan(ipsiHPC_cat) | isnan(contraHPC_cat));
ipsiV1_cat = ipsiV1_cat(validIdx);
contraV1_cat = contraV1_cat(validIdx);
ipsiHPC_cat = ipsiHPC_cat(validIdx);
contraHPC_cat = contraHPC_cat(validIdx);
ipsiRipples_cat = ipsiRipples_cat(validIdx);
contraRipples_cat = contraRipples_cat(validIdx);
UP_DOWN_index = UP_DOWN_index(validIdx);

%     ipsi_ripple_power(ipsiRipples_cat) = ripples_info.ipsi_first_ripples_power_UP(UP_DOWN_index(ipsiRipples_cat))';
%     contra_ripple_power(contraRipples_cat) = ripples_info.contra_first_ripples_power_UP(UP_DOWN_index(contraRipples_cat))';
%     ipsi_ripple_lag(ipsiRipples_cat) = abs(ripples_info.ipsi_first_ripples_lag_UP(UP_DOWN_index(ipsiRipples_cat))');
%     contra_ripple_lag(contraRipples_cat) = abs(ripples_info.contra_first_ripples_lag_UP(UP_DOWN_index(contraRipples_cat))');


ipsi_ripple_power = nan(length(UP_DOWN_index),1);
contra_ripple_power = nan(length(UP_DOWN_index),1);
ipsi_ripple_duration = nan(length(UP_DOWN_index),1);
contra_ripple_duration = nan(length(UP_DOWN_index),1);
% ripple_lag = nan(length(DOWN_UP_index),1);
ipsi_ripple_lag = nan(length(UP_DOWN_index),1);
contra_ripple_lag = nan(length(UP_DOWN_index),1);
ipsi_time_from_last_ripple = nan(length(UP_DOWN_index),1);
contra_time_from_last_ripple = nan(length(UP_DOWN_index),1);
ipsi_Delta_peaks_zscore_UD = nan(length(UP_DOWN_index),1);
contra_Delta_peaks_zscore_UD = nan(length(UP_DOWN_index),1);


ipsi_ripple_power(ipsiRipples_cat) = UP_DOWN_info.ipsi_last_ripples_power_UP(UP_DOWN_index(ipsiRipples_cat))';
contra_ripple_power(contraRipples_cat) = UP_DOWN_info.contra_last_ripples_power_UP(UP_DOWN_index(contraRipples_cat))';
ipsi_ripple_duration(ipsiRipples_cat) = UP_DOWN_info.ipsi_last_ripples_duration_UP(UP_DOWN_index(ipsiRipples_cat))';
contra_ripple_duration(contraRipples_cat) = UP_DOWN_info.contra_last_ripples_duration_UP(UP_DOWN_index(contraRipples_cat))';
ipsi_ripple_lag(ipsiRipples_cat) = abs(UP_DOWN_info.ipsi_last_ripples_lag_diff_UP(UP_DOWN_index(ipsiRipples_cat))');
contra_ripple_lag(contraRipples_cat) = abs(UP_DOWN_info.contra_last_ripples_lag_diff_UP(UP_DOWN_index(contraRipples_cat))');
ipsi_time_from_last_ripple(ipsiRipples_cat) = abs(UP_DOWN_info.ipsi_time_from_last_ripples_UP(UP_DOWN_index(ipsiRipples_cat))');
contra_time_from_last_ripple(contraRipples_cat) = abs(UP_DOWN_info.contra_time_from_last_ripples_UP(UP_DOWN_index(contraRipples_cat))');


ipsi_Delta_peaks_zscore_UD = (UP_DOWN_info.SWpeakmag_UD(UP_DOWN_index)');
contra_Delta_peaks_zscore_UD = (UP_DOWN_info.contra_Delta_peaks_zscore_UD(UP_DOWN_index)');
DOWN_duration = UP_DOWN_info.next_DOWN_duration(UP_DOWN_index);
ipsi_spindles = double(sum(UP_DOWN_info.ipsi_spindles_DOWN(UP_DOWN_index,:),2,'omitnan')>0);
contra_spindles = double(sum(UP_DOWN_info.contra_spindles_DOWN(UP_DOWN_index,:),2,'omitnan')>0);


ripples_index = find(ipsiRipples_cat==1 | contraRipples_cat==1);

if isempty(output)
    % [b, logl, H, stats] = coxphfit(UP_DOWN_lag(index), ipsi_time_from_last_ripple(index), 'Strata', subject_id(index));

    clear output
    %%%%%%%% All other models
    parfor iBoot = 1:nBoot
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index = randsample(s, ripples_index, size(ripples_index,1), true);
        % index = randsample(s, size(UP_DOWN_index,2), size(UP_DOWN_index,2), true);
        tbl = table(normalize(ipsi_ripple_duration(index)),normalize(contra_ripple_duration(index)),normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), ...
            zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)), ...
            zscore(ipsiV1_cat(index)), zscore(contraV1_cat(index)), ...
            normalize(ipsi_Delta_peaks_zscore_UD(index))',  normalize(contra_Delta_peaks_zscore_UD(index))', ...
            normalize(DOWN_duration(index)),ipsi_spindles(index),contra_spindles(index), ...
            zscore(UP_DOWN_lag(index)), categorical(subject_id(index)), ...
            'VariableNames', {'ipsiDuration','contraDuration','ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC','ipsiV1','contraV1','ipsiDelta','contraDelta', ...
            'downDuration','ipsiSpindle','contraSpindle', ...
            'UP_DOWN_lag','subjectID'});


        % Model list
        modelList = {
            'ipsiV1 ~ ipsiPower + contraPower + ipsiPower*contraPower + (1|subjectID)';
            'ipsiV1 ~ ipsiPower + (1|subjectID)';
            'ipsiV1 ~ contraPower + (1|subjectID)';

            'ipsiV1 ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)';
            'ipsiV1 ~ ipsiHPC + (1|subjectID)';
            'ipsiV1 ~ contraHPC + (1|subjectID)';
            'ipsiV1 ~ ipsiLag + (1|subjectID)';
            'ipsiV1 ~ ipsiDuration + (1|subjectID)';


            'UP_DOWN_lag ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)';
            'UP_DOWN_lag ~ ipsiHPC + (1|subjectID)';
            'UP_DOWN_lag ~ contraHPC + (1|subjectID)';

            'UP_DOWN_lag ~ ipsiPower + contraPower + ipsiPower*contraPower + (1|subjectID)';
            'UP_DOWN_lag ~ ipsiPower + (1|subjectID)';
            'UP_DOWN_lag ~ contraPower + (1|subjectID)';
            'UP_DOWN_lag ~ ipsiLag + (1|subjectID)';
            'UP_DOWN_lag ~ ipsiDuration + (1|subjectID)';


            'ipsiDelta ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)';
            'ipsiDelta ~ ipsiHPC + (1|subjectID)';
            'ipsiDelta ~ contraHPC + (1|subjectID)';

            'ipsiDelta ~ ipsiPower + contraPower + ipsiPower*contraPower + (1|subjectID)';
            'ipsiDelta ~ ipsiPower + (1|subjectID)';
            'ipsiDelta ~ contraPower + (1|subjectID)';

            'ipsiDelta ~ ipsiLag + (1|subjectID)';
            'ipsiDelta ~ ipsiDuration + (1|subjectID)';


            %             'ipsiDuration ~ ipsiPower + contraPower + ipsiPower*contraPower + (1|subjectID)'
            %             'ipsiDuration ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)'
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



    %%%% Ripple occurance
    parfor iBoot = 1:nBoot
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index = randsample(s, size(UP_DOWN_index,2), size(UP_DOWN_index,2), true);
        %         index = intersect(index,ripples_index)
        tbl = table(normalize(ipsi_ripple_duration(index)),normalize(contra_ripple_duration(index)),normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), ...
            zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)), ...
            zscore(ipsiV1_cat(index)), zscore(contraV1_cat(index)), ...
            normalize(ipsi_Delta_peaks_zscore_UD(index))',  normalize(contra_Delta_peaks_zscore_UD(index))', ...
            normalize(DOWN_duration(index)),ipsi_spindles(index),contra_spindles(index), ...
            zscore(UP_DOWN_lag(index)), categorical(subject_id(index)), ...
            'VariableNames', {'ipsiDuration','contraDuration','ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC','ipsiV1','contraV1','ipsiDelta','contraDelta', ...
            'downDuration','ipsiSpindle','contraSpindle', ...
            'UP_DOWN_lag','subjectID'});

        % Define model list (only one for now)
        modelList = {
            'ipsiV1 ~ ipsiOccurance  + (1|subjectID)';
            'ipsiV1 ~ contraOccurance  + (1|subjectID)';
            'ipsiV1 ~ ipsiOccurance+ contraOccurance  + (1|subjectID)';
            'ipsiV1 ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)'

            'UP_DOWN_lag ~ ipsiOccurance  + (1|subjectID)';
            'UP_DOWN_lag ~ contraOccurance  + (1|subjectID)';
            'UP_DOWN_lag ~ ipsiOccurance+ contraOccurance  + (1|subjectID)';

            'ipsiDelta ~ ipsiOccurance + (1|subjectID)';
            'ipsiDelta ~ contraOccurance + (1|subjectID)';
            'ipsiDelta ~ ipsiOccurance + contraOccurance + ipsiOccurance*contraOccurance + (1|subjectID)';
            };

        nModels = numel(modelList);

        % Preallocate local outputs
        local_b = cell(nModels,1);
        local_R2 = zeros(nModels,1);
        local_tstat = cell(nModels,1);
        local_pval = cell(nModels,1);
        local_model = modelList; % store model formulas
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
        output2(iBoot).b = local_b;
        output2(iBoot).R2 = local_R2;
        output2(iBoot).t_stat = local_tstat;
        output2(iBoot).pval = local_pval;
        output2(iBoot).model = local_model;
        output2(iBoot).variable = local_variable;
        output2(iBoot).type = repmat({'All'}, nModels, 1);
        toc
    end



    nBoot = numel(output2); % number of bootstraps

    % Get field names automatically
    fieldNames = fieldnames(output);

    % Initialize combined structure
    combined_output = struct();

    for iBoot = 1:nBoot
        for f = 1:numel(fieldNames)
            thisField = fieldNames{f};

            % Combine the fields
            combined_output(iBoot).(thisField) = [output(iBoot).(thisField); output2(iBoot).(thisField)];
        end
    end

    output = combined_output; clear combined_output output2
end



if plot_option == 1


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------------------ Extract bootstrapped stats ------------------
    for nBoot = 1:1000
        ipsi_t(nBoot)      = output(nBoot).t_stat{10};
        contra_t(nBoot)    = output(nBoot).t_stat{11};
        pval_ipsi(nBoot)   = output(nBoot).pval{10};
        pval_contra(nBoot) = output(nBoot).pval{11};
        R2_ipsi(nBoot)     = output(nBoot).R2(10);
        R2_contra(nBoot)   = output(nBoot).R2(11);
    end

    % ------------------ Plot setup ------------------
    fig = figure('Color', 'w');
    fig.Position = [350 59 800 930/3*2];
    fig.Name = 'UP-DOWN synchrony predicted by HPC excitation';
    customColors = [0,90,50; 74,20,134]/256;

    % ------------------ IPSI ------------------
    nexttile
    x = ipsiHPC_cat(ripples_index);x = x(:);
    y = UP_DOWN_lag(ripples_index);y = y(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.05)
    xlabel('ipsi HPC excitation'); ylabel('UP DOWN transition ipsi-contra lags')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12); hold on

    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);

    % Plot visual regression line
    plot(x_fit, y_fit, 'k', 'LineWidth', 2)

    % Annotate with bootstrap medians
    medR2 = prctile(R2_ipsi,50); medP = prctile(pval_ipsi,50);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    xlim([-0.5 15])


    % ------------------ CONTRA ------------------
    nexttile
    x = contraHPC_cat(ripples_index); x = x(:);
    y = UP_DOWN_lag(ripples_index);y = y(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.05)
    xlabel('contra HPC excitation'); ylabel('UP DOWN transition ipsi-contra lags')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12); hold on

    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);

    % Plot visual regression line
    plot(x_fit, y_fit, 'k', 'LineWidth', 2)

    % Annotate with bootstrap medians
    medR2 = prctile(R2_contra,50); medP = prctile(pval_contra,50);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    xlim([-0.5 15])


    % ------------------ Bar Plot: Bootstrapped t-statistics ------------------
    mean_ipsi   = mean(ipsi_t);
    mean_contra = mean(contra_t);
    ci_ipsi     = prctile(ipsi_t, [2.5 97.5]);
    ci_contra   = prctile(contra_t, [2.5 97.5]);

    barData  = [mean_ipsi, mean_contra];
    lowerCI  = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI  = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    nexttile
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
    end
    xlim([0.5 2.5]); xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped t-statistic'); title('HPC excitation')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------------ Extract data ------------------
    for nBoot = 1:1000
        % Ripple power
        ipsi_t_power(nBoot) = output(nBoot).t_stat{13};
        contra_t_power(nBoot) = output(nBoot).t_stat{14};
        R2_ipsi_power(nBoot) = output(nBoot).R2(13);
        R2_contra_power(nBoot) = output(nBoot).R2(14);
        pval_ipsi_power(nBoot) = output(nBoot).pval{13};
        pval_contra_power(nBoot) = output(nBoot).pval{14};

        % Ripple occurrence
        ipsi_t_occ(nBoot) = output(nBoot).t_stat{29};
        contra_t_occ(nBoot) = output(nBoot).t_stat{30};
    end

    customColors = [0,90,50; 74,20,134]/256;
    x_pos = [1 2];

    fig = figure('Color','w');
    fig.Position = [350 59 800 930];
    fig.Name = 'UP DOWN lag predicted by ripples';


    % ------------------ Scatter: ripple power vs lag ------------------
    nexttile
    x = ipsi_ripple_power(:); y = UP_DOWN_lag(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2)
    xlabel('ipsi ripple power'); ylabel('UP DOWN transition ipsi-contra lags')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on

    % Simple regression for visualization only
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);

    % Median R² and p from bootstrap
    medR2 = prctile(R2_ipsi_power,50); medP = prctile(pval_ipsi_power,50)
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)

    nexttile
    x = contra_ripple_power(:); y = UP_DOWN_lag(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.2)
    xlabel('contra ripple power'); ylabel('UP DOWN transition ipsi-contra lags')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on

    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);

    medR2 = prctile(R2_contra_power,50); medP = prctile(pval_contra_power,50)

    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)


    % ------------------ Histogram: lag with vs without ripples ------------------
    nexttile
    histogram(UP_DOWN_lag(ipsiRipples_cat==1), 0:0.005:0.15, 'Normalization','probability'); hold on
    histogram(UP_DOWN_lag(ipsiRipples_cat==0), 0:0.005:0.15, 'Normalization','probability');
    [~, pVal_ks_ipsi] = kstest2(UP_DOWN_lag(ipsiRipples_cat==1), UP_DOWN_lag(ipsiRipples_cat==0));
    legend('With ipsi ripple', 'Without ipsi ripple', 'Box','off')
    xlabel('UP DOWN transition ipsi-contra lags'); ylabel('Proportion');
    text(0.1, 0.08, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    nexttile
    histogram(UP_DOWN_lag(contraRipples_cat==1), 0:0.005:0.15, 'Normalization','probability'); hold on
    histogram(UP_DOWN_lag(contraRipples_cat==0), 0:0.005:0.15, 'Normalization','probability');
    [~, pVal_ks_contra] = kstest2(UP_DOWN_lag(contraRipples_cat==1), UP_DOWN_lag(contraRipples_cat==0));
    legend('With contra ripple', 'Without contra ripple', 'Box','off')
    xlabel('UP DOWN transition ipsi-contra lags'); ylabel('Proportion');
    text(0.1, 0.08, sprintf('p = %.2e', pVal_ks_contra), 'FontSize', 10, 'Color', 'k')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)


    % ------------------ Bar plot: t-statistics (ripple power) ------------------
    mean_ipsi = mean(ipsi_t_power); ci_ipsi = prctile(ipsi_t_power, [2.5 97.5]);
    mean_contra = mean(contra_t_power); ci_contra = prctile(contra_t_power, [2.5 97.5]);
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
    end
    xticks(x_pos); xticklabels({'ipsi','contra'});
    ylabel('Bootstrapped t-statistic'); title('Ripple power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    % ------------------ Bar plot: t-statistics (ripple occurrence) ------------------
    mean_ipsi = mean(ipsi_t_occ); ci_ipsi = prctile(ipsi_t_occ, [2.5 97.5]);
    mean_contra = mean(contra_t_occ); ci_contra = prctile(contra_t_occ, [2.5 97.5]);
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
    end
    xticks(x_pos); xticklabels({'ipsi','contra'});
    ylabel('Bootstrapped t-statistic'); title('Ripple occurrence')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------------------ Extract bootstrapped stats ------------------
    for nBoot = 1:1000
        ipsi_t(nBoot)   = output(nBoot).t_stat{5};
        contra_t(nBoot) = output(nBoot).t_stat{6};
        pval_ipsi(nBoot)   = output(nBoot).pval{5};
        pval_contra(nBoot) = output(nBoot).pval{6};
        R2_ipsi(nBoot)     = output(nBoot).R2(5);
        R2_contra(nBoot)   = output(nBoot).R2(6);
    end

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
        sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    xlim([-2 15])
    ylim([-2 10])

    % ------------------ Scatter: contra ------------------
    nexttile
    x = contraHPC_cat(ripples_index);x = x(:);
    y = ipsiV1_cat(ripples_index);y = y(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.05)
    xlabel('contra HPC excitation'); ylabel('ipsi V1 depression')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12); hold on

    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);

    medR2 = prctile(R2_contra,50); medP = prctile(pval_contra,50);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    xlim([-2 15])
    ylim([-2 10])

    % ------------------ Bar Plot: Bootstrapped t-statistics ------------------
    mean_ipsi   = mean(ipsi_t);
    mean_contra = mean(contra_t);
    ci_ipsi     = prctile(ipsi_t, [2.5 97.5]);
    ci_contra   = prctile(contra_t, [2.5 97.5]);

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
    end
    xlim([0.5 2.5]); xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped t-statistic'); title('HPC excitation')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off


%%%%%%%%%%%%%%%%
    % ------------------ Extract bootstrapped stats ------------------
    for nBoot = 1:1000
        % Ripple power
        ipsi_t_power(nBoot)   = output(nBoot).t_stat{2};
        contra_t_power(nBoot) = output(nBoot).t_stat{3};
        R2_ipsi_power(nBoot)  = output(nBoot).R2(2);
        R2_contra_power(nBoot)= output(nBoot).R2(3);
        pval_ipsi_power(nBoot)= output(nBoot).pval{2};
        pval_contra_power(nBoot)= output(nBoot).pval{3};

        % Ripple occurrence
        ipsi_t_occ(nBoot)     = output(nBoot).t_stat{25};
        contra_t_occ(nBoot)   = output(nBoot).t_stat{26};
    end

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
    medR2 = prctile(R2_ipsi_power,50); medP = prctile(pval_ipsi_power,50);
    % medR2 = median(R2_ipsi_power); medP = median(pval_ipsi_power);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    ylim([-1 10])


    nexttile
    x = contra_ripple_power(:); y = ipsiV1_cat(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.2)
    xlabel('contra ripple power'); ylabel('ipsi V1 depression')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);
    medR2 = prctile(R2_contra_power,50); medP = prctile(pval_contra_power,50);
    % medR2 = median(R2_contra_power); medP = median(pval_contra_power);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)
    ylim([-1 10])
    % ------------------ Histogram: V1 depression with vs without ripples ------------------

    % --- IPSI
    nexttile
    histogram(ipsiV1_cat(ipsiRipples_cat == 1), -1:0.1:6, 'Normalization', 'probability'); hold on
    histogram(ipsiV1_cat(ipsiRipples_cat == 0), -1:0.1:6, 'Normalization', 'probability');
    [~, pVal_ks_ipsi] = kstest2(ipsiV1_cat(ipsiRipples_cat == 1), ipsiV1_cat(ipsiRipples_cat == 0));
    legend('With ipsi ripple', 'Without ipsi ripple', 'Box','off')
    xlabel('ipsi V1 depression'); ylabel('Proportion');
    text(6, 0.04, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize', 10, 'Color', 'k')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

    % --- CONTRA
    nexttile
    histogram(ipsiV1_cat(contraRipples_cat == 1), -1:0.1:6, 'Normalization', 'probability'); hold on
    histogram(ipsiV1_cat(contraRipples_cat == 0), -1:0.1:6, 'Normalization', 'probability');
    [~, pVal_ks_contra] = kstest2(ipsiV1_cat(contraRipples_cat == 1), ipsiV1_cat(contraRipples_cat == 0));
    legend('With contra ripple', 'Without contra ripple', 'Box','off')
    xlabel('ipsi V1 depression'); ylabel('Proportion');
    text(6, 0.04, sprintf('p = %.2e', pVal_ks_contra), 'FontSize', 10, 'Color', 'k')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)


    % ------------------ Bar plot: t-statistics (ripple power) ------------------
    mean_ipsi   = mean(ipsi_t_power); ci_ipsi = prctile(ipsi_t_power, [2.5 97.5]);
    mean_contra = mean(contra_t_power); ci_contra = prctile(contra_t_power, [2.5 97.5]);
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
    end
    xticks(x_pos); xticklabels({'ipsi','contra'});
    ylabel('Bootstrapped t-statistic'); title('Ripple power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    % ------------------ Bar plot: t-statistics (ripple occurrence) ------------------
    mean_ipsi   = mean(ipsi_t_occ); ci_ipsi = prctile(ipsi_t_occ, [2.5 97.5]);
    mean_contra = mean(contra_t_occ); ci_contra = prctile(contra_t_occ, [2.5 97.5]);
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
    end
    xticks(x_pos); xticklabels({'ipsi','contra'});
    ylabel('Bootstrapped t-statistic'); title('Ripple occurrence')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% V1 delta power predicted by HPC excitation

    % Extract bootstrapped stats
    for nBoot = 1:1000
        ipsi_t(nBoot)      = output(nBoot).t_stat{18};
        contra_t(nBoot)    = output(nBoot).t_stat{19};
        pval_ipsi(nBoot)   = output(nBoot).pval{18};
        pval_contra(nBoot) = output(nBoot).pval{19};
        R2_ipsi(nBoot)     = output(nBoot).R2(18);
        R2_contra(nBoot)   = output(nBoot).R2(19);
    end

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
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)

    % CONTRA
    nexttile
    x = contraHPC_cat(ripples_index);x = x(:);
    y = ipsi_Delta_peaks_zscore_UD(ripples_index);y = y(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.05)
    xlabel('contra HPC excitation'); ylabel('ipsi V1 delta power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    valid = ~isnan(x) & ~isnan(y);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(polyfit(x(valid), y(valid), 1), x_fit);
    medR2 = prctile(R2_contra,50); medP = prctile(pval_contra,50);
    % medR2 = median(R2_contra); medP = median(pval_contra);
    col = 'r'; if medP > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.2f\np = %.2e', medR2, medP), ...
        'FontSize', 10, 'Color', col)

    % T-stat bar
    nexttile
    barData = [mean(ipsi_t), mean(contra_t)];
    ci_ipsi = prctile(ipsi_t, [2.5 97.5]);
    ci_contra = prctile(contra_t, [2.5 97.5]);
    lowerCI = [barData(1) - ci_ipsi(1), barData(2) - ci_contra(1)];
    upperCI = [barData(1) - ci_ipsi(2), barData(2) - ci_contra(2)];
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5)
    end
    xticks(x_pos); xticklabels({'ipsi','contra'})
    ylabel('Bootstrapped t-statistic'); title('HPC excitation')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% V1 delta power predicted by ripple features

    % Extract stats
    for nBoot = 1:1000
        % Ripple power
        ipsi_t_power(nBoot)   = output(nBoot).t_stat{21};
        contra_t_power(nBoot) = output(nBoot).t_stat{22};
        R2_ipsi_power(nBoot)  = output(nBoot).R2(21);
        R2_contra_power(nBoot)= output(nBoot).R2(22);
        pval_ipsi_power(nBoot)= output(nBoot).pval{21};
        pval_contra_power(nBoot)= output(nBoot).pval{22};

        % Ripple occurrence
        ipsi_t_occ(nBoot)     = output(nBoot).t_stat{32};
        contra_t_occ(nBoot)   = output(nBoot).t_stat{33};
    end

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
    R2 = prctile(R2_ipsi_power, 50); p = prctile(pval_ipsi_power, 50);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.2f\np = %.2e', R2, p), ...
        'FontSize', 10, 'Color', col)

    nexttile
    x = contra_ripple_power(:); y = ipsi_Delta_peaks_zscore_UD(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.2)
    xlabel('contra ripple power'); ylabel('ipsi V1 delta power')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); hold on
    valid = ~isnan(x) & ~isnan(y);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(polyfit(x(valid), y(valid), 1), x_fit);
    R2 = prctile(R2_contra_power, 50); p = prctile(pval_contra_power, 50);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2*max(x_fit), 0.8*max(y_fit), sprintf('R^2 = %.2f\np = %.2e', R2, p), ...
        'FontSize', 10, 'Color', col)

    % Histogram: ripple occurrence
    nexttile
    histogram(ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 1), 'BinWidth', 0.2, 'Normalization','probability'); hold on
    histogram(ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 0), 'BinWidth', 0.2, 'Normalization','probability');
    [~, pVal_ks_ipsi] = kstest2(ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 1), ...
        ipsi_Delta_peaks_zscore_UD(ipsiRipples_cat == 0));
    legend('With ipsi ripple','Without','Box','off')
    xlabel('ipsi V1 delta power'); ylabel('Proportion');
    text(min(ipsi_Delta_peaks_zscore_UD), 0.08, sprintf('p = %.2e', pVal_ks_ipsi), 'FontSize',10)
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    nexttile
    histogram(ipsi_Delta_peaks_zscore_UD(contraRipples_cat == 1), 'BinWidth', 0.2, 'Normalization','probability'); hold on
    histogram(ipsi_Delta_peaks_zscore_UD(contraRipples_cat == 0), 'BinWidth', 0.2, 'Normalization','probability');
    [~, pVal_ks_contra] = kstest2(ipsi_Delta_peaks_zscore_UD(contraRipples_cat == 1), ...
        ipsi_Delta_peaks_zscore_UD(contraRipples_cat == 0));
    legend('With contra ripple','Without','Box','off')
    xlabel('ipsi V1 delta power'); ylabel('Proportion');
    text(min(ipsi_Delta_peaks_zscore_UD), 0.08, sprintf('p = %.2e', pVal_ks_contra), 'FontSize',10)
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off


    % Bar: t-stats
    nexttile
    barData = [mean(ipsi_t_power), mean(contra_t_power)];
    ci1 = prctile(ipsi_t_power, [2.5 97.5]);
    ci2 = prctile(contra_t_power, [2.5 97.5]);
    lowerCI = [barData(1) - ci1(1), barData(2) - ci2(1)];
    upperCI = [barData(1) - ci1(2), barData(2) - ci2(2)];
    x_pos = [1 2];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle','none','linewidth',1.5)
    end
    xticks(x_pos); xticklabels({'ipsi','contra'}); title('Ripple power')
    ylabel('Bootstrapped t-statistic')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    nexttile
    barData = [mean(ipsi_t_occ), mean(contra_t_occ)];
    ci1 = prctile(ipsi_t_occ, [2.5 97.5]);
    ci2 = prctile(contra_t_occ, [2.5 97.5]);
    lowerCI = [barData(1) - ci1(1), barData(2) - ci2(1)];
    upperCI = [barData(1) - ci1(2), barData(2) - ci2(2)];
    for i = 1:2
        hold on
        bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', ...
            'CData', customColors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
        errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), ...
            'k', 'linestyle','none','linewidth',1.5)
    end
    xticks(x_pos); xticklabels({'ipsi','contra'}); title('Ripple occurrence')
    ylabel('Bootstrapped t-statistic')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% ipsi ripple lag predicts V1 depression, delta power, and UP-DOWN lag

    % Load custom color scheme
    customColors = [0,90,50; 74,20,134]/256;

    % Extract stats
    for nBoot = 1:1000
        % V1 depression (model 7)
        t_V1depr(nBoot)   = output(nBoot).t_stat{7};
        R2_V1depr(nBoot)  = output(nBoot).R2(7);
        pval_V1depr(nBoot)= output(nBoot).pval{7};

        % V1 delta power (model 23)
        t_V1delta(nBoot)   = output(nBoot).t_stat{23};
        R2_V1delta(nBoot)  = output(nBoot).R2(23);
        pval_V1delta(nBoot)= output(nBoot).pval{23};

        % UP-DOWN lag (model 23 again)
        t_UPDOWN(nBoot)   = output(nBoot).t_stat{15};
        R2_UPDOWN(nBoot)  = output(nBoot).R2(15);
        pval_UPDOWN(nBoot)= output(nBoot).pval{15};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------- Plot 1: V1 depression -----------------------
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3];
    fig.Name = 'V1 depression vs ripple lag';

    % Scatter
    nexttile
    x = ipsi_ripple_lag(:); y = ipsiV1_cat(:);
    scatter(x, y, 20, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2); hold on
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);
    R2 = prctile(R2_V1depr, 50); p = prctile(pval_V1depr, 50);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    xlabel('ipsi ripple lag'); ylabel('V1 depression')
    text(0.2*max(x_fit), 0.8*max(y(valid)), sprintf('R^2 = %.2f\np = %.2e', R2, p), ...
        'FontSize', 10, 'Color', col)
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    % T-stat bar
    nexttile
    barData = mean(t_V1depr);
    ci = prctile(t_V1depr, [2.5 97.5]);
    bar(1, barData, 0.4, 'FaceColor', customColors(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    errorbar(1, barData, barData - ci(1), barData - ci(2), 'k', 'LineWidth', 1.5)
    xlim([0.5 1.5]); xticks(1); xticklabels({'V1 depression'})
    ylabel('Bootstrapped t-statistic')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------- Plot 2: V1 delta power -----------------------
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3];
    fig.Name = 'V1 delta power vs ripple lag';

    % Scatter
    nexttile
    x = ipsi_ripple_lag(:); y = ipsi_Delta_peaks_zscore_UD(:);
    scatter(x, y, 20, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2); hold on
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);
    R2 = prctile(R2_V1delta, 50); p = prctile(pval_V1delta, 50);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    xlabel('ipsi ripple lag'); ylabel('V1 delta power')
    text(0.2*max(x_fit), 0.8*max(y(valid)), sprintf('R^2 = %.2f\np = %.2e', R2, p), ...
        'FontSize', 10, 'Color', col)
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    % T-stat bar
    nexttile
    barData = mean(t_V1delta);
    ci = prctile(t_V1delta, [2.5 97.5]);
    bar(1, barData, 0.4, 'FaceColor', customColors(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    errorbar(1, barData, barData - ci(1), barData - ci(2), 'k', 'LineWidth', 1.5)
    xlim([0.5 1.5]); xticks(1); xticklabels({'V1 delta power'})
    ylabel('Bootstrapped t-statistic')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------- Plot 3: UP→DOWN lag -------------------------
    fig = figure('Color','w');
    fig.Position = [350 59 800 930/3];
    fig.Name = 'UP-DOWN lag vs ripple lag';

    % Scatter
    nexttile
    x = ipsi_ripple_lag(:); y = UP_DOWN_lag(:);
    scatter(x, y, 20, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.2); hold on
    valid = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(valid), y(valid), 1);
    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(coeffs, x_fit);
    R2 = prctile(R2_UPDOWN, 50); p = prctile(pval_UPDOWN, 50);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    xlabel('ipsi ripple lag'); ylabel('UP-DOWN transition lag')
    text(0.2*max(x_fit), 0.8*max(y(valid)), sprintf('R^2 = %.2f\np = %.2e', R2, p), ...
        'FontSize', 10, 'Color', col)
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off

    % T-stat bar
    nexttile
    barData = mean(t_UPDOWN);
    ci = prctile(t_UPDOWN, [2.5 97.5]);
    bar(1, barData, 0.4, 'FaceColor', customColors(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    errorbar(1, barData, barData - ci(1), barData - ci(2), 'k', 'LineWidth', 1.5)
    xlim([0.5 1.5]); xticks(1); xticklabels({'UP→DOWN lag'})
    ylabel('Bootstrapped t-statistic')
    set(gca,'TickDir','out','Box','off','Color','none','FontSize',12); grid off



    if exist('D:\corticohippocampal_replay')>0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay')>0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (150ms windows)'),[],'ContentType','image')


% 
% colour_lines = []
% colour_lines{1} = [
%     0.6289, 0.8477, 0.6055;
%     0.5617, 0.8172, 0.5672;
%     0.4945, 0.7867, 0.5289;
%     0.4273, 0.7562, 0.4906;
%     0.3602, 0.7257, 0.4523;
%     0.3056, 0.6654, 0.4066;
%     0.2510, 0.6051, 0.3609;
%     0.1965, 0.5448, 0.3152;
%     0.1419, 0.4844, 0.2695;
%     0.0874, 0.4241, 0.2238;
% ];
% 
% colour_lines{2} = [
%     0.7344, 0.7383, 0.8594;
%     0.6914, 0.6969, 0.8359;
%     0.6484, 0.6555, 0.8125;
%     0.6055, 0.6141, 0.7891;
%     0.5625, 0.5727, 0.7656;
%     0.5156, 0.5098, 0.7344;
%     0.4688, 0.4470, 0.7031;
%     0.4219, 0.3841, 0.6719;
%     0.3750, 0.3213, 0.6406;
%     0.3281, 0.2584, 0.6094
% ];
colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 

end
