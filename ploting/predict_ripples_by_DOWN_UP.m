function output = predict_ripples_by_DOWN_UP_synchrony(ipsi_V1_MUA,contra_V1_MUA,ipsi_HPC_MUA,contra_HPC_MUA,ipsi_ripples,contra_ripples,UP_DOWN_info,varargin)
% SETTINGS
p = inputParser;
addParameter(p, 'nPerm', 1000);
addParameter(p, 'nBoot', 1000);
addParameter(p, 'output', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
% addParameter(p, 'DOWN_UP_lag', []);
% addParameter(p, 'DOWN_UP_index', []);
% addParameter(p, 'ripples_info', []);
% addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
output = p.Results.output;
subject_id = p.Results.subject_id;
% DOWN_UP_lag = p.Results.DOWN_UP_lag;
% DOWN_UP_index = p.Results.DOWN_UP_index;
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
%%%%%% Windows
key_window_idx1 = find(timeVec < 0 & timeVec >= -0.05);
key_window_idx2 = find(timeVec > 0 & timeVec <= 0.1);

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
key_window_idx2 = find(timeVec_prob >= 0 & timeVec_prob <= 0.1);

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
    parfor iBoot = 1:nBoot
        tic
        s = RandStream('mrg32k3a','Seed',iBoot); % Set random seed for resampling
        index = randsample(s, size(DOWN_UP_index,2), size(DOWN_UP_index,2), true);
        %         index = intersect(index,ripples_index)
        tbl = table(normalize(ipsi_ripple_duration(index)),normalize(contra_ripple_duration(index)),normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), normalize(ipsiHPC_cat(index)),normalize(contraHPC_cat(index)), ...
            normalize(ipsi_Delta_peaks_zscore_DU(index))',  normalize(contra_Delta_peaks_zscore_DU(index))', ...
            normalize(DOWN_duration(index)),ipsi_spindles(index),contra_spindles(index), ...
            categorical(subject_id(index)), ...
            'VariableNames', {'ipsiDuration','contraDuration','ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC' ...
            'ipsiDelta','contraDelta','downDuration','ipsiSpindle','contraSpindle','subjectID'});

        % Define model list (only one for now)
        modelList = {

        'ipsiOccurance ~ ipsiDelta + (1|subjectID)';
        'contraOccurance ~ ipsiDelta + (1|subjectID)';

        'ipsiOccurance ~ downDuration + (1|subjectID)';
        'contraOccurance ~ downDuration + (1|subjectID)';

        'ipsiOccurance ~ ipsiSpindle + (1|subjectID)';
        'contraOccurance ~ ipsiSpindle + (1|subjectID)';

        'ipsiOccurance ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID)';
        'contraOccurance ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID)';

        };

        nModels = numel(modelList);

        % Preallocate local outputs
        local_b = cell(nModels,1);
        local_R2 = zeros(nModels,1);
        local_R2_McFadden = zeros(nModels,1);
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
        output2(iBoot).b = local_b;
        output2(iBoot).R2 = local_R2;
        % output2(iBoot).local_R2_McFadden = local_R2_McFadden;
        output2(iBoot).t_stat = local_tstat;
        output2(iBoot).pval = local_pval;
        output2(iBoot).model = local_model;
        output2(iBoot).variable = local_variable;
        output2(iBoot).type = repmat({'All'}, nModels, 1);
        toc
    end


    %%%%%%%% All other models
    parfor iBoot = 1:nBoot
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        % index = randsample(s, ripples_index, size(ripples_index,1), true);
        index = randsample(s, size(DOWN_UP_index,2), size(DOWN_UP_index,2), true);
        tbl = table(normalize(ipsi_ripple_duration(index)),normalize(contra_ripple_duration(index)),normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), normalize(ipsiHPC_cat(index)),normalize(contraHPC_cat(index)), ...
            normalize(ipsi_Delta_peaks_zscore_DU(index))',  normalize(contra_Delta_peaks_zscore_DU(index))', ...
            normalize(DOWN_duration(index)),ipsi_spindles(index),contra_spindles(index), ...
            categorical(subject_id(index)), ...
            'VariableNames', {'ipsiDuration','contraDuration','ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC' ...
            'ipsiDelta','contraDelta','downDuration','ipsiSpindle','contraSpindle','subjectID'});


        % Model list
        modelList = {

        'ipsiHPC ~ ipsiDelta + (1|subjectID)';
        'contraHPC ~ ipsiDelta + (1|subjectID)';

        'ipsiHPC ~ downDuration + (1|subjectID)';
        'contraHPC ~ downDuration + (1|subjectID)';

        'ipsiHPC ~ ipsiSpindle + (1|subjectID)';
        'contraHPC ~ ipsiSpindle + (1|subjectID)';

        'ipsiPower ~ ipsiDelta + (1|subjectID)';
        'contraPower ~ ipsiDelta + (1|subjectID)';

        'ipsiPower ~ downDuration  +(1|subjectID)';
        'contraPower ~ downDuration + (1|subjectID)';

        'ipsiPower ~ ipsiSpindle + (1|subjectID)';
        'contraPower ~ ipsiSpindle + (1|subjectID)';

        'ipsiPower ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID)';
        'contraPower ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID)';

        'ipsiHPC ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID)';
        'contraHPC ~ ipsiDelta + downDuration + ipsiSpindle + (1|subjectID) + (1|subjectID)';

        %             'ipsiHPC ~ DOWN_UP_lag + (1|subjectID)';
        %             'contraHPC ~ DOWN_UP_lag + (1|subjectID)';

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
    %% 6. PLOTTING


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract bootstrapped stats
    % Extract bootstrapped stats
    for nBoot = 1:1000
        ipsi_t(nBoot)   = output(nBoot).t_stat{10};
        contra_t(nBoot) = output(nBoot).t_stat{11};
        pval_ipsi(nBoot)   = output(nBoot).pval{10};
        pval_contra(nBoot) = output(nBoot).pval{11};
        R2_ipsi(nBoot)     = output(nBoot).R2(10);
        R2_contra(nBoot)   = output(nBoot).R2(11);
    end

    % Plot setup
    fig = figure('Color', 'w');
    fig.Position = [350 59 800 600];
    fig.Name = 'UP→DOWN lag predicted by HPC excitation';
    customColors = [0,90,50; 74,20,134]/256;

    % --- IPSI
    nexttile
    x = ipsiHPC_cat(:);
    y = UP_DOWN_lag(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(1,:), 'MarkerFaceAlpha', 0.05)
    xlabel('ipsi HPC excitation'); ylabel('UP→DOWN transition lag')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12); hold on

    % Linear fit
    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'k', 'LineWidth', 2)

    % Mixed-effects model
    tbl = table(x(validIdx), y(validIdx), 'VariableNames', {'excitation','lag'});
    lme = fitlme(tbl, 'lag ~ excitation');
    R2 = lme.Rsquared.Ordinary;
    p = lme.Coefficients.pValue(2);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', R2, p), 'FontSize', 10, 'Color', col)

    % --- CONTRA
    nexttile
    x = contraHPC_cat(:);
    y = UP_DOWN_lag(:);
    scatter(x, y, 'filled', 'MarkerFaceColor', customColors(2,:), 'MarkerFaceAlpha', 0.05)
    xlabel('contra HPC excitation'); ylabel('UP→DOWN transition lag')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12); hold on

    validIdx = ~isnan(x) & ~isnan(y);
    coeffs = polyfit(x(validIdx), y(validIdx), 1);
    x_fit = linspace(min(x(validIdx)), max(x(validIdx)), 100);
    y_fit = polyval(coeffs, x_fit);

    tbl = table(x(validIdx), y(validIdx), 'VariableNames', {'excitation','lag'});
    lme = fitlme(tbl, 'lag ~ excitation');
    R2 = lme.Rsquared.Ordinary;
    p = lme.Coefficients.pValue(2);
    col = 'r'; if p > 0.05, col = 'k'; end
    plot(x_fit, y_fit, 'Color', col, 'LineWidth', 2)
    text(0.2 * max(x_fit), 0.8 * max(y_fit), ...
        sprintf('R^2 = %.2f\np = %.2e', R2, p), 'FontSize', 10, 'Color', col)


    % Compute mean and CI
    mean_ipsi   = mean(ipsi_t);
    mean_contra = mean(contra_t);
    ci_ipsi     = prctile(ipsi_t, [2.5 97.5]);
    ci_contra   = prctile(contra_t, [2.5 97.5]);

    barData  = [mean_ipsi, mean_contra];
    lowerCI  = [mean_ipsi - ci_ipsi(1), mean_contra - ci_contra(1)];
    upperCI  = [mean_ipsi - ci_ipsi(2), mean_contra - ci_contra(2)];

    % Plot
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
    ylabel('Bootstrapped t-statistic'); title('t-stat: HPC → UP→DOWN lag')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)
    grid off







    if exist('D:\corticohippocampal_replay')>0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay')>0
        analysis_folder = 'P:\corticohippocampal_replay';
    end
    save_all_figures(fullfile(analysis_folder,'V1-HPC bilateral interaction','mixed effect regression (full windows)'),[],'ContentType','image')


end
end

