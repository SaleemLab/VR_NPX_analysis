function output = predict_UP_DOWN_V1_MUA_by_ripples(ipsi_V1_MUA,contra_V1_MUA,ipsi_HPC_MUA,contra_HPC_MUA,ipsi_ripples,contra_ripples,ripples_info,varargin)

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
plot_option = p.Results.plot_option;
% ripples_info = p.Results.ripples_info;
% ripples_power = p.Results.ripples_power;

if isempty(output)
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
    key_window_idx2 = find(timeVec < 0 & timeVec >= -0.05);
    key_window_idx1 = find(timeVec > 0 & timeVec <= 0.05);

%     ipsiV1_key = max(ipsi_V1_MUA(:, key_window_idx2)')'-min(ipsi_V1_MUA(:, key_window_idx1)')' ;
%     contraV1_key = max(contra_V1_MUA(:, key_window_idx2)')' - min(contra_V1_MUA(:, key_window_idx1)')';
%     ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')'-min(ipsi_HPC_MUA(:, key_window_idx1)')' ;
%     contraHPC_key = max(contra_HPC_MUA(:, key_window_idx2)')' - min(contra_HPC_MUA(:, key_window_idx1)')';

    ipsiV1_key = mean(ipsi_V1_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_V1_MUA(:, key_window_idx1)',"omitnan")' ;
    contraV1_key = mean(contra_V1_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_V1_MUA(:, key_window_idx1)',"omitnan")';
    ipsiHPC_key = mean(ipsi_HPC_MUA(:, key_window_idx2)',"omitnan")'-mean(ipsi_HPC_MUA(:, key_window_idx1)',"omitnan")' ;
    contraHPC_key = mean(contra_HPC_MUA(:, key_window_idx2)',"omitnan")' - mean(contra_HPC_MUA(:, key_window_idx1)',"omitnan")';
    %     combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

    %%%% HPC ripples
    key_window_idx2= find(timeVec_prob < 0 & timeVec_prob >= -0.05);
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
%     UP_DOWN_index = UP_DOWN_index(validIdx);

    %     ipsi_ripple_power(ipsiRipples_cat) = ripples_info.ipsi_first_ripples_power_UP(UP_DOWN_index(ipsiRipples_cat))';
    %     contra_ripple_power(contraRipples_cat) = ripples_info.contra_first_ripples_power_UP(UP_DOWN_index(contraRipples_cat))';
    %     ipsi_ripple_lag(ipsiRipples_cat) = abs(ripples_info.ipsi_first_ripples_lag_UP(UP_DOWN_index(ipsiRipples_cat))');
    %     contra_ripple_lag(contraRipples_cat) = abs(ripples_info.contra_first_ripples_lag_UP(UP_DOWN_index(contraRipples_cat))');


    ipsi_ripple_power = nan(length(ipsiV1_cat),1);
    contra_ripple_power = nan(length(ipsiV1_cat),1);
    % ripple_lag = nan(length(DOWN_UP_index),1);
    ipsi_ripple_lag = nan(length(ipsiV1_cat),1);
    contra_ripple_lag = nan(length(ipsiV1_cat),1);
    ipsi_time_from_last_ripple = nan(length(ipsiV1_cat),1);
    contra_time_from_last_ripple = nan(length(ipsiV1_cat),1);

    ipsi_ripple_power(ipsiRipples_cat) = ripples_info.ipsi_last_ripples_power_UP(ipsiRipples_cat)';
    contra_ripple_power(contraRipples_cat) = ripples_info.contra_last_ripples_power_UP(contraRipples_cat)';
    ipsi_ripple_lag(ipsiRipples_cat) = abs(ripples_info.ipsi_last_ripples_lag_UP(ipsiRipples_cat)');
    contra_ripple_lag(contraRipples_cat) = abs(ripples_info.contra_last_ripples_lag_UP(contraRipples_cat)');
    ipsi_time_from_last_ripple(ipsiRipples_cat) = abs(ripples_info.ipsi_time_from_last_ripples_UP(ipsiRipples_cat)');
    contra_time_from_last_ripple(contraRipples_cat) = abs(ripples_info.contra_time_from_last_ripples_UP(contraRipples_cat)');

    ripples_index = find(ipsiRipples_cat==1 | contraRipples_cat==1);
    
    
% [b, logl, H, stats] = coxphfit(UP_DOWN_lag(index), ipsi_time_from_last_ripple(index), 'Strata', subject_id(index));

    clear output
    %%%%%%%% All other models
    parfor iBoot = 1:1000
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index = randsample(s, ripples_index, size(ripples_index,1), true);

        tbl = table(normalize(ipsi_time_from_last_ripple(index)),normalize(contra_time_from_last_ripple(index)),normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), ...
            zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)), ...
            zscore(ipsiV1_cat(index)), zscore(contraV1_cat(index)), ...
            categorical(subject_id(index)), ...
            'VariableNames', {'ipsiDuration','contraDuration','ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC', ...
            'ipsiV1','contraV1','subjectID'});

        % Model list
        modelList = {
            'ipsiV1 ~ ipsiPower + contraPower + ipsiPower*contraPower + (1|subjectID)';
            'ipsiV1 ~ ipsiPower + (1|subjectID)';
            'ipsiV1 ~ contraPower + (1|subjectID)';

            'ipsiV1 ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)';
            'ipsiV1 ~ ipsiHPC + (1|subjectID)';
            'ipsiV1 ~ contraHPC + (1|subjectID)';

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
    parfor iBoot = 1:1000
        tic
        s = RandStream('philox4x32_10', 'Seed', iBoot);
        index = randsample(s, size(ipsiV1_cat,1), size(ipsiV1_cat,1), true);
%         index = intersect(index,ripples_index)
        tbl = table(normalize(ipsi_ripple_lag(index)), normalize(contra_ripple_lag(index)), ...
            normalize(ipsi_ripple_power(index)), normalize(contra_ripple_power(index)), ...
            ipsiRipples_cat(index), contraRipples_cat(index), ...
            zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)), ...
            zscore(ipsiV1_cat(index)), zscore(contraV1_cat(index)), ...
            categorical(subject_id(index)), ...
            'VariableNames', {'ipsiLag','contraLag','ipsiPower','contraPower', ...
            'ipsiOccurance','contraOccurance','ipsiHPC','contraHPC', ...
            'ipsiV1','contraV1','subjectID'});

        % Define model list (only one for now)
        modelList = {
            'ipsiV1 ~ ipsiOccurance  + (1|subjectID)';
            'ipsiV1 ~ contraOccurance  + (1|subjectID)';
            'ipsiV1 ~ ipsiOccurance+ contraOccurance  + (1|subjectID)';
            'ipsiV1 ~ ipsiHPC + contraHPC + ipsiHPC*contraHPC + (1|subjectID)'
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
lag_thresholds = prctile(UP_DOWN_lag,[0:10:100]);

for ngroup = 1:length(lag_thresholds)-1

   mean_ipsi_power(ngroup) = nanmean(ipsi_ripple_power(UP_DOWN_lag>lag_thresholds(ngroup)&UP_DOWN_lag<lag_thresholds(ngroup+1)));
end

fig = figure('Color','w');
fig.Position = [350 59 800 620];
fig.Name ='HPC MUA predicted by DOWN UP synchrony';


customColors = [0,90,50;74,20,134; 228,42,168 ]/256; % dark purple, dark green, magenta

nexttile
scatter(abs(DOWN_UP_lag),ipsi_ripple_power,'filled','MarkerFaceColor',customColors(1,:),'MarkerFaceAlpha',0.3)
ylim([5 30])
xlabel('DOWN UP transition ipsi-contra lags')
ylabel('ipsi ripple peak power (z)')

set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

nexttile
scatter(abs(DOWN_UP_lag),contra_ripple_power,'filled','MarkerFaceColor',customColors(2,:),'MarkerFaceAlpha',0.3)
ylim([5 30])
xlabel('DOWN UP transition ipsi-contra lags')
ylabel('contra ripple peak power (z)')
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)

nexttile
scatter(ipsi_ripple_power,contra_ripple_power,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.1)
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
xlabel('ipsi ripple peak power (z)')
ylabel('contra ripple peak power (z)')
xlim([5 30])
ylim([5 30])
% xlim([-2 12])
% ylim([-2 12])


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
% barColors = [0,90,50;65,171,93;74,20,134;128,125,186]/256;
barColors = [0,90,50;74,20,134]/256;
% Define custom x-positions
x_pos = [1, 2, 4, 5]; % spacing

%     b = bar(1:4, barData, 'FaceColor', 'flat');
for i = 1:2
    hold on
    bar(x_pos(i), barData(i), 0.4, 'FaceColor', 'flat', 'CData', barColors(i,:), 'EdgeColor', 'none');
    errorbar(x_pos(i), barData(i), lowerCI(i), upperCI(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
end

% Customize plot
xlim([0.5 2.5])
xticks(x_pos)
xticklabels({'ipsi HPC','ipsi HPC shuffled','contra HPC','contra HPC shuffled'})
ylabel('Bootstrapped t-statistic')
title('t-statistic of HPC excitation')
%     legend({'Data','Shuffled'}, 'Box','off')
set(gca,'TickDir','out','Box','off','Color','none','FontSize',12)
grid off
hold off


% 
% % 
% % colour_lines = []
% % colour_lines{1} = [
% %     0.6289, 0.8477, 0.6055;
% %     0.5617, 0.8172, 0.5672;
% %     0.4945, 0.7867, 0.5289;
% %     0.4273, 0.7562, 0.4906;
% %     0.3602, 0.7257, 0.4523;
% %     0.3056, 0.6654, 0.4066;
% %     0.2510, 0.6051, 0.3609;
% %     0.1965, 0.5448, 0.3152;
% %     0.1419, 0.4844, 0.2695;
% %     0.0874, 0.4241, 0.2238;
% % ];
% % 
% % colour_lines{2} = [
% %     0.7344, 0.7383, 0.8594;
% %     0.6914, 0.6969, 0.8359;
% %     0.6484, 0.6555, 0.8125;
% %     0.6055, 0.6141, 0.7891;
% %     0.5625, 0.5727, 0.7656;
% %     0.5156, 0.5098, 0.7344;
% %     0.4688, 0.4470, 0.7031;
% %     0.4219, 0.3841, 0.6719;
% %     0.3750, 0.3213, 0.6406;
% %     0.3281, 0.2584, 0.6094
% % ];
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 
% 
% % ipsiRipples_cat+contraRipples_cat > 0
% 
% figure
% nexttile
% amp_threshold = prctile(ipsiV1_cat(ipsiRipples_cat == 1),[0:20:100]);
% 
% % amp_threshold = [0:1:5];
% % amp_threshold = prctile(ipsiV1_cat(ipsiRipples_cat == 0),[0:20:100]);
% % amp_threshold = prctile(ipsiV1_cat,[0:20:100]);
% 
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(ipsi_V1_MUA(ipsiRipples_cat == 1&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi V1 MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     plot(timeVec,nanmean(ipsi_HPC_MUA(ipsiRipples_cat == 1&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi HPC MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec_prob,nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% legend({'Top 0-20%','Top 20-40%','Top 40-60%','Top 60-80%','Top 80-100%','Shuffled'})
% title('ipsi HPC Ripples')
% 
% 
% nexttile
% amp_threshold = prctile(ipsiV1_cat(contraRipples_cat == 1),[0:20:100]);
% % amp_threshold = prctile(ipsiV1_cat,[0:20:100]);
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(contra_V1_MUA(contraRipples_cat == 1&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra V1 MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     plot(timeVec,nanmean(contra_HPC_MUA(contraRipples_cat == 1&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra HPC MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec_prob,nanmean(contra_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% legend({'Top 0-20%','Top 20-40%','Top 40-60%','Top 60-80%','Top 80-100%','Shuffled'})
% title('contra HPC Ripples')
% sgtitle('ipsi UP DOWN magnitude (ripples)')
% 
% % 
% 
% figure
% nexttile
% amp_threshold = prctile(ipsiV1_cat(ipsiRipples_cat == 0),[0:20:100]);
% 
% % amp_threshold = [0:1:5];
% % amp_threshold = prctile(ipsiV1_cat(ipsiRipples_cat == 0),[0:20:100]);
% % amp_threshold = prctile(ipsiV1_cat,[0:20:100]);
% 
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(ipsi_V1_MUA(ipsiRipples_cat == 0&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi V1 MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     plot(timeVec,nanmean(ipsi_HPC_MUA(ipsiRipples_cat == 0&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi HPC MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec_prob,nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% legend({'Top 0-20%','Top 20-40%','Top 40-60%','Top 60-80%','Top 80-100%','Shuffled'})
% title('ipsi HPC Ripples')
% 
% 
% nexttile
% % amp_threshold = prctile(ipsiV1_cat(contraRipples_cat == 1),[0:20:100]);
% % amp_threshold = prctile(ipsiV1_cat,[0:20:100]);
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(contra_V1_MUA(contraRipples_cat == 0&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra V1 MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     plot(timeVec,nanmean(contra_HPC_MUA(contraRipples_cat == 0&ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra HPC MUA')
% 
% nexttile
% for ngroup = 1:length(amp_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec_prob,nanmean(contra_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% legend({'Top 0-20%','Top 20-40%','Top 40-60%','Top 60-80%','Top 80-100%','Shuffled'})
% title('contra HPC Ripples')
% sgtitle('ipsi UP DOWN magnitude without ripples')

% 
% 
% figure
% nexttile
% [~,index]=sort(ipsiV1_cat);imagesc(ipsi_V1_MUA(index,:));colorbar;clim([-1 1])
% nexttile
% [~,index]=sort(ipsiV1_cat);imagesc(movmedian(ipsi_HPC_MUA(index,:),10));colorbar;clim([-0.3 0.5])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% nexttile
% lag_threshold = prctile(UP_DOWN_lag,[0:20:100]);
% 
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(ipsi_V1_MUA(UP_DOWN_lag>lag_threshold(ngroup)&UP_DOWN_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi V1 MUA')
% 
% nexttile
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     plot(timeVec,nanmean(ipsi_HPC_MUA(UP_DOWN_lag>lag_threshold(ngroup)&UP_DOWN_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi HPC MUA')
% 
% nexttile
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec_prob,nanmean(ipsi_ripples(UP_DOWN_lag>lag_threshold(ngroup)&UP_DOWN_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('ipsi HPC Ripples')
% 
% 
% 
% 
% 
% nexttile
% % lag_threshold = prctile(UP_DOWN_lag,[0:10:100]);
% 
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(contra_V1_MUA(UP_DOWN_lag>lag_threshold(ngroup)&UP_DOWN_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra V1 MUA')
% 
% nexttile
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     plot(timeVec,nanmean(contra_HPC_MUA(UP_DOWN_lag>lag_threshold(ngroup)&UP_DOWN_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra HPC MUA')
% 
% % [p, tbl, stats] = kruskalwallis(data, group, 'off');  % 'off' suppresses plot
% 
% 
% nexttile
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec_prob,nanmean(contra_ripples(UP_DOWN_lag>lag_threshold(ngroup)&UP_DOWN_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra HPC Ripples')
% % legend('')
% sgtitle('ipsi-contra UP DOWN lag')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
for n = 1:6
    amp_threshold = prctile(ipsiV1_cat(ipsiRipples_cat == 1 & subject_id==n),[0:20:100]);

    nexttile
    for ngroup = 1:5
        hold on
%         plot(nanmean( ipsi_ripples(subject_id==n &ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
        plot(nanmean(ipsi_HPC_MUA(ipsiRipples_cat == 1 &ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
        % plot(nanmean(contra_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
        % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
    end
end

figure
colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
for n = 1:6
    amp_threshold = prctile(ipsiV1_cat(ipsiRipples_cat == 1 & subject_id==n),[0:20:100]);

    nexttile
    for ngroup = 1:5
        hold on
        plot(nanmean( ipsi_ripples(subject_id==n &ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%         plot(nanmean(ipsi_HPC_MUA(ipsiRipples_cat == 1 &ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
        % plot(nanmean(contra_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
        % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
    end
end
end
