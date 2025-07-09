function outputStruct = predict_ripples_by_DOWN_UP_V1_MUA(ipsi_V1_MUA,contra_V1_MUA,ipsi_HPC_MUA,contra_HPC_MUA,ipsi_ripples,contra_ripples,ripples_info,varargin)
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
addParameter(p, 'outputStruct', []);
addParameter(p, 'plot_option', 0);
addParameter(p, 'subject_id', 0);
addParameter(p, 'DOWN_UP_lag', []);
addParameter(p, 'DOWN_UP_index', []);
% addParameter(p, 'ripples_info', []);
% addParameter(p, 'ripples_power', []);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
nBoot = p.Results.nBoot;
outputStruct = p.Results.outputStruct;
subject_id = p.Results.subject_id;
DOWN_UP_lag = p.Results.DOWN_UP_lag;
DOWN_UP_index = p.Results.DOWN_UP_index;
% ripples_info = p.Results.ripples_info;
% ripples_power = p.Results.ripples_power;

ipsi_ripple_power = nan(length(DOWN_UP_index),1);
contra_ripple_power = nan(length(DOWN_UP_index),1);
% ripple_lag = nan(length(DOWN_UP_index),1);
ipsi_ripple_lag = nan(length(DOWN_UP_index),1);
contra_ripple_lag = nan(length(DOWN_UP_index),1);

if isempty(outputStruct)
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
    key_window_idx1 = find(timeVec < 0 & timeVec >= -0.05);
    key_window_idx2 = find(timeVec >= 0 & timeVec <= 0.05);

    ipsiV1_key = max(ipsi_V1_MUA(:, key_window_idx2)')'-min(ipsi_V1_MUA(:, key_window_idx1)')' ;
    contraV1_key = max(contra_V1_MUA(:, key_window_idx2)')' - min(contra_V1_MUA(:, key_window_idx1)')';
    ipsiHPC_key = max(ipsi_HPC_MUA(:, key_window_idx2)')'-min(ipsi_HPC_MUA(:, key_window_idx1)')' ;
    contraHPC_key = max(contra_HPC_MUA(:, key_window_idx2)')' - min(contra_HPC_MUA(:, key_window_idx1)')';
    %     combinedV1_key = mean((ipsi_V1_MUA_norm(:, key_window_idx)+contra_V1_MUA_norm(:, key_window_idx)),2,'omitnan');

    %%%% HPC ripples
    key_window_idx1 = find(timeVec_prob < 0 & timeVec_prob >= -0.05);
    key_window_idx2 = find(timeVec_prob >= 0 & timeVec_prob <= 0.2);

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

    ipsi_ripple_power(ipsiRipples_cat) = ripples_info.ipsi_first_ripples_power_UP(DOWN_UP_index(ipsiRipples_cat))';
    contra_ripple_power(contraRipples_cat) = ripples_info.contra_first_ripples_power_UP(DOWN_UP_index(contraRipples_cat))';
    ipsi_ripple_lag(ipsiRipples_cat) = abs(ripples_info.ipsi_first_ripples_lag_UP(DOWN_UP_index(ipsiRipples_cat))');
    contra_ripple_lag(contraRipples_cat) = abs(ripples_info.contra_first_ripples_lag_UP(DOWN_UP_index(contraRipples_cat))');



   
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
    tbl = table(ipsi_ripple_power,contra_ripple_power,ipsiHPC_cat, contraHPC_cat,ipsiV1_cat,contraV1_cat,zscore(DOWN_UP_lag), categorical(subject_id), 'VariableNames', {'ipsiPower','contraPower','ipsiOccurance','contraOccurance','ipsiV1','contraV1','DOWN_UP_lag','subjectID'});

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
%         index = intersect(index,ripples_index);
        tbl = table(normalize(ipsi_ripple_lag(index)),normalize(contra_ripple_lag(index)),normalize(ipsi_ripple_power(index)),normalize(contra_ripple_power(index)),ipsiRipples_cat(index), contraRipples_cat(index),...
            zscore(ipsiHPC_cat(index)), zscore(contraHPC_cat(index)),zscore(ipsiV1_cat(index)),zscore(contraV1_cat(index)),zscore(DOWN_UP_lag(index)), categorical(subject_id(index)),...
            'VariableNames', {'ipsiLag','contraLag','ipsiPower','contraPower','ipsiOccurance','contraOccurance','ipsiHPC','contraHPC','ipsiV1','contraV1','DOWN_UP_lag','subjectID'});

        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1*DOWN_UP_lag + contraV1*DOWN_UP_lag + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1 + contraV1 + DOWN_UP_lag + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ ipsiV1  + (1|subjectID)');
        %         ipsi_glme = fitglme(tbl, 'ipsiHPC ~ contraV1  + (1|subjectID)');

ipsi_glme = fitglme(tbl, 'ipsiOccurance ~ DOWN_UP_lag  + (1|subjectID)','Distribution', 'Binomial', 'Link', 'logit');
ipsi_glme = fitglme(tbl, 'ipsiOccurance ~ ipsiV1 + contraV1 + ipsiV1 * contraV1  + (1|subjectID)','Distribution', 'Binomial', 'Link', 'logit');
ipsi_glme = fitglme(tbl, 'contraOccurance ~ ipsiV1  + (1|subjectID)','Distribution', 'Binomial', 'Link', 'logit');
ipsi_glme = fitglme(tbl, 'contraOccurance ~ contraV1  + (1|subjectID)','Distribution', 'Binomial', 'Link', 'logit');

ipsi_glme = fitglme(tbl, 'ipsiOccurance ~ ipsiV1 * contraV1   + (1|subjectID)','Distribution', 'Binomial', 'Link', 'logit');
ipsi_glme = fitglme(tbl, 'contraOccurance ~ ipsiV1 * contraV1   + (1|subjectID)','Distribution', 'Binomial', 'Link', 'logit');



ipsi_glme = fitglme(tbl, 'ipsiPower ~ DOWN_UP_lag  + (1|subjectID)');

ipsi_glme = fitglme(tbl, 'ipsiPower ~ ipsiV1 + contraV1 + ipsiV1 * contraV1  + (1|subjectID)');

        %         ipsi_glme.Rsquared.Ordinary
        ipsi_R2(iBoot) = ipsi_glme.Rsquared.Ordinary;
        ipsi_t_stat(iBoot) = ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
        ipsi_pval(iBoot) = ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'DOWN_UP_lag'));
        %         ipsi_glme.Coefficients.tStat(strcmp(ipsi_glme.Coefficients.Name,'ipsiV1'));
        %         ipsi_glme.Coefficients.pValue(strcmp(ipsi_glme.Coefficients.Name,'contraV1'));

        %         contra_glme = fitglme(tbl, 'contraHPC ~ ipsiV1*DOWN_UP_lag + contraV1*DOWN_UP_lag + (1|subjectID)');
        %         contra_glme = fitglme(tbl, 'contraHPC ~ ipsiV1 + contraV1 + DOWN_UP_lag + (1|subjectID)');
        contra_glme = fitglme(tbl, 'contraPower ~ DOWN_UP_lag  + (1|subjectID)');
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







end

lag_thresholds = prctile(DOWN_UP_lag,[0:10:100]);

for ngroup = 1:length(lag_thresholds)-1

   mean_ipsi_power(ngroup) = nanmean(ipsi_ripple_power(DOWN_UP_lag>lag_thresholds(ngroup)&DOWN_UP_lag<lag_thresholds(ngroup+1)));
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



% 
% 
% 
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
%     0.0874, 0.4241, 0.2238
% ];
% colour_lines{1} = [161,217,155;116,196,118;65,171,93;35,139,69;0,90,50]/256;% 5 green for 
% colour_lines{2} = [188,189,220;158,154,200;128,125,186;106,81,163;74,20,134]/256;% 5 purple for 
% 
% nexttile
% for ngroup = 1:10
%     hold on
%     % plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     % plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% 
% 
% for n = 1:6
%     amp_threshold = prctile(ipsiV1_cat(subject_id==n),[0:20:100]);
% 
%     nexttile
%     for ngroup = 1:5
%         hold on
%         plot(nanmean( ipsi_ripples(subject_id==n &ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%         % plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%         % plot(nanmean(contra_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%         % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
%     end
% end
% 
% 
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
% sgtitle('ipsi DOWN UP magnitude (ripples)')
% 
% 
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
% sgtitle('ipsi DOWN UP magnitude without ripples')
% 


% 
% 
% figure
% nexttile
% lag_threshold = prctile(DOWN_UP_lag,[0:20:100]);
% 
% for ngroup = 1:length(lag_threshold)-1
%     hold on
%     %     plot(nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     plot(timeVec,nanmean(ipsi_V1_MUA(DOWN_UP_lag>lag_threshold(ngroup)&DOWN_UP_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
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
%     plot(timeVec,nanmean(ipsi_HPC_MUA(DOWN_UP_lag>lag_threshold(ngroup)&DOWN_UP_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
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
%     plot(timeVec_prob,nanmean(ipsi_ripples(DOWN_UP_lag>lag_threshold(ngroup)&DOWN_UP_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
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
%     plot(timeVec,nanmean(contra_V1_MUA(DOWN_UP_lag>lag_threshold(ngroup)&DOWN_UP_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
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
%     plot(timeVec,nanmean(contra_HPC_MUA(DOWN_UP_lag>lag_threshold(ngroup)&DOWN_UP_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
%     %     plot(nanmean(ipsi_V1_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
%     %     plot(nanmean(ipsi_ripples(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:)),'Color',colour_lines{1}(ngroup,:));
% 
% %     mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
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
%     plot(timeVec_prob,nanmean(contra_ripples(DOWN_UP_lag>lag_threshold(ngroup)&DOWN_UP_lag<lag_threshold(ngroup+1),:)),'Color',colour_lines{2}(ngroup,:));
% 
%     % mean_PSTH(ngroup,:) = nanmean(ipsi_HPC_MUA(ipsiV1_cat>amp_threshold(ngroup)&ipsiV1_cat<amp_threshold(ngroup+1),:));
% 
% end
% title('contra HPC Ripples')
% % legend('')
% sgtitle('ipsi-contra UP DOWN lag')
% 

end


