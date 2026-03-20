function Bayesian_log_odds_RUN_validation = Bayesian_log_odds_ROC_test(decoded_Bayesian_RUN,decoded_Bayesian_RUN_shuffled,options,varargin)

% Default values
p = inputParser;
addParameter(p,'bin_selection_option',1,@isnumeric)
addParameter(p,'shuffle_option',1,@isnumeric)
addParameter(p,'plot_option',1,@isnumeric)
addParameter(p,'event_type','RUN',@ischar)
addParameter(p,'time_bin_size',0.1,@isnumeric)
addParameter(p,'good_bins',[],@isnumeric)


% assign parameters (either defaults or given)
parse(p,varargin{:});
plot_option = p.Results.plot_option;
shuffle_option = p.Results.shuffle_option;
event_type = p.Results.event_type;
time_bin_size= p.Results.time_bin_size;
good_bins= p.Results.good_bins;
positions = decoded_Bayesian_RUN.position_index;


if isempty(good_bins) == 1
    for pos = 1:max(decoded_Bayesian_RUN.position_index)
        labels=decoded_Bayesian_RUN.track_id(positions==pos);
        log_odds=decoded_Bayesian_RUN.position_log_odds(positions==pos);
        %     log_odds_shuffled
        %     log_odds_shuffled(i,:)
        % -------------------- STEP 6: ROC + Bootstrap --------------------
        fpr = 0:0.01:1;
        TPR_real = zeros(1000, length(fpr));
        AUC_real = zeros(1,1000);
        parfor i = 1:1000
            s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
            idx = randsample(s,length(labels), length(labels), true);
            %         [~, ~, ~, AUC_real(i)] = perfcurve(labels(idx), log_odds(idx),1,'XVals', fpr);
            [x,y,~,AUC_real(i)] = perfcurve(labels(idx), log_odds(idx),1, 'XVals', fpr);
            [C,ia,ic] = unique(x);
            y = y(ia);
            x = x(ia);

            TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
            TPR_real(i,:) = TPR_interp;
        end
        Bayesian_log_odds_RUN_validation.AUC_real(pos,:) = AUC_real;

        TPR_shuf = zeros(1000, length(fpr));
        AUC_shuf = zeros(1,1000);
        for i = 1:1000
            log_odds_shuffled = decoded_Bayesian_RUN_shuffled{i}.position_log_odds(positions==pos);
            s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
            idx = randsample(s,length(labels), length(labels), true);
            [x, y, ~, AUC_shuf(i)] = perfcurve(labels, log_odds_shuffled, 1, 'XVals', fpr);

            %         [x,y,~,AUC_shuf(i)] = perfcurve(labels(idx), log_odds(:),1, 'XVals', fpr);
            [C,ia,ic] = unique(x);
            y = y(ia);
            x = x(ia);

            TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
            TPR_shuf(i,:) = TPR_interp;
        end
        Bayesian_log_odds_RUN_validation.AUC_cell_ID_shuffle(pos,:) = AUC_shuf;


        TPR_shuf = zeros(1000, length(fpr));
        AUC_shuf = zeros(1,1000);
        for i = 1:1000
            log_odds_shuffled = decoded_Bayesian_RUN_shuffled{i}.position_log_odds(positions==pos);
            s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
            idx = randsample(s,length(labels), length(labels), true);
            %         [x, y, ~, AUC_shuf(i)] = perfcurve(labels, log_odds_shuffled, 1, 'XVals', fpr);

            [x,y,~,AUC_shuf(i)] = perfcurve(labels(idx), log_odds(:),1, 'XVals', fpr);
            [C,ia,ic] = unique(x);
            y = y(ia);
            x = x(ia);

            TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
            TPR_shuf(i,:) = TPR_interp;
        end
        Bayesian_log_odds_RUN_validation.AUC_shuffle(pos,:) = AUC_shuf;


        %     % -------------------- STEP 7: Plot ROC --------------------
        %     if plot_option
        %         title_text = sprintf('%s %s %s Bayesian log odds ROC validation %.0f ms bins', ...
        %             options.SUBJECT, options.SESSION, event_type, time_bin_size*1000);
        %         fig = figure('Name', title_text, 'Position', [200 100 475 880]); hold on;
        %
        %         nexttile
        %         CI_real = prctile(TPR_real, [0.5 99.5]);
        %         CI_shuf = prctile(TPR_shuf, [0.5 99.5]);
        %         hold on
        %         plot([0 1],[0 1],'--k')
        %         PLOT(1) = fill([fpr fliplr(fpr)], [CI_real(2,:) fliplr(CI_real(1,:))], ...
        %             [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        %         PLOT(2) = fill([fpr fliplr(fpr)], [CI_shuf(2,:) fliplr(CI_shuf(1,:))], ...
        %             'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        %         plot(fpr, mean(TPR_real,1), 'Color', [231,41,138]/256);
        %         plot(fpr, mean(TPR_shuf,1), 'k', 'LineWidth', 1.5);
        %         xlabel('True positive rate')
        %         ylabel('False Positive rate')
        %         %     xline(0.5, '--k'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
        %         title(title_text); legend(PLOT(1:2),{ 'Real', 'Shuffle'},'box', 'off');
        %         set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
        %
        %         nexttile
        %
        %         AUC = [AUC_real;AUC_shuf];
        %         bar_colors = [231,41,138; 0, 0, 0]/255;
        %         x_pos=[1 2];
        %         for i = 1:2
        %
        %             mean_AUC = mean(AUC(i,:),'omitnan');
        %             auc_errors = [mean_AUC-prctile(AUC(i,:), [0.5]),...  % lower error
        %                 prctile(AUC(i,:), [99.5])-mean_AUC];   % upper error
        %             hold on
        %             bar(x_pos(i), mean_AUC, 0.4, 'FaceAlpha',0.5, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
        %             errorbar(x_pos(i), mean_AUC, auc_errors(1), auc_errors(2), 'k', 'linestyle', 'none', 'linewidth', 1);
        %         end
        %
        %
        %         xticks([1 2]);
        %         xticklabels({'Real', 'Shuffled'});
        %         ylabel('AUC');
        %         ylim([0 1]);
        %         title('Mean AUC – Real vs Shuffled');
        %         set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
        %     end

    end

%     Bayesian_log_odds_RUN_validation.good_bins = find(prctile(Bayesian_log_odds_RUN_validation.AUC_real',[2.5]) > prctile(Bayesian_log_odds_RUN_validation.AUC_cell_ID_shuffle',[97.5]));
    Bayesian_log_odds_RUN_validation.good_bins = find(prctile(Bayesian_log_odds_RUN_validation.AUC_real',[2.5]) > prctile(Bayesian_log_odds_RUN_validation.AUC_shuffle',[97.5])...
        & prctile(Bayesian_log_odds_RUN_validation.AUC_real',[2.5])>0.6);
%     Bayesian_log_odds_RUN_validation.good_bins = good_idx;
else
     Bayesian_log_odds_RUN_validation.good_bins = good_bins;
end


is_good = ismember(decoded_Bayesian_RUN.position_index, Bayesian_log_odds_RUN_validation.good_bins);
good_idx = find(is_good);

% pos = [];
labels=decoded_Bayesian_RUN.track_id(good_idx);
log_odds=decoded_Bayesian_RUN.position_log_odds(good_idx);

fpr = 0:0.01:1;
TPR_real = zeros(1000, length(fpr));
AUC_real = zeros(1,1000);
parfor i = 1:1000
    s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
    idx = randsample(s,length(labels), length(labels), true);
    %         [~, ~, ~, AUC_real(i)] = perfcurve(labels(idx), log_odds(idx),1,'XVals', fpr);
    [x,y,~,AUC_real(i)] = perfcurve(labels(idx), log_odds(idx),1, 'XVals', fpr);
    [C,ia,ic] = unique(x);
    y = y(ia);
    x = x(ia);

    TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
    TPR_real(i,:) = TPR_interp;
end
Bayesian_log_odds_RUN_validation.AUC_real_good_bins = AUC_real;

TPR_shuf = zeros(1000, length(fpr));
AUC_shuf = zeros(1,1000);
for i = 1:1000
    log_odds_shuffled = decoded_Bayesian_RUN_shuffled{i}.position_log_odds(good_idx);
    s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
    idx = randsample(s,length(labels), length(labels), true);
            [x, y, ~, AUC_shuf(i)] = perfcurve(labels, log_odds_shuffled, 1, 'XVals', fpr);

%     [x,y,~,AUC_shuf(i)] = perfcurve(labels(idx), log_odds(:),1, 'XVals', fpr);
    [C,ia,ic] = unique(x);
    y = y(ia);
    x = x(ia);

    TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
    TPR_shuf(i,:) = TPR_interp;
end
Bayesian_log_odds_RUN_validation.AUC_cell_ID_shuffle_good_bins = AUC_shuf;

TPR_shuf = zeros(1000, length(fpr));
AUC_shuf = zeros(1,1000);
for i = 1:1000
%     log_odds_shuffled = decoded_Bayesian_RUN_shuffled{i}.position_log_odds();
    s = RandStream('mrg32k3a','Seed',i); % Set random seed for resampling
    idx = randsample(s,length(labels), length(labels), true);
    %         [x, y, ~, AUC_shuf(i)] = perfcurve(labels, log_odds_shuffled, 1, 'XVals', fpr);

    [x,y,~,AUC_shuf(i)] = perfcurve(labels(idx), log_odds(:),1, 'XVals', fpr);
    [C,ia,ic] = unique(x);
    y = y(ia);
    x = x(ia);

    TPR_interp = interp1(x, y, fpr, 'linear', 'extrap');  % force to match fpr
    TPR_shuf(i,:) = TPR_interp;
end
Bayesian_log_odds_RUN_validation.AUC_shuffle_good_bins = AUC_shuf;





% -------------------- STEP 3: Plot ROC --------------------
if plot_option
    title_text = sprintf('%s %s %s Bayesian log odds ROC validation %.0f ms bins', ...
        options.SUBJECT, options.SESSION, event_type, time_bin_size*1000);
    fig = figure('Name', title_text, 'Position', [200 100 475 880]); hold on;

    nexttile
    CI_real = prctile(TPR_real, [0.5 99.5]);
    CI_shuf = prctile(TPR_shuf, [0.5 99.5]);
    hold on
    plot([0 1],[0 1],'--k')
    PLOT(1) = fill([fpr fliplr(fpr)], [CI_real(2,:) fliplr(CI_real(1,:))], ...
        [231,41,138]/256, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    PLOT(2) = fill([fpr fliplr(fpr)], [CI_shuf(2,:) fliplr(CI_shuf(1,:))], ...
        'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(fpr, mean(TPR_real,1), 'Color', [231,41,138]/256);
    plot(fpr, mean(TPR_shuf,1), 'k', 'LineWidth', 1.5);
    xlabel('True positive rate')
    ylabel('False Positive rate')
    %     xline(0.5, '--k'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title(title_text); legend(PLOT(1:2),{ 'Real', 'Shuffle'},'box', 'off');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

    nexttile

    AUC = [AUC_real;AUC_shuf];
    bar_colors = [231,41,138; 0, 0, 0]/255;
    x_pos=[1 2];
    for i = 1:2

        mean_AUC = mean(AUC(i,:),'omitnan');
        auc_errors = [mean_AUC-prctile(AUC(i,:), [0.5]),...  % lower error
            prctile(AUC(i,:), [99.5])-mean_AUC];   % upper error
        hold on
        bar(x_pos(i), mean_AUC, 0.4, 'FaceAlpha',0.5, 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
        errorbar(x_pos(i), mean_AUC, auc_errors(1), auc_errors(2), 'k', 'linestyle', 'none', 'linewidth', 1);
    end


    xticks([1 2]);
    xticklabels({'Real', 'Shuffled'});
    ylabel('AUC');
    ylim([0 1]);
    title('Mean AUC – Real vs Shuffled');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
end

