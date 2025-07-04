function boot_output = plot_UP_survival_probability(ripple_feature, timeToTransition, ripple_counts,varargin)
    % Inputs:
    % ripple_feature{1,2}      - Cell with ipsi and contra ripple features (e.g., power, lag, etc.)
    % timeToTransition{1,2}    - Cell with time from last ripple to UP-DOWN transition (ipsi, contra)
    % ripple_counts{1,2}       - Cell with ripple counts (ipsi, contra)


    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'event_option', []);
    addParameter(p, 'count_option', []);
    addParameter(p, 'title_name', [], @ischar);
    addParameter(p, 'timebin', 0.015);
    addParameter(p, 'subject_id', []);


    parse(p, varargin{:});
    title_name = p.Results.title_name;
    event_option = p.Results.event_option;
    count_option = p.Results.count_option;
    timebin = p.Results.timebin;
    subject_id = p.Results.subject_id;

    binEdges = 0:timebin:0.3;
    binCenters = binEdges(1:end-1) + diff(binEdges)/2;

    % Colors: first row for low, second row for high
    colour_lines{1} = [116,196,118; 0,90,50] / 256;      % ipsi
    colour_lines{2} = [158,154,200; 74,20,134] / 256;    % contra

    group_name = {'ipsi','contra'};

    for hemi = 1:2  % 1 = ipsi, 2 = contra

        feature = ripple_feature{hemi};
        ttt = timeToTransition{hemi};
        counts = ripple_counts{hemi};

        % Get valid entries
        if isempty(count_option)
            valid_idx = intersect(find(~isnan(feature)), find(counts > 0));
        elseif count_option == 1
            valid_idx = intersect(find(~isnan(feature)), find(counts == 1));
        elseif count_option == 2
            valid_idx = intersect(find(~isnan(feature)), find(counts > 1));
        end
        feature_valid = feature(valid_idx)';
        ttt_valid = ttt(valid_idx);
        subject_used =  subject_id(valid_idx);

        if isempty(event_option)
            % Compute thresholds
            low_thresh = prctile(feature_valid, 25);
            high_thresh = prctile(feature_valid, 75);
            % else
            %
        end

        cox_time = ttt_valid;
        cox_event = ones(size(cox_time));  % All UPs transition
        feature_used = feature_valid;

        if ~isempty(subject_used)
            % Bootstrap Cox Regression (1000x)
            nBoot = 1000;
            boot_b = nan(nBoot, 1);
            boot_p = nan(nBoot, 1);
            parfor iBoot = 1:nBoot
                s = RandStream('mrg32k3a','Seed',iBoot);
                boot_idx = datasample(s, 1:length(cox_time), length(cox_time));
                feat_sample = feature_used(boot_idx);
                time_sample = cox_time(boot_idx);
                subj_sample = subject_used(boot_idx);
                try
                    [b_tmp, ~, ~, stats_tmp] = coxphfit(feat_sample, time_sample, ...
                        'Censoring', zeros(size(time_sample)), ...
                        'Strata', subj_sample);
                    boot_b(iBoot) = b_tmp;
                    boot_p(iBoot) = stats_tmp.p;
                catch
                    boot_b(iBoot) = NaN;
                    boot_p(iBoot) = NaN;
                end
            end

            % Clean and save
            boot_b = boot_b(~isnan(boot_b));
            boot_p = boot_p(~isnan(boot_p));
            boot_output(hemi).b = boot_b;
            boot_output(hemi).p = boot_p;
            boot_output(hemi).name = group_name{hemi};

            % Store summary for bar plot
            barData(hemi) = median(boot_b);
            lowerCI(hemi) = median(boot_b) - prctile(boot_b, 2.5);
            upperCI(hemi) = prctile(boot_b, 97.5) - median(boot_b);
            p50(hemi) = prctile(boot_p, 50);
        end
    end

    % Initialize outputs
    boot_output = struct('b', [], 'p', [], 'name', []);
    barColors = [0,90,50;74,20,134]/256;
    % barData = nan(1,2); lowerCI = nan(1,2); upperCI = nan(1,2); p50 = nan(1,2);

    fig = figure
    fig.Position = [350 59 800/2 930/3];

    if ~isempty(title_name)
        fig.Name = [title_name,' COX regression'];
    end

    x_pos = [1 2];
    for hemi = 1:2
        hold on
        bar(x_pos(hemi), barData(hemi), 0.4, 'FaceColor', 'flat', ...
            'CData', barColors(hemi,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        errorbar(x_pos(hemi), barData(hemi), lowerCI(hemi), upperCI(hemi), ...
            'k', 'linestyle', 'none', 'linewidth', 1.5);
        text(x_pos(hemi), barData(hemi) + upperCI(hemi) + 0.0001, ...
            sprintf('p_{50%%} = %.3e', p50(hemi)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end

    xlim([0.5 2.5])
    xticks(x_pos)
    xticklabels({'ipsi','contra'})
    ylabel('Cox coefficient (b)')
    title('Bootstrap Cox regression by hemisphere')
    set(gca, 'TickDir', 'out', 'Box', 'off', 'Color', 'none', 'FontSize', 12)




    fig = figure
    fig.Position = [700 120 790 650];

    if ~isempty(title_name)
        fig.Name = title_name;
    end


    for hemi = 1:2  % 1 = ipsi, 2 = contra
        subplot_offset = (hemi - 1) * 2;
        feature = ripple_feature{hemi};
        ttt = timeToTransition{hemi};
        counts = ripple_counts{hemi};

        % Get valid entries
        if isempty(count_option)
            valid_idx = intersect(find(~isnan(feature)), find(counts > 0));
        elseif count_option == 1
            valid_idx = intersect(find(~isnan(feature)), find(counts == 1));
        elseif count_option == 2
            valid_idx = intersect(find(~isnan(feature)), find(counts > 1));
        end
        feature_valid = feature(valid_idx)';
        ttt_valid = ttt(valid_idx);

        if isempty(event_option)
            % Compute thresholds
            low_thresh = prctile(feature_valid, 25);
            high_thresh = prctile(feature_valid, 75);
        % else
        % 
        end


        for n = 1:2  % 1 = low, 2 = high

            if isempty(event_option)
                if n == 1
                    event_id = find(feature_valid < low_thresh);
                else
                    event_id = find(feature_valid > high_thresh);
                end
            else
                % event_id = event_index{n};
            end



            % Bootstrapped survival curves
            y_boot = [];
            for iBoot = 1:1000
                s = RandStream('mrg32k3a','Seed',iBoot);
                sampled_times = ttt_valid(datasample(s, event_id, length(event_id)));
                [temp_y, temp_x] = ecdf(sampled_times, 'Function', 'survivor');
                temp_x(isnan(temp_y)) = []; temp_y(isnan(temp_y)) = [];
                [temp_x, idx] = unique(temp_x);
                temp_y = temp_y(idx);
                y_boot(iBoot,:) = interp1(temp_x, temp_y, binCenters, 'previous', 'extrap')';
            end

            % Select color for this hemi/quartile
            col = colour_lines{hemi}(n,:);

            % Plot survival
            x = binCenters';
            y = mean(y_boot)';
            LCI = prctile(y_boot, 2.5)';
            UCI = prctile(y_boot, 97.5)';
            subplot(2,2,1 + subplot_offset)
            hold on;
            plot(x, y, 'Color', col);
            patch([x' fliplr(x')], [UCI' fliplr(LCI')], col, ...
                'FaceAlpha', 0.3, 'LineStyle', 'none');
            ylabel('Survival probability of UP')
            xlabel('Time (s)')
            xlim([0 0.3])
            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)

            % Plot derivative
            subplot(2,2,2 + subplot_offset)
            hold on;
            y_diff = mean(diff(y_boot')')';
            LCI = prctile(diff(y_boot')', 2.5)';
            UCI = prctile(diff(y_boot')', 97.5)';
            y_diff = [0; y_diff];
            LCI = [0; LCI];
            UCI = [0; UCI];
            plot(x, y_diff, 'Color', col);
            P(n) = patch([x' fliplr(x')], [UCI' fliplr(LCI')], col, ...
                'FaceAlpha', 0.3, 'LineStyle', 'none');
            ylabel('Slope of survival probability');
            xlabel('Time (s)')
            xlim([0 0.3])
            set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 12)
        end

        legend(P(1:2),{'low','high'},'box','off','Location','southeast')
    end
end
