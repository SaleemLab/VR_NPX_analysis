function [normalised_raw_matrix PPvector shuffled_globalRemap_PPvector shuffled_rateRemap_PPvector] = plot_place_cell_map_correlation(clusters,place_fields_all,place_fields_even,place_fields_odd,position,lap_times,options)
% Function to plot spatial cell activity
% Plot 1 is each cell's spatial ratemap for each lap
% Plot 2 is each cell's spatial tuning stability by plotting place field
% calculated by using all laps, odd no laps and even no laps.
ROOTPATH = options.ROOTPATH;

if isfield(options,'region')
    region = options.region;
else
    region = [];
end

if isfield(options,'probe_combined') && options.probe_combined == 1
    nprobe = 'combined';
else
    nprobe = options.probe_no;
    probe_hemisphere = options.probe_hemisphere;
    probe_hemisphere_text = {'left','right'}
end

% probe_hemisphere_text = {'left','right'}

x_bins_width = 10;


%Create linearised position matrix
index_track = [];
sorted_cells= [];

for type = 1:3

    if type == 3
        place_fields = place_fields_all;
        %         sgtitle('Whole')
    elseif type == 1
        place_fields = place_fields_odd;
        %         sgtitle('Odd laps')
    elseif type == 2
        place_fields = place_fields_even;
        %         sgtitle('Even laps')
    end

    for track_id = 1:2
        %                 max_FR = max([place_fields.track(1).raw_peak(place_fields_all.good_place_cells_LIBERAL);...
        %                     place_fields.track(1).raw_peak(place_fields_all.good_place_cells_LIBERAL)])';
        max_FR = place_fields.track(track_id).raw_peak(place_fields_all.good_place_cells_LIBERAL)';
        raw_matrix = cat(1,place_fields.track(track_id).raw{place_fields_all.good_place_cells_LIBERAL});
        normalised_raw_matrix{type}{track_id} = bsxfun(@rdivide, raw_matrix, max_FR);
        normalised_raw_matrix{type}{track_id}(isnan(normalised_raw_matrix{type}{track_id})) = 0;
        [~,index_track{type}(track_id,:)] = max(normalised_raw_matrix{type}{track_id},[],2);
        %     unsorted_cells(track_id,:) =
        [~,sorted_cells{type}(track_id,:)] = sort(index_track{type}(track_id,:));
        %     ordered_matrix = normalised_raw_matrix(new_order,:);
    end
end

laps_pairs = [1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2;1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2];
track_pairs = [1 1;1 1;1 2;1 2;1 1;1 1;1 2;1 2;2 1;2 1;2 2;2 2;2 1;2 1;2 2;2 2];
laps_type_text = {'odd','even','all'};
%plot heat map position
fig = figure
fig.Position = [500 70 1300 930];

if isempty(region)
    if ischar(nprobe)
        sgtitle(sprintf('%s %s Even vs Odd laps place cell maps probes %s',options.SUBJECT,options.SESSION,nprobe))
        fig.Name = sprintf('%s %s Even vs Odd laps place cell maps probes %s',options.SUBJECT,options.SESSION,nprobe);
    else
        sgtitle(sprintf('%s %s Even vs Odd laps place cell maps probe %s',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s Even vs Odd laps place cell maps probe %s',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere});
    end
else
    if ischar(nprobe)
        sgtitle(sprintf('%s %s Even vs Odd laps %s cell maps probes %s',options.SUBJECT,options.SESSION,region,nprobe))
        fig.Name = sprintf('%s %s Even vs Odd laps %s cell maps probes %s',options.SUBJECT,options.SESSION,region,nprobe);
    else
        sgtitle(sprintf('%s %s Even vs Odd laps %s cell maps probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}))
       fig.Name = sprintf('%s %s Even vs Odd laps %s cell maps probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});
    end
end

for nplot = 1:16
    subplot(4,4,nplot)
    ordered_matrix =  normalised_raw_matrix{laps_pairs(nplot,1)}{track_pairs(nplot,1)}(sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:),:);
    imagesc(ordered_matrix);
    xlabel('Position');
    ylabel('Cell ID');
    yt=sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:);
    set(gca,'ytick',1:length(sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:)));
    set(gca,'yticklabel',yt);
    axis xy;
    h = colorbar;
    set(get(h,'label'),'string','Normalised firing rate');
    clim([0.3 1])

    xticks([30 50 70 90 110 140]/10)
    xticklabels([30 50 70 90 110 140])
    xline(10,'r')

    title(sprintf('%s laps T%i sorted by %s laps T%i',...
        laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
end
set(gca,"TickDir","out",'box', 'off','Color','none')

%         title(sprintf('%s %s Even vs Odd laps place cell maps probe %i',options.SUBJECT,options.SESSION,nprobe))


% Calculate cell population vector for each track comparison by correlating each position bin between tracks
shuffled_globalRemap_PPvector.population_vector = [];
shuffled_globalRemap_PPvector.pval = [];
shuffled_rateRemap_PPvector.population_vector =[];
shuffled_rateRemap_PPvector.pval = [];
laps_pairs = [1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2;1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2];
track_pairs = [1 1;1 1;1 2;1 2;1 1;1 1;1 2;1 2;2 1;2 1;2 2;2 2;2 1;2 1;2 2;2 2];
PPvector = [];

% Original raw rate map population vector correlation
for i = 1 : length(laps_pairs)
    for j = 1 : size( normalised_raw_matrix{1}{track_id},2) % for each position bin
        % Correlation to global remapping shuffle
        [Grho,Gpval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), normalised_raw_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j)); %corr between position bins across cells
        PPvector.population_vector(j,i) = Grho;
        PPvector.pval(j,i) = Gpval;
    end
end

for nshuffle = 1 : 1000
    cellID_shuffled_matrix = []; rate_remap_matrix= [];

    % Create shuffled matrix for global remapping: shuffle cell ID & one for rate remapping: multiply each ratemap with a random num within 0.1 to 1.5
    for track_id = 1 : length(place_fields_all.track)
        s = RandStream('mrg32k3a','Seed',nshuffle*1000); % Set random seed for resampling
        cellID_shuffled_matrix{1}{track_id} = normalised_raw_matrix{1}{track_id}(randperm(s,size(normalised_raw_matrix{1}{track_id},1)),:);
        s = RandStream('mrg32k3a','Seed',nshuffle*1000+1); % Set random seed for resampling
        cellID_shuffled_matrix{2}{track_id} = normalised_raw_matrix{2}{track_id}(randperm(s,size(normalised_raw_matrix{2}{track_id},1)),:);

        s = RandStream('mrg32k3a','Seed',nshuffle*1000+2); % Set random seed for resampling
        rate_change = 0.1 + (1.5-0.1) .* rand(s,size(normalised_raw_matrix{1}{track_id},1),1); %create a vector with random values from 0.1 to 1.5
        rate_remap_matrix{1}{track_id}  = normalised_raw_matrix{1}{track_id}.*rate_change; % multiply ratemaps by the correspondin random num

        s = RandStream('mrg32k3a','Seed',nshuffle*1000+3); % Set random seed for resampling
        rate_change = 0.1 + (1.5-0.1) .* rand(s,size(normalised_raw_matrix{2}{track_id},1),1); %create a vector with random values from 0.1 to 1.5
        rate_remap_matrix{2}{track_id}  = normalised_raw_matrix{2}{track_id}.*rate_change; % multiply ratemaps by the correspondin random num
    end

    curr_size = size(shuffled_rateRemap_PPvector.population_vector,1);

    for i = 1 : length(laps_pairs)
        for j = 1 : size( normalised_raw_matrix{1}{track_id},2) % for each position bin
            % Correlation to global remapping shuffle
            [Grho,Gpval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), cellID_shuffled_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j)); %corr between position bins across cells
            shuffled_globalRemap_PPvector.population_vector(curr_size+j,i) = Grho; %each column is a comparison, each row a position bin
            shuffled_globalRemap_PPvector.pval(curr_size+j,i) = Gpval;
            % Correlation to rate remapping shuffle
            [Rrho,Rpval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), rate_remap_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j));
            shuffled_rateRemap_PPvector.population_vector(curr_size+j,i) = Rrho; %each column is a comparison, each row a position bin
            shuffled_rateRemap_PPvector.pval(curr_size+j,i) = Rpval;
        end
    end
end

%         % Run Kruskal-Wallis
%         all_sig_diff_idx = []; sig_diff_idx= [];
%         [pv3,~,stats3] = kruskalwallis(combined,[],'off');
%         if pv3 < 0.05
%             [all_c,~,~,~] = multcompare(stats3,'ctype','dunn-sidak','Display','off'); % if anova pval is < 0.05, run multiple comparisons
%             protocol(p).all_multiple_comparisons = all_c;
%             all_sig_diff_idx = find(all_c(:,6)<0.05);
%         end

fig = figure
fig.Position = [500 70 1300 930];
if isempty(region)
    if ischar(nprobe)
        sgtitle(sprintf('%s %s place cell PV correlation cumulative frequency comparision probes %s.pdf',options.SUBJECT,options.SESSION,nprobe))
        fig.Name = sprintf('%s %s place cell PV correlation cumulative frequency comparision probes %s.pdf',options.SUBJECT,options.SESSION,nprobe);
    else
        sgtitle(sprintf('%s %s place cell PV correlation cumulative frequency comparision probe %s.fig',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s place cell PV correlation cumulative frequency comparision probe %s.fig',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere});
    end
else
    if ischar(nprobe)
        sgtitle(sprintf('%s %s %s cell PV correlation cumulative frequency comparision probes %s',options.SUBJECT,options.SESSION,region,nprobe));
        fig.Name = sprintf('%s %s %s cell PV correlation cumulative frequency comparision probes %s',options.SUBJECT,options.SESSION,region,nprobe);
    else
        sgtitle(sprintf('%s %s %s cell PV correlation cumulative frequency comparision probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s %s cell PV correlation cumulative frequency comparision probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});
    end
end

for nplot = 1:16
    if nplot == 1| nplot == 6 | nplot == 11| nplot == 16
        continue
    end
    subplot(4,4,nplot)
    h(1) = cdfplot(PPvector.population_vector(:,nplot));
    hold on
    h(2) = cdfplot(shuffled_globalRemap_PPvector.population_vector(:,nplot));
    set( h(1), 'LineStyle', '-', 'Color', 'r');
    set( h(2), 'LineStyle', '-.', 'Color', 'k');
    xlabel('PV correlation')
    ylabel('Cumulative frequency')
    legend('Original','Cell ID shuffled','Color','none')
    title(sprintf('%s laps T%i sorted by %s laps T%i',...
        laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
    set(gca,"TickDir","out",'box', 'off','Color','none')
end


% Place cell stability and remapping
fig = figure
fig.Position = [500 70 1300 930];
if isempty(region)
    if ischar(nprobe)
        sgtitle(sprintf('%s %s PV correlation comparision for place cell stability and remapping probes %s',options.SUBJECT,options.SESSION,nprobe))
        fig.Name = sprintf('%s %s PV correlation comparision for place cell stability and remapping probes %s',options.SUBJECT,options.SESSION,nprobe);
    else
        sgtitle(sprintf('%s %s PV correlation comparision for place cell stability and remapping probe %s',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s PV correlation comparision for place cell stability and remapping probe %s',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere});
    end
else
    if ischar(nprobe)
        sgtitle(sprintf('%s %s PV correlation comparision for place cell stability and remapping probes %s',options.SUBJECT,options.SESSION,region,nprobe));
        fig.Name = sprintf('%s %s PV correlation comparision for place cell stability and remapping probes %s',options.SUBJECT,options.SESSION,region,nprobe);
    else
        sgtitle(sprintf('%s %s PV correlation comparision for %s cell stability and remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}));
        fig.Name = sprintf('%s %s PV correlation comparision for %s cell stability and remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});
    end
end

for nplot = 1:16
    if nplot == 1| nplot == 6 | nplot == 11| nplot == 16
        continue
    end
    subplot(4,4,nplot)
    hold on
    x = 1:14;
    CI_shuffle = prctile(reshape(shuffled_globalRemap_PPvector.population_vector(:,nplot),14,1000)',[1 99]);
    plot(x, CI_shuffle(2,:), 'k--', 'LineWidth', 1);
    plot(x, CI_shuffle(1,:), 'k--', 'LineWidth', 1);
    x2 = [x, fliplr(x)];
    inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
    h(2) = fill(x2, inBetween, 'k','FaceAlpha',0.2);

    h(1) = plot(x,PPvector.population_vector(:,nplot),'r')

    xticks([30 50 70 90 110 140]/10)
    xticklabels([30 50 70 90 110 140])
    xline(10,'r')

    ylabel('PV Correlation')
    xlabel('Position bin')

    title(sprintf('%s laps T%i VS %s laps T%i',...
        laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
    set(gca,"TickDir","out",'box', 'off','Color','none')
end
legend([h(1) h(2)],{'Original','Cell ID shuffled'},'Color','none')
box off


% Place cell stability and remapping
fig = figure
fig.Position = [500 70 1300 930];
if isempty(region)
    if ischar(nprobe)
        sgtitle(sprintf('%s %s PV correlation comparision for place cell stability and rate remapping probes %s',options.SUBJECT,options.SESSION,nprobe))
        fig.Name = sprintf('%s %s PV correlation comparision for place cell stability and rate remapping probes %s',options.SUBJECT,options.SESSION,nprobe);
    else
        sgtitle(sprintf('%s %s PV correlation comparision for place cell stability and rate remapping probe %s',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s PV correlation comparision for place cell stability and rate remapping probe %s',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere});
    end
else
    if ischar(nprobe)
        sgtitle(sprintf('%s %s PV correlation comparision for place cell stability and rate remapping probes %s',options.SUBJECT,options.SESSION,region,nprobe));
        fig.Name = sprintf('%s %s PV correlation comparision for place cell stability and rate remapping probes %s',options.SUBJECT,options.SESSION,region,nprobe);
    else
        sgtitle(sprintf('%s %s PV correlation comparision for %s cell stability and rate remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}));
        fig.Name = sprintf('%s %s PV correlation comparision for %s cell stability and rate remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});
    end
end
for nplot = 1:16
    if nplot == 1| nplot == 6 | nplot == 11| nplot == 16
        continue
    end
    subplot(4,4,nplot)
    hold on
    x = 1:14;
    CI_shuffle = prctile(reshape(shuffled_rateRemap_PPvector.population_vector(:,nplot),14,1000)',[1 99]);
    plot(x, CI_shuffle(2,:), 'k--', 'LineWidth', 1);
    plot(x, CI_shuffle(1,:), 'k--', 'LineWidth', 1);
    x2 = [x, fliplr(x)];
    inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
    h(2) = fill(x2, inBetween, 'k','FaceAlpha',0.2);

    h(1) = plot(x,PPvector.population_vector(:,nplot),'r')

    ylabel('PV Correlation')
    xlabel('Position bin')
    xticks([30 50 70 90 110 140]/10)
    xticklabels([30 50 70 90 110 140])
    xline(10,'r')

    title(sprintf('%s laps T%i VS %s laps T%i',...
        laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
    set(gca,"TickDir","out",'box', 'off','Color','none')
end
legend([h(1) h(2)],{'Original','Cell ID shuffled'},'Color','none')
box off

% figure
% plot(V1_clusters.probe(1).MUA_tvec,V1_clusters.probe(1).MUA_zscore);
% hold on
% plot(V1_clusters.probe(2).MUA_tvec,-V1_clusters.probe(2).MUA_zscore)