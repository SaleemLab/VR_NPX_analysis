function [normalised_raw_matrix PPvector shuffled_globalRemap_PPvector shuffled_rateRemap_PPvector] = plot_place_cell_map_correlation(place_fields_all,good_cell_index,Task_info,Behaviour,options)

% Function to plot spatial cell activity
% Plot 1 is each cell's spatial ratemap for each lap
% Plot 2 is each cell's spatial tuning stability by plotting place field
% calculated by using all laps, odd no laps and even no laps.
if isfield(options,'region')
    region = options.region;
else
    region = [];
end

if isfield(options,'probe_combined') && options.probe_combined == 1
    nprobe = 'combined';
elseif isfield(options,'probe_MEC') || isfield(options,'probe_V1')
    nprobe = options.probe_id+1;
    probe_hemisphere = 1;
    probe_hemisphere_text = options.text;
else
    nprobe = options.probe_id+1;
    probe_hemisphere = options.probe_hemisphere;
    probe_hemisphere_text = {'left','right'};
end



% probe_hemisphere_text = {'left','right'}
% place_fields_all = calculate_spatial_cells(clusters,Task_info,Behaviour,[0 140],10);

%Create linearised position matrix
index_track = [];
sorted_cells= [];
normalised_raw_matrix = [];

spatial_cell_index = unique([find(place_fields_all(1).peak_percentile>0.95 & place_fields_all(1).odd_even_stability>0.95)...
    find(place_fields_all(2).peak_percentile>0.95 & place_fields_all(2).odd_even_stability>0.95)]);
% spatial_cell_index = unique([find(place_fields_all(1).peak_percentile>0.99 & place_fields_all(1).odd_even_stability>0.99)...
%     find(place_fields_all(2).peak_percentile>0.99 & place_fields_all(2).odd_even_stability>0.99)]);
% spatial_cell_index = unique([find(place_fields_all(1).odd_even_stability>0.95)...
%     find(place_fields_all(2).odd_even_stability>0.95)]);
good_cell_index = intersect(spatial_cell_index,find(ismember(place_fields_all(1).cluster_id,good_cell_index)));
% good_cell_index = find(intersect(place_fields_all(1).cluster_id,good_cell_index));
if isempty(good_cell_index)
    disp('No spatial cells')
else
    for type = 1:3
        ratemap_matrix = [];
        for track_id = 1:2
            ratemap_matrix{track_id} = [place_fields_all(track_id).raw{good_cell_index}];
            ratemap_matrix{track_id} = reshape(ratemap_matrix{track_id},size(place_fields_all(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
        end
        %         ratemap_matrix = [place_fields_all(track_id).raw{place_fields(1).all_good_cells_LIBERAL}];
        %         ratemap_matrix = reshape(ratemap_matrix,size(place_fields_all(track_id).raw{1},1),[],length(place_fields(1).all_good_cells_LIBERAL));%laps X position bins X cells

        if type == 3 %even laps
            average_maps{1}= squeeze(mean(ratemap_matrix{1},1));
            average_maps{2} = squeeze(mean(ratemap_matrix{2},1));
        elseif type == 1 %odd laps
            average_maps{1}= squeeze(mean(ratemap_matrix{1}(1:2:end,:,:),1));
            average_maps{2} = squeeze(mean(ratemap_matrix{2}(1:2:end,:,:),1));

        elseif type == 2 %even laps
            average_maps{1}= squeeze(mean(ratemap_matrix{1}(2:2:end,:,:),1));
            average_maps{2} = squeeze(mean(ratemap_matrix{2}(2:2:end,:,:),1));
        end

        temp_map = normalize([average_maps{1}; average_maps{2}],'range');
        average_maps{1} = temp_map(1:size(average_maps{1},1),:);
        average_maps{2} = temp_map(1+size(average_maps{1},1):size(average_maps{1},1)+size(average_maps{2},1),:);

        average_maps{1}(isnan(average_maps{1}))=0;
        average_maps{2}(isnan(average_maps{2}))=0;

        for track_id = 1:2
            normalised_raw_matrix{type}{track_id} = average_maps{track_id}';

            [~,index_track{type}(track_id,:)] = max(normalised_raw_matrix{type}{track_id},[],2);
            %     unsorted_cells(track_id,:) =
            [~,sorted_cells{type}(track_id,:)] = sort(index_track{type}(track_id,:));
            %     ordered_matrix = normalised_raw_matrix(new_order,:);
        end
    end
    %
    % ncount= 1
    % figure
    % for ncell = 1:size( average_maps{1},2)
    %
    %     if ncount == 26
    %         figure
    %         ncount = 1
    %     end
    %     subplot(5,5,ncount)
    %
    %     plot(average_maps{1}(:,ncell))
    %     hold on
    %     plot(average_maps{2}(:,ncell))
    %
    %     ncount = 1 + ncount;
    % end

    laps_pairs = [1 2;1 2;1 2;1 2];
    track_pairs = [1 1;1 2;2 1;2 2];
    laps_type_text = {'odd','even','all'};
    %plot heat map position
    fig = figure
    fig.Position = [500 70 650 490];


    if ischar(nprobe)
        sgtitle(sprintf('%s %s %s cell map remapping probes %s',options.SUBJECT,options.SESSION,region,nprobe))
        fig.Name = sprintf('%s %s %s cell map remapping probes %s',options.SUBJECT,options.SESSION,region,nprobe);
    else
        sgtitle(sprintf('%s %s %s cell map remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}))
        fig.Name = sprintf('%s %s %s cell map remapping probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});
    end


    for nplot = 1:4
        subplot(2,2,nplot)
        ordered_matrix =  normalised_raw_matrix{laps_pairs(nplot,1)}{track_pairs(nplot,1)}(sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:),:);
        imagesc(ordered_matrix);
        xlabel('Position');
        ylabel('Cell ID');
        yt=sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:);
        set(gca,'ytick',1:length(sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:)));
        set(gca,'yticklabel',yt);
        axis xy;
        h = colorbar;
        colormap(flip(gray))
        set(get(h,'label'),'string','Normalised firing rate across tracks');
        clim([0.3 1])

        xticks([30 50 70 90 110 140]/(place_fields_all(1).x_bin_centres(2)-place_fields_all(1).x_bin_centres(1)))
        xticklabels([30 50 70 90 110 140])
        xline(100/(place_fields_all(1).x_bin_centres(2)-place_fields_all(1).x_bin_centres(1)),'r')

        title(sprintf('%s laps T%i sorted by %s laps T%i',...
            laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end


    %         title(sprintf('%s %s Even vs Odd laps place cell maps probe %i',options.SUBJECT,options.SESSION,nprobe))


    laps_pairs = [1 2;2 1;1 2;2 1];
    track_pairs = [1 1;1 1;2 2;2 2];
    laps_type_text = {'odd','even','all'};
    %plot heat map position
    fig = figure
    fig.Position = [500 70 650 490];


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

    for nplot = 1:4
        subplot(2,2,nplot)
        ordered_matrix =  normalize(normalised_raw_matrix{laps_pairs(nplot,1)}{track_pairs(nplot,1)}(sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:),:),2,'range');
        imagesc(ordered_matrix);
        xlabel('Position');
        ylabel('Cell ID');
        yt=sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:);
        set(gca,'ytick',1:length(sorted_cells{laps_pairs(nplot,2)}(track_pairs(nplot,2),:)));
        set(gca,'yticklabel',yt);
        axis xy;
        h = colorbar;
        colormap(flip(gray))
        set(get(h,'label'),'string','Normalised firing rate within track');
        clim([0.3 1])

        xticks([30 50 70 90 110 140]/(place_fields_all(1).x_bin_centres(2)-place_fields_all(1).x_bin_centres(1)))
        xticklabels([30 50 70 90 110 140])
        xline(100/(place_fields_all(1).x_bin_centres(2)-place_fields_all(1).x_bin_centres(1)),'r')

        title(sprintf('%s laps T%i sorted by %s laps T%i',...
            laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end

    %         title(sprintf('%s %s Even vs Odd laps place cell maps probe %i',options.SUBJECT,options.SESSION,nprobe))


    % Calculate cell population vector for each track comparison by correlating each position bin between tracks
    shuffled_globalRemap_PPvector.population_vector = [];
    shuffled_globalRemap_PPvector.pval = [];
    shuffled_rateRemap_PPvector.population_vector =[];
    shuffled_rateRemap_PPvector.pval = [];
    laps_pairs = [1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2;1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2];
    track_pairs = [1 1;1 1;1 2;1 2;1 1;1 1;1 2;1 2;2 1;2 1;2 2;2 2;2 1;2 1;2 2;2 2];
    PPvector = [];

    tic
    % Original raw rate map population vector correlation
    for i = 1 : length(laps_pairs)
        for j = 1 : size( normalised_raw_matrix{1}{track_id},2) % for each position bin
            % Correlation to global remapping shuffle
            [Grho,Gpval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), normalised_raw_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j)); %corr between position bins across cells
            PPvector.population_vector(j,i) = Grho;
            PPvector.pval(j,i) = Gpval;
        end
    end
    Gpv= cell(1,1000);
    Gpval= cell(1,1000);
    Rpv= cell(1,1000);
    Rpval= cell(1,1000);

    parfor nshuffle = 1 : 1000
        cellID_shuffled_matrix = []; rate_remap_matrix= [];

        % Create shuffled matrix for global remapping: shuffle cell ID & one for rate remapping: multiply each ratemap with a random num within 0.1 to 1.5
        for track_id = 1 : length(place_fields_all)
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

        %     curr_size = size(shuffled_rateRemap_PPvector.population_vector,1);
        %
        %     for i = 1 : length(laps_pairs)
        %         for j = 1 : size( normalised_raw_matrix{1}{track_id},2) % for each position bin
        %             % Correlation to global remapping shuffle
        %             [Grho,Gpval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), cellID_shuffled_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j)); %corr between position bins across cells
        %             shuffled_globalRemap_PPvector.population_vector(curr_size+j,i) = Grho; %each column is a comparison, each row a position bin
        %             shuffled_globalRemap_PPvector.pval(curr_size+j,i) = Gpval;
        %             % Correlation to rate remapping shuffle
        %             [Rrho,Rpval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), rate_remap_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j));
        %             shuffled_rateRemap_PPvector.population_vector(curr_size+j,i) = Rrho; %each column is a comparison, each row a position bin
        %             shuffled_rateRemap_PPvector.pval(curr_size+j,i) = Rpval;
        %         end
        %     end

        for i = 1 : length(laps_pairs)
            for j = 1 : size( normalised_raw_matrix{1}{track_id},2) % for each position bin
                % Correlation to global remapping shuffle
                [Grho,pval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), cellID_shuffled_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j)); %corr between position bins across cells
                Gpv{nshuffle}(j,i) = Grho; %each column is a comparison, each row a position bin
                Gpval{nshuffle}(j,i)  = pval;
                % Correlation to rate remapping shuffle
                [Rrho,pval] = corr(normalised_raw_matrix{laps_pairs(i,1)}{track_pairs(i,1)}(:,j), rate_remap_matrix{laps_pairs(i,2)}{track_pairs(i,2)}(:,j));
                Rpv{nshuffle}(j,i) = Rrho; %each column is a comparison, each row a position bin
                Rpval{nshuffle}(j,i)  = pval;
            end
        end
    end

    shuffled_globalRemap_PPvector.population_vector = [Gpv{:}];
    shuffled_globalRemap_PPvector.population_vector = reshape(shuffled_globalRemap_PPvector.population_vector,length(place_fields_all(1).x_bin_centres),length(laps_pairs),1000);

    shuffled_globalRemap_PPvector.pval = [Gpval{:}];
    shuffled_globalRemap_PPvector.pval = reshape(shuffled_globalRemap_PPvector.pval,length(place_fields_all(1).x_bin_centres),length(laps_pairs),1000);

    shuffled_rateRemap_PPvector.population_vector = [Rpv{:}];
    shuffled_rateRemap_PPvector.population_vector = reshape(shuffled_rateRemap_PPvector.population_vector,length(place_fields_all(1).x_bin_centres),length(laps_pairs),1000);

    shuffled_rateRemap_PPvector.pval = [Rpval{:}];
    shuffled_rateRemap_PPvector.pval = reshape(shuffled_rateRemap_PPvector.pval,length(place_fields_all(1).x_bin_centres),length(laps_pairs),1000);

    toc

    %         % Run Kruskal-Wallis
    %         all_sig_diff_idx = []; sig_diff_idx= [];
    %         [pv3,~,stats3] = kruskalwallis(combined,[],'off');
    %         if pv3 < 0.05
    %             [all_c,~,~,~] = multcompare(stats3,'ctype','dunn-sidak','Display','off'); % if anova pval is < 0.05, run multiple comparisons
    %             protocol(p).all_multiple_comparisons = all_c;
    %             all_sig_diff_idx = find(all_c(:,6)<0.05);
    %         end

    % fig = figure
    % fig.Position = [500 70 1300 930];
    % if isempty(region)
    %     if ischar(nprobe)
    %         sgtitle(sprintf('%s %s place cell PV correlation cumulative frequency comparision probes %s.pdf',options.SUBJECT,options.SESSION,nprobe))
    %         fig.Name = sprintf('%s %s place cell PV correlation cumulative frequency comparision probes %s.pdf',options.SUBJECT,options.SESSION,nprobe);
    %     else
    %         sgtitle(sprintf('%s %s place cell PV correlation cumulative frequency comparision probe %s.fig',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere}))
    %         fig.Name = sprintf('%s %s place cell PV correlation cumulative frequency comparision probe %s.fig',options.SUBJECT,options.SESSION,probe_hemisphere_text{probe_hemisphere});
    %     end
    % else
    %     if ischar(nprobe)
    %         sgtitle(sprintf('%s %s %s cell PV correlation cumulative frequency comparision probes %s',options.SUBJECT,options.SESSION,region,nprobe));
    %         fig.Name = sprintf('%s %s %s cell PV correlation cumulative frequency comparision probes %s',options.SUBJECT,options.SESSION,region,nprobe);
    %     else
    %         sgtitle(sprintf('%s %s %s cell PV correlation cumulative frequency comparision probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere}))
    %         fig.Name = sprintf('%s %s %s cell PV correlation cumulative frequency comparision probe %s',options.SUBJECT,options.SESSION,region,probe_hemisphere_text{probe_hemisphere});
    %     end
    % end
    %
    % for nplot = 1:16
    %     if nplot == 1| nplot == 6 | nplot == 11| nplot == 16
    %         continue
    %     end
    %     subplot(4,4,nplot)
    %     h(1) = cdfplot(PPvector.population_vector(:,nplot));
    %     hold on
    %     h(2) = cdfplot(shuffled_globalRemap_PPvector.population_vector(:,nplot));
    %     set( h(1), 'LineStyle', '-', 'Color', 'r');
    %     set( h(2), 'LineStyle', '-.', 'Color', 'k');
    %     xlabel('PV correlation')
    %     ylabel('Cumulative frequency')
    %     legend('Original','Cell ID shuffled','Color','none')
    %     title(sprintf('%s laps T%i sorted by %s laps T%i',...
    %         laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
    %     set(gca,"TickDir","out",'box', 'off','Color','none')
    % end


    % Place cell stability and remapping
    laps_pairs = [1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2;1 1; 1 2;1 1;1 2;2 1;2 2;2 1;2 2];
    track_pairs = [1 1;1 1;1 2;1 2;1 1;1 1;1 2;1 2;2 1;2 1;2 2;2 2;2 1;2 1;2 2;2 2];
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
        x = 1:length(place_fields_all(1).x_bin_centres);
        CI_shuffle = prctile(squeeze(shuffled_globalRemap_PPvector.population_vector(:,nplot,:))',[1 99]);
        %     CI_shuffle = prctile(reshape(shuffled_globalRemap_PPvector.population_vector(:,nplot),length(place_fields_all(1).x_bin_centres),1000)',[1 99]);
        plot(x, CI_shuffle(2,:), 'k--', 'LineWidth', 1);
        plot(x, CI_shuffle(1,:), 'k--', 'LineWidth', 1);
        x2 = [x, fliplr(x)];
        inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
        h(2) = fill(x2, inBetween, 'k','FaceAlpha',0.2);

        h(1) = plot(x,PPvector.population_vector(:,nplot),'r')

        xticks([30 50 70 90 110 140]/mean(diff(place_fields_all(1).x_bin_centres)))
        xticklabels([30 50 70 90 110 140])
        xline(100/mean(diff(place_fields_all(1).x_bin_centres)),'r')

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
        x = 1:length(place_fields_all(1).x_bin_centres);
        %     CI_shuffle = prctile(reshape(shuffled_rateRemap_PPvector.population_vector(:,nplot),length(place_fields_all(1).x_bin_centres),1000)',[1 99]);
        CI_shuffle = prctile(squeeze(shuffled_rateRemap_PPvector.population_vector(:,nplot,:))',[1 99]);
        plot(x, CI_shuffle(2,:), 'k--', 'LineWidth', 1);
        plot(x, CI_shuffle(1,:), 'k--', 'LineWidth', 1);
        x2 = [x, fliplr(x)];
        inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
        h(2) = fill(x2, inBetween, 'k','FaceAlpha',0.2);

        h(1) = plot(x,PPvector.population_vector(:,nplot),'r')

        ylabel('PV Correlation')
        xlabel('Position bin')
        xticks([30 50 70 90 110 140]/mean(diff(place_fields_all(1).x_bin_centres)))
        xticklabels([30 50 70 90 110 140])
        xline(100/mean(diff(place_fields_all(1).x_bin_centres)),'r')

        title(sprintf('%s laps T%i VS %s laps T%i',...
            laps_type_text{laps_pairs(nplot,1)},track_pairs(nplot,1),laps_type_text{laps_pairs(nplot,2)},track_pairs(nplot,2)))
        set(gca,"TickDir","out",'box', 'off','Color','none')
    end
    legend([h(1) h(2)],{'Original','Cell ID shuffled'},'Color','none')
    box off
end


%% lap by lap PV correlation (odd/even laps representation)

%
% clear single_lap_PV_correlation
%
% for type = 1:2
%     for track_id = 1:2
%         ratemap_matrix = [place_fields_all(track_id).raw{good_cell_index}];
%         ratemap_matrix = reshape(ratemap_matrix,size(place_fields_all(track_id).raw{1},1),[],length(good_cell_index));%laps X position bins X cells
%
%
%         if type == 2
%             % even lap template
%             average_maps = normalize(squeeze(mean(ratemap_matrix(2:2:end,:,:),1)),'range');
%             average_maps(isnan(average_maps))=0;
%             selected_laps = 1:2:sum(Task_info.track_ID_all == track_id);
%
%         else
%             % odd lap template
%             average_maps = normalize(squeeze(mean(ratemap_matrix(1:2:end,:,:),1)),'range');
%             average_maps(isnan(average_maps))=0;
%             selected_laps = 2:2:sum(Task_info.track_ID_all == track_id);
%         end
%         count = 1;
%         for nlap = selected_laps
%             lap_PV = normalize(squeeze(ratemap_matrix(nlap,:,:)),1,'range');
%             for j = 1:length(place_fields_all(1).x_bin_centres)
%                 [Rrho,Rpval] = corr(lap_PV(j,:)', average_maps(j,:)');
%                 single_lap_PV_correlation(track_id).population_vector{type}(count,j) =Rrho;
%                 single_lap_PV_correlation(track_id).pval{type}(count,j) =Rpval;
%             end
%             count = count+1;
%         end
%         % odd lap template
%     end
% end
