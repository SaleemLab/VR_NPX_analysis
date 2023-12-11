%% Main Code for place cell analysis and Bayesian decoding
% Go to main_NPX_data_preprocessing for preprocessing before running this
% script

addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
%     ROOTPATH = 'X:\ibn-vision';
    ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive
%     ROOTPATH = '/research';
end


%% Load Spike train data for place cell and ripple analysis

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD

        raw_LFP = [];
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            column = 1;
            session_info(n).probe(nprobe).task_type = stimulus_name{n};
            options = session_info(n).probe(nprobe);
            options.importMode = 'LF';
            options.probe_no = nprobe; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            if nprobe ~= 1
                session_info.probe(1).importMode = 'KS';

                [~, imecMeta, ~, ~] = extract_NPX_channel_config(session_info.probe(1),column);
                [raw_LFP{nprobe} tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column,'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp);
            else
                [raw_LFP{nprobe} tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column);
            end

            if ~isempty(tvec)
                LFP_tvec = tvec;
            else
                LFP_tvec = [];
            end
        end

        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            options = session_info(n).probe(nprobe);
            options.importMode = 'KS';
            options.gFileNum = gFileNum(n);
            % Load all spike data sorted according to the channel position

            [SUA.probe{nprobe} chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'group','by channel','cell_exporer','on');
            % for n = 1:length(SUA)
            %     if ~isempty(SUA(n).spike_times)
            %         spike_count(n) = length(SUA(n).spike_times);
            %     else
            %         spike_count(n) = 0;
            %     end
            %
            % end
            %
            % hold on
            % plot(spike_count/max(spike_count),chan_config.Ks_ycoord','Color','g')
            % ylim([0 4000])


            % L4 spike data (roughly 100 micron or 10 channels)
            [L4_clusters.probe(nprobe) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{nprobe}.L4_channel-5 best_channels{nprobe}.L4_channel+5],'group','by region','cell_exporer','on');
            % [L4_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L4_channel-10 best_channels.L4_channel],'group','by region');
            % L5 spike data (roughly 200 micron or 20 channels))
            [L5_clusters.probe(nprobe) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{nprobe}.L5_channel-10 best_channels{nprobe}.L5_channel+10],'group','by region','cell_exporer','on');

            % all V1 spike data (roughly 300 micron or 30 channels))
            [superficial_clusters.probe(nprobe) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{nprobe}.first_in_brain_channel-30 best_channels{nprobe}.first_in_brain_channel],'group','by region','cell_exporer','on');
            % [superficial_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.first_in_brain_channel-15 best_channels.first_in_brain_channel],'group','by region');

            % all V1 spike data
            [V1_clusters.probe(nprobe) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{nprobe}.L5_channel-10 best_channels{nprobe}.first_in_brain_channel],'group','by region','cell_exporer','on');
            V1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters.probe(nprobe),[]);

            % CA1 spike data (roughly 200 micron or 20 channels))
            [CA1_clusters.probe(nprobe) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{nprobe}.CA1_channel-10 best_channels{nprobe}.CA1_channel+10],'group','by region','cell_exporer','on');
            CA1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,CA1_clusters.probe(nprobe),[]);
            
            % 'Whole' HPC from 100 micron above CA1 cell layer to 1000 micron
            % below
            [HPC_clusters.probe(nprobe) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{nprobe}.CA1_channel-100 best_channels{nprobe}.CA1_channel+10],'group','by region','cell_exporer','on');
            HPC_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,HPC_clusters.probe(nprobe),[]);
        end

        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load extracted_laps.mat
           

        for nprobe = 1:length(session_info(n).probe)
            CA1_place_fields_even.probe(nprobe)  = calculate_place_fields_masa_NPX(x_bins_width,position,CA1_clusters.probe(nprobe),'even laps');
            CA1_place_fields_odd.probe(nprobe) = calculate_place_fields_masa_NPX(x_bins_width,position,CA1_clusters.probe(nprobe),'odd laps');
            %             CA1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,CA1_clusters.probe(nprobe),[]);
            %             HPC_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,HPC_clusters.probe(nprobe),[]);
            %             V1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters.probe(nprobe),[]);
        end

        if contains(stimulus_name{n},'RUN')
            save extracted_HPC_place_fields HPC_place_fields
            save extracted_CA1_place_fields CA1_place_fields
            save extracted_CA1_place_fields_even CA1_place_fields_even
            save extracted_CA1_place_fields_odd CA1_place_fields_odd
            save extracted_V1_place_fields V1_place_fields
        end

        save(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'HPC_clusters')
        save(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'CA1_clusters')
        save(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters')





        for nprobe = 1:length(session_info(n).probe)
            clusters = HPC_clusters;
                        clusters = CA1_clusters;
            clusters = V1_clusters;

            x_bins_width = 10;
            
            place_fields_even = CA1_place_fields_even.probe(nprobe);
            place_fields_odd = CA1_place_fields_odd.probe(nprobe);
            place_fields_all = CA1_place_fields.probe(nprobe);
            
            place_fields_lap = [];
            for track_id = 1:2
                for nlap = 1:length(lap_times(track_id).start)
                    field_this_lap = get_lap_place_fields_masa(x_bins_width,position,place_fields_all.probe(nprobe),...
                        clusters.probe(nprobe),track_id,lap_times(track_id).start(nlap),lap_times(track_id).end(nlap));

                    place_fields_lap{track_id}{nlap} = field_this_lap.raw;
                end
            end



            for track_id = 1:2
                figure
                for ncell = 1:length(place_fields.probe(nprobe).all_cells)
                    this_cell_place_field = [];
                    for nlap = 1:2:length(lap_times(track_id).start)
                        this_cell_place_field(nlap,:) = place_fields_lap{track_id}{nlap}{ncell};
                    end

                    subplot(ceil(sqrt(length(place_fields.probe(nprobe).all_cells))),ceil(sqrt(length(place_fields.probe(nprobe).all_cells))),ncell)
                    imagesc(this_cell_place_field./max(this_cell_place_field')')
%                     if sum(ncell ==
%                         title(ncell,'Color','r')
%                     else
                        title(ncell,'Color','k')
%                     end
                end
                sgtitle(sprintf('Track %i',track_id))
            end

        end

        %         x_bins_width = 10;
        %
        %         place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,CA1_clusters.probe(nprobe),'even laps')
        %         place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,CA1_clusters.probe(nprobe),'odd laps')

        
        for ncell = place_fields_all.good_place_cells_LIBERAL
            figure
            subplot(2,2,1)
            plot(1:x_bins_width:140,place_fields_all.track(1).raw{ncell}/max(place_fields_all.track(1).raw{ncell})); hold on; plot(1:x_bins_width:140,place_fields_all.track(2).raw{ncell}/max(place_fields_all.track(2).raw{ncell}))
            ylabel('Normalised FR')
            xlabel('Position bins')
            title('all trials')
            subplot(2,2,2)
            plot(1:x_bins_width:140,place_fields_odd.track(1).raw{ncell}/max(place_fields_odd.track(1).raw{ncell})); hold on; plot(1:x_bins_width:140,place_fields_odd.track(2).raw{ncell}/max(place_fields_odd.track(2).raw{ncell}))
            title('odd trials')
            subplot(2,2,3)
            plot(1:x_bins_width:140,place_fields_even.track(1).raw{ncell}/max(place_fields_even.track(1).raw{ncell})); hold on; plot(1:x_bins_width:140,place_fields_even.track(2).raw{ncell}/max(place_fields_even.track(2).raw{ncell}))
            title('even trials')
        end


        for test = 1:3
            figure;

            if test == 1
                place_fields = place_fields_all;
                sgtitle('Whole')
            elseif test == 2
                place_fields = place_fields_odd;
                sgtitle('Odd laps')
            elseif test == 3
                place_fields = place_fields_even;
                sgtitle('Even laps')
            end

            c=1;
            for kk=1:length(place_fields_all.track)
                for j=1:length(place_fields_all.track)
                    y_vector=[];
                    matrix=[];
                    normalized_matrix=[];
                    for ii=1:length(place_fields_all.track(j).sorted_good_cells_LIBERAL)
                        %plot sorted
                        %             matrix=[];
                        %             normalized_matrix=[];
                        matrix(ii,:)=place_fields.track(kk).raw{place_fields_all.track(j).sorted_good_cells_LIBERAL(ii)};
                        normalized_matrix(ii,:)=(matrix(ii,:)-min(matrix(ii,:)))/(max(matrix(ii,:))-min(matrix(ii,:)));
                        subplot(length(place_fields.track),length(place_fields.track),c)
                        plfield_row= normalized_matrix(ii,:)+(1.5*ii-1);
                        plot(1:length(plfield_row),plfield_row,'k'); hold on;
                        xx = [1:length(plfield_row), fliplr(1:length(plfield_row))];
                        inBetween = [(1.5*ii-1)*ones(size(plfield_row)), fliplr(plfield_row)];
                        fill(xx, inBetween,[139,0,0]/255);
                        y_vector= [y_vector, 1.5*ii-1];
                    end
                    xlim([0 size(normalized_matrix,2)+2]);
                    ylim([0 max(y_vector)+1.2]);
                    yt=place_fields_all.track(j).sorted_good_cells_LIBERAL;
                    set(gca,'ytick',y_vector);
                    set(gca,'yticklabel',yt);
                    ylabel('Unit ID');
                    xlabel('sorted linearized position (bins)');
                    c=c+1;
                    title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(kk)]}]);
                end
            end
        end



        %Create linearised position matrix
        index_track = [];
        sorted_cells= [];

        for type = 1:3

            if type == 1
                place_fields = place_fields_all;
                sgtitle('Whole')
            elseif type == 2
                place_fields = place_fields_odd;
                sgtitle('Odd laps')
            elseif type == 3
                place_fields = place_fields_even;
                sgtitle('Even laps')
            end

            for track_id = 1:2
                raw_matrix = cat(1,place_fields.track(track_id).raw{place_fields_all.good_place_cells_LIBERAL});
                normalised_raw_matrix{type}{track_id} = bsxfun(@rdivide, raw_matrix, (place_fields.track(track_id).raw_peak(place_fields_all.good_place_cells_LIBERAL))');
                [~,index_track{type}(track_id,:)] = max(normalised_raw_matrix{type}{track_id},[],2);
                %     unsorted_cells(track_id,:) =
                [~,sorted_cells{type}(track_id,:)] = sort(index_track{type}(track_id,:));
                %     ordered_matrix = normalised_raw_matrix(new_order,:);
            end
        end
        
        laps_pairs = [2 2; 2 3;2 2;2 3;3 2;3 3;3 2;3 3;2 2; 2 3;2 2;2 3;3 2;3 3;3 2;3 3];
        track_pairs = [1 1;1 1;2 1;2 1;1 1;1 1;2 1;2 1;2 1;2 1;2 2;2 2;2 1;2 1;2 2;2 2];
        laps_type_text = {'all','odd','even'};
        %plot heat map position
        fig = figure
        fig.Position = [500 70 1300 930];
        for n = 1:16
            subplot(4,4,n)
            ordered_matrix =  normalised_raw_matrix{laps_pairs(n,1)}{track_pairs(n,1)}(sorted_cells{laps_pairs(n,2)}(track_pairs(n,2),:),:);
            imagesc(ordered_matrix);
            colorbar;
            xlabel('Position');
            ylabel('Cell ID');
            axis xy;
            h = colorbar;
            set(get(h,'label'),'string','Normalised firing rate');
            title(sprintf('%s laps T%i sorted by %s laps T%i',...
                laps_type_text{laps_pairs(n,1)},track_pairs(n,1),laps_type_text{laps_pairs(n,2)},track_pairs(n,2)))
        end
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
%         title(sprintf('%s %s Even vs Odd laps place cell maps probe %i',options.SUBJECT,options.SESSION,nprobe))
        filename = sprintf('%s %s Even vs Odd laps place cell maps probe %i.pdf',options.SUBJECT,options.SESSION,nprobe)
        saveas(gcf,filename)
        filename = sprintf('%s %s Even vs Odd laps place cell maps probe %i.fig',options.SUBJECT,options.SESSION,nprobe)
        saveas(gcf,filename)

        % bayesian decoding
        time_on_T2 = position.t(~isnan(position.linear(2).linear));
        bayesian_spike_count = create_spike_count_masa(place_fields,clusters,time_on_T2(2),time_on_T2(end),[])
        estimated_position = bayesian_decoding(place_fields,bayesian_spike_count,[])

        % Ripple detection
        channel_to_use = find(sorted_config.Channel == best_channels{nprobe}.CA1_channel);
        [replay,reactivations] = detect_candidate_events_masa(tvec,raw_LFP(channel_to_use,:),...
            CA1_clusters.MUA_zscore,[CA1_clusters.spike_id CA1_clusters.spike_times],peripherals,zscore_min,zscore_max,options)
        save extracted_candidate_events replay reactivations
    end
end

 % bayesian decoding
time_on_T2 = position.t(~isnan(position.linear(2).linear));
bayesian_spike_count = create_spike_count_masa(place_fields,clusters,time_on_T2(2),time_on_T2(end),[])
 estimated_position = bayesian_decoding(place_fields,bayesian_spike_count,[])

figure;
subplot(1,2,1)
plot(position.t,estimated_position(2).discrete_position,'k')
hold on
scatter(estimated_position(1).run_time_centered,estimated_position(1).peak_position,5,'r','filled')

subplot(1,2,2)
plot(position.t,estimated_position(2).discrete_position,'k')
hold on
scatter(estimated_position(2).run_time_centered,estimated_position(2).peak_position,'b','filled')
imagesc(estimated_position(2).run_time_centered,[1:12:120],estimated_position(2).run);



place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'even laps')
place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'odd laps')


V1_track_1_cells =  cell_id(find(place_fields.track(1).mean_rate_track > 1 & place_fields.track(2).mean_rate_track < 0.5));
V1_track_2_cells =  cell_id(find(place_fields.track(1).mean_rate_track < 0.5 & place_fields.track(2).mean_rate_track > 1));
% cell_selected = find((place_fields.track(1).mean_rate_track  - place_fields.track(2).mean_rate_track)>5)

save track_cells V1_track_1_cells V1_track_2_cells CA1_track_1_cells CA1_track_2_cells

V1_track_1_cells =  find(place_fields.track(1).mean_rate_track > 1 & place_fields.track(2).mean_rate_track < 0.5);
figure
subplot(2,2,1)

bar(place_fields.track(1).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate')
ylim([0 20])
title('Mean Firing Rate on Track 1 (when moving)')

subplot(2,2,2)
bar(place_fields.track(2).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate')
ylim([0 20])
title('Mean Firing Rate on Track 2 (when moving)')

subplot(2,2,3)
bar(place_fields.track(1).mean_rate_track(V1_track_1_cells) - place_fields.track(2).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate Difference')
ylim([0 20])
title('Firing rate difference (when moving)')
sgtitle('Track 1 ''preferred'' V1 cells')


% CA1
clusters.spike_times = CA1_spike_times(:,2);
clusters.spike_id = CA1_spike_times(:,1);
cell_id = unique(CA1_spike_times(:,1));
for cell = 1:length(cell_id)
    clusters.spike_id(find(CA1_spike_times(:,1) == cell_id(cell))) = cell;
end

 place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[])
 save extracted_place_fields_CA1 place_fields
 place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'even laps')
  place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'odd laps')

% track_1_cells =  cell_id(find(place_fields.track(1).mean_rate_track > 3 & place_fields.track(2).mean_rate_track < 2);
CA1_track_1_cells = find((place_fields.track(1).mean_rate_track  - place_fields.track(2).mean_rate_track)>0.5);
CA1_track_2_cells = find((place_fields.track(1).mean_rate_track  - place_fields.track(2).mean_rate_track)<-0.5);
pyramidal_cells = find(place_fields.mean_rate < 2);
CA1_track_1_cells = cell_id(intersect(CA1_track_1_cells,pyramidal_cells));
CA1_track_2_cells = intersect(CA1_track_2_cells,pyramidal_cells);


figure
subplot(2,2,1)
bar(place_fields.track(1).mean_rate_track(CA1_track_1_cells))
ylim([0 20])
title('Mean Firing Rate on Track 1 (when moving)')

subplot(2,2,2)
bar(place_fields.track(2).mean_rate_track(CA1_track_1_cells))
ylim([0 20])
title('Mean Firing Rate on Track 2 (when moving)')

subplot(2,2,3)
bar(place_fields.track(1).mean_rate_track(CA1_track_1_cells) - place_fields.track(2).mean_rate_track(CA1_track_1_cells))
ylim([0 20])
title('Firing rate difference (when moving)')
sgtitle('Track 1 ''preferred'' CA1 cells')



 place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'even laps')
  place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'odd laps')

for test = 1:3
        figure;
       
        if test == 1
            place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[]);
            sgtitle('Whole')
        elseif test == 2

            place_fields = place_fields_odd;
            sgtitle('Odd')
        elseif test == 3
            place_fields = place_fields_even;
            sgtitle('Even')
        end

        c=1;
    for kk=1:length(place_fields.track)
        for j=1:length(place_fields.track)
            y_vector=[];
            matrix=[];
            normalized_matrix=[];
            for ii=1:length(place_fields.track(j).sorted_good_cells)
                %plot sorted
                %             matrix=[];
                %             normalized_matrix=[];
                matrix(ii,:)=place_fields.track(kk).raw{place_fields.track(j).sorted_good_cells(ii)};
                normalized_matrix(ii,:)=(matrix(ii,:)-min(matrix(ii,:)))/(max(matrix(ii,:))-min(matrix(ii,:)));
                subplot(length(place_fields.track),length(place_fields.track),c)
                plfield_row= normalized_matrix(ii,:)+(1.5*ii-1);
                plot(1:length(plfield_row),plfield_row,'k'); hold on;
                xx = [1:length(plfield_row), fliplr(1:length(plfield_row))];
                inBetween = [(1.5*ii-1)*ones(size(plfield_row)), fliplr(plfield_row)];
                fill(xx, inBetween,[139,0,0]/255);
                y_vector= [y_vector, 1.5*ii-1];
            end
            xlim([0 size(normalized_matrix,2)+2]);
            ylim([0 max(y_vector)+1.2]);
            yt=place_fields.track(j).sorted_good_cells;
            set(gca,'ytick',y_vector);
            set(gca,'yticklabel',yt);
            ylabel('Unit ID');
            xlabel('sorted linearized position (bins)');
            c=c+1;
            title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(kk)]}]);
        end
    end
end


%Create linearised position matrix
index_track = [];
 sorted_cells= [];

 for type = 1:3
     if type == 1
         place_fields = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[]);
     elseif type == 2
         place_fields = place_fields_odd;
     elseif type == 3
            place_fields = place_fields_even;
     end

     for track_id = 1:2
         raw_matrix = cat(1,place_fields.track(track_id).raw{:});
         normalised_raw_matrix{type}{track_id} = bsxfun(@rdivide, raw_matrix, (place_fields.track(track_id).raw_peak)');
         [~,index_track{type}(track_id,:)] = max(normalised_raw_matrix{type}{track_id},[],2);
         %     unsorted_cells(track_id,:) =
         [~,sorted_cells{type}(track_id,:)] = sort(index_track{type}(track_id,:));
         %     ordered_matrix = normalised_raw_matrix(new_order,:);
     end
 end

 %plot heat map position
 subplot(2,2,1)
 ordered_matrix =  normalised_raw_matrix{2}{1}(sorted_cells{3}(1,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

  subplot(2,2,2)
 ordered_matrix =  normalised_raw_matrix{2}{2}(sorted_cells{3}(2,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

  subplot(2,2,3)
 ordered_matrix =  normalised_raw_matrix{2}{1}(sorted_cells{2}(1,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

  subplot(2,2,4)
 ordered_matrix =  normalised_raw_matrix{3}{1}(sorted_cells{3}(1,:),:);
 imagesc(ordered_matrix);
 colorbar;
 xlabel('Position');
 ylabel('Cell ID');
 axis xy;
 h = colorbar;
 set(get(h,'label'),'string','Normalised firing rate');

 %% Biasing track selective neruons in CA1 by V1
events = ripples;
events = reactivations;
CA1_track1_SWR_spikes = [];
V1_track1_SWR_spikes = [];
V1_track2_SWR_spikes = [];
active_cells_SWR = [];


 for event = 1:length(events.onset)
     onset = events.onset(event)-0.2;
     offset = events.offset(event);
     spike_index = find(V1_spike_times(:,2)>onset & V1_spike_times(:,2)<offset);
     for cell = 1:length(V1_track_1_cells)
         V1_track1_SWR_spikes(cell,event) = sum(find(V1_spike_times(spike_index,1) == V1_track_1_cells(cell)));
     end
     
     for cell = 1:length(V1_track_2_cells)
         V1_track2_SWR_spikes(cell,event) = sum(find(V1_spike_times(spike_index,1) == V1_track_2_cells(cell)));
     end

     onset = events.onset(event);
     offset = events.offset(event);
     spike_index = find(CA1_spike_times(:,2)>onset & CA1_spike_times(:,2)<offset);
     for cell = 1:length(CA1_track_1_cells)
         CA1_track1_SWR_spikes(cell,event) = sum(find(CA1_spike_times(spike_index,1) == CA1_track_1_cells(cell)));
     end

     for cell = 1:length(CA1_track_2_cells)
         CA1_track2_SWR_spikes(cell,event) = sum(find(CA1_spike_times(spike_index,1) == CA1_track_2_cells(cell)));
     end

%      active_cells_SWR(1,event) =  sum(V1_track1_SWR_spikes(:,event) > 0);
%      active_cells_SWR(2,event) =  sum(V1_track2_SWR_spikes(:,event) > 0) ;
%      active_cells_SWR(3,event) =  sum(CA1_track1_SWR_spikes(:,event) > 0);

          active_cells_SWR(1,event) =  sum(V1_track1_SWR_spikes(:,event) > 0) >=3;
     active_cells_SWR(2,event) =   sum(V1_track2_SWR_spikes(:,event) > 0) >=1;
     active_cells_SWR(3,event) =   sum(CA1_track1_SWR_spikes(:,event) > 0) >=1;
     active_cells_SWR(4,event) =   sum(CA1_track2_SWR_spikes(:,event) >= 0) >=0;
 end

find(reactivations.ripple_peak>3)


figure
% p1 = plot(events.onset(find(reactivations.ripple_peak>3)),cumsum(active_cells_SWR(1,find(reactivations.ripple_peak>3))...
%     /max(cumsum(active_cells_SWR(1,find(reactivations.ripple_peak>3))))),'r')
% subplot(2,2,1)
p1 = plot(events.onset,cumsum(active_cells_SWR(1,:))/max(cumsum(active_cells_SWR(1,:))),'r')
hold on
% p2 = plot(events.onset,cumsum(active_cells_SWR(4,:))/max(cumsum(active_cells_SWR(4,:))),'k')
p2 = plot(events.onset,cumsum(active_cells_SWR(3,:))/max(cumsum(active_cells_SWR(3,:))),'b')
p3 = plot(events.onset,cumsum(active_cells_SWR(4,:))/max(cumsum(active_cells_SWR(4,:))),'k')
plot(position.t,position.linear(1).linear/100,'k')
% plot(position.t,-position.linear(2).linear/100,'k')
% s1 = scatter(MousePos.stimuli_onset(MousePos.stimuli_track==1),0.5*ones(1,sum(MousePos.stimuli_track==1)),'r');
% hold on
% s2 = scatter(MousePos.stimuli_onset(MousePos.stimuli_track==2),-0.5*ones(1,sum(MousePos.stimuli_track==2)),'b');
% legend([p s1 s2],{'Cumulative CA1 population bursting events with activation of track 1 preferred V1 neurons ','Track 1 stimuli','Track 2 stimuli'})
% legend([p1 p2 s1 s2],{'Events with Track 1 V1 neurons co-activation','All ripples','Track 1 stimuli','Track 2 stimuli'})
legend([p1 p2 p3],{'Events with Track 1 V1 neurons co-activation','Events with Track 1 CA1 cells activation','All Ripples'})
xlabel('time (s)')
% ylabel('Cumulative CA1 population bursting events events with activation of track 1 preferred V1 neurons (more than 4 neurons)')
ylabel('Cumulative ripple events')
% scatter(MousePos.sglxTime(locs),400*ones(1,length(locs)))





 [rho,pval] = corr(active_cells_SWR(1,find(active_cells_SWR(1,:)>1 & active_cells_SWR(2,:)==0))'...
     ,active_cells_SWR(3,find(active_cells_SWR(1,:)>1 & active_cells_SWR(2,:)==0))','Type','Spearman')

 [rho,pval] = corr(active_cells_SWR(1,find(active_cells_SWR(2,:)>1 & active_cells_SWR(1,:)==0))'...
     ,active_cells_SWR(3,find(active_cells_SWR(1,:)>1 & active_cells_SWR(2,:)==0))','Type','Spearman')


[rho,pval] = corr(active_cells_SWR(1,:)',active_cells_SWR(3,:)','Type','Spearman')
[rho,pval] = corr(active_cells_SWR(2,:)',active_cells_SWR(3,:)','Type','Spearman')

[rho,pval] = corr(sum(V1_track1_SWR_spikes,1)',sum(CA1_track1_SWR_spikes,1)','Type','Spearman')
[rho,pval] = corr(sum(V1_track2_SWR_spikes,1)',sum(CA1_track1_SWR_spikes,1)','Type','Spearman')


figure
subplot(2,2,1)
 hold on
 arrayfun(@(x) scatter(active_cells_SWR(1,x)',active_cells_SWR(3,x)',86,'k','filled','o'),1:size(active_cells_SWR,2))
subplot(2,2,2)
 hold on
 arrayfun(@(x) scatter(active_cells_SWR(2,x)',active_cells_SWR(3,x)',86,'k','filled','o'),1:size(active_cells_SWR,2))

% figure
% subplot(2,2,1)
%  hold on
%  arrayfun(@(x) scatter(sum(V1_track1_SWR_spikes(:,x))',sum(CA1_track1_SWR_spikes(:,x))',30,'k','filled','o'),1:size(active_cells_SWR,2))
% subplot(2,2,2)
%  hold on
%  arrayfun(@(x) scatter(sum(V1_track2_SWR_spikes(:,x))',sum(CA1_track1_SWR_spikes(:,x))',30,'k','filled','o'),1:size(active_cells_SWR,2))


 mdl2 = fitlm(active_cells_SWR(2,:)',active_cells_SWR(3,:)');
 mdl1 = fitlm(active_cells_SWR(1,:)',active_cells_SWR(3,:)');
 [pval,F_stat,~] = coefTest(mdl1);
  [pval,F_stat,~] = coefTest(mdl2);
R2 = mdl1.Rsquared.Adjusted;
R2 = mdl2.Rsquared.Adjusted;

 x =[min(awake_rate) max(awake_rate)];
 b = mdl.Coefficients.Estimate';
 y_est = polyval(fliplr(b),x);
 plot(x,y_est,':','Color','k','LineWidth',3)
 xlabel('Rate of awake replay rate (log2)')
 ylabel('Number of awake replay number (log2)')

V1_track1_SWR_spikes(V1_track1_SWR_spikes > 0)

CA1_track_1_cells
V1_track_1_cells