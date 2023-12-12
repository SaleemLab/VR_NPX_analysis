
%% Load Spike train data for place cell analysis
clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        column = 1;


        LFP_tvec = [];
        
        x_bins_width = 10; % bin for decoding
        
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            options = session_info(n).probe(nprobe);
            options.importMode = 'KS';
            options.gFileNum = gFileNum(n);
            probe_no = session_info(n).probe(nprobe).probe_id + 1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            options.probe_no = probe_no;
            options.ROOTPATH = ROOTPATH;
            % Load all spike data sorted according to the channel position

            [SUA.probe{probe_no} chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'group','by channel','cell_exporer','on');


            if ~isempty(best_channels{probe_no}.L4_channel)
                % L4 spike data (roughly 100 micron or 10 channels)
                [L4_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.L4_channel-5 best_channels{probe_no}.L4_channel+5],'group','by region','cell_exporer','on');
            end

            if ~isempty(best_channels{probe_no}.L5_channel)
                % [L4_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L4_channel-10 best_channels.L4_channel],'group','by region');
                % L5 spike data (roughly 200 micron or 20 channels))
                if isempty(best_channels{probe_no}.L4_channel) | find(chan_config_KS.Channel ==best_channels{probe_no}.L5_channel)+10 < find(chan_config_KS.Channel ==best_channels{probe_no}.L4_channel)-6   
                    [L5_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.L5_channel-10 best_channels{probe_no}.L5_channel+10],'group','by region','cell_exporer','on');
                else
                    [L5_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.L5_channel-10 best_channels{probe_no}.L4_channel-6],'group','by region','cell_exporer','on');
                end

                
            end

            % all V1 spike data (roughly 300 micron or 30 channels))
            [superficial_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.first_in_brain_channel-30 best_channels{probe_no}.first_in_brain_channel],'group','by region','cell_exporer','on');
            % [superficial_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.first_in_brain_channel-15 best_channels.first_in_brain_channel],'group','by region');

            % all V1 spike data
            if best_channels{probe_no}.first_in_brain_channel-100 > best_channels{probe_no}.CA1_channel+10
                [V1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.first_in_brain_channel-100 best_channels{probe_no}.first_in_brain_channel],'group','by region','cell_exporer','on');
            else
                [V1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.CA1_channel+11 best_channels{probe_no}.first_in_brain_channel],'group','by region','cell_exporer','on');
                
            end
            V1_place_fields.probe(probe_no) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters.probe(probe_no),[]);
            % CA1 spike data (roughly 200 micron or 20 channels))
            if ~isempty(best_channels{probe_no}.CA1_channel)
                [CA1_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.CA1_channel-10 best_channels{probe_no}.CA1_channel+10],'group','by region','cell_exporer','on');
                CA1_place_fields.probe(probe_no) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,CA1_clusters.probe(probe_no),[]);
            end
            % 'Whole' HPC from 100 micron above CA1 cell layer to 1000 micron
            % below
            [HPC_clusters.probe(probe_no) chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels{probe_no}.CA1_channel-100 best_channels{probe_no}.CA1_channel+10],'group','by region','cell_exporer','on');
            HPC_place_fields.probe(probe_no) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,HPC_clusters.probe(probe_no),[]);
        end
        
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load extracted_laps.mat

        save(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'HPC_clusters')
        save(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'CA1_clusters')
        save(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters')
        
        save(sprintf('extracted_superficial_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'superficial_clusters')
        if ~isempty(best_channels{probe_no}.L4_channel)
            save(sprintf('extracted_L4_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'L4_clusters')
        end
        if ~isempty(best_channels{probe_no}.L5_channel)
            save(sprintf('extracted_L5_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')),'L5_clusters')
        end

        clear HPC_clusters CA1_clusters V1_clusters superficial_clusters L5_clusters L4_clusters

        if contains(stimulus_name{n},'RUN')
            save extracted_HPC_place_fields HPC_place_fields
            save extracted_CA1_place_fields CA1_place_fields
            %             save extracted_CA1_place_fields_even CA1_place_fields_even
            %             save extracted_CA1_place_fields_odd CA1_place_fields_odd
            save extracted_V1_place_fields V1_place_fields
        end

       if length(session_info(n).probe)>1
            options = session_info(n).probe(1);
            V1_clusters_combined = combine_clusters_from_multiple_probes(V1_clusters);

            V1_place_fields_combined = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters_combined,[]);

            save('extracted_V1_place_fields_combined.mat','V1_place_fields_combined');
            save(sprintf('extracted_V1_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters_combined');
            close all
       end


        if length(session_info(n).probe) > 1
            clusters_combined = combine_clusters_from_multiple_probes(HPC_clusters);

            HPC_place_fields_combined = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,clusters_combined,[]);

        
            %         % Bayesian decoding with 10 folder cross validation
            %         for nprobe = 1:length(session_info(n).probe)
            %             estimated_position_lap_CV = bayesian_decoding_lap_cross_validation(HPC_clusters.probe(nprobe),HPC_place_fields.probe(nprobe),position,lap_times)
            %         end

%             estimated_position_lap_CV_HPC_combined = bayesian_decoding_RUN_lap_cross_validation(clusters_combined,HPC_place_fields_combined,position,lap_times)

            HPC_clusters_combined = clusters_combined;
            save('extracted_HPC_clusters_combined',HPC_clusters_combined)
            save(sprintf('extracted_HPC_place_fields_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')),'HPC_place_fields_combined')
            
        end

        clear  HPC_clusters HPC_place_fields_combined HPC_place_fields V1_place_fields CA1_place_fields
    end

    %         for track_id = 1:length(place_fields_all.track)
    %             bayesian_spike_count = create_spike_count_masa(place_fields_all,CA1_clusters.probe(nprobe),...
    %                 lap_times(track_id).start(1),lap_times(track_id).end(end),[]);
    %             estimated_position = bayesian_decoding(place_fields_all,bayesian_spike_count,position,[]);
    %         end

end


%% place cell and visual cell summary

clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_types_all = {'RUN','POST'};
Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;
x_bins_width = 10;
colour_line = {'r','b'};
multi_colours{1} = {[215,48,39]/256,[244,109,67]/256,[253,174,97]/256};
multi_colours{2} =      {[69,117,180]/256,[116,173,209]/256,[171,217,233]/256};


for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(session_info)
        continue
    end

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)

        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load extracted_laps
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        column = 1;
        load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('decoded_replay_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))

        load(sprintf('probability_ratio_original%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('probability_ratio_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        
        load(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load('extracted_HPC_place_fields_combined.mat')
        load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load('extracted_V1_place_fields.mat')
        load('extracted_CA1_place_fields.mat')
        load('extracted_HPC_place_fields.mat')



        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.gFileNum = gFileNum(n);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            options.ROOTPATH = ROOTPATH;


            %             CA1_place_fields_even.probe(nprobe)  = calculate_place_fields_masa_NPX(x_bins_width,position,CA1_clusters.probe(nprobe),'even laps');
            %             CA1_place_fields_odd.probe(nprobe) = calculate_place_fields_masa_NPX(x_bins_width,position,CA1_clusters.probe(nprobe),'odd laps');

            %             CA1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,CA1_clusters.probe(nprobe),[]);
            %             HPC_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,HPC_clusters.probe(nprobe),[]);
            %             V1_place_fields.probe(nprobe) = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters.probe(nprobe),[]);


            V1_place_fields_even.probe(nprobe)  = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters.probe(nprobe),'even laps');
            V1_place_fields_odd.probe(nprobe) = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters.probe(nprobe),'odd laps');
            HPC_place_fields_combined_even  = calculate_place_fields_masa_NPX(x_bins_width,position,HPC_clusters_combined,'even laps');
            HPC_place_fields_combined_odd = calculate_place_fields_masa_NPX(x_bins_width,position,HPC_clusters_combined,'odd laps');

            options.probe_combined = 1;
            options.region = 'HPC';
            place_fields_lap = plot_spatial_cell_tuning(HPC_clusters_combined,HPC_place_fields_combined,HPC_place_fields_combined_even,...
                HPC_place_fields_combined_odd,position,lap_times,options);
            save_all_figures((fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells')),[])
            options = rmfield(options,'probe_combined');

            options.region = 'V1';
            spatial_tuning_lap = plot_spatial_cell_tuning(V1_clusters,V1_place_fields.probe(probe_no),V1_place_fields_even.probe(probe_no),...
                V1_place_fields_odd.probe(probe_no),position,lap_times,options);
            save_all_figures((fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells')),[])
            close all
            
           
            HPC_all_cells{nsession}{options.probe_hemisphere} = unique(HPC_clusters.probe(probe_no).id_conversion(:,2));
            HPC_spatial_cells{nsession}{options.probe_hemisphere}{1} = HPC_clusters.probe(probe_no).id_conversion(HPC_place_fields.probe(probe_no).track(1).good_cells_LIBERAL,2)';
            HPC_spatial_cells{nsession}{options.probe_hemisphere}{2} = HPC_clusters.probe(probe_no).id_conversion(HPC_place_fields.probe(probe_no).track(2).good_cells_LIBERAL,2)';
            HPC_all_cells_FR_bias{nsession}{options.probe_hemisphere} = (HPC_place_fields.probe(probe_no).track(1).mean_rate_track-HPC_place_fields.probe(probe_no).track(2).mean_rate_track)...
                ./(HPC_place_fields.probe(probe_no).track(1).mean_rate_track + HPC_place_fields.probe(probe_no).track(2).mean_rate_track);

            V1_all_cells{nsession}{options.probe_hemisphere} = unique(V1_clusters.probe(probe_no).id_conversion(:,2));
            V1_spatial_cells{nsession}{options.probe_hemisphere}{1} = V1_clusters.probe(probe_no).id_conversion(V1_place_fields.probe(probe_no).track(1).good_cells_LIBERAL,2);
            V1_spatial_cells{nsession}{options.probe_hemisphere}{2} = V1_clusters.probe(probe_no).id_conversion(V1_place_fields.probe(probe_no).track(2).good_cells_LIBERAL,2);
            V1_all_cells_FR_bias{nsession}{options.probe_hemisphere} = (V1_place_fields.probe(probe_no).track(1).mean_rate_track-V1_place_fields.probe(probe_no).track(2).mean_rate_track)...
                ./(V1_place_fields.probe(probe_no).track(1).mean_rate_track + V1_place_fields.probe(probe_no).track(2).mean_rate_track);


            rmfield(options,'probe_combined')
            options.region = 'V1';
            x_bins_width = 5;
            V1_spatial_tuning_lap = lap_spatial_cell_tuning(V1_clusters,V1_place_fields.probe(probe_no),V1_place_fields_even.probe(probe_no),...
                V1_place_fields_odd.probe(probe_no),position,lap_times,x_bins_width,options);

%             for ncell = 1:length(V1_place_fields.probe(probe_no).all_cells)
%                 V1_spatial_tuning_lap{track_id}{nlap}
%             end

            for ncell = 1:length(V1_place_fields.probe(probe_no).track(track_id).good_cells_LIBERAL)
                V1_place_fields.probe(probe_no).track(track_id).raw{ncell}

                max_FR = max([V1_place_fields.probe(probe_no).track(1).raw{ncell} V1_place_fields.probe(probe_no).track(2).raw{ncell}]);
                plot(1:x_bins_width:140,V1_place_fields.probe(probe_no).track(1).raw{ncell}/max_FR); hold on; plot(1:x_bins_width:140,V1_place_fields.probe(probe_no).track(2).raw{ncell}/max_FR)
                title('even laps')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                xticks([30 50 70 90 110 140]) % 1-5 landmarks and end
                xline(100,'r')


                for track_id = 1:2
                    for ncell = V1_place_fields.probe(probe_no).all_cells
                        this_cell_spatial_tuning = [];
                        for nlap = 1:length(lap_times(track_id).start)
                            this_cell_spatial_tuning(nlap,:) = normalize(V1_spatial_tuning_lap{track_id}{nlap}{ncell},'range');

                            if sum(this_cell_spatial_tuning(nlap,8:11)) + sum(this_cell_spatial_tuning(nlap,16:19)) == 0
                                V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell}(nlap) = nan;
                            elseif max(this_cell_spatial_tuning(nlap,8:11)) < max(this_cell_spatial_tuning(nlap,16:19))
                                V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell}(nlap) = max(this_cell_spatial_tuning(nlap,8:11))/max(this_cell_spatial_tuning(nlap,16:19));
                            else
                                V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell}(nlap) = max(this_cell_spatial_tuning(nlap,16:19))/max(this_cell_spatial_tuning(nlap,8:11));
                            end


                        end
                        V1_spatial_modulation_session{track_id}(ncell) = mean(V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell}(isfinite(V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell})));
                    end
                    V1_spatial_modulation_session{track_id}(isnan(V1_spatial_modulation_session{track_id})) == 1;

                end


                for track_id = 1:2
                    for ncell = V1_place_fields.probe(probe_no).all_cells
                        this_cell_spatial_tuning = V1_place_fields.probe(probe_no).track(track_id).raw{ncell};

                        if sum(this_cell_spatial_tuning(4:5)) + sum(this_cell_spatial_tuning(8:9)) == 0
                            V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell} = nan;

                        elseif max(this_cell_spatial_tuning(4:5)) < max(this_cell_spatial_tuning(8:9))
                            V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell} = max(this_cell_spatial_tuning(4:5))/max(this_cell_spatial_tuning(8:9));

                        else
                            V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell} = max(this_cell_spatial_tuning(8:9))/max(this_cell_spatial_tuning(4:5));

                        end

                        V1_spatial_modulation_session{track_id}(ncell) = mean(V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell}(isfinite(V1_spatial_modulation{nsession}{probe_hemisphere}{track_id}{ncell})));
                    end

                    V1_spatial_modulation_session{track_id}(isnan(V1_spatial_modulation_session{track_id})) == 1;

                end
                histogram(V1_spatial_modulation_session{1},30,'FaceColor','r','FaceAlpha',0.3,'Normalization','cdf');hold on
                histogram(V1_spatial_modulation_session{2},30,'FaceColor','b','FaceAlpha',0.3,'Normalization','cdf')

            end

        end

    end
end

cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
save HPC_all_cells HPC_all_cells
save HPC_spatial_cells HPC_spatial_cells
save V1_all_cells V1_all_cells
save V1_spatial_cells V1_spatial_cells
save V1_all_cells_FR_bias V1_all_cells_FR_bias
save HPC_all_cells_FR_bias HPC_all_cells_FR_bias

cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
load HPC_all_cells
load HPC_spatial_cells
load V1_all_cells
load V1_spatial_cells
load V1_all_cells_FR_bias
load HPC_all_cells_FR_bias
hemisphere_text = {'left','right'};
c = 1;
close all
for nsession = 1:10
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'));
%     load('extracted_V1_place_fields.mat')
%     load('extracted_CA1_place_fields.mat')
%     load('extracted_HPC_place_fields.mat')
    for nprobe = 1:2

        if nprobe > length(session_info(n).probe)
            continue
        end
        option = session_info(n).probe(nprobe);
        probe_hemisphere = option.probe_hemisphere;



        if isempty(HPC_all_cells{nsession}{probe_hemisphere})
            continue
        else
            if probe_hemisphere == 1
                fig(1) = figure(1);
                fig(1).Position = [500 100 1200 900];
                fig(1).Name = 'HPC left probe cell info';
            else
                fig(2) = figure(2);
                fig(2).Position = [500 100 1200 900];
                fig(2).Name = 'HPC Right probe cell info';
            end


            Track_1_cells = setdiff(HPC_spatial_cells{nsession}{probe_hemisphere}{1},HPC_spatial_cells{nsession}{probe_hemisphere}{2});
            Track_2_cells = setdiff(HPC_spatial_cells{nsession}{probe_hemisphere}{2},HPC_spatial_cells{nsession}{probe_hemisphere}{1});
            common_cells = intersect(HPC_spatial_cells{nsession}{probe_hemisphere}{2},HPC_spatial_cells{nsession}{probe_hemisphere}{1});
            other_cells = setdiff(HPC_all_cells{nsession}{probe_hemisphere},unique([HPC_spatial_cells{nsession}{probe_hemisphere}{1} HPC_spatial_cells{nsession}{probe_hemisphere}{2}]));

            newLabels = [];
            data = [length(Track_1_cells) length(common_cells) length(Track_2_cells) length(other_cells)];
            %              Labels = ["T1 - ","Common - ","T2 - ",'Other - '];
            actual_data = data;
            for i=1:length(data)
                %                  newLabels = [newLabels {sprintf('%s %i%', Labels{i}, actual_data(i))}];
                newLabels = [newLabels {sprintf('%i%',actual_data(i))}];
            end

            Explode = [1 0 1 0];
            subplot(5,2,nsession)
            h = pie(data, Explode, newLabels);
            view([90 90])
            legend("T1","Common","T2",'Other','color','none')
            patchHand = findobj(h, 'Type', 'Patch'); 
            patchHand(1).FaceColor = 'red';
            patchHand(2).FaceColor = 'magenta';
            patchHand(3).FaceColor = 'blue';
            patchHand(4).FaceColor = 'yellow';
            for npie = 1:4
            patchHand(npie).FaceAlpha = 0.5;
            end
            title(sprintf('HPC cells session %i %s',nsession,hemisphere_text{probe_hemisphere}));
            view([90 90])
            legend("T1","Common","T2",'Other','color','none')
            set(gca,"TickDir","out",'box', 'off','Color','none')

            if probe_hemisphere == 1
                fig(3) = figure(3);
                fig(3).Position = [500 100 1200 900];
                fig(3).Name = 'V1 Left probe cell info';
            else
                fig(4) = figure(4);
                fig(4).Position = [500 100 1200 900];
                fig(4).Name = 'V1 Right probe cell info';
            end

            Track_1_cells = setdiff(V1_spatial_cells{nsession}{probe_hemisphere}{1},V1_spatial_cells{nsession}{probe_hemisphere}{2});
            Track_2_cells = setdiff(V1_spatial_cells{nsession}{probe_hemisphere}{2},V1_spatial_cells{nsession}{probe_hemisphere}{1});
            common_cells = intersect(V1_spatial_cells{nsession}{probe_hemisphere}{2},V1_spatial_cells{nsession}{probe_hemisphere}{1});
            other_cells = setdiff(V1_all_cells{nsession}{probe_hemisphere},unique([V1_spatial_cells{nsession}{probe_hemisphere}{1}; V1_spatial_cells{nsession}{probe_hemisphere}{2}]));

            newLabels = [];
            data = [length(Track_1_cells) length(common_cells) length(Track_2_cells) length(other_cells)];
            %              Labels = ["T1 - ","Common - ","T2 - ",'Other - '];
            actual_data = data;
            for i=1:length(data)
                %                  newLabels = [newLabels {sprintf('%s %i%', Labels{i}, actual_data(i))}];
                newLabels = [newLabels {sprintf('%i%',actual_data(i))}];
            end
            Explode = [1 0 1 0];
            subplot(5,2,nsession)
            h = pie(data, Explode, newLabels);
            view([90 90])
            legend("T1","Common","T2",'Other','color','none')
            title(sprintf('V1 cells session %i %s',nsession,hemisphere_text{probe_hemisphere}));
            set(gca,"TickDir","out",'box', 'off','Color','none')
            patchHand = findobj(h, 'Type', 'Patch'); 
            patchHand(1).FaceColor = 'red';
            patchHand(2).FaceColor = 'magenta';
            patchHand(3).FaceColor = 'blue';
            patchHand(4).FaceColor = 'yellow';
            for npie = 1:4
            patchHand(npie).FaceAlpha = 0.5;
            end

            if probe_hemisphere == 1
                fig(5) = figure(5);
                fig(5).Position = [500 100 1200 900];
                fig(5).Name = 'V1 Left probe Normalised firing rate difference (T1-T2)';
            else
                fig(6) = figure(6);
                fig(6).Position = [500 100 1200 900];
                fig(6).Name = 'V1 Right probe Normalised firing rate difference (T1-T2)';
            end

            subplot(5,2,nsession)
            [~,index1,~] = intersect(V1_all_cells{nsession}{probe_hemisphere},Track_1_cells);
            bar(1:length(index1),V1_all_cells_FR_bias{nsession}{probe_hemisphere}(index1),'FaceColor','r','FaceAlpha',0.5)
            hold on
            [~,index2,~] = intersect(V1_all_cells{nsession}{probe_hemisphere},Track_2_cells);
            bar_count = 1+length(index1);
            bar(bar_count:bar_count+length(index2)-1,V1_all_cells_FR_bias{nsession}{probe_hemisphere}(index2),'FaceColor','b','FaceAlpha',0.5)
            [~,index3,~] = intersect(V1_all_cells{nsession}{probe_hemisphere},common_cells);
            bar_count = 1+length(index1)+length(index2);
            bar(bar_count:bar_count+length(index3)-1,V1_all_cells_FR_bias{nsession}{probe_hemisphere}(index3),'FaceColor','m','FaceAlpha',0.5)
            [~,index4,~] = intersect(V1_all_cells{nsession}{probe_hemisphere},other_cells);
            bar_count = 1+length(index1)+length(index2)+length(index3);
            bar(bar_count:bar_count+length(index4)-1,...
                V1_all_cells_FR_bias{nsession}{probe_hemisphere}(index4),'FaceColor','y','FaceAlpha',0.5)
            xticks(1:4:length(V1_all_cells{nsession}{probe_hemisphere}))
            set(gca,"TickDir","out",'box', 'off','Color','none')
            title(sprintf('V1 cells session %i %s',nsession,hemisphere_text{probe_hemisphere}));
            ylabel('Normalized T1-T2 FR')
        end
    end
end

cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])


%% just decoding assuming place field calculation and cluster extraction is done.
clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
x_bins_width = 10;

for nsession =1:10
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load('extracted_laps.mat')
        column = 1;
        options = session_info(n).probe(1);
        options.ROOTPATH = ROOTPATH;
        if exist('spatial cells') == 0
            mkdir('spatial cells')
        end


        if length(session_info(n).probe) > 1
           %             estimated_position_lap_CV_HPC_combined = bayesian_decoding_RUN_lap_cross_validation(clusters_combined,HPC_place_fields_combined,position,lap_times)
            load('extracted_HPC_place_fields_combined')
            load(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            HPC_place_fields_combined_odd  = calculate_place_fields_masa_NPX(x_bins_width,position,HPC_clusters_combined,'even laps');
            HPC_place_fields_combined_even = calculate_place_fields_masa_NPX(x_bins_width,position,HPC_clusters_combined,'odd laps');

            options.probe_combined = 1;
            options.region = 'HPC';
            place_fields_lap = plot_spatial_cell_tuning(HPC_clusters_combined,HPC_place_fields_combined,HPC_place_fields_combined_even,HPC_place_fields_combined_odd,...
                position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])


            [normalised_raw_matrix,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = plot_place_cell_map_correlation(HPC_clusters_combined,HPC_place_fields_combined,HPC_place_fields_combined_even,HPC_place_fields_combined_odd,...
                position,lap_times,options);
            save_all_figures((fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells')),[])
            options = rmfield(options,'probe_combined');
            close all
            [probability_ratio_RUN_lap estimated_position_lap_CV_HPC_combined] = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters_combined,HPC_place_fields_combined,position,lap_times);
            save('estimated_position_lap_CV_HPC_combined.mat','estimated_position_lap_CV_HPC_combined')
            save('probability_ratio_RUN_lap.mat','probability_ratio_RUN_lap')
            estimated_position_lap_CV = estimated_position_lap_CV_HPC_combined;

        else

            load('extracted_HPC_place_fields')
            load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))

            HPC_place_fields_odd.probe(1)  = calculate_place_fields_masa_NPX(x_bins_width,position,HPC_clusters.probe(1),'even laps');
            HPC_place_fields_even.probe(1) = calculate_place_fields_masa_NPX(x_bins_width,position,HPC_clusters.probe(1),'odd laps');

            options.region = 'HPC';
            probe_no = session_info(n).probe(1).probe_id + 1;
            options.probe_no = probe_no;

            place_fields_lap = plot_spatial_cell_tuning(HPC_clusters,HPC_place_fields.probe(1),HPC_place_fields_even.probe(1),...
                HPC_place_fields_odd.probe(1),position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])


            [normalised_raw_matrix,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = plot_place_cell_map_correlation(HPC_clusters,HPC_place_fields.probe(1),HPC_place_fields_even.probe(1),...
                HPC_place_fields_odd.probe(1),position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

            [probability_ratio_RUN_lap estimated_position_lap_CV_HPC]  = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters.probe(1),HPC_place_fields.probe(1),position,lap_times)
            save('estimated_position_lap_CV_HPC.mat','estimated_position_lap_CV_HPC')
            save('probability_ratio_RUN_lap.mat','probability_ratio_RUN_lap')

            estimated_position_lap_CV = estimated_position_lap_CV_HPC;
        end


        pcount = 1;
        nfigure = 1;
        for track_id = 1:length(lap_times)
            for lap_id = lap_times(track_id).completeLaps_id(6:2:20)
                %             for lap_id = lap_times(track_id).completeLaps_id(6:2:40)
                if pcount == 17
                    nfigure = nfigure + 1;
                    pcount = 1;
                end

                fig = figure(nfigure)
                fig.Position = [300 150 945 800];
                fig.Name = (sprintf('%s %s CV HPC Bayesian decoding visualisation probe combined',options.SUBJECT,options.SESSION));
                subplot(4,4,pcount)
                if ~isempty(estimated_position_lap_CV(track_id).lap(lap_id))
                    imagesc([estimated_position_lap_CV(track_id).lap(lap_id).track(1).run; estimated_position_lap_CV(track_id).lap(lap_id).track(2).run])
                    hold on
                    plot(estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_actual_position/10,'r')
                    plot(estimated_position_lap_CV(track_id).lap(lap_id).track(2).run_actual_position/10 + 15,'b')
                    yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
                    yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
                    yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])
                    run_time_edges = estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_time_edges;

                    xticks(linspace(1,length(run_time_edges),5))
                    xticklabels(linspace(run_time_edges(1),run_time_edges(end),5))
                end
                set(gca,"TickDir","out",'box', 'off','Color','none')
                title(sprintf('Lap %i',lap_id))
                colorbar
                colormap(flip(bone))
                pcount = pcount + 1;
            end
        end
        sgtitle(sprintf('%s %s CV HPC Bayesian decoding visualisation probe combined',options.SUBJECT,options.SESSION))
        save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

        clear estimated_position_lap_CV estimated_position_lap_CV_HPC HPC_clusters HPC_place_fields_combined HPC_place_fields V1_place_fields CA1_place_fields
    end
end


%% just decoding for V1
clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
Hemisphere = {'Left','Right'};
x_bins_width = 10;

for nsession =1:10
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load('extracted_laps.mat')
        column = 1;

        load('extracted_V1_place_fields')
        load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))

        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no;
         
            options.region = 'V1';

            V1_place_fields_even.probe(probe_no)  = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters.probe(probe_no),'even laps');
            V1_place_fields_odd.probe(probe_no) = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters.probe(probe_no),'odd laps');

            plot_spatial_cell_tuning(V1_clusters,V1_place_fields.probe(probe_no),V1_place_fields_even.probe(probe_no),...
                V1_place_fields_odd.probe(probe_no),position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

            plot_place_cell_map_correlation(V1_clusters,V1_place_fields.probe(probe_no),V1_place_fields_even.probe(probe_no),...
                V1_place_fields_odd.probe(probe_no),position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

            [probability_ratio_RUN_lap_V1{probe_no} estimated_position_lap_CV_V1(probe_no).track]  = bayesian_decoding_RUN_lap_cross_validation(V1_clusters.probe(probe_no),V1_place_fields.probe(probe_no),position,lap_times)
            
            estimated_position_lap_CV = estimated_position_lap_CV_V1(probe_no).track;
            pcount = 1;
            nfigure = 1;
            for track_id = 1:length(lap_times)
                for lap_id = lap_times(track_id).completeLaps_id(6:2:20)
                    %             for lap_id = lap_times(track_id).completeLaps_id(6:2:40)
                    if pcount == 17
                        nfigure = nfigure + 1;
                        pcount = 1;
                    end

                    fig = figure(nfigure)
                    fig.Position = [300 150 945 800];
                    fig.Name = sprintf('%s %s CV V1 Bayesian decoding visualisation probe %s',options.SUBJECT,options.SESSION,Hemisphere{options.probe_hemisphere})
                    
                    subplot(4,4,pcount)
                    if ~isempty(estimated_position_lap_CV(track_id).lap(lap_id))
                        imagesc([estimated_position_lap_CV(track_id).lap(lap_id).track(1).run; estimated_position_lap_CV(track_id).lap(lap_id).track(2).run])
                        hold on
                        plot(estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_actual_position/10,'r')
                        plot(estimated_position_lap_CV(track_id).lap(lap_id).track(2).run_actual_position/10 + 15,'b')
                        yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
                        yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
                        yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])
                        run_time_edges = estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_time_edges;

                        xticks(linspace(1,length(run_time_edges),5))
                        xticklabels(linspace(run_time_edges(1),run_time_edges(end),5))
                    end
                    set(gca,"TickDir","out",'box', 'off','Color','none')
                    title(sprintf('Lap %i',lap_id))
                    colorbar
                    colormap(flip(bone))
                    pcount = pcount + 1;
                end
            end
            sgtitle(sprintf('%s %s CV V1 Bayesian decoding visualisation probe %s',options.SUBJECT,options.SESSION,Hemisphere{options.probe_hemisphere}))
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

        end

%         [probability_ratio_RUN_lap estimated_position_lap_CV_HPC]  = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters.probe(1),HPC_place_fields.probe(1),position,lap_times)
        save('estimated_position_lap_CV_V1.mat','estimated_position_lap_CV_V1')
        save('probability_ratio_RUN_lap_V1.mat','probability_ratio_RUN_lap_V1')
    end
end


%% Decoding for V1 combined
clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
Hemisphere = {'Left','Right'};
x_bins_width = 10;


for nsession =1:10
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load('extracted_laps.mat')
        column = 1;

        load('extracted_V1_place_fields')
        load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        if length(session_info(n).probe)>1
            options = session_info(n).probe(1);
            V1_clusters_combined = combine_clusters_from_multiple_probes(V1_clusters);

%             V1_place_fields_combined = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters_combined,[]);
% %             V1_place_fields_combined_odd  = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters_combined,'even laps');
% %             V1_place_fields_combined_even = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters_combined,'odd laps');
% % 
% %             options.probe_combined = 1;
% %             options.ROOTPATH = ROOTPATH;
% %             options.region = 'V1';
% %             place_fields_lap = plot_spatial_cell_tuning(V1_clusters_combined,V1_place_fields_combined,V1_place_fields_combined_even,...
% %                 V1_place_fields_combined_odd,position,lap_times,options);
% %             save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])
% % 
% %             [normalised_raw_matrix,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = plot_place_cell_map_correlation(V1_clusters,V1_place_fields_combined,V1_place_fields_combined_even,...
% %                 V1_place_fields_combined_odd,position,lap_times,options);
% %             save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])
% %             options = rmfield(options,'probe_combined');

%             save('extracted_V1_place_fields_combined.mat','V1_place_fields_combined');
            save(sprintf('extracted_V1_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters_combined');
            close all
        end
    end
end

for nsession =1:10
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
        load('extracted_laps.mat')
        column = 1;

        load('extracted_V1_place_fields')
        load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))

        if length(session_info(n).probe)>1
            options = session_info(n).probe(1);
            V1_clusters_combined = combine_clusters_from_multiple_probes(V1_clusters);

            V1_place_fields_combined = calculate_place_fields_masa_NPX_against_shuffle(x_bins_width,position,V1_clusters_combined,[]);
            V1_place_fields_combined_odd  = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters_combined,'even laps');
            V1_place_fields_combined_even = calculate_place_fields_masa_NPX(x_bins_width,position,V1_clusters_combined,'odd laps');

            options.probe_combined = 1;
            options.ROOTPATH = ROOTPATH;
            options.region = 'V1';
            place_fields_lap = plot_spatial_cell_tuning(V1_clusters_combined,V1_place_fields_combined,V1_place_fields_combined_even,...
                V1_place_fields_combined_odd,position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

            [normalised_raw_matrix,PPvector,shuffled_globalRemap_PPvector,shuffled_rateRemap_PPvector] = plot_place_cell_map_correlation(V1_clusters,V1_place_fields_combined,V1_place_fields_combined_even,...
                V1_place_fields_combined_odd,position,lap_times,options);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])
            options = rmfield(options,'probe_combined');

            save('extracted_V1_place_fields_combined.mat','V1_place_fields_combined');
            save(sprintf('extracted_V1_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')),'V1_clusters_combined');
            close all

            [probability_ratio_RUN_lap_V1_combined estimated_position_lap_CV_V1_combined.track]  = bayesian_decoding_RUN_lap_cross_validation(V1_clusters_combined,V1_place_fields_combined,position,lap_times)
            save('estimated_position_lap_CV_V1_combined.mat','estimated_position_lap_CV_V1_combined')
            save('probability_ratio_RUN_lap_V1_combined.mat','probability_ratio_RUN_lap_V1_combined')

            estimated_position_lap_CV = estimated_position_lap_CV_V1_combined.track;
            pcount = 1;
            nfigure = 1;
            for track_id = 1:length(lap_times)
                for lap_id = lap_times(track_id).completeLaps_id(6:2:20)
                    %             for lap_id = lap_times(track_id).completeLaps_id(6:2:40)
                    if pcount == 17
                        nfigure = nfigure + 1;
                        pcount = 1;
                    end

                    fig = figure(nfigure)
                    fig.Position = [300 150 945 800];
                    fig.Name = sprintf('%s %s CV V1 Bayesian decoding visualisation probe combined',options.SUBJECT,options.SESSION);
                    subplot(4,4,pcount)
                    if ~isempty(estimated_position_lap_CV(track_id).lap(lap_id))
                        imagesc([estimated_position_lap_CV(track_id).lap(lap_id).track(1).run; estimated_position_lap_CV(track_id).lap(lap_id).track(2).run])
                        colormap(flip(bone))
                        hold on
                        plot(estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_actual_position/10,'r')
                        plot(estimated_position_lap_CV(track_id).lap(lap_id).track(2).run_actual_position/10 + 15,'b')
                        yticks([30 50 70 90 110 140 170 190 210 230 250 280]/10)
                        yline(14.5,'LineWidth',2,'Color','k','DisplayName','Track 2')
                        yticklabels([30 50 70 90 110 140 30 50 70 90 110 140])
                        run_time_edges = estimated_position_lap_CV(track_id).lap(lap_id).track(1).run_time_edges;

                        xticks(linspace(1,length(run_time_edges),5))
                        xticklabels(linspace(run_time_edges(1),run_time_edges(end),5))
                    end
                    set(gca,"TickDir","out",'box', 'off','Color','none')
                    title(sprintf('Lap %i',lap_id))
                    colorbar
                    colormap(flip(bone))
                    pcount = pcount + 1;
                end
            end
            sgtitle(sprintf('%s %s CV V1 Bayesian decoding visualisation probe combined',options.SUBJECT,options.SESSION))
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

        else
            disp(sprintf('Only one porbe for session %i',nsession))
        end

     
%         [probability_ratio_RUN_lap estimated_position_lap_CV_HPC]  = bayesian_decoding_RUN_lap_cross_validation(HPC_clusters.probe(1),HPC_place_fields.probe(1),position,lap_times)

    end
end

% %% RUN log odds ROC curve
% clear all
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% % addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% % addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% %
% rmpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% rmpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'))
% rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% 
% SUBJECTS = {'M23017','M23028','M23029'};
% experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% % Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% % Stimulus_type = 'POST'; % extract LFP during RUN
% ROOTPATH = 'Z:\ibn-vision';
% c = 1;
% 
% 
% % fontsize(fig, 14, "points")
% for nsession =1:length(experiment_info)
%     tic
%     session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
% 
%     if isempty(session_info)
%         continue
%     end
% 
%     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
%         cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
%         load best_channels
%         load extracted_PSD
%         column = 1;
% 
%         load('probability_ratio_RUN_lap.mat')
%         % probability_ratio_RUN_lap{1}{track_id}
% 
%         z_log_odds = [];
%         track_label = [];
% 
%         % RUN lap log odds (for track 2 use 1/track 2 bias )
%         % Because decoding now speed thresholded (nan for wrong track)
%         for track_id = 1:length(probability_ratio_RUN_lap{1})
%             for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
%                 data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
% 
%                 if track_id == 1
%                     data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
%                 else
%                     data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
%                 end
%                
%                 for nshuffle = 1:1000
%                     if track_id == 1
%                         T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
%                     else
%                         T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
%                     end
%                 end
% 
%                 shuffled_data = log(T1_T2_ratio_shuffled);
%                 z_log_odds{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
%                 track_label{track_id}(nlap) = track_id;
% 
%                 HPC_bayesian_bias{track_id}(nlap) = sum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(1).run_bias)/...
%                     (sum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(1).run_bias)+sum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(2).run_bias)) ;
% 
%             end
%         end
% 
% 
%         fig(1) = figure(1);
%         fig(1).Position = [500 100 1200 900];
%         fig(1).Name = 'Place cell lap log odds distribution during RUN';
% 
%         subplot(2,5,nsession)
%         hold on
%         scatter(track_label{1}  .* (rand(1,length(z_log_odds{1}))),z_log_odds{1},'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')
%         hold on
%         scatter(track_label{2}./2 .* (2+rand(1,length(z_log_odds{2}))),z_log_odds{2},'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
%         xticks([0.5 2.5])
%         xticklabels({'Track 1','Track 2'})
%         xlabel('Track ID')
%         ylabel('zscored log odds')
%         set(gca,"TickDir","out",'box', 'off','Color','none')
%         sgtitle('Place cell lap log odds distribution')
% 
%         z_log_odds = [z_log_odds{1} z_log_odds{2}];
%         track_label = [track_label{1} track_label{2}];
% 
%         FPR = [];
%         TPR = [];
%         AUC = [];
%         z_log_odds_resampled = [];
%         track_label_resampled = [];
% 
%         for nboot = 1:1000
%             s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling
% 
%             index = datasample(s,1:length(z_log_odds),length(z_log_odds));
%             z_log_odds_resampled(nboot,:) = z_log_odds(index);
%             track_label_resampled(nboot,:) = track_label(index)-1;
%             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1,'NBoot',1);
%             %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);
% 
%             FPR = X;
%             TPR(nboot,:) = Y(:,1);
%             AUC(nboot) = A(1);
%         end
% 
%         fig(2) = figure(2);
%         fig(2).Position = [500 100 1200 900];
%         fig(2).Name = 'Place cell map log odds two track discrimination ROC during RUN';
%         % ROC#
%         subplot(2,5,nsession)
%         hold on
%         x = FPR';
%         CI_shuffle = prctile(TPR,[2.5 97.5]);
%         plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
%         plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
%         x2 = [x, fliplr(x)];
%         inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
%         h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);
% 
%         h(1) = plot([0 1],[0 1],'k--')
%         set(gca,"TickDir","out",'box', 'off','Color','none')
%         legend([h(2) h(1)],{'Real','chance'})
%         title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
%         sgtitle('Place cell map log odds two track discrimination ROC during RUN')
%     end
% 
% 
%     toc
% end
% cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
% 
% save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])
% 
% % 
% % saveas(gcf,'Place cell map log odds two track discrimination ROC during RUN.pdf')
% % saveas(gcf,'Place cell map log odds two track discrimination ROC during RUN.fig')
% 
% 
% %% RUN log odds ROC curve V1
% clear all
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% % addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% % addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% % 
% rmpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% rmpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'))
% rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% 
% SUBJECTS = {'M23017','M23028','M23029'};
% experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% % Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% % Stimulus_type = 'POST'; % extract LFP during RUN
% ROOTPATH = 'Z:\ibn-vision';
% c = 1;
% 
% 
% 
% % fontsize(fig, 14, "points")
% for nsession =1:length(experiment_info)
%     tic
%     session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
% 
%     if isempty(session_info)
%         continue
%     end
% 
%     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
%         cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
%         load best_channels
%         load extracted_PSD
%         column = 1;
% 
% 
% 
%         load('probability_ratio_RUN_lap_V1.mat')
%         % probability_ratio_RUN_lap{1}{track_id}
%         for nprobe = 1:length(session_info(n).probe)
%             options = session_info(n).probe(nprobe);
%             options.ROOTPATH = ROOTPATH;
%             probe_no = session_info(n).probe(nprobe).probe_id + 1;
%             options.probe_no = probe_no;
% 
%             probability_ratio_RUN_lap = probability_ratio_RUN_lap_V1{probe_no};
%             z_log_odds = [];
%             track_label = [];
% 
%             for track_id = 1:length(probability_ratio_RUN_lap{1})
%                 for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
%                     data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
% 
%                     if track_id == 1
%                         data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
%                     else
%                         data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
%                     end
% 
%                     for nshuffle = 1:1000
%                         if track_id == 1
%                             T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
%                         else
%                             T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
%                         end
%                     end
% 
%                     shuffled_data = log(T1_T2_ratio_shuffled);
%                     z_log_odds{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
%                     track_label{track_id}(nlap) = track_id;
%                 end
%             end
% 
%             if session_info(n).probe(nprobe).probe_hemisphere == 1
%                 fig(1) = figure(1);
%                 fig(1).Position = [500 100 1200 900];
%                 fig(1).Name = 'lap log odds ditribution in V1 for left probe';
%                 subplot(2,5,nsession)
%                 hold on
%                 scatter(track_label{1}  .* (rand(1,length(z_log_odds{1}))),z_log_odds{1},'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')
%                 hold on
%                 scatter(track_label{2}./2 .* (2+rand(1,length(z_log_odds{2}))),z_log_odds{2},'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
%                 xticks([0.5 2.5])
%                 xticklabels({'Track 1','Track 2'})
%                 xlabel('Track ID')
%                 ylabel('zscored log odds')
%                 set(gca,"TickDir","out",'box', 'off','Color','none')
%                 sgtitle('lap log odds ditribution in V1 for left probe')
%             elseif  session_info(n).probe(nprobe).probe_hemisphere == 2
%                 fig(2) = figure(2);
%                 fig(2).Position = [500 100 1200 900];
%                 fig(2).Name = 'lap log odds ditribution in V1 for right probe';
%                 subplot(2,5,nsession)
%                 hold on
%                 scatter(track_label{1}  .* (rand(1,length(z_log_odds{1}))),z_log_odds{1},'MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha','0.1')
%                 hold on
%                 scatter(track_label{2}./2 .* (2+rand(1,length(z_log_odds{2}))),z_log_odds{2},'MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha','0.1')
%                 xticks([0.5 2.5])
%                 xticklabels({'Track 1','Track 2'})
%                 xlabel('Track ID')
%                 ylabel('zscored log odds')
%                 set(gca,"TickDir","out",'box', 'off','Color','none')
%                 sgtitle('lap log odds ditribution in V1 for right probe')
%             end
% 
% 
% 
%             z_log_odds = [z_log_odds{1} z_log_odds{2}];
%             track_label = [track_label{1} track_label{2}];
% 
%             FPR = [];
%             TPR = [];
%             AUC = [];
%             z_log_odds_resampled = [];
%             track_label_resampled = [];
% 
%             for nboot = 1:1000
%                 s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling
% 
%                 index = datasample(s,1:length(z_log_odds),length(z_log_odds));
%                 z_log_odds_resampled(nboot,:) = z_log_odds(index);
%                 track_label_resampled(nboot,:) = track_label(index)-1;
%                 [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1,'NBoot',1);
%                 %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);
% 
%                 FPR = X;
%                 TPR(nboot,:) = Y(:,1);
%                 AUC(nboot) = A(1);
%             end
% 
%             if session_info(n).probe(nprobe).probe_hemisphere == 1
%                 fig(3) = figure(3);
%                 fig(3).Position = [500 100 1200 900];
%                 fig(3).Name = 'lap log odds ROC two track discrimination in V1 for left probe';
%                 % ROC#
%                 subplot(2,5,nsession)
%                 hold on
%                 x = FPR';
%                 CI_shuffle = prctile(TPR,[2.5 97.5]);
%                 plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
%                 plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
%                 x2 = [x, fliplr(x)];
%                 inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
%                 h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);
% 
%                 h(1) = plot([0 1],[0 1],'k--')
%                 set(gca,"TickDir","out",'box', 'off','Color','none')
%                 legend([h(2) h(1)],{'Real','chance'})
%                 title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
%                 sgtitle('lap log odds ROC two track discrimination in V1 for left probe')
%             elseif session_info(n).probe(nprobe).probe_hemisphere == 2
%                 fig(4) = figure(4);
%                 fig(4).Position = [500 100 1200 900];
%                 fig(4).Name = 'lap log odds ROC two track discrimination in V1 for right probe';
%                 % ROC#
%                 subplot(2,5,nsession)
%                 hold on
%                 x = FPR';
%                 CI_shuffle = prctile(TPR,[2.5 97.5]);
%                 plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
%                 plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
%                 x2 = [x, fliplr(x)];
%                 inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
%                 h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);
% 
%                 h(1) = plot([0 1],[0 1],'k--')
%                 set(gca,"TickDir","out",'box', 'off','Color','none')
%                 legend([h(2) h(1)],{'Real','chance'})
%                 title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
%                 sgtitle('lap log odds ROC two track discrimination in V1 for right probe')
%             end
%         end
%     end
%     toc
% end
% cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
% save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])
% 
% saveas(gcf,'V1 cell lap log odds distribution during RUN.pdf')
% saveas(gcf,'V1 cell lap log odds distribution during RUN.fig')
% 
% % sgtitle('Place cell map log odds two track discrimination ROC during RUN')
% % saveas(gcf,'Place cell map log odds two track discrimination ROC during RUN.pdf')
% % saveas(gcf,'Place cell map log odds two track discrimination ROC during RUN.fig')

%% RUN lap log odds HPC and V1 interaction
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% 
rmpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
rmpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'))
rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_types_all = {'RUN','POST'};
Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;


% fontsize(fig, 14, "points")
for nsession =1:length(experiment_info)
    tic
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(session_info)
        continue
    end

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
                load best_channels
                load extracted_PSD
        load extracted_laps
        column = 1;


        load('probability_ratio_RUN_lap.mat')
        load('probability_ratio_RUN_lap_V1.mat')
        load estimated_position_lap_CV_V1
%         estimated_position_lap_CV_V1_temp = estimated_position_lap_CV_V1;
%         clear estimated_position_lap_CV_V1
% 
%         estimated_position_lap_CV_V1(1).track = estimated_position_lap_CV_V1_temp(1).probe;
%         if length(session_info(n).probe)>1
%         estimated_position_lap_CV_V1(2).track = estimated_position_lap_CV_V1_temp(2).probe;
%         end
%         save estimated_position_lap_CV_V1 estimated_position_lap_CV_V1
% 
%     end
% end


        if length(session_info(n).probe)>1
            load estimated_position_lap_CV_HPC_combined
            estimated_position_lap_CV_HPC = estimated_position_lap_CV_HPC_combined;
        else
            load estimated_position_lap_CV_HPC
        end

        z_log_odds = [];
        track_label = [];
        HPC_bayesian_bias = [];
        V1_z_log_odds= [];
        V1_track_label= [];
        V1_bayesian_bias = [];

        % RUN lap log odds (for track 2 use 1/track 2 bias )
        % Because decoding now speed thresholded (nan for wrong track)
        for track_id = 1:length(probability_ratio_RUN_lap{1})
            for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
                data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));

                if track_id == 1
                    data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                else
                    data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                end

                for nshuffle = 1:1000
                    if track_id == 1
                        T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                    else
                        T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                    end
                end

                shuffled_data = log(T1_T2_ratio_shuffled);
                z_log_odds{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
                track_label{track_id}(nlap) = track_id;

                HPC_bayesian_bias{track_id}(nlap) = nansum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(1).run_bias)/...
                    (nansum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(1).run_bias)+nansum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(2).run_bias)) ;

            end
        end


        % probability_ratio_RUN_lap{1}{track_id}
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
            options.probe_no = probe_no;

            probability_ratio_RUN_lap = probability_ratio_RUN_lap_V1{probe_no};

            V1_z_log_odds{probe_hemisphere} = [];
            V1_track_label{probe_hemisphere} = [];
            V1_bayesian_bias{probe_hemisphere} = [];

            for track_id = 1:length(probability_ratio_RUN_lap{1})
                for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
                    data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));

                    if track_id == 1
                        data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                    else
                        data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                    end

                    for nshuffle = 1:1000
                        if track_id == 1
                            T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                        else
                            T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                        end
                    end

                    shuffled_data = log(T1_T2_ratio_shuffled);
                    V1_z_log_odds{probe_hemisphere}{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
                    V1_track_label{probe_hemisphere}{track_id}(nlap) = track_id;


                    V1_bayesian_bias{probe_hemisphere}{track_id}(nlap) = nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(1).run_bias)/...
                        (nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(1).run_bias)+nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(2).run_bias)) ;

                end
            end

        end


        z_log_odds = [z_log_odds{1}, z_log_odds{2}];
        HPC_bayesian_bias = [HPC_bayesian_bias{1} HPC_bayesian_bias{2}];
        if ~isempty(V1_z_log_odds{1})
            V1_z_log_odds{1} = [V1_z_log_odds{1}{1}, V1_z_log_odds{1}{2}];
            V1_bayesian_bias{1} = [V1_bayesian_bias{1}{1}, V1_bayesian_bias{1}{2}];
        end

        if ~isempty(V1_z_log_odds{2})
            V1_z_log_odds{2} = [V1_z_log_odds{2}{1}, V1_z_log_odds{2}{2}];
            V1_bayesian_bias{2} = [V1_bayesian_bias{2}{1}, V1_bayesian_bias{2}{2}];
        end





        %         subplot(2,3,2)
        %         bar(V1_bayesian_bias{2}(sorted_id),'b','FaceAlpha',0.3,'EdgeColor','none');hold on;
        %
        %         bar(HPC_bayesian_bias(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        %         for nlap = 1:length(track_orders)
        %             if track_orders(nlap) == 1
        %                 scatter(nlap,track_orders(nlap) ,'r')
        %             elseif track_orders(nlap) == 2
        %                 scatter(nlap,track_orders(nlap) -2.2,'b')
        %             end
        %         end
        
        [~,sorted_id] = sort([lap_times(1).start  lap_times(2).start]);
        track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];
        track_orders = track_orders(sorted_id);

        fig(nsession) = figure(nsession);
        fig(nsession).Position = [500 100 1200 900];
        fig(nsession).Name = sprintf('%s %s RUN lap bayesian bias and log odds',options.SUBJECT,options.SESSION)
        colour_lines = {'b','r'};
        clear h s
        sgtitle(sprintf('%s %s RUN lap bayesian bias and log odds',options.SUBJECT,options.SESSION))
        subplot(2,2,1)
        for nprobe = 1:length(session_info(n).probe)
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            h(nprobe) = plot(V1_bayesian_bias{probe_hemisphere}(sorted_id),colour_lines{probe_hemisphere});hold on;
        end
        %         plot(HPC_bayesian_bias(sorted_id),'k');hold on;

        h(3) = bar(HPC_bayesian_bias(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        for nlap = 1:length(track_orders)
            if track_orders(nlap) == 1
                s(1) = scatter(nlap,track_orders(nlap) ,3,'r','filled','MarkerFaceAlpha',1)
            elseif track_orders(nlap) == 2
                s(2) = scatter(nlap,track_orders(nlap) - 2.1,3,'b','filled','MarkerFaceAlpha',1)
            end
        end
        ylabel('Bayesian Bias')
        xlabel('lap id')
        yline(0.5,'--')
        set(gca,"TickDir","out",'box', 'off','Color','none')

        %         legend([h(1:3),s(1),s(2)],{'V1 Left','V1 Right','HPC','Track 1','Track 2'})



        subplot(2,2,2)
        for nprobe = 1:length(session_info(n).probe)
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            h(probe_hemisphere) = plot(V1_z_log_odds{probe_hemisphere}(sorted_id),colour_lines{probe_hemisphere});hold on;
        end
        %         plot(HPC_bayesian_bias(sorted_id),'k');hold on;

        h(3) = bar(z_log_odds(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        for nlap = 1:length(track_orders)
            if track_orders(nlap) == 1
                s(1) = scatter(nlap,track_orders(nlap) +2.5,3,'r','filled','MarkerFaceAlpha',1)
            elseif track_orders(nlap) == 2
                s(2) = scatter(nlap,track_orders(nlap) -5.5,3,'b','filled','MarkerFaceAlpha',1)
            end
        end
        ylabel('Log odds (z)')
        xlabel('lap id')
        yline(0,'--')
        set(gca,"TickDir","out",'box', 'off','Color','none')
        
        if isempty(V1_z_log_odds{2})
            legend([h(1),h(3),s(1),s(2)],{'V1 Left','HPC','Track 1','Track 2'})
        elseif isempty(V1_z_log_odds{1})
            legend([h(2:3),s(1),s(2)],{'V1 Right','HPC','Track 1','Track 2'})
        else
            legend([h(1:3),s(1),s(2)],{'V1 Left','V1 Right','HPC','Track 1','Track 2'})
        end


        track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];

        subplot(2,4,5)
        if ~isempty(V1_bayesian_bias{1})
            scatter(HPC_bayesian_bias(track_orders == 1), V1_bayesian_bias{1}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(HPC_bayesian_bias(track_orders == 2), V1_bayesian_bias{1}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(HPC_bayesian_bias',V1_bayesian_bias{1});
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias') max(HPC_bayesian_bias')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Left Bayesian bias')
            xlabel('HPC Bayesian bias')
            yline(0.5,'--')
            xline(0.5,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,6)
        if ~isempty(V1_bayesian_bias{2})
            scatter(HPC_bayesian_bias(track_orders == 1), V1_bayesian_bias{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(HPC_bayesian_bias(track_orders == 2), V1_bayesian_bias{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(HPC_bayesian_bias',V1_bayesian_bias{2});
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias') max(HPC_bayesian_bias')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Right Bayesian bias')
            xlabel('HPC Bayesian bias')
            yline(0.5,'--')
            xline(0.5,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,7)
        if ~isempty(V1_z_log_odds{1})
            scatter(z_log_odds(track_orders == 1), V1_z_log_odds{1}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(z_log_odds(track_orders == 2), V1_z_log_odds{1}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(z_log_odds',V1_z_log_odds{1});
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds') max(z_log_odds')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Left Log odds (z)')
            xlabel('HPC log odds (z)')
            yline(0,'--')
            xline(0,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,8)
        if ~isempty(V1_z_log_odds{2})
            scatter(z_log_odds(track_orders == 1), V1_z_log_odds{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(z_log_odds(track_orders == 2), V1_z_log_odds{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(z_log_odds',V1_z_log_odds{2});
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds') max(z_log_odds')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                %         xlim([-30 120])

                %     title(sprintf('Session %i',s),'Color','red')
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                %         xlim([-30 120])
                %     title(sprintf('Session %i',s))
                %         text(gca,.7,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontSize',12,'FontName','Arial');
                text(gca,0.5,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end
            ylabel('V1 Right Log odds (z)')
            xlabel('HPC log odds (z)')
            yline(0,'--')
            xline(0,'--')
            %             track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];
            %             scatter3(z_log_odds(track_orders == 1), V1_z_log_odds{1}(track_orders == 1),V1_z_log_odds{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.5)
            %             hold on
            %             scatter3(z_log_odds(track_orders == 2), V1_z_log_odds{1}(track_orders == 2),V1_z_log_odds{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.5)
            %             xlabel('HPC log odds')
            %             ylabel('left V1 log odds')
            %             zlabel('right V1 log odds')
            %             set(gca,"TickDir","out",'box', 'off','Color','none')
            %             legend([h(1:3)],{'V1 Left','V1 Right','HPC'})

        end
        set(gca,"TickDir","out",'box', 'off','Color','none')
        fontsize(gcf,14,"points")


        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            data = V1_bayesian_bias{probe_hemisphere};
            

            FPR = [];
            TPR = [];
            AUC = [];
            data_resampled = [];
            track_label_resampled = [];

            for nboot = 1:1000
                s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

                index = datasample(s,1:length(data),length(data));
                data_resampled(nboot,:) = data(index);
                track_label_resampled(nboot,:) = track_orders(index)-1;
                [X,Y,T,A] = perfcurve(track_orders(index),data(index),1,'XVals',0:0.05:1,'NBoot',1);
                %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

                FPR = X;
                TPR(nboot,:) = Y(:,1);
                AUC(nboot) = A(1);
            end

            if session_info(n).probe(nprobe).probe_hemisphere == 1
                fig(11) = figure(11);
                fig(11).Position = [500 100 1200 900];
                fig(11).Name = 'lap Bayesian Bias ROC two track discrimination in V1 for left probe';
                % ROC#
                subplot(2,5,nsession)
                hold on
                x = FPR';
                CI_shuffle = prctile(TPR,[2.5 97.5]);
                plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
                plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
                x2 = [x, fliplr(x)];
                inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
                h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

                h(1) = plot([0 1],[0 1],'k--')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                legend([h(2) h(1)],{'Real','chance'})
                title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
                sgtitle('lap Bayesian Bias ROC two track discrimination in V1 for left probe')
                        fontsize(gcf,14,"points")
            elseif session_info(n).probe(nprobe).probe_hemisphere == 2
                fig(12) = figure(12);
                fig(12).Position = [500 100 1200 900];
                fig(12).Name = 'lap Bayesian Bias ROC two track discrimination in V1 for right probe';
                % ROC#
                subplot(2,5,nsession)
                hold on
                x = FPR';
                CI_shuffle = prctile(TPR,[2.5 97.5]);
                plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
                plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
                x2 = [x, fliplr(x)];
                inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
                h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

                h(1) = plot([0 1],[0 1],'k--')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                legend([h(2) h(1)],{'Real','chance'})
                title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
                sgtitle('lap Bayesian Bias ROC two track discrimination in V1 for right probe')
                        fontsize(gcf,14,"points")
            end
        end

        data = HPC_bayesian_bias;

        FPR = [];
        TPR = [];
        AUC = [];
        data_resampled = [];
        track_label_resampled = [];

        for nboot = 1:1000
            s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

            index = datasample(s,1:length(data),length(data));
            data_resampled(nboot,:) = data(index);
            track_label_resampled(nboot,:) = track_orders(index)-1;
            [X,Y,T,A] = perfcurve(track_orders(index),data(index),1,'XVals',0:0.05:1,'NBoot',1);
            %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

            FPR = X;
            TPR(nboot,:) = Y(:,1);
            AUC(nboot) = A(1);
        end

        fig(13) = figure(13);
        fig(13).Position = [500 100 1200 900];
        fig(13).Name = 'lap Bayesian Bias ROC two track discrimination in HPC';

        % ROC#
        subplot(2,5,nsession)
        hold on
        x = FPR';
        CI_shuffle = prctile(TPR,[2.5 97.5]);
        plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
        plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
        x2 = [x, fliplr(x)];
        inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
        h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

        h(1) = plot([0 1],[0 1],'k--')
        set(gca,"TickDir","out",'box', 'off','Color','none')
        legend([h(2) h(1)],{'Real','chance'})
        title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
        sgtitle('lap Bayesian Bias ROC two track discrimination in HPC')
        fontsize(gcf,14,"points")

    end
    toc
end

cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])



% sgtitle('Place cell map log odds two track discrimination ROC during RUN')
% saveas(gcf,'Place cell map log odds two track discrimination ROC during RUN.pdf')
% saveas(gcf,'Place cell map log odds two track discrimination ROC during RUN.fig')

%% RUN lap log odds HPC and V1 interaction (with regression within track)
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% 
rmpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
rmpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'))
rmpath(genpath('Z:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_types_all = {'RUN','POST'};
Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;


% fontsize(fig, 14, "points")
for nsession =1:length(experiment_info)
    tic
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(session_info)
        continue
    end

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
                load best_channels
                load extracted_PSD
        load extracted_laps
        column = 1;


        load('probability_ratio_RUN_lap.mat')
        load('probability_ratio_RUN_lap_V1.mat')
        load estimated_position_lap_CV_V1
%         estimated_position_lap_CV_V1_temp = estimated_position_lap_CV_V1;
%         clear estimated_position_lap_CV_V1
% 
%         estimated_position_lap_CV_V1(1).track = estimated_position_lap_CV_V1_temp(1).probe;
%         if length(session_info(n).probe)>1
%         estimated_position_lap_CV_V1(2).track = estimated_position_lap_CV_V1_temp(2).probe;
%         end
%         save estimated_position_lap_CV_V1 estimated_position_lap_CV_V1
% 
%     end
% end


        if length(session_info(n).probe)>1
            load estimated_position_lap_CV_HPC_combined
            estimated_position_lap_CV_HPC = estimated_position_lap_CV_HPC_combined;
        else
            load estimated_position_lap_CV_HPC
        end

        z_log_odds = [];
        track_label = [];
        HPC_bayesian_bias = [];
        V1_z_log_odds= [];
        V1_track_label= [];
        V1_bayesian_bias = [];

        % RUN lap log odds (for track 2 use 1/track 2 bias )
        % Because decoding now speed thresholded (nan for wrong track)
        for track_id = 1:length(probability_ratio_RUN_lap{1})
            for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
                data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));

                if track_id == 1
                    data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                else
                    data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                end

                for nshuffle = 1:1000
                    if track_id == 1
                        T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                    else
                        T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                    end
                end

                shuffled_data = log(T1_T2_ratio_shuffled);
                z_log_odds{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
                track_label{track_id}(nlap) = track_id;

                HPC_bayesian_bias{track_id}(nlap) = nansum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(1).run_bias)/...
                    (nansum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(1).run_bias)+nansum(estimated_position_lap_CV_HPC(track_id).lap(nlap).track(2).run_bias)) ;

            end
        end


        % probability_ratio_RUN_lap{1}{track_id}
        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere ;
            options.probe_no = probe_no;

            probability_ratio_RUN_lap = probability_ratio_RUN_lap_V1{probe_no};

            V1_z_log_odds{probe_hemisphere} = [];
            V1_track_label{probe_hemisphere} = [];
            V1_bayesian_bias{probe_hemisphere} = [];

            for track_id = 1:length(probability_ratio_RUN_lap{1})
                for nlap = 1:length(probability_ratio_RUN_lap{1}{track_id}{1})
                    data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));

                    if track_id == 1
                        data = log(cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                    else
                        data = log(1/cell2mat(probability_ratio_RUN_lap{1}{track_id}{track_id}(nlap)));
                    end

                    for nshuffle = 1:1000
                        if track_id == 1
                            T1_T2_ratio_shuffled(nshuffle) = cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                        else
                            T1_T2_ratio_shuffled(nshuffle) = 1/cell2mat(probability_ratio_RUN_lap{2}{nshuffle}{track_id}{track_id}(nlap));
                        end
                    end

                    shuffled_data = log(T1_T2_ratio_shuffled);
                    V1_z_log_odds{probe_hemisphere}{track_id}(nlap) = (data-mean(shuffled_data))/std(shuffled_data);
                    V1_track_label{probe_hemisphere}{track_id}(nlap) = track_id;


                    V1_bayesian_bias{probe_hemisphere}{track_id}(nlap) = nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(1).run_bias)/...
                        (nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(1).run_bias)+nansum(estimated_position_lap_CV_V1(probe_no).track(track_id).lap(nlap).track(2).run_bias)) ;

                end
            end

        end


        z_log_odds = [z_log_odds{1}, z_log_odds{2}];
        HPC_bayesian_bias = [HPC_bayesian_bias{1} HPC_bayesian_bias{2}];
        if ~isempty(V1_z_log_odds{1})
            V1_z_log_odds{1} = [V1_z_log_odds{1}{1}, V1_z_log_odds{1}{2}];
            V1_bayesian_bias{1} = [V1_bayesian_bias{1}{1}, V1_bayesian_bias{1}{2}];
        end

        if ~isempty(V1_z_log_odds{2})
            V1_z_log_odds{2} = [V1_z_log_odds{2}{1}, V1_z_log_odds{2}{2}];
            V1_bayesian_bias{2} = [V1_bayesian_bias{2}{1}, V1_bayesian_bias{2}{2}];
        end





        %         subplot(2,3,2)
        %         bar(V1_bayesian_bias{2}(sorted_id),'b','FaceAlpha',0.3,'EdgeColor','none');hold on;
        %
        %         bar(HPC_bayesian_bias(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        %         for nlap = 1:length(track_orders)
        %             if track_orders(nlap) == 1
        %                 scatter(nlap,track_orders(nlap) ,'r')
        %             elseif track_orders(nlap) == 2
        %                 scatter(nlap,track_orders(nlap) -2.2,'b')
        %             end
        %         end
        
        [~,sorted_id] = sort([lap_times(1).start  lap_times(2).start]);
        track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];
        track_orders = track_orders(sorted_id);

        fig(nsession) = figure(nsession);
        fig(nsession).Position = [500 100 1200 900];
        fig(nsession).Name = sprintf('%s %s RUN lap bayesian bias and log odds',options.SUBJECT,options.SESSION)
        colour_lines = {'b','r'};
        clear h s
        sgtitle(sprintf('%s %s RUN lap bayesian bias and log odds',options.SUBJECT,options.SESSION))
        subplot(2,2,1)
        for nprobe = 1:length(session_info(n).probe)
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            h(nprobe) = plot(V1_bayesian_bias{probe_hemisphere}(sorted_id),colour_lines{probe_hemisphere});hold on;
        end
        %         plot(HPC_bayesian_bias(sorted_id),'k');hold on;

        h(3) = bar(HPC_bayesian_bias(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        for nlap = 1:length(track_orders)
            if track_orders(nlap) == 1
                s(1) = scatter(nlap,track_orders(nlap) ,3,'r','filled','MarkerFaceAlpha',1)
            elseif track_orders(nlap) == 2
                s(2) = scatter(nlap,track_orders(nlap) - 2.1,3,'b','filled','MarkerFaceAlpha',1)
            end
        end
        ylabel('Bayesian Bias')
        xlabel('lap id')
        yline(0.5,'--')
        set(gca,"TickDir","out",'box', 'off','Color','none')

        %         legend([h(1:3),s(1),s(2)],{'V1 Left','V1 Right','HPC','Track 1','Track 2'})



        subplot(2,2,2)
        for nprobe = 1:length(session_info(n).probe)
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            h(probe_hemisphere) = plot(V1_z_log_odds{probe_hemisphere}(sorted_id),colour_lines{probe_hemisphere});hold on;
        end
        %         plot(HPC_bayesian_bias(sorted_id),'k');hold on;

        h(3) = bar(z_log_odds(sorted_id),'k','FaceAlpha',0.3,'EdgeColor','none');hold on;
        for nlap = 1:length(track_orders)
            if track_orders(nlap) == 1
                s(1) = scatter(nlap,track_orders(nlap) +2.5,3,'r','filled','MarkerFaceAlpha',1)
            elseif track_orders(nlap) == 2
                s(2) = scatter(nlap,track_orders(nlap) -5.5,3,'b','filled','MarkerFaceAlpha',1)
            end
        end
        ylabel('Log odds (z)')
        xlabel('lap id')
        yline(0,'--')
        set(gca,"TickDir","out",'box', 'off','Color','none')
        
        if isempty(V1_z_log_odds{2})
            legend([h(1),h(3),s(1),s(2)],{'V1 Left','HPC','Track 1','Track 2'})
        elseif isempty(V1_z_log_odds{1})
            legend([h(2:3),s(1),s(2)],{'V1 Right','HPC','Track 1','Track 2'})
        else
            legend([h(1:3),s(1),s(2)],{'V1 Left','V1 Right','HPC','Track 1','Track 2'})
        end


        track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];

        subplot(2,4,5)
        if ~isempty(V1_bayesian_bias{1})
            scatter(HPC_bayesian_bias(track_orders == 1), V1_bayesian_bias{1}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(HPC_bayesian_bias(track_orders == 2), V1_bayesian_bias{1}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(HPC_bayesian_bias(track_orders == 1)',V1_bayesian_bias{1}(track_orders == 1));
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias(track_orders == 1)') max(HPC_bayesian_bias(track_orders == 1)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            mdl = fitlm(HPC_bayesian_bias(track_orders == 2)',V1_bayesian_bias{1}(track_orders == 2));
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias(track_orders == 2)') max(HPC_bayesian_bias(track_orders == 2)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'b:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','blue');
            else
                plot(x,y_est,'k:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end


            ylabel('V1 Left Bayesian bias')
            xlabel('HPC Bayesian bias')
            yline(0.5,'--')
            xline(0.5,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,6)
        if ~isempty(V1_bayesian_bias{2})
            scatter(HPC_bayesian_bias(track_orders == 1), V1_bayesian_bias{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(HPC_bayesian_bias(track_orders == 2), V1_bayesian_bias{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(HPC_bayesian_bias(track_orders == 1)',V1_bayesian_bias{2}(track_orders == 1));
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias(track_orders == 1)') max(HPC_bayesian_bias(track_orders == 1)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            mdl = fitlm(HPC_bayesian_bias(track_orders == 2)',V1_bayesian_bias{2}(track_orders == 2));
            [pval,~,~] = coefTest(mdl);
            x =[min(HPC_bayesian_bias(track_orders == 2)') max(HPC_bayesian_bias(track_orders == 2)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);


            if pval <= 0.05
                plot(x,y_est,'b:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','blue');
            else
                plot(x,y_est,'k:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            ylabel('V1 Right Bayesian bias')
            xlabel('HPC Bayesian bias')
            yline(0.5,'--')
            xline(0.5,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,7)
        if ~isempty(V1_z_log_odds{1})
            scatter(z_log_odds(track_orders == 1), V1_z_log_odds{1}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(z_log_odds(track_orders == 2), V1_z_log_odds{1}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(z_log_odds(track_orders == 1)',V1_z_log_odds{1}(track_orders == 1));
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds(track_orders == 1)') max(z_log_odds(track_orders == 1)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            mdl = fitlm(z_log_odds(track_orders == 2)',V1_z_log_odds{1}(track_orders == 2));
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds(track_orders == 2)') max(z_log_odds(track_orders == 2)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);


            if pval <= 0.05
                plot(x,y_est,'b:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','blue');
            else
                plot(x,y_est,'k:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            ylabel('V1 Left Log odds (z)')
            xlabel('HPC log odds (z)')
            yline(0,'--')
            xline(0,'--')
        end
        set(gca,"TickDir","out",'box', 'off','Color','none')

        subplot(2,4,8)
        if ~isempty(V1_z_log_odds{2})
            scatter(z_log_odds(track_orders == 1), V1_z_log_odds{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.2)
            hold on
            scatter(z_log_odds(track_orders == 2), V1_z_log_odds{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.2)

            mdl = fitlm(z_log_odds(track_orders == 1)',V1_z_log_odds{2}(track_orders == 1));
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds(track_orders == 1)') max(z_log_odds(track_orders == 1)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'r:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','red');
            else
                plot(x,y_est,'k:')
                text(gca,0.8,0.5,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            mdl = fitlm(z_log_odds(track_orders == 2)',V1_z_log_odds{2}(track_orders == 2));
            [pval,~,~] = coefTest(mdl);
            x =[min(z_log_odds(track_orders == 2)') max(z_log_odds(track_orders == 2)')];
            b = mdl.Coefficients.Estimate';
            y_est = polyval(fliplr(b),x);

            if pval <= 0.05
                plot(x,y_est,'b:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','blue');
            else
                plot(x,y_est,'k:')
                text(gca,0.3,0.2,['p = ' num2str(pval,2)],'Units','Normalized','FontName','Arial','Color','black');
            end

            ylabel('V1 Right Log odds (z)')
            xlabel('HPC log odds (z)')
            yline(0,'--')
            xline(0,'--')
            %             track_orders = [ones(1,length(lap_times(1).start))  2*ones(1,length(lap_times(2).start))];
            %             scatter3(z_log_odds(track_orders == 1), V1_z_log_odds{1}(track_orders == 1),V1_z_log_odds{2}(track_orders == 1),'r','filled','MarkerFaceAlpha',0.5)
            %             hold on
            %             scatter3(z_log_odds(track_orders == 2), V1_z_log_odds{1}(track_orders == 2),V1_z_log_odds{2}(track_orders == 2),'b','filled','MarkerFaceAlpha',0.5)
            %             xlabel('HPC log odds')
            %             ylabel('left V1 log odds')
            %             zlabel('right V1 log odds')
            %             set(gca,"TickDir","out",'box', 'off','Color','none')
            %             legend([h(1:3)],{'V1 Left','V1 Right','HPC'})

        end
        set(gca,"TickDir","out",'box', 'off','Color','none')
        fontsize(gcf,14,"points")


        for nprobe = 1:length(session_info(n).probe)
            options = session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no;
            probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere

            data = V1_bayesian_bias{probe_hemisphere};
            

            FPR = [];
            TPR = [];
            AUC = [];
            data_resampled = [];
            track_label_resampled = [];

            for nboot = 1:1000
                s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

                index = datasample(s,1:length(data),length(data));
                data_resampled(nboot,:) = data(index);
                track_label_resampled(nboot,:) = track_orders(index)-1;
                [X,Y,T,A] = perfcurve(track_orders(index),data(index),1,'XVals',0:0.05:1,'NBoot',1);
                %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

                FPR = X;
                TPR(nboot,:) = Y(:,1);
                AUC(nboot) = A(1);
            end

            if session_info(n).probe(nprobe).probe_hemisphere == 1
                fig(11) = figure(11);
                fig(11).Position = [500 100 1200 900];
                fig(11).Name = 'lap Bayesian Bias ROC two track discrimination in V1 for left probe';
                % ROC#
                subplot(2,5,nsession)
                hold on
                x = FPR';
                CI_shuffle = prctile(TPR,[2.5 97.5]);
                plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
                plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
                x2 = [x, fliplr(x)];
                inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
                h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

                h(1) = plot([0 1],[0 1],'k--')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                legend([h(2) h(1)],{'Real','chance'})
                title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
                sgtitle('lap Bayesian Bias ROC two track discrimination in V1 for left probe')
                        fontsize(gcf,14,"points")
            elseif session_info(n).probe(nprobe).probe_hemisphere == 2
                fig(12) = figure(12);
                fig(12).Position = [500 100 1200 900];
                fig(12).Name = 'lap Bayesian Bias ROC two track discrimination in V1 for right probe';
                % ROC#
                subplot(2,5,nsession)
                hold on
                x = FPR';
                CI_shuffle = prctile(TPR,[2.5 97.5]);
                plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
                plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
                x2 = [x, fliplr(x)];
                inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
                h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

                h(1) = plot([0 1],[0 1],'k--')
                set(gca,"TickDir","out",'box', 'off','Color','none')
                legend([h(2) h(1)],{'Real','chance'})
                title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
                sgtitle('lap Bayesian Bias ROC two track discrimination in V1 for right probe')
                        fontsize(gcf,14,"points")
            end
        end

        data = HPC_bayesian_bias;

        FPR = [];
        TPR = [];
        AUC = [];
        data_resampled = [];
        track_label_resampled = [];

        for nboot = 1:1000
            s = RandStream('mrg32k3a','Seed',nboot); % Set random seed for resampling

            index = datasample(s,1:length(data),length(data));
            data_resampled(nboot,:) = data(index);
            track_label_resampled(nboot,:) = track_orders(index)-1;
            [X,Y,T,A] = perfcurve(track_orders(index),data(index),1,'XVals',0:0.05:1,'NBoot',1);
            %             [X,Y,T,A] = perfcurve(track_label(index),z_log_odds(index),1,'XVals',0:0.05:1);

            FPR = X;
            TPR(nboot,:) = Y(:,1);
            AUC(nboot) = A(1);
        end

        fig(13) = figure(13);
        fig(13).Position = [500 100 1200 900];
        fig(13).Name = 'lap Bayesian Bias ROC two track discrimination in HPC';

        % ROC#
        subplot(2,5,nsession)
        hold on
        x = FPR';
        CI_shuffle = prctile(TPR,[2.5 97.5]);
        plot(x, CI_shuffle(2,:), 'r--', 'LineWidth', 1);
        plot(x, CI_shuffle(1,:), 'r--', 'LineWidth', 1);
        x2 = [x, fliplr(x)];
        inBetween = [CI_shuffle(1,:), fliplr(CI_shuffle(2,:))];
        h(2) = fill(x2, inBetween, 'r','FaceAlpha',0.2);

        h(1) = plot([0 1],[0 1],'k--')
        set(gca,"TickDir","out",'box', 'off','Color','none')
        legend([h(2) h(1)],{'Real','chance'})
        title(sprintf('Session %i AUC %.2f',nsession,mean(AUC)))
        sgtitle('lap Bayesian Bias ROC two track discrimination in HPC')
        fontsize(gcf,14,"points")

    end
    toc
end

cd('Z:\ibn-vision\USERS\Masa\V1_HPC_project')
save_all_figures('Z:\ibn-vision\USERS\Masa\V1_HPC_project',[])

% 
% %% Firing rate bias
% clear all
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
% 
% 
% SUBJECTS = {'M23017','M23028','M23029'};
% experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% % Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% % Stimulus_type = 'POST'; % extract LFP during RUN
% ROOTPATH = 'Z:\ibn-vision';
% c = 1;
% x_bins_width = 10;
% 
% for nsession =1:length(experiment_info)
%     session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
% 
%     if isempty(session_info)
%         continue
%     end
% 
%     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
% 
%         cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
%         for nprobe = 1:length(session_info(n).probe)
%             options = session_info(n).probe(nprobe);
%             probe_no = options.probe_id + 1;
%             probe_hemisphere = options.probe_hemisphere;
% 
% 
%             %         load best_channels
%             %         load extracted_PSD
%             load extracted_laps
%             load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             column = 1;
%             load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             load(sprintf('decoded_replay_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
% 
%             load(sprintf('probability_ratio_original%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             load(sprintf('probability_ratio_global_remapped%s.mat',erase(stimulus_name{n},'Masa2tracks')))
% 
%             load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             load(sprintf('extracted_HPC_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')))
%             load('extracted_V1_place_fields.mat')
%             load('extracted_CA1_place_fields.mat')
%             load('extracted_HPC_place_fields.mat')
%             temp.probe(probe_no)  = calculate_FR_immobility_masa_NPX(x_bins_width,position,HPC_clusters.probe(probe_no),[]);
%             
%             immobility_FR.probe(probe_no)
%         end
%     end
% end

