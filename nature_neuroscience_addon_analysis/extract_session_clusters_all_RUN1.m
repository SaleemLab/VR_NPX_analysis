function extract_session_clusters_all_RUN1()
    addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
    addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))

    clear all
    SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
    option = 'bilateral';
    experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
    experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);

    % For RUN1, the stimulus type is Masa2tracks or similar
    Stimulus_type = 'Masa2tracks'; 

    session_clusters_all = struct();
    session_count = 0;

    for nsession =1:length(experiment_info)
        % First check if there is Masa2tracks or Sleep, we just need to load RUN1
        stimulus_name = experiment_info(nsession).StimulusName;
        
        session_info = experiment_info(nsession).session(1);
        options = session_info.probe(1);

        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        if isempty(DIR1)
            fprintf('No session_clusters_RUN1.mat found for session %d\n', nsession);
            continue
        end

        session_count = session_count + 1;
        
        SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info.probe(1).SUBJECT},option);
        iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info.probe(1).SESSION));

        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        session_clusters_RUN = session_clusters;
        clear session_clusters;

        % From cell structure back to spike times and spike id
        session_clusters_RUN.spike_id_all = vertcat(session_clusters_RUN.spike_id{:});
        session_clusters_RUN.spike_times_all = vertcat(session_clusters_RUN.spike_times{:});
        [session_clusters_RUN.spike_times_all,index] = sort(session_clusters_RUN.spike_times_all);
        session_clusters_RUN.spike_id_all = session_clusters_RUN.spike_id_all(index);

        spatial_cell_id = session_clusters_RUN.cluster_id((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
            | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));
        spatial_cell_index = find((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
            | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));
        
        spike_index = ismember(session_clusters_RUN.spike_id_all, spatial_cell_id);
        spike_times = session_clusters_RUN.spike_times_all(spike_index);
        spike_id = session_clusters_RUN.spike_id_all(spike_index);
        
        [spike_times, sort_idx] = sort(spike_times);
        spike_id = spike_id(sort_idx);

        mean_activity_track1 = nan(length(spatial_cell_id), 1);
        mean_activity_track2 = nan(length(spatial_cell_id), 1);

        for i = 1:length(spatial_cell_id)
            % Track 1
            data1 = session_clusters_RUN.spatial_response_raw{spatial_cell_index(i), 1};
            if ~isempty(data1)
                mean_activity_track1(i) = mean(data1(:), 'omitnan');
            end

            % Track 2
            data2 = session_clusters_RUN.spatial_response_raw{spatial_cell_index(i), 2};
            if ~isempty(data2)
                mean_activity_track2(i) = mean(data2(:), 'omitnan');
            end
        end

        session_clusters_all.subject(session_count,:) = options.SUBJECT;
        session_clusters_all.session(session_count,:) = options.SESSION;
        session_clusters_all.session_day(session_count) = iDate;

        session_clusters_all.spike_id{nsession} = spike_id;
        session_clusters_all.spike_times{nsession} = spike_times;
        session_clusters_all.spatial_cell_id{nsession} = spatial_cell_id;

        session_clusters_all.peak_channel{nsession} = session_clusters_RUN.peak_channel(spatial_cell_index);
        session_clusters_all.peak_depth{nsession} = session_clusters_RUN.peak_depth(spatial_cell_index);
        session_clusters_all.peak_channel_waveforms{nsession} = session_clusters_RUN.peak_channel_waveforms(spatial_cell_index,:);
        session_clusters_all.cell_type{nsession} = session_clusters_RUN.cell_type(spatial_cell_index);  % double
        session_clusters_all.receptive_field{nsession} = session_clusters_RUN.receptive_field(spatial_cell_index);  % cell
        session_clusters_all.region{nsession} = session_clusters_RUN.region(spatial_cell_index);  % string

        session_clusters_all.spatial_response{nsession} = session_clusters_RUN.spatial_response(spatial_cell_index,:);
        session_clusters_all.spatial_response_raw{nsession} = session_clusters_RUN.spatial_response_raw(spatial_cell_index,:);
        session_clusters_all.mean_FR{nsession} = [mean_activity_track1 mean_activity_track2];

        session_clusters_all.SMI_peak_RUN1{nsession} = session_clusters_RUN.SMI_peak_RUN1(spatial_cell_index,:);
        session_clusters_all.SMI_area_RUN1{nsession} = session_clusters_RUN.SMI_area_RUN1(spatial_cell_index,:);
        session_clusters_all.SMI_peak_RUN2{nsession} = session_clusters_RUN.SMI_peak_RUN2(spatial_cell_index,:);
        session_clusters_all.SMI_area_RUN2{nsession} = session_clusters_RUN.SMI_area_RUN2(spatial_cell_index,:);

        session_clusters_all.peak_percentile_RUN1{nsession} = session_clusters_RUN.peak_percentile_RUN1(spatial_cell_index,:);
        session_clusters_all.odd_even_stability_RUN1{nsession} = session_clusters_RUN.odd_even_stability_RUN1(spatial_cell_index,:);
        session_clusters_all.first_second_stability_RUN1{nsession} = session_clusters_RUN.first_second_stability_RUN1(spatial_cell_index,:);
        session_clusters_all.peak_percentile_RUN2{nsession} = session_clusters_RUN.peak_percentile_RUN2(spatial_cell_index,:);
        session_clusters_all.odd_even_stability_RUN2{nsession} = session_clusters_RUN.odd_even_stability_RUN2(spatial_cell_index,:);
        session_clusters_all.first_second_stability_RUN2{nsession} = session_clusters_RUN.first_second_stability_RUN2(spatial_cell_index,:);

        session_clusters_all.landmark_tuned_status{nsession} = session_clusters_RUN.landmark_tuned_status(spatial_cell_index,:);
        session_clusters_all.landmark_centres{nsession} = session_clusters_RUN.landmark_centres(spatial_cell_index,:);
        session_clusters_all.landmark_peaks_all{nsession} = session_clusters_RUN.landmark_peaks_all(spatial_cell_index,:);
        session_clusters_all.landmark_areas_all{nsession} = session_clusters_RUN.landmark_areas_all(spatial_cell_index,:);
        session_clusters_all.interneurons_all{nsession} = session_clusters_RUN.interneurons_all(spatial_cell_index,:);

        session_clusters_all.lap_session_ID{nsession} = session_clusters_RUN.lap_seesion_ID;
        session_clusters_all.pass_landmarks_laps{nsession} = session_clusters_RUN.pass_landmarks_laps;
        session_clusters_all.spatial_response_all_2cm{nsession} = session_clusters_RUN.spatial_response_all_2cm(spatial_cell_index,:,:);

        % Save behavior variables (time vector, track position, track id)
        if isfield(session_clusters_RUN, 'tvec')
            session_clusters_all.tvec{nsession} = session_clusters_RUN.tvec;
        end
        if isfield(session_clusters_RUN, 'position')
            session_clusters_all.position{nsession} = session_clusters_RUN.position;
        end
        if isfield(session_clusters_RUN, 'track_ID_all')
            session_clusters_all.track_ID_all{nsession} = session_clusters_RUN.track_ID_all;
        end
        if isfield(session_clusters_RUN, 'speed')
            session_clusters_all.speed{nsession} = session_clusters_RUN.speed;
        end
    end

    if exist('D:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'D:\corticohippocampal_replay';
    elseif exist('P:\corticohippocampal_replay', 'dir') > 0
        analysis_folder = 'P:\corticohippocampal_replay';
    end

    save(fullfile(analysis_folder,'session_clusters_all_RUN1.mat'),'session_clusters_all','-v7.3')
    fprintf('Saved session_clusters_all_RUN1.mat successfully!\n');
end
