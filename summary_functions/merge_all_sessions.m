%%merge all sessions
addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))

mouse = {'M23031','M23032','M23034','M23037','M23038'};
date = {
    ['20230711';'20230712';'20230713';'20230714'];
    ['20230718';'20230719';'20230720';'20230721';'20230722'];
    ['20230804';'20230805';'20230806';'20230807'];
    ['20230810';'20230811';'20230812';'20230813'];
    ['20230816';'20230817']
    };
session_count =0;
    
clusters_all = struct;
clusters_all.region = [];
clusters_all.cluster_id = [];
clusters_all.cluster_spike_id = [];
clusters_all.spike_times = [];
clusters_all.raw = [];
clusters_all.session_count = [];
clusters_all.skaggs_info = [];
clusters_all.spatial_xcorr = [];
clusters_all.dwell_map = {};
clusters_all.within_track_corr = {};
clusters_all.within_track_xcorr = {};
clusters_all.first_second_corr = [];
clusters_all.odd_even_corr = [];
clusters_all.first_second_corr_shuffled = [];
clusters_all.first_second_stability = [];
clusters_all.odd_even_corr_shuffled = [];
clusters_all.odd_even_stability = [];
clusters_all.raw_peak = [];
clusters_all.peak_percentile = [];
clusters_all.skaggs_percentile = [];
clusters_all.across_tracks_correlation = {};
clusters_all.t1_t2_corr = [];
clusters_all.t1_t2_corr_shuffled = [];
clusters_all.t1_t2_remapping = [];
clusters_all.explained_variance = [];
clusters_all.session_label = [];
clusters_all.peak_channel = [];
clusters_all.peak_depth = [];
clusters_all.peak_channel_waveforms = [];
clusters_all.cell_type = [];


for iMouse = 1:5
    for iDate = 1:size(date{iMouse,1},1)
        session_count = session_count+1
        base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
        session_folder = fullfile(base_folder,mouse{iMouse},'analysis',date{iMouse}(iDate,:));

        DIR = dir(session_folder);
        stimulus_list = {DIR(3:6).name};

        stimulus = 'Track';

        %
        load(fullfile(session_folder,stimulus,'extracted_behaviour.mat'))
        load(fullfile(session_folder,stimulus,'extracted_peripherals.mat'))

        load(fullfile(session_folder,stimulus,'merged_clusters.mat'))
        load(fullfile(session_folder,stimulus,'session_info.mat'))
        load(fullfile(session_folder,'best_channels.mat'))
        load(fullfile(session_folder,stimulus,'extracted_place_fields.mat'))
        load(fullfile(session_folder,stimulus,'extracted_task_info.mat'))
        base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
        ANALYSIS_DATAPATH = session_info(1).probe(1).ANALYSIS_DATAPATH;
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH, 'Z:\ibn-vision\DATA\SUBJECTS', base_folder);
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH,'\', '/');

        if size(session_info.probe,2) == 1
            session_info.probe(2).SUBJECT = [];

        end
        if isempty(session_info.probe(2).SUBJECT)
            session_info.probe(1).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
            options = session_info.probe(1);

            V1_channels = determine_region_channels(best_channels{1},options,'region','V1','group','by probe');
            HPC_channels = determine_region_channels(best_channels{1},options,'region','HPC','group','by probe');
            %add 10000 to cluster ids so V1 probe always has 10000 added

            %combine the merged clusters
            clusters_combined = merged_clusters;
            clusters_combined = select_clusters(clusters_combined,metric_param);
            clusters_combined.region = strings(length(clusters_combined.cluster_id),1);
            clusters_combined.region(:) = 'n.a';
            clusters_combined.region(find(ismember(clusters_combined.peak_channel,V1_channels))) = 'V1';
            clusters_combined.region(find(ismember(clusters_combined.peak_channel,HPC_channels))) = 'HPC';
            [unique_cluster_ids,first_index] = unique(clusters_combined.merged_cluster_id);
            single_cluster_spike_id = cell(length(unique_cluster_ids),1);
            single_cluster_spike_times = cell(length(unique_cluster_ids),1);
            clear metric_param
            metric_param.merged_cluster_id = @(x) ismember(find(x), first_index);
            clusters_combined = select_clusters(clusters_combined,metric_param);
            clusters_combined.merged_spike_id = clusters_combined.merged_spike_id + iDate* 10^6 + str2double(mouse{iMouse}(2:end))*10^8 + 10^5;
            clusters_combined.merged_cluster_id = clusters_combined.merged_cluster_id + iDate* 10^6 + str2double(mouse{iMouse}(2:end))*10^8 + 10^5;
            [unique_cluster_ids,first_index] = unique(clusters_combined.merged_cluster_id);
            for ncell =1:length(unique_cluster_ids)
                single_cluster_spike_id{ncell} = clusters_combined.merged_spike_id(clusters_combined.merged_spike_id == unique_cluster_ids(ncell));
                single_cluster_spike_times{ncell} = clusters_combined.spike_times(clusters_combined.merged_spike_id == unique_cluster_ids(ncell));
            end
            
            clusters_all.region = [clusters_all.region;clusters_combined.region];
            
            clusters_all.cluster_id = [clusters_all.cluster_id;clusters_combined.merged_cluster_id];
            clusters_all.cluster_spike_id = [clusters_all.cluster_spike_id;single_cluster_spike_id];
            clusters_all.spike_times = [clusters_all.spike_times;single_cluster_spike_times];
            cluster_fields = fieldnames(clusters_combined);

            task_info_fields = fieldnames(Task_info);
            for iF = 1:length(task_info_fields)
                if session_count == 1
                    clusters_all.(task_info_fields{iF}) = {};
                end
                clusters_all.(task_info_fields{iF}){session_count} = Task_info.(task_info_fields{iF});
            end
            behaviour_fields = fieldnames(Behaviour);
            for iF = 1:length(behaviour_fields)
                if session_count == 1
                    clusters_all.(behaviour_fields{iF}) = {};
                end
                if ~(contains(behaviour_fields{iF},'eye_coordinates'))
                    clusters_all.(behaviour_fields{iF}){session_count} = Behaviour.(behaviour_fields{iF});
                    elseif contains(behaviour_fields{iF},'face_motion_SVD')
                        clusters_all.eye_coordinates{session_count} = Behaviour.face_motion_SVD(:,10);
                end

            end
            if size(place_fields,2) >1
                clusters_all.raw = [ clusters_all.raw;[place_fields(1).raw,place_fields(2).raw]];
                clusters_all.session_count = [clusters_all.session_count;repmat(session_count,[size(unique_cluster_ids)])];
                clusters_all.skaggs_info = [clusters_all.skaggs_info; [place_fields(1).skaggs_info,place_fields(2).skaggs_info]];
                clusters_all.spatial_xcorr = [clusters_all.spatial_xcorr;[place_fields(1).spatial_xcorr,place_fields(2).spatial_xcorr]];
                clusters_all.dwell_map{session_count,1} = place_fields(1).dwell_map;
                clusters_all.dwell_map{session_count,2} = place_fields(2).dwell_map;
                clusters_all.within_track_corr{session_count,1} = place_fields(1).within_track_corr;
                clusters_all.within_track_corr{session_count,2} = place_fields(2).within_track_corr;
                clusters_all.within_track_xcorr{session_count,1} = place_fields(1).within_track_xcorr;
                clusters_all.within_track_xcorr{session_count,2} = place_fields(2).within_track_xcorr;
                clusters_all.first_second_corr = [clusters_all.first_second_corr;[place_fields(1).first_second_corr',place_fields(2).first_second_corr']];
                clusters_all.odd_even_corr = [clusters_all.odd_even_corr;[place_fields(1).odd_even_corr',place_fields(2).odd_even_corr']];
                first_second_corr_shuffled_tmp = nan([size(place_fields(1).first_second_corr_shuffled),2]);
                first_second_corr_shuffled_tmp(:,:,1) = place_fields(1).first_second_corr_shuffled;
                first_second_corr_shuffled_tmp(:,:,2) = place_fields(2).first_second_corr_shuffled;
                clusters_all.first_second_corr_shuffled = [clusters_all.first_second_corr_shuffled;first_second_corr_shuffled_tmp];
                clusters_all.first_second_stability = [clusters_all.first_second_stability;[place_fields(1).first_second_stability',place_fields(2).first_second_stability']];
                odd_even_corr_shuffled_tmp = nan([size(place_fields(1).odd_even_corr_shuffled),2]);
                odd_even_corr_shuffled_tmp(:,:,1) = place_fields(1).odd_even_corr_shuffled;
                odd_even_corr_shuffled_tmp(:,:,2) = place_fields(2).odd_even_corr_shuffled;
                clusters_all.odd_even_corr_shuffled = [clusters_all.odd_even_corr_shuffled;odd_even_corr_shuffled_tmp];
                clusters_all.odd_even_stability = [clusters_all.odd_even_stability;[place_fields(1).odd_even_stability',place_fields(2).odd_even_stability']];
                clusters_all.raw_peak = [clusters_all.raw_peak;[place_fields(1).raw_peak',place_fields(2).raw_peak']];
                clusters_all.peak_percentile = [clusters_all.peak_percentile;[place_fields(1).peak_percentile',place_fields(2).peak_percentile']];
                clusters_all.skaggs_percentile = [clusters_all.skaggs_percentile;[place_fields(1).skaggs_percentile',place_fields(2).skaggs_percentile']];
                clusters_all.across_tracks_correlation{session_count,1} = place_fields(1).across_tracks_correlation;
                clusters_all.across_tracks_correlation{session_count,2} = place_fields(2).across_tracks_correlation;
                clusters_all.t1_t2_corr = [clusters_all.t1_t2_corr;[place_fields(1).t1_t2_corr',place_fields(2).t1_t2_corr']];
                t1_t2_corr_shuffled_tmp = nan([size(place_fields(1).t1_t2_corr_shuffled),2]);
                t1_t2_corr_shuffled_tmp(:,:,1) = place_fields(1).t1_t2_corr_shuffled;
                t1_t2_corr_shuffled_tmp(:,:,2) = place_fields(2).t1_t2_corr_shuffled;
                clusters_all.t1_t2_corr_shuffled = [clusters_all.t1_t2_corr_shuffled;t1_t2_corr_shuffled_tmp];
                clusters_all.t1_t2_remapping = [clusters_all.t1_t2_remapping;[place_fields(1).t1_t2_remapping',place_fields(2).t1_t2_remapping']];
                explained_variance_tmp = nan([size(place_fields(1).explained_variance),2]);
                explained_variance_tmp(:,:,1) = place_fields(1).explained_variance;
                explained_variance_tmp(:,:,2) = place_fields(2).explained_variance;
                clusters_all.explained_variance = [clusters_all.explained_variance;explained_variance_tmp];
            else
                % first session - preset the fields
                clusters_all.raw = [place_fields(1).raw,cell(size(place_fields(1).raw))];
                clusters_all.session_count = [repmat(session_count,[size(unique_cluster_ids)])];
                clusters_all.skaggs_info = [place_fields(1).skaggs_info,nan(size(place_fields(1).skaggs_info))];
                clusters_all.spatial_xcorr = [place_fields(1).spatial_xcorr,cell(size(place_fields(1).spatial_xcorr))];
                clusters_all.dwell_map{session_count,1} = place_fields(1).dwell_map;
                clusters_all.dwell_map{session_count,2} = [];
                clusters_all.within_track_corr{session_count,1} = place_fields(1).within_track_corr;
                clusters_all.within_track_corr{session_count,2} = [];
                clusters_all.within_track_xcorr{session_count,1} = place_fields(1).within_track_xcorr;
                clusters_all.within_track_xcorr{session_count,2} = [];
                clusters_all.first_second_corr = [place_fields(1).first_second_corr',nan(size(place_fields(1).first_second_corr'))];
                clusters_all.odd_even_corr = [place_fields(1).odd_even_corr',nan(size(place_fields(1).odd_even_corr'))];
                first_second_corr_shuffled_tmp = nan([size(place_fields(1).first_second_corr_shuffled),2]);
                first_second_corr_shuffled_tmp(:,:,1) = place_fields(1).first_second_corr_shuffled;
                clusters_all.first_second_corr_shuffled =[first_second_corr_shuffled_tmp];
                clusters_all.first_second_stability = [place_fields(1).first_second_stability',nan(size(place_fields(1).first_second_stability'))];
                odd_even_corr_shuffled_tmp = nan([size(place_fields(1).odd_even_corr_shuffled),2]);
                odd_even_corr_shuffled_tmp(:,:,1) = place_fields(1).odd_even_corr_shuffled;
                clusters_all.odd_even_corr_shuffled = [odd_even_corr_shuffled_tmp];
                clusters_all.odd_even_stability = [place_fields(1).odd_even_stability',nan(size(place_fields(1).odd_even_stability'))];
                clusters_all.raw_peak = [place_fields(1).raw_peak',nan(size(place_fields(1).raw_peak'))];
                clusters_all.peak_percentile = [place_fields(1).peak_percentile',nan(size(place_fields(1).peak_percentile'))];
                clusters_all.skaggs_percentile = [place_fields(1).skaggs_percentile',nan(size(place_fields(1).skaggs_percentile'))];
                clusters_all.across_tracks_correlation{session_count,1} = place_fields(1).across_tracks_correlation;
                clusters_all.across_tracks_correlation{session_count,2} = [];
                clusters_all.t1_t2_corr = nan(length(unique_cluster_ids),2);
                t1_t2_corr_shuffled_tmp = nan([length(unique_cluster_ids),1000,2]);
                clusters_all.t1_t2_corr_shuffled = [t1_t2_corr_shuffled_tmp];
                clusters_all.t1_t2_remapping = nan(length(unique_cluster_ids),2);
                explained_variance_tmp = nan([size(place_fields(1).explained_variance),2]);
                explained_variance_tmp(:,:,1) = place_fields(1).explained_variance;
                clusters_all.explained_variance = [explained_variance_tmp];
            end


        
        for iF = 2:28
            if session_count ==1
                clusters_all.(cluster_fields{iF}) = clusters_combined.(cluster_fields{iF});
            else
                clusters_all.(cluster_fields{iF}) = [clusters_all.(cluster_fields{iF});clusters_combined.(cluster_fields{iF})];
            end
        end
        clusters_all.peak_channel = [clusters_all.peak_channel;clusters_combined.peak_channel];
        clusters_all.peak_depth = [clusters_all.peak_depth;clusters_combined.peak_depth];
        clusters_all.peak_channel_waveforms = [clusters_all.peak_channel_waveforms;clusters_combined.peak_channel_waveforms];
        clusters_all.cell_type = [clusters_all.cell_type;clusters_combined.cell_type];
        clusters_all.session_label = [clusters_all.session_label;repmat([mouse{iMouse},date{iMouse}(iDate,:)],[size(unique_cluster_ids)])];

            

        else
            session_info.probe(1).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
            session_info.probe(2).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
            options = session_info.probe(2);
            V1_channels = determine_region_channels(best_channels{2},options,'region','V1','group','by probe');
            HPC_channels = determine_region_channels(best_channels{2},options,'region','HPC','group','by probe');
            merged_clusters(2).region = strings(length(merged_clusters(2).cluster_id),1);
            merged_clusters(2).region(:) = 'n.a';
            merged_clusters(2).region(find(ismember(merged_clusters(2).peak_channel,V1_channels))) = 'V1';
            merged_clusters(2).region(find(ismember(merged_clusters(2).peak_channel,HPC_channels))) = 'HPC';
            options = session_info.probe(1);
            MEC_channels = determine_region_channels(best_channels{1},options,'region','MEC_entry','group','by probe');
            HVA_channels = determine_region_channels(best_channels{1},options,'region','HVA','group','by probe');
            merged_clusters(1).region = strings(length(merged_clusters(1).cluster_id),1);
            merged_clusters(1).region(:) = 'n.a';
            merged_clusters(1).region(find(ismember(merged_clusters(1).peak_channel,MEC_channels))) = 'MEC';
            merged_clusters(1).region(find(ismember(merged_clusters(1).peak_channel,HVA_channels))) = 'HVA';
            %combine the merged clusters
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
            [unique_cluster_ids,first_index] = unique(clusters_combined.merged_cluster_id);
            single_cluster_spike_id = cell(length(unique_cluster_ids),1);
            single_cluster_spike_times = cell(length(unique_cluster_ids),1);
            clear metric_param
            metric_param.merged_cluster_id = @(x) ismember(find(x), first_index);
            clusters_combined = select_clusters(clusters_combined,metric_param);
            clusters_combined.merged_spike_id = clusters_combined.merged_spike_id + iDate* 10^6 + str2double(mouse{iMouse}(2:end))*10^8;
            clusters_combined.merged_cluster_id = clusters_combined.merged_cluster_id + iDate* 10^6 + str2double(mouse{iMouse}(2:end))*10^8;
            [unique_cluster_ids,first_index] = unique(clusters_combined.merged_cluster_id);
            for ncell =1:length(unique_cluster_ids)
                single_cluster_spike_id{ncell} = clusters_combined.merged_spike_id(clusters_combined.merged_spike_id == unique_cluster_ids(ncell));
                single_cluster_spike_times{ncell} = clusters_combined.spike_times(clusters_combined.merged_spike_id == unique_cluster_ids(ncell));
            end
            clusters_all.region = [clusters_all.region;clusters_combined.region];
            

            
            clusters_all.cluster_id = [clusters_all.cluster_id;clusters_combined.merged_cluster_id];
            clusters_all.cluster_spike_id = [clusters_all.cluster_spike_id;single_cluster_spike_id];
            clusters_all.spike_times = [clusters_all.spike_times;single_cluster_spike_times];
            cluster_fields = fieldnames(clusters_combined);

             task_info_fields = fieldnames(Task_info);
            for iF = 1:length(task_info_fields)
                if session_count == 1
                    clusters_all.(task_info_fields{iF}) = {};
                end
                clusters_all.(task_info_fields{iF}){session_count} = Task_info.(task_info_fields{iF});
            end

            behaviour_fields = fieldnames(Behaviour);
            for iF = 1:length(behaviour_fields)
                if ~(contains(behaviour_fields{iF},'eye_coordinates'))
                clusters_all.(behaviour_fields{iF}){session_count} = Behaviour.(behaviour_fields{iF});
                elseif contains(behaviour_fields{iF},'face_motion_SVD')
                    clusters_all.eye_coordinates{session_count} = Behaviour.face_motion_SVD(:,10);
                end
            end
            
            clusters_all.raw = [ clusters_all.raw;[place_fields(1).raw,place_fields(2).raw]];
            clusters_all.session_count = [clusters_all.session_count;repmat(session_count,[size(unique_cluster_ids)])];
            clusters_all.skaggs_info = [clusters_all.skaggs_info; [place_fields(1).skaggs_info,place_fields(2).skaggs_info]];
            clusters_all.spatial_xcorr = [clusters_all.spatial_xcorr;[place_fields(1).spatial_xcorr,place_fields(2).spatial_xcorr]];
            clusters_all.dwell_map{session_count,1} = place_fields(1).dwell_map;
            clusters_all.dwell_map{session_count,2} = place_fields(2).dwell_map;
            clusters_all.within_track_corr{session_count,1} = place_fields(1).within_track_corr;
            clusters_all.within_track_corr{session_count,2} = place_fields(2).within_track_corr;
            clusters_all.within_track_xcorr{session_count,1} = place_fields(1).within_track_xcorr;
            clusters_all.within_track_xcorr{session_count,2} = place_fields(2).within_track_xcorr;
            clusters_all.first_second_corr = [clusters_all.first_second_corr;[place_fields(1).first_second_corr',place_fields(2).first_second_corr']];
            clusters_all.odd_even_corr = [clusters_all.odd_even_corr;[place_fields(1).odd_even_corr',place_fields(2).odd_even_corr']];
            first_second_corr_shuffled_tmp = nan([size(place_fields(1).first_second_corr_shuffled),2]);
            first_second_corr_shuffled_tmp(:,:,1) = place_fields(1).first_second_corr_shuffled;
            first_second_corr_shuffled_tmp(:,:,2) = place_fields(2).first_second_corr_shuffled;
            clusters_all.first_second_corr_shuffled = [clusters_all.first_second_corr_shuffled;first_second_corr_shuffled_tmp];
            clusters_all.first_second_stability = [clusters_all.first_second_stability;[place_fields(1).first_second_stability',place_fields(2).first_second_stability']];
            odd_even_corr_shuffled_tmp = nan([size(place_fields(1).odd_even_corr_shuffled),2]);
            odd_even_corr_shuffled_tmp(:,:,1) = place_fields(1).odd_even_corr_shuffled;
            odd_even_corr_shuffled_tmp(:,:,2) = place_fields(2).odd_even_corr_shuffled;
            clusters_all.odd_even_corr_shuffled = [clusters_all.odd_even_corr_shuffled;odd_even_corr_shuffled_tmp];
            clusters_all.odd_even_stability = [clusters_all.odd_even_stability;[place_fields(1).odd_even_stability',place_fields(2).odd_even_stability']];
            clusters_all.raw_peak = [clusters_all.raw_peak;[place_fields(1).raw_peak',place_fields(2).raw_peak']];
            clusters_all.peak_percentile = [clusters_all.peak_percentile;[place_fields(1).peak_percentile',place_fields(2).peak_percentile']];
            clusters_all.skaggs_percentile = [clusters_all.skaggs_percentile;[place_fields(1).skaggs_percentile',place_fields(2).skaggs_percentile']];
            clusters_all.across_tracks_correlation{session_count,1} = place_fields(1).across_tracks_correlation;
            clusters_all.across_tracks_correlation{session_count,2} = place_fields(2).across_tracks_correlation;
            clusters_all.t1_t2_corr = [clusters_all.t1_t2_corr;[place_fields(1).t1_t2_corr',place_fields(2).t1_t2_corr']];
            t1_t2_corr_shuffled_tmp = nan([size(place_fields(1).t1_t2_corr_shuffled),2]);
            t1_t2_corr_shuffled_tmp(:,:,1) = place_fields(1).t1_t2_corr_shuffled;
            t1_t2_corr_shuffled_tmp(:,:,2) = place_fields(2).t1_t2_corr_shuffled;
            clusters_all.t1_t2_corr_shuffled = [clusters_all.t1_t2_corr_shuffled;t1_t2_corr_shuffled_tmp];
            clusters_all.t1_t2_remapping = [clusters_all.t1_t2_remapping;[place_fields(1).t1_t2_remapping',place_fields(2).t1_t2_remapping']];
            explained_variance_tmp = nan([size(place_fields(1).explained_variance),2]);
            explained_variance_tmp(:,:,1) = place_fields(1).explained_variance;
            explained_variance_tmp(:,:,2) = place_fields(2).explained_variance;
            clusters_all.explained_variance = [clusters_all.explained_variance;explained_variance_tmp];

        
        for iF = 2:28
            try
            clusters_all.(cluster_fields{iF}) = [clusters_all.(cluster_fields{iF});clusters_combined.(cluster_fields{iF})];
            catch
            end
        end
        clusters_all.peak_channel = [clusters_all.peak_channel;clusters_combined.peak_channel];
        clusters_all.peak_depth = [clusters_all.peak_depth;clusters_combined.peak_depth];
        clusters_all.peak_channel_waveforms = [clusters_all.peak_channel_waveforms;clusters_combined.peak_channel_waveforms];
        clusters_all.cell_type = [clusters_all.cell_type;clusters_combined.cell_type];
        clusters_all.session_label = [clusters_all.session_label;repmat([mouse{iMouse},date{iMouse}(iDate,:)],[size(unique_cluster_ids)])];
            %add temporal smoothed lap responses to the place fields


%{
             spatial_response = spatial_responses(clusters_combined,Behaviour,Task_info,'within',0);
            spatial_response_norm = spatial_responses(clusters_combined,Behaviour,Task_info,'within',1);
            spatial_response_extended = spatial_responses(clusters_combined,Behaviour,Task_info,'extension',0);
            spatial_response_extended_norm = spatial_responses(clusters_combined,Behaviour,Task_info,'extension',1);
            

            place_fields_all_t1{session_count}.spatial_response = spatial_response{1};
            place_fields_all_t2{session_count}.spatial_response = spatial_response{2};
            place_fields_all_t1{session_count}.spatial_response_extended = spatial_response_extended{1};
            place_fields_all_t2{session_count}.spatial_response_extended = spatial_response_extended{2};
            place_fields_all_t1{session_count}.spatial_response_norm = spatial_response_norm{1};
            place_fields_all_t2{session_count}.spatial_response_norm = spatial_response_norm{2};
            place_fields_all_t1{session_count}.spatial_response_extended_norm = spatial_response_extended_norm{1};
            place_fields_all_t2{session_count}.spatial_response_extended_norm = spatial_response_extended_norm{2}; 
%}





        end

    end
end




save(fullfile(base_folder,'clusters_all'),"clusters_all",'-v7.3')
%% pre-compute the spatial response
% for V1 and HVA, a delay factor is added to the spike times and it's 0.06s now

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
V1_index = (clusters_all.region == "V1");
HVA_index =  clusters_all.region =='HVA';
spatially_tuned_neurons = clusters_all.odd_even_stability >=0.95 ...
    & clusters_all.peak_percentile >=0.95;
overall_cluster_index = V1_index & spatially_tuned_neurons(:,1) | HVA_index;
delay = 0.08;
spatial_response = cell(sum(overall_cluster_index),2);
spatial_response_extended = cell(sum(overall_cluster_index),2);
cluster_counter = 0;
for iS = 1:19
    clusters_this_session = clusters_all.session_count == iS & overall_cluster_index;
    no_clusters_this_session = sum(clusters_this_session);
    
    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    spike_times_this_session = clusters_all.spike_times(clusters_this_session);
    parfor iC = 1:no_clusters_this_session
        spike_times = spike_times_this_session{iC};
        spatial_response(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within',delay);
        spatial_response_extended(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay);
    end
    cluster_counter = cluster_counter+no_clusters_this_session;
end
save(fullfile(base_folder, 'spatial_responses_V1_HVA.mat'),'spatial_response_extended','spatial_response','-v7.3')


% for V1 and HVA, a delay factor is added to the spike times and it's 0.06s now
%for MEC no delay is added
% select clusters you want to plot for population analysis e.g. only plot neurons from V1
MEC_index = (clusters_all.region == "MEC");

spatially_tuned_neurons = clusters_all.odd_even_stability >=0.95 ...
    & clusters_all.peak_percentile >=0.95;
overall_cluster_index = MEC_index & spatially_tuned_neurons(:,1) ;
delay = 0;
spatial_response = cell(sum(overall_cluster_index),2);
spatial_response_extended = cell(sum(overall_cluster_index),2);
cluster_counter = 0;
for iS = 1:19
    clusters_this_session = clusters_all.session_count == iS & overall_cluster_index;
    no_clusters_this_session = sum(clusters_this_session);
    
    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    spike_times_this_session = clusters_all.spike_times(clusters_this_session);
    parfor iC = 1:no_clusters_this_session
        spike_times = spike_times_this_session{iC};
        spatial_response(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within',delay);
        spatial_response_extended(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay);
    end
    cluster_counter = cluster_counter+no_clusters_this_session;

end
save(fullfile(base_folder, 'spatial_responses_MEC.mat'),'spatial_response_extended','spatial_response','-v7.3')

