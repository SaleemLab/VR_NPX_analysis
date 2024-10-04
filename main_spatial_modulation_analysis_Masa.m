%% Main place cell and V1 spatial modulation analysis code
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

%% setting metrics to screen good clusters
clear all
params = create_cluster_selection_params('sorting_option','masa');

%% organize clusters from each session
session_count = 0;
% SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS={'M24016','M24017','M24018'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([6 9 14 21 22 27 35 38 40]);
Stimulus_type = 'RUN1'; % has to be RUN1 or RUN2 
% Stimulus_type = 'RUN1'; % has to be RUN1 or RUN2 
base_folder='Z:\ibn-vision\DATA\SUBJECTS';
% for iSub = 1:length(SUBJECTS)
%     load(fullfile(base_folder,SUBJECTS{iSub},'analysis','experiment_info.mat'))
%     Stimulus_type = 'Track';
for nsession = 1:length(experiment_info)
    session_clusters = struct();
    iDate = nsession;
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    session_folder = fullfile(base_folder,session_info.probe(1).SUBJECT,'analysis',session_info.probe(1).SESSION);

    for n = 1:length(session_info) % should just be 1 stimulus type normally (RUN1)
        % change/update analysis folder depending on the pc used
        ANALYSIS_DATAPATH = session_info(n).probe(1).ANALYSIS_DATAPATH;
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH, 'Z:\ibn-vision\DATA\SUBJECTS', base_folder);
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH,'\', '/');
        session_info(n).probe(1).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
        session_info(n).probe(2).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;

        options = session_info(n).probe(1);
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        end

        clusters = clusters_ks3;
        sorting_option = 'masa';
        load(fullfile(options.ANALYSIS_DATAPATH,'..','best_channels.mat'));
        receptive_field_file =fullfile(session_folder,'SparseNoise_fullscreen','receptiveFields_ks3.mat');
        if exist(receptive_field_file,"file")
            load(receptive_field_file);
            receptive_field=RF;
        else
            receptive_field_file =fullfile(session_folder,'SparseNoise','receptiveFields_ks3.mat');
            if exist(receptive_field_file,"file")
                load(receptive_field_file);
                receptive_field=RF;
            else
                receptive_field = [];
            end
        end

        for nprobe = 1:length(clusters)
            if isempty(receptive_field)
                clusters(nprobe).receptive_field = cell(size(clusters(nprobe).cluster_id));
            else
                clusters(nprobe).receptive_field = receptive_field{nprobe};
            end
        end

        task_info_fields = fieldnames(Task_info);
        for iField =1:length(task_info_fields)
            session_clusters.(task_info_fields{iField})= {Task_info.(task_info_fields{iField})};
        end
        % add behaviour info for each session
        behaviour_fields = fieldnames(Behaviour);
        for iF = 1:length(behaviour_fields)
            if ~(contains(behaviour_fields{iF},'eye_coordinates'))
                session_clusters.(behaviour_fields{iF}) = {Behaviour.(behaviour_fields{iF})};
            elseif contains(behaviour_fields{iF},'face_motion_SVD')
                session_clusters.eye_coordinates = {Behaviour.face_motion_SVD(:,10)};
            end
        end

        if isfield(clusters,'sorter')
            clusters = rmfield(clusters,'sorter');
        end

        if isfield(clusters,'probe_id')
            clusters = rmfield(clusters,'probe_id');
        end

        for nprobe = 1:length(clusters)
            clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
            options = session_info.probe(nprobe);
            if clusters(nprobe).probe_hemisphere == 1
                clusters(nprobe).region(:) = 'n.a_L';
%                                 V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by shank');
%                                 HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by shank');

                kClusters=kmeans(clusters(nprobe).peak_depth,2);
                if mean(clusters(nprobe).peak_depth(kClusters==1))>mean(clusters(nprobe).peak_depth(kClusters==2))
                    % if mean ocation of cluster one is above cluster two, it is
                    % Cortex.
                    % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                    V1_cell_id = find(kClusters==1);
                    HPC_cell_id = find(kClusters==2);
                else
                    % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                    V1_cell_id = find(kClusters==2);
                    HPC_cell_id = find(kClusters==1);
                end
                clusters(nprobe).region(V1_cell_id) = 'V1_L';
                clusters(nprobe).region(HPC_cell_id) = 'HPC_L';
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1_L';
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HPC_channels))) = 'HPC_L';

            elseif clusters(nprobe).probe_hemisphere == 2
                clusters(nprobe).region(:) = 'n.a_R';
%                 V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by shank');
%                 HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by shank');
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,MEC_channels))) = 'V1_R';
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HVA_channels))) = 'HPC_R';


                kClusters=kmeans(clusters(nprobe).peak_depth,2);
                if mean(clusters(nprobe).peak_depth(kClusters==1))>mean(clusters(nprobe).peak_depth(kClusters==2))
                    % if mean ocation of cluster one is above cluster two, it is
                    % Cortex.
                    % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==1);
                    V1_cell_id = find(kClusters==1);
                    HPC_cell_id = find(kClusters==2);
                else
                    % V1_cell_id = find(RF.probe(nprobe).shank == nshank & kClusters==2);
                    V1_cell_id = find(kClusters==2);
                    HPC_cell_id = find(kClusters==1);
                end
                clusters(nprobe).region(V1_cell_id) = 'V1_R';
                clusters(nprobe).region(HPC_cell_id) = 'HPC_R';
            end
        end
        field_names = fieldnames(clusters);

        for nprobe = 1:length(clusters)
            options = session_info.probe(nprobe);
            clusters_probe = select_clusters(clusters(nprobe),params); %only look at good clusters
            single_cluster_spike_id = cell(size(clusters_probe.cluster_id));
            single_cluster_spike_times = cell(size(clusters_probe.cluster_id));
            unique_spike_id = clusters_probe.spike_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
            unique_cluster_id = clusters_probe.cluster_id + nprobe*10^4 + iDate* 10^6 + str2double(options.SUBJECT(2:end))*10^8;
            [unique_cluster_ids,first_index] = unique(unique_cluster_id);
            for ncell =1:length(unique_cluster_ids)
                single_cluster_spike_id{ncell,1} = unique_spike_id(unique_spike_id == unique_cluster_ids(ncell));
                single_cluster_spike_times{ncell,1} = clusters_probe.spike_times(unique_spike_id == unique_cluster_ids(ncell));
            end
            clusters_probe.cluster_id = unique_cluster_id(first_index);
            clusters_probe.spike_id = single_cluster_spike_id;
            clusters_probe.spike_times = single_cluster_spike_times;

            for iField =1:length(field_names)
                if nprobe == 1
                    session_clusters.(field_names{iField}) = [clusters_probe.(field_names{iField})];

                else
                    session_clusters.(field_names{iField}) = [session_clusters.(field_names{iField});clusters_probe.(field_names{iField})];

                end
            end
        end
        clear clusters clusters_ks3 clusters_probe
        tvec = Behaviour.tvec;
        track_ID =Behaviour.track_ID;
        position = Behaviour.position;
        speed = Behaviour.speed;
        track_ID_all = Task_info.track_ID_all;
        start_time_all = Task_info.start_time_all;
        end_time_all = Task_info.end_time_all;
        no_clusters_this_session = length(session_clusters.cluster_id);
        spatial_response = cell(no_clusters_this_session,2);
        delay = 0;
        within_track_corr = cell(no_clusters_this_session,2);
        across_track_corr = cell(no_clusters_this_session,2);
        first_second_stability = zeros(no_clusters_this_session,2);
        odd_even_stability = zeros(no_clusters_this_session,2);
        peak_percentile = zeros(no_clusters_this_session,2);
        reliablity = zeros(no_clusters_this_session,2);
        odd_even_corr_shuffled = cell(no_clusters_this_session,2);
        first_second_corr_shuffled = cell(no_clusters_this_session,2);


        clear Behaviour Task_info

        all_spike_times = vertcat(session_clusters.spike_times{:});
        cluster_id = session_clusters.cluster_id;
        all_spike_id = vertcat(session_clusters.spike_id{:});
        spatial_response = fast_spatial_response(all_spike_id,cluster_id,all_spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all);

        nbatch = 10;
        % Calculate the indices to split the clusters into batches
        clusters_index_batch = round(linspace(0, no_clusters_this_session, nbatch + 1));

        position_shuffled = cell(1000,1);   %shuffled positions for shuffles

        parfor nshuffle = 1:1000
            %     tic
            disp(sprintf('Shuffle %i',nshuffle))
            position_shuffled{nshuffle} = position;

            for nlap = 1:length(start_time_all)
                s = RandStream('mrg32k3a','Seed',nlap+nshuffle*1000); % Set random seed for resampling

                position_shuffled{nshuffle}(tvec  >= start_time_all(nlap) & tvec <= end_time_all(nlap)) =...
                    circshift(position(tvec  >= start_time_all(nlap) & tvec <= end_time_all(nlap))...
                    ,randi(s,length(tvec  >= start_time_all(nlap) & tvec <= end_time_all(nlap))));
            end

        end
        for iBatch = 1:nbatch
            disp(['batch: ',num2str(iBatch)])
            all_spike_times_batch = vertcat(session_clusters.spike_times{clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1)});
            cluster_id_batch = session_clusters.cluster_id(clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1));
            all_spike_id_batch = vertcat(session_clusters.spike_id{clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1)});
            spatial_response_shuffled = cell(clusters_index_batch(iBatch+1)-clusters_index_batch(iBatch),2,1000); %pre-set shuffled response profiles
            parfor nshuffle = 1: 1000
                disp(sprintf('Shuffle %i',nshuffle))
                spatial_response_shuffled(:,:,nshuffle) = fast_spatial_response(all_spike_id_batch,cluster_id_batch,all_spike_times_batch,tvec...
                    ,position_shuffled{nshuffle},speed,track_ID_all,start_time_all,end_time_all);
            end
            clear all_spike_times_batch cluster_id_batch all_spike_id_batch

            for iCluster = clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1)
                if max(track_ID)==2
                    across_track_corr{iCluster} = corr(normalize(spatial_response{iCluster,1}','range'),...
                        normalize(spatial_response{iCluster,2}','range')); % lap by lap correlation;
                end

                tic
                % Stability
                for track_id = 1:max(track_ID)
                    start_times = start_time_all(track_ID_all == track_id);

                    lap_correlation = corr(normalize(spatial_response{iCluster,track_id}','range'),...
                        normalize(spatial_response{iCluster,track_id}','range')); % lap by lap correlation

                    first_second_corr = mean(mean(lap_correlation(1:2:end,2:2:end),'omitnan'),'omitnan'); % First half vs Second half
                    odd_even_corr = mean(mean(lap_correlation(1:round(length(start_times)/2),round(length(start_times)/2)+1:end),'omitnan'),'omitnan'); % Odd laps vs even laps

                    within_track_corr{iCluster,track_id} = lap_correlation;

                    normalised_ratemap = normalize(spatial_response{iCluster,track_id}','range');

                    peak_shuffled = zeros(1000,1);
                    odd_even_corr_shuffled_temp = zeros(1000,1);
                    first_second_corr_shuffled_temp = zeros(1000,1);
                    sr_shuffled_temp = spatial_response_shuffled(iCluster-clusters_index_batch(iBatch),track_id,:);
                    parfor nshuffle = 1:1000
                        % Range between min and max FR
                        peak_shuffled(nshuffle,1) = max(mean(sr_shuffled_temp{1,1,nshuffle}))-min(mean(sr_shuffled_temp{1,1,nshuffle}));

                        lap_correlation_shuffled = corr(normalised_ratemap,...
                            normalize(sr_shuffled_temp{1,1,nshuffle}','range')); % original X shuffled correlation

                        odd_even_corr_shuffled_temp(nshuffle) = mean(mean(lap_correlation_shuffled(1:2:end,2:2:end),'omitnan'),'omitnan'); % odd lap from
                        first_second_corr_shuffled_temp(nshuffle) = mean(mean(lap_correlation_shuffled(1:round(length(start_times)/2),round(length(start_times)/2)+1:end),'omitnan'),'omitnan');
                    end
                    clear lap_correlation_shuffled
                    first_second_stability(iCluster,track_id) = sum(first_second_corr > first_second_corr_shuffled_temp)/length(first_second_corr_shuffled_temp);
                    first_second_corr_shuffled{iCluster,track_id} = first_second_corr_shuffled_temp;

                    odd_even_stability(iCluster,track_id) = sum(odd_even_corr > odd_even_corr_shuffled_temp)/length(odd_even_corr_shuffled_temp);
                    odd_even_corr_shuffled{iCluster,track_id} = odd_even_corr_shuffled_temp;



                    peak_percentile(iCluster,track_id) = sum(max(mean(spatial_response{iCluster,track_id}))-min(mean(spatial_response{iCluster,track_id}))>...
                        peak_shuffled)/length(peak_shuffled); % max-min FR relative to max-min of shuffled data


                end
                toc
            end
        end
        clear spatial_response_shuffled
        if max(track_ID)==2
            session_clusters.across_track_corr = across_track_corr;
        end
        session_clusters.within_track_corr = within_track_corr;
        session_clusters.first_second_stability = first_second_stability;
        session_clusters.odd_even_stability = odd_even_stability;
        session_clusters.peak_percentile = peak_percentile;
        session_clusters.reliablity = reliablity;
        session_clusters.odd_even_corr_shuffled = odd_even_corr_shuffled;
        session_clusters.first_second_corr_shuffled = first_second_corr_shuffled;
        session_clusters.spatial_response = spatial_response;
        
        % add task info for each session
        if contains(stimulus_name{n},'Masa2tracks')
            save(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'session_clusters');
        end
    end
end
%% Peri event raster plots

for nsession =[1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'PSD');

    plot_raster_both_track(clusters_combined.spike_times,clusters_combined.spike_id,Task_info,Behaviour,[5 1],[0 140],2,...
        'unit_depth',clusters_combined.peak_depth(ia),'unit_region',clusters_combined.region(ia),'unit_id',C);
end
%% Add ripple and replay response for each cluster


% end


%% merge all sessions into clusters_all

base_folder = '/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data';
SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};
option = 'V1-MEC';
session_count = 0;
clusters_all = struct();
for iSub = 1:5
    load(fullfile(base_folder,SUBJECTS{iSub},'analysis','experiment_info.mat'))
    Stimulus_type = 'Track';

    for nsession = 1:length(experiment_info)
        session_count = session_count + 1;

        iDate = nsession;
        session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        session_folder = fullfile(base_folder,SUBJECTS{iSub},'analysis',session_info.probe(1).SESSION);
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        load(fullfile(session_folder,stimulus_name{1},'extracted_behaviour.mat'));
        load(fullfile(session_folder,stimulus_name{1},'extracted_task_info.mat'));
        load(fullfile(session_folder,"session_clusters.mat"))
        session_clusters = rmfield(session_clusters,'timevec');
        session_clusters.session_label = repmat([SUBJECTS{iSub},session_info.probe(1).SESSION],[size(session_clusters.cluster_id)]);

        % put all clusters in one struct
        field_names = fieldnames(session_clusters);
        for iField =1:length(field_names)
            if session_count == 1
                clusters_all.(field_names{iField}) = session_clusters.(field_names{iField});
            else
                try
                    clusters_all.(field_names{iField}) = [clusters_all.(field_names{iField});session_clusters.(field_names{iField})];
                catch
                end
            end
        end

        
        session_count_clusters = repmat(session_count,size(session_clusters.cluster_id));
        if session_count == 1
            clusters_all.session_count = session_count_clusters;
        else
            clusters_all.session_count = [clusters_all.session_count;session_count_clusters];
        end
    end
end
% clusters_all.region(contains(clusters_all.region,'MEC_entry')) = 'MEC';



 %% Add/modify particular attributes to the clusters, e.g. modifying regions or adding layers
% 
% params.amplitude_cutoff = @(x) x<=0.01; %0.01 if strict and removes lots of units
% params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
% params.num_negative_peaks = @(x) x<=1;
% params.num_positive_peaks = @(x) x<=2;
% params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0002;
% params.amplitude_cv_median = @(x) x<=0.7;
% params.amplitude_cv_range = @(x) x<=0.7;
% params.firing_rate = @(x) x>= 0.05;
% params.presence_ratio = @(x) x>= 0.5;
% params.sd_ratio = @(x) x<=4;
% params.snr = @(x) x>= 2;
% params.half_width = @(x) x<= 0.0003;
% 
% 
% load(fullfile(base_folder,'clusters_all_ks3'))
% 
% SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};
% option = 'V1-MEC';
% session_count = 0;
% for iSub = 1:5
%     load(fullfile(base_folder,SUBJECTS{iSub},'analysis','experiment_info.mat'))
%     Stimulus_type = 'Track';
% 
%     for nsession = 1:length(experiment_info)
%         session_count = session_count + 1;
%         load(fullfile(session_folder,"session_clusters.mat"))
%         iDate = nsession;
%         session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%         session_folder = fullfile(base_folder,SUBJECTS{iSub},'analysis',session_info.probe(1).SESSION);
%         stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%         ANALYSIS_DATAPATH = session_info(1).probe(1).ANALYSIS_DATAPATH;
%         ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH, 'Z:\ibn-vision\DATA\SUBJECTS', base_folder);
%         ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH,'\', '/');
%         session_info.probe(1).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
%         session_info.probe(2).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
%         load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
%         clusters = clusters_ks3;
%         options = session_info.probe(1);
% 
% 
%         for nprobe = 1:length(clusters)
%             clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
%             clusters(nprobe).region(:) = 'n.a';
%             options = session_info.probe(nprobe);
%             if options.probe_id == options.probe_V1
%                 V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
%                 HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HPC_channels))) = 'HPC';
% 
%             elseif options.probe_id == options.probe_MEC
%                 %             V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
%                 MEC_channels = determine_region_channels(best_channels{nprobe},options,'region','MEC_entry','group','by probe');
%                 HVA_channels = determine_region_channels(best_channels{nprobe},options,'region','HVA','group','by probe');
%                 %             clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,MEC_channels))) = 'MEC';
%                 clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HVA_channels))) = 'HVA';
%             end
%         end
% 
% 
%         for nprobe = 1:length(clusters)
%             options = session_info.probe(nprobe);
%             clusters_probe = select_clusters(clusters(nprobe),params); %only look at good clusters
%             single_cluster_spike_id = cell(size(clusters_probe.cluster_id));
%             single_cluster_spike_times = cell(size(clusters_probe.cluster_id));
%             field_names = {'region'};
%             for iField =1:length(field_names)
%                 if nprobe == 1
%                     session_clusters.(field_names{iField}) = [clusters_probe.(field_names{iField})];
%                 else
%                     session_clusters.(field_names{iField}) = [session_clusters.(field_names{iField});clusters_probe.(field_names{iField})];
%                 end
%             end
%         end
% 
%         save(fullfile(session_folder,'session_clusters.mat'),'session_clusters')
%         if session_count == 1
%             clusters_all.(field_names{1}) = session_clusters.(field_names{1});
%         else
%             clusters_all.(field_names{1}) = [clusters_all.(field_names{1}); session_clusters.(field_names{1})];
%         end
%     end
% end
% save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')
%% decoding and reliability of each cluster

addpath(genpath('/Users/atom/Documents/GitHub/Saleem-Lab-Code/'))
clusters_all.reliability = nan(size(clusters_all.cluster_id,1),2);
V1_index = (clusters_all.region == "V1");
HVA_index =  clusters_all.region =='HVA';
overall_cluster_index = V1_index|HVA_index;
for iS = 1:19
    
    clusters_this_session = clusters_all.session_count == iS & overall_cluster_index;

    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    w = gausswin(11);
    w = w / sum(w);
    speed(isnan(speed)) = 0;
    speed= filtfilt(w,1,speed')';
    tvec_edges = [tvec(1)-1/(60*2) tvec+1/(60*2)];
    spike_times_this_session = clusters_all.spike_times(clusters_this_session);
    cluster_responses = zeros(length(spike_times_this_session),length(tvec));
    parfor iC = 1:length(spike_times_this_session)
        spike_times = spike_times_this_session{iC,1};
        spike_speed =  interp1(tvec,speed,spike_times,'nearest');
        spike_times_run = spike_times(spike_speed >= 1);
        cluster_responses(iC,:) = histcounts(spike_times_run,tvec_edges);
        % cluster_responses(iC,:) = filtfilt(w,1,cluster_responses(iC,:)')';
    end
    t1 = track_ID_all == 1;
    t1_start = start_time_all(t1);
    t1_end = end_time_all(t1);
    t1_time = logical(zeros(size(tvec)));
    for iT = 1:sum(t1)
        temp_time = tvec >= t1_start(iT) & tvec <= t1_end(iT);
        t1_time = t1_time | temp_time;
    end
    t1_pos = zeros(size(tvec));
    t1_pos(t1_time) = position(t1_time);
    t1_pos(isnan(t1_pos)) = 0;

    t2 = track_ID_all == 2;
    t2_start = start_time_all(t2);
    t2_end = end_time_all(t2);
    t2_time = logical(zeros(size(tvec)));
    for iT = 1:sum(t2)
        temp_time = tvec >= t2_start(iT) & tvec <= t2_end(iT);
        t2_time = t2_time | temp_time;
    end
    t2_pos = zeros(size(tvec));
    t2_pos(t2_time) = position(t2_time);
    t2_pos(isnan(t2_pos)) = 0;
    % Create an instance of bayesDecoder
    decoder_t1 = bayesDecoder;
    decoder_t1.numBins = 140;
    decoder_t1.fixedSmth = 1;
    decoder_t1.variable = [1 140];
    % Define your variable X and firing Y
    t = speed>= 1 & t1_time & t1_pos>=1;
    X = t1_pos(t)'; % replace with your data
    Y = cluster_responses(:,t)'; % replace with your data
    [decoder_t1, prediction, X, Posterior, nPosterior] = decoder_t1.trainDecoder(X, Y);

    mean_EV_t1 = mean(decoder_t1.model.EV);

    clusters_all.reliability(clusters_this_session,1) = mean_EV_t1';
    if iS > 1
        decoder_t2 = bayesDecoder;
        decoder_t2.numBins = 140;
        decoder_t2.fixedSmth = 1;
        decoder_t2.variable = [1 140];
        % Define your variable X and firing Y
        t = speed>= 1 & t2_time & t2_pos>=1;
        X = t2_pos(t)'; % replace with your data
        Y = cluster_responses(:,t)'; % replace with your data
        [decoder_t2, prediction, X, Posterior, nPosterior] = decoder_t2.trainDecoder(X, Y);

        mean_EV_t2 = mean(decoder_t2.model.EV);

        clusters_all.reliability(clusters_this_session,2) = mean_EV_t2';
    end
end

MEC_index = (clusters_all.region == "MEC");
overall_cluster_index = MEC_index;
for iS = 5:19
    % MEC_index = (clusters_all.region == "MEC");
    % overall_cluster_index = MEC_index ;


    clusters_this_session = clusters_all.session_count == iS & overall_cluster_index;

    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    w = gausswin(11);
    w = w / sum(w);
    speed(isnan(speed)) = 0;
    speed= filtfilt(w,1,speed')';
    tvec_edges = [tvec(1)-1/(60*2) tvec+1/(60*2)];
    spike_times_this_session = clusters_all.spike_times(clusters_this_session);
    cluster_responses = zeros(length(spike_times_this_session),length(tvec));
    parfor iC = 1:length(spike_times_this_session)
        spike_times = spike_times_this_session{iC,1};
        spike_speed =  interp1(tvec,speed,spike_times,'nearest');
        spike_times_run = spike_times(spike_speed >= 1);
        cluster_responses(iC,:) = histcounts(spike_times_run,tvec_edges);
        % cluster_responses(iC,:) = filtfilt(w,1,cluster_responses(iC,:)')';
    end
    t1 = track_ID_all == 1;
    t1_start = start_time_all(t1);
    t1_end = end_time_all(t1);
    t1_time = logical(zeros(size(tvec)));
    for iT = 1:sum(t1)
        temp_time = tvec >= t1_start(iT) & tvec <= t1_end(iT);
        t1_time = t1_time | temp_time;
    end
    t1_pos = zeros(size(tvec));
    t1_pos(t1_time) = position(t1_time);
    t1_pos(isnan(t1_pos)) = 0;

    t2 = track_ID_all == 2;
    t2_start = start_time_all(t2);
    t2_end = end_time_all(t2);
    t2_time = logical(zeros(size(tvec)));
    for iT = 1:sum(t2)
        temp_time = tvec >= t2_start(iT) & tvec <= t2_end(iT);
        t2_time = t2_time | temp_time;
    end
    t2_pos = zeros(size(tvec));
    t2_pos(t2_time) = position(t2_time);
    t2_pos(isnan(t2_pos)) = 0;
    % Create an instance of bayesDecoder
    decoder_t1 = bayesDecoder;
    decoder_t1.numBins = 140;
    decoder_t1.fixedSmth = 1;
    decoder_t1.variable = [1 140];
    % Define your variable X and firing Y
    t = speed>= 1 & t1_time & t1_pos>=1;
    X = t1_pos(t)'; % replace with your data
    Y = cluster_responses(:,t)'; % replace with your data
    [decoder_t1, prediction, X, Posterior, nPosterior] = decoder_t1.trainDecoder(X, Y);

    mean_EV_t1 = mean(decoder_t1.model.EV);

    clusters_all.reliability(clusters_this_session,1) = mean_EV_t1';

    decoder_t2 = bayesDecoder;
    decoder_t2.numBins = 140;
    decoder_t2.fixedSmth = 1;
    decoder_t2.variable = [1 140];
    % Define your variable X and firing Y
    t = speed>= 1 & t2_time & t2_pos>=1;
    X = t2_pos(t)'; % replace with your data
    Y = cluster_responses(:,t)'; % replace with your data
    [decoder_t2, prediction, X, Posterior, nPosterior] = decoder_t2.trainDecoder(X, Y);

    mean_EV_t2 = mean(decoder_t2.model.EV);

    clusters_all.reliability(clusters_this_session,2) = mean_EV_t2';
end
HPC_index = (clusters_all.region == "HPC");
overall_cluster_index = HPC_index;
for iS = 5:19
    % MEC_index = (clusters_all.region == "MEC");
    % overall_cluster_index = MEC_index ;


    clusters_this_session = clusters_all.session_count == iS & overall_cluster_index;

    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    w = gausswin(11);
    w = w / sum(w);
    speed(isnan(speed)) = 0;
    speed= filtfilt(w,1,speed')';
    tvec_edges = [tvec(1)-1/(60*2) tvec+1/(60*2)];
    spike_times_this_session = clusters_all.spike_times(clusters_this_session);
    cluster_responses = zeros(length(spike_times_this_session),length(tvec));
    parfor iC = 1:length(spike_times_this_session)
        spike_times = spike_times_this_session{iC,1};
        spike_speed =  interp1(tvec,speed,spike_times,'nearest');
        spike_times_run = spike_times(spike_speed >= 1);
        cluster_responses(iC,:) = histcounts(spike_times_run,tvec_edges);
        % cluster_responses(iC,:) = filtfilt(w,1,cluster_responses(iC,:)')';
    end
    t1 = track_ID_all == 1;
    t1_start = start_time_all(t1);
    t1_end = end_time_all(t1);
    t1_time = logical(zeros(size(tvec)));
    for iT = 1:sum(t1)
        temp_time = tvec >= t1_start(iT) & tvec <= t1_end(iT);
        t1_time = t1_time | temp_time;
    end
    t1_pos = zeros(size(tvec));
    t1_pos(t1_time) = position(t1_time);
    t1_pos(isnan(t1_pos)) = 0;

    t2 = track_ID_all == 2;
    t2_start = start_time_all(t2);
    t2_end = end_time_all(t2);
    t2_time = logical(zeros(size(tvec)));
    for iT = 1:sum(t2)
        temp_time = tvec >= t2_start(iT) & tvec <= t2_end(iT);
        t2_time = t2_time | temp_time;
    end
    t2_pos = zeros(size(tvec));
    t2_pos(t2_time) = position(t2_time);
    t2_pos(isnan(t2_pos)) = 0;
    % Create an instance of bayesDecoder
    decoder_t1 = bayesDecoder;
    decoder_t1.numBins = 140;
    decoder_t1.fixedSmth = 1;
    decoder_t1.variable = [1 140];
    % Define your variable X and firing Y
    t = speed>= 1 & t1_time & t1_pos>=1;
    X = t1_pos(t)'; % replace with your data
    Y = cluster_responses(:,t)'; % replace with your data
    [decoder_t1, prediction, X, Posterior, nPosterior] = decoder_t1.trainDecoder(X, Y);

    mean_EV_t1 = mean(decoder_t1.model.EV);

    clusters_all.reliability(clusters_this_session,1) = mean_EV_t1';

    decoder_t2 = bayesDecoder;
    decoder_t2.numBins = 140;
    decoder_t2.fixedSmth = 1;
    decoder_t2.variable = [1 140];
    % Define your variable X and firing Y
    t = speed>= 1 & t2_time & t2_pos>=1;
    X = t2_pos(t)'; % replace with your data
    Y = cluster_responses(:,t)'; % replace with your data
    [decoder_t2, prediction, X, Posterior, nPosterior] = decoder_t2.trainDecoder(X, Y);

    mean_EV_t2 = mean(decoder_t2.model.EV);

    clusters_all.reliability(clusters_this_session,2) = mean_EV_t2';
end


save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')

%% pre-compute the spatial response
% for V1 and HVA, a delay factor is added to the spike times and it's 0.06s now

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
V1_index = (clusters_all.region == "V1");
HVA_index =  clusters_all.region =='HVA';
spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = (V1_index|HVA_index) & spatially_tuned_neurons;
delay = 0.4;
spatial_response = cell(sum(overall_cluster_index),2);
spatial_response_extended = cell(sum(overall_cluster_index),2);
spatial_response_all = cell(sum(overall_cluster_index),1);
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
            position,speed,track_ID_all,start_time_all,end_time_all,'within',delay,'temporal');
        spatial_response_extended(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay,'temporal');
        spatial_response_all(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within_all',delay,'temporal');
        spatial_response_extended_all(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension_all',delay,'temporal');

    end
    cluster_counter = cluster_counter+no_clusters_this_session;
end
save(fullfile(base_folder, 'spatial_responses_V1_HVA_ks3.mat'),'spatial_response_extended','spatial_response','spatial_response_all','spatial_response_extended_all','-v7.3')


% for V1 and HVA, a delay factor is added to the spike times and it's 0.06s now
%for MEC no delay is added
% select clusters you want to plot for population analysis e.g. only plot neurons from V1
MEC_index = (clusters_all.region == "MEC");
spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = MEC_index & spatially_tuned_neurons ;
delay = 0;
spatial_response = cell(sum(overall_cluster_index),2);
spatial_response_extended = cell(sum(overall_cluster_index),2);
spatial_response_all = cell(sum(overall_cluster_index),1);
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
            position,speed,track_ID_all,start_time_all,end_time_all,'within',delay,'spatial');
        spatial_response_extended(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay,'spatial');
        spatial_response_all(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within_all',delay,'spatial');
        spatial_response_extended_all(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension_all',delay,'spatial');
    end
    cluster_counter = cluster_counter+no_clusters_this_session;

end
save(fullfile(base_folder, 'spatial_responses_MEC_ks3.mat'),'spatial_response_extended','spatial_response','spatial_response_all','spatial_response_extended_all','-v7.3')

HPC_index = (clusters_all.region == "HPC");
spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = HPC_index & spatially_tuned_neurons ;
delay = 0;
spatial_response = cell(sum(overall_cluster_index),2);
spatial_response_extended = cell(sum(overall_cluster_index),2);
spatial_response_all = cell(sum(overall_cluster_index),1);
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
            position,speed,track_ID_all,start_time_all,end_time_all,'within',delay,'spatial');
        spatial_response_extended(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension',delay,'spatial');
        spatial_response_all(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within_all',delay,'spatial');
        spatial_response_extended_all(cluster_counter+iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'extension_all',delay,'spatial');
    end
    cluster_counter = cluster_counter+no_clusters_this_session;

end
save(fullfile(base_folder, 'spatial_responses_HPC_ks3.mat'),'spatial_response_extended','spatial_response','spatial_response_all','spatial_response_extended_all','-v7.3')


%% k-mean clustering for remapping - mainly replicating Giocomo lab paper
addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
% Load clusters_all_all from base_folder
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all_ks3.mat'));
load(fullfile(base_folder,'spatial_responses_MEC_ks3.mat'))
% select clusters you want to plot for population analysis e.g. only plot neurons from V1
MEC_index = (clusters_all.region == "MEC");

spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = MEC_index & spatially_tuned_neurons;
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);
lap_idx_all = cell(15,1);
session_count = 0;
fig1 = figure;
fig2 = figure;
fig3 = figure;
%only look at clusters pass the criteria - i.e. the ones in spatial tuning
for iS = 5:19
    session_count= session_count + 1;
    MEC_clusters = clusters_all.session_count(overall_cluster_index) == iS;
    %convert the spatial response to a matrix
    no_clusters_this_session = sum(MEC_clusters);

    rewarded_laps_index = clusters_all.selected_laps_idx{iS};
    sr_extended_all_session = cellfun(@(x) x(rewarded_laps_index,:), spatial_response_extended_all(MEC_clusters), 'UniformOutput', false);
    no_lap = length(sr_extended_all_session{1,1});
    max_length = 200; % Fixed length for each 1D array
    % Initialize the 3D matrix with NaN values

    spatial_response_matrix = nan(no_clusters_this_session,no_lap,max_length);
    % Fill the 3D matrix
    for i = 1:no_clusters_this_session
        for j = 1:no_lap
            current_array_length = min(length(sr_extended_all_session{i}{j}), max_length);
            spatial_response_matrix(i,j,1:current_array_length) = sr_extended_all_session{i}{j}(1:current_array_length);
        end
    end

    spatial_response_matrix(isnan(spatial_response_matrix)) = 0;
    w = gausswin(9);
    w = w / sum(w);
    smooth_spatial_response_matrix = permute(filtfilt(w,1,permute(spatial_response_matrix,[3,1,2])),[2,3,1]);
    norm_sr_matrix = normalization_cutoff(smooth_spatial_response_matrix,[10,90]);
    all_lap_spatial_response = reshape(permute(smooth_spatial_response_matrix,[2,1,3]),[],no_clusters_this_session*max_length);
    t1_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 1 ;
    t2_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 2 ;
    all_lap_spatial_response(isnan(all_lap_spatial_response)) = 0;
    t1_lap_spatial_response = all_lap_spatial_response(t1_index,:);
    t2_lap_spatial_response = all_lap_spatial_response(t2_index,:);
    no_k = 10;
    no_laps = sum(rewarded_laps_index);
    all_lap_idx = zeros(no_laps,no_k);
    t1_lap_idx = zeros(sum(t1_index),no_k);
    t2_lap_idx = zeros(sum(t2_index),no_k);
    all_sse = zeros(1,10);
    t1_sse = zeros(1,10);
    t2_sse = zeros(1,10);
    for iK = 1:10
        [all_idx, ~,all_sumd] = kmeans(all_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        [t1_idx, ~,t1_sumd] = kmeans(t1_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        [t2_idx, ~,t2_sumd] = kmeans(t2_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        all_lap_idx(:,iK) = all_idx;
        t1_lap_idx(:,iK) = t1_idx;
        t2_lap_idx(:,iK) = t2_idx;
        all_sse(iK) = sum(all_sumd);
        t1_sse(iK) = sum(t1_sumd);
        t2_sse(iK) = sum(t2_sumd);
    end

    figure(fig1);
    subplot(4,4,session_count)
    plot(1:10,all_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS)])
    figure(fig2);
    subplot(4,4,session_count)
    plot(1:10,t1_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS),' Track 1'])
    figure(fig3);
    subplot(4,4,session_count)
    plot(1:10,t2_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS),' Track 2'])

    lap_idx_all{session_count} = all_lap_idx;
    lap_idx_t1{session_count} = t1_lap_idx;
    lap_idx_t2{session_count} = t2_lap_idx;
end
MEC_lap_idx_all = cell(19,1);
MEC_lap_idx_all(5:19) = lap_idx_all;
MEC_lap_idx_t1 = cell(19,1);
MEC_lap_idx_t1(5:19) = lap_idx_t1;
MEC_lap_idx_t2 = cell(19,1);
MEC_lap_idx_t2(5:19) = lap_idx_t2';
clusters_all.MEC_lap_index_all = MEC_lap_idx_all;
clusters_all.MEC_lap_index_t1 = MEC_lap_idx_t1;
clusters_all.MEC_lap_index_t2 = MEC_lap_idx_t2;
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')


%% use k-means and correlation to remove early laps on all clusters
% found multiple sessions that contain early laps having different firing
% rates from later laps - could be something happended in the recording.
% Find these laps and remove them in some analysis.

addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
% Load clusters_all_all from base_folder
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all_ks3.mat'));

%use all the clusters

delay = 0;
lap_idx_all = cell(15,1);
session_count = 0;
fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;
fig5 = figure;
fig6 = figure; % for correlation
%only look at clusters pass the criteria - i.e. the ones in spatial tuning
spike_times_non_empty = cellfun(@(x) ~isempty(x),clusters_all.spike_times,'UniformOutput',false);
for iS = 5:19
    session_count= session_count + 1;
    
    session_clusters = clusters_all.session_count == iS  ...
    & clusters_all.peak_percentile(:,1) >= 0.60 & clusters_all.peak_percentile(:,2) >= 0.60;
    %convert the spatial response to a matrix
    no_clusters_this_session = sum(session_clusters);
    rewarded_laps_index = zeros(length(clusters_all.track_ID_all{iS}),1);
    rewarded_laps_index(clusters_all.rewarded_lap_id{iS}) = 1;
    rewarded_laps_index = logical(rewarded_laps_index);
    
    max_length = 140; % Fixed length for each 1D array
    % Initialize the 3D matrix with NaN values

    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    spike_times_this_session = clusters_all.spike_times(session_clusters);
    spatial_response_all = cell(no_clusters_this_session,1);
    parfor iC = 1:no_clusters_this_session
        spike_times = spike_times_this_session{iC};
        spatial_response_all(iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within_all',delay,'spatial');
    end
    spatial_response_matrix = zeros(no_clusters_this_session,length(track_ID_all),max_length);
    parfor iC = 1:no_clusters_this_session
        spatial_response_matrix(iC,:,:) = spatial_response_all{iC,1};
    
    end
    spatial_response_matrix = spatial_response_matrix(:,rewarded_laps_index,:);
    spatial_response_matrix(isnan(spatial_response_matrix)) = 0;
    w = gausswin(9);
    w = w / sum(w);
    smooth_spatial_response_matrix = permute(filtfilt(w,1,permute(spatial_response_matrix,[3,1,2])),[2,3,1]);
    norm_sr_matrix = normalization_cutoff(smooth_spatial_response_matrix,[10,90]);
    all_lap_spatial_response = reshape(permute(smooth_spatial_response_matrix,[2,1,3]),[],no_clusters_this_session*max_length);
    t1_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 1 ;
    t2_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 2 ;
    all_lap_spatial_response(isnan(all_lap_spatial_response)) = 0;
    t1_lap_spatial_response = all_lap_spatial_response(t1_index,:);
    t2_lap_spatial_response = all_lap_spatial_response(t2_index,:);
    no_k = 10;
    no_laps = sum(rewarded_laps_index);
    all_lap_idx = zeros(no_laps,no_k);
    t1_lap_idx = zeros(sum(t1_index),no_k);
    t2_lap_idx = zeros(sum(t2_index),no_k);
    all_sse = zeros(1,10);
    t1_sse = zeros(1,10);
    t2_sse = zeros(1,10);
    for iK = 1:10
        [all_idx, ~,all_sumd] = kmeans(all_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        [t1_idx, ~,t1_sumd] = kmeans(t1_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        [t2_idx, ~,t2_sumd] = kmeans(t2_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        all_lap_idx(:,iK) = all_idx;
        t1_lap_idx(:,iK) = t1_idx;
        t2_lap_idx(:,iK) = t2_idx;
        all_sse(iK) = sum(all_sumd);
        t1_sse(iK) = sum(t1_sumd);
        t2_sse(iK) = sum(t2_sumd);
    end

    figure(fig1);
    subplot(4,4,session_count)
    plot(1:10,all_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS)])
    figure(fig2);
    subplot(4,4,session_count)
    plot(1:10,t1_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS),' Track 1'])
    figure(fig3);
    subplot(4,4,session_count)
    plot(1:10,t2_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS),' Track 2'])

    figure(fig4);
    subplot(4,4,session_count)

    [R,P] = corrcoef(all_lap_spatial_response');
    imagesc(R)
    colorbar
    clim([0,0.9])
    colormap(flip(gray(256)))
    title(['Session ',num2str(iS),' All Laps'])
    color = hsv(2);
    label_index = all_lap_idx(:,2);
        for i = 1:length(label_index)
            text(0.5,i,num2str(label_index(i)),'Color',color(label_index(i),:),'FontSize',8)
        end

    figure(fig5);
    subplot(4,4,session_count)
    [R,P] = corrcoef(t1_lap_spatial_response');
    imagesc(R)
    colorbar
    clim([0,0.9])
    colormap(flip(gray(256)))
    title(['Session ',num2str(iS),' Track 1'])
    color = hsv(2);
    label_index = t1_lap_idx(:,2);
        for i = 1:length(label_index)
            text(0.5,i,num2str(label_index(i)),'Color',color(label_index(i),:),'FontSize',8)
        end

    figure(fig6);
    subplot(4,4,session_count)
    [R,P] = corrcoef(t2_lap_spatial_response');
    imagesc(R)
    colorbar
    clim([0,0.9])
    colormap(flip(gray(256)))
    title(['Session ',num2str(iS),' Track 2'])
    color = hsv(2);
    label_index = t2_lap_idx(:,2);
        for i = 1:length(label_index)
            text(0.5,i,num2str(label_index(i)),'Color',color(label_index(i),:),'FontSize',8)
        end

    lap_idx_all{session_count} = all_lap_idx;
    lap_idx_t1{session_count} = t1_lap_idx;
    lap_idx_t2{session_count} = t2_lap_idx; 
end
early_lap_idx_all = cell(19,1);
early_lap_idx_all(5:19) = lap_idx_all;
early_lap_idx_t1 = cell(19,1);
early_lap_idx_t1(5:19) = lap_idx_t1;
early_lap_idx_t2 = cell(19,1);
early_lap_idx_t2(5:19) = lap_idx_t2;
clusters_all.early_lap_index_all = early_lap_idx_all;
clusters_all.early_lap_index_t1 = early_lap_idx_t1;
clusters_all.early_lap_index_t2 = early_lap_idx_t2;
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')
mkdir(fullfile(base_folder,'correlation_map'))
saveas(fig4,fullfile(base_folder,'all_laps_all_clusters'),'png')
saveas(fig5,fullfile(base_folder,'track1_all_clusters'),'png')
saveas(fig6,fullfile(base_folder,'track2_all_clusters'),'png')




%% k-mean clustering for remapping - mainly replicating Giocomo lab paper
addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
% Load clusters_all_all from base_folder
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all_ks3.mat'));
load(fullfile(base_folder,'spatial_responses_MEC_ks3.mat'))
% select clusters you want to plot for population analysis e.g. only plot neurons from V1
MEC_index = (clusters_all.region == "MEC");

spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = MEC_index & spatially_tuned_neurons(:,1);
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);

lap_idx_all = cell(15,1);
session_count = 0;
fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;
%only look at clusters pass the criteria - i.e. the ones in spatial tuning
for iS = 5:19
    session_count= session_count + 1;
    MEC_clusters = clusters_all.session_count(overall_cluster_index) == iS;
    %convert the spatial response to a matrix
    no_clusters_this_session = sum(MEC_clusters);
    rewarded_laps_index = clusters_all.selected_laps_idx{iS};
    sr_extended_all_session = cellfun(@(x) x(rewarded_laps_index,:), spatial_response_extended_all(MEC_clusters), 'UniformOutput', false);
    no_lap = length(sr_extended_all_session{1,1});
    max_length = 200; % Fixed length for each 1D array
    % Initialize the 3D matrix with NaN values

    spatial_response_matrix = nan(no_clusters_this_session,no_lap,max_length);
    % Fill the 3D matrix
    for i = 1:no_clusters_this_session
        for j = 1:no_lap
            current_array_length = min(length(sr_extended_all_session{i}{j}), max_length);
            spatial_response_matrix(i,j,1:current_array_length) = sr_extended_all_session{i}{j}(1:current_array_length);
        end
    end

    spatial_response_matrix(isnan(spatial_response_matrix)) = 0;
    w = gausswin(9);
    w = w / sum(w);
    smooth_spatial_response_matrix = permute(filtfilt(w,1,permute(spatial_response_matrix,[3,1,2])),[2,3,1]);
    norm_sr_matrix = normalization_cutoff(smooth_spatial_response_matrix,[0,100]);
    all_lap_spatial_response = reshape(permute(norm_sr_matrix,[2,1,3]),[],no_clusters_this_session*max_length);
    t1_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 1 ;
    t2_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 2 ;
    all_lap_spatial_response(isnan(all_lap_spatial_response)) = 0;
    t1_lap_spatial_response = all_lap_spatial_response(t1_index,:);
    t2_lap_spatial_response = all_lap_spatial_response(t2_index,:);
    no_k = 10;
    no_laps = sum(rewarded_laps_index);
    all_lap_idx = zeros(no_laps,no_k);
    t1_lap_idx = zeros(sum(t1_index),no_k);
    t2_lap_idx = zeros(sum(t2_index),no_k);
    all_sse = zeros(1,10);
    t1_sse = zeros(1,10);
    t2_sse = zeros(1,10);
    for iK = 1:10
        [all_idx, ~,all_sumd] = kmeans(all_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        [t1_idx, ~,t1_sumd] = kmeans(t1_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        [t2_idx, ~,t2_sumd] = kmeans(t2_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        all_lap_idx(:,iK) = all_idx;
        t1_lap_idx(:,iK) = t1_idx;
        t2_lap_idx(:,iK) = t2_idx;
        all_sse(iK) = sum(all_sumd);
        t1_sse(iK) = sum(t1_sumd);
        t2_sse(iK) = sum(t2_sumd);
    end

    figure(fig1);
    subplot(4,4,session_count)
    plot(1:10,all_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS)])
    figure(fig2);
    subplot(4,4,session_count)
    plot(1:10,t1_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS),' Track 1'])
    figure(fig3);
    subplot(4,4,session_count)
    plot(1:10,t2_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS),' Track 2'])
    figure(fig4);
    subplot(4,5,session_count)

    [R,P] = corrcoef(all_lap_spatial_response');
    imagesc(R)
    colorbar
    clim([0,0.75])
    colormap(flip(gray(256)))
    title(['Session ',num2str(iS),' All Laps'])
    color = hsv(2);
    label_index = all_lap_idx(:,2);
        for i = 1:length(label_index)
            text(0.5,i,num2str(label_index(i)),'Color',color(label_index(i),:),'FontSize',8)
        end


    lap_idx_all{session_count} = all_lap_idx;
    lap_idx_t1{session_count} = t1_lap_idx;
    lap_idx_t2{session_count} = t2_lap_idx;
end
MEC_lap_idx_all = cell(19,1);
MEC_lap_idx_all(5:19) = lap_idx_all;
MEC_lap_idx_t1 = cell(19,1);
MEC_lap_idx_t1(5:19) = lap_idx_t1;
MEC_lap_idx_t2 = cell(19,1);
MEC_lap_idx_t2(5:19) = lap_idx_t2';
clusters_all.MEC_lap_index_all = MEC_lap_idx_all;
clusters_all.MEC_lap_index_t1 = MEC_lap_idx_t1;
clusters_all.MEC_lap_index_t2 = MEC_lap_idx_t2;
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')


%% use k-means and correlation to remove early laps on all clusters
% found multiple sessions that contain early laps having different firing
% rates from later laps - could be something happended in the recording.
% Find these laps and remove them in some analysis.

addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
% Load clusters_all_all from base_folder
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all_ks3.mat'));

%use all the clusters

delay = 0;
lap_idx_all = cell(19,1);
session_count = 0;
fig1 = figure;
fig2 = figure;
spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005)& (clusters_all.region == "MEC");
%only look at clusters pass the criteria - i.e. the ones in spatial tuning
spike_times_non_empty = cellfun(@(x) ~isempty(x),clusters_all.spike_times,'UniformOutput',false);
for iS = 5:19
    session_count= session_count + 1;
    
    session_clusters = clusters_all.session_count == iS & spatially_tuned_neurons;
    %convert the spatial response to a matrix
    no_clusters_this_session = sum(session_clusters);
    rewarded_laps_index = zeros(length(clusters_all.track_ID_all{iS}),1);
    rewarded_laps_index(clusters_all.rewarded_lap_id{iS}) = 1;
    rewarded_laps_index = logical(rewarded_laps_index);
    
    max_length = 140; % Fixed length for each 1D array
    % Initialize the 3D matrix with NaN values

    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    spike_times_this_session = clusters_all.spike_times(session_clusters);
    spatial_response_all = cell(no_clusters_this_session,1);
    parfor iC = 1:no_clusters_this_session
        spike_times = spike_times_this_session{iC};
        spatial_response_all(iC,:) = cluster_spatial_responses(spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all,'within_all',delay,'spatial');
    end
    spatial_response_matrix = zeros(no_clusters_this_session,length(track_ID_all),max_length);
    parfor iC = 1:no_clusters_this_session
        spatial_response_matrix(iC,:,:) = spatial_response_all{iC,1};
    
    end
    spatial_response_matrix = spatial_response_matrix(:,rewarded_laps_index,:);
    spatial_response_matrix(isnan(spatial_response_matrix)) = 0;
    w = gausswin(9);
    w = w / sum(w);
    smooth_spatial_response_matrix = permute(filtfilt(w,1,permute(spatial_response_matrix,[3,1,2])),[2,3,1]);
    norm_sr_matrix = normalization_cutoff(smooth_spatial_response_matrix,[10,90]);
    all_lap_spatial_response = reshape(permute(norm_sr_matrix,[2,1,3]),[],no_clusters_this_session*max_length);
    t1_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 1 ;
    t2_index = clusters_all.track_ID_all{iS}(rewarded_laps_index) == 2 ;
    all_lap_spatial_response(isnan(all_lap_spatial_response)) = 0;
    t1_lap_spatial_response = all_lap_spatial_response(t1_index,:);
    t2_lap_spatial_response = all_lap_spatial_response(t2_index,:);
    no_k = 10;
    no_laps = sum(rewarded_laps_index);
    all_lap_idx = zeros(no_laps,no_k);
    t1_lap_idx = zeros(sum(t1_index),no_k);
    t2_lap_idx = zeros(sum(t2_index),no_k);
    all_sse = zeros(1,10);
    t1_sse = zeros(1,10);
    t2_sse = zeros(1,10);
    for iK = 1:10
        [all_idx, ~,all_sumd] = kmeans(all_lap_spatial_response,iK,'Replicates',10,'MaxIter',1000);
        all_lap_idx(:,iK) = all_idx;
        all_sse(iK) = sum(all_sumd);
    end

    figure(fig1);
    subplot(4,5,session_count)
    plot(1:10,all_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS)])

    figure(fig2);
    subplot(4,5,session_count)

    [R,P] = corrcoef(all_lap_spatial_response');
    imagesc(R)
    colorbar
    clim([0,0.75])
    colormap(flip(gray(256)))
    title(['Session ',num2str(iS),' All Laps'])
    color = hsv(2);
    label_index = all_lap_idx(:,2);
        for i = 1:length(label_index)
            text(0.5,i,num2str(label_index(i)),'Color',color(label_index(i),:),'FontSize',8)
        end


    lap_idx_all{session_count} = all_lap_idx;

end

early_lap_idx_all = lap_idx_all;
clusters_all.early_lap_index_all = early_lap_idx_all;

base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')
saveas(fig4,fullfile(base_folder,'correlation_map','all_laps_all_clusters'),'png')
saveas(fig5,fullfile(base_folder,'correlation_map','track1_all_clusters'),'png')
saveas(fig6,fullfile(base_folder,'correlation_map','track2_all_clusters'),'png')
%% use k-means and correlation to remove early laps on all clusters
% found multiple sessions that contain early laps having different firing
% rates from later laps - could be something happended in the recording.
% Find these laps and remove them in some analysis.

addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
% Load clusters_all_all from base_folder
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all_ks3.mat'));

%use all the clusters

delay = 0;
time_all_idx = cell(15,1);
session_count = 0;
fig1 = figure;
fig2 = figure;
fig3 =figure;
time_bin_size = 5;
%only look at clusters pass the criteria - i.e. the ones in spatial tuning
spike_times_non_empty = cellfun(@(x) ~isempty(x),clusters_all.spike_times,'UniformOutput',false);
for iS = 5:19
    session_count= session_count + 1;
    
    session_clusters = clusters_all.session_count == iS ;
    %convert the spatial response to a matrix
    no_clusters_this_session = sum(session_clusters);
    rewarded_laps_index = zeros(length(clusters_all.track_ID_all{iS}),1);
    rewarded_laps_index(clusters_all.rewarded_lap_id{iS}) = 1;
    rewarded_laps_index = logical(rewarded_laps_index);

    tvec = clusters_all.tvec{iS};
    start_time_all = clusters_all.start_time_all{iS};
    end_time_all = clusters_all.end_time_all{iS};
    track_ID_all = clusters_all.track_ID_all{iS};
    position = clusters_all.position{iS};
    speed = clusters_all.speed{iS};
    spike_times_this_session = clusters_all.spike_times(session_clusters);
    time_edges = [0:time_bin_size:floor(max(tvec))];
    cluster_responses = zeros(no_clusters_this_session,length(time_edges)-1);
    parfor iC = 1:no_clusters_this_session
        spike_times = spike_times_this_session{iC};
        cluster_responses(iC,:) = histcounts(spike_times,time_edges)';
    end
    norm_cluster_responses = normalize(cluster_responses,2);
    no_k = 4;
    no_time_bin = size(norm_cluster_responses,2);
    time_bin_idx = zeros(no_time_bin,no_k);
    time_bin_sse = zeros(1,no_k);

    for iK = 1:no_k
        [time_idx, ~,time_sumd] = kmeans(norm_cluster_responses',iK,'Replicates',10,'MaxIter',1000);
        time_bin_idx(:,iK) = time_idx;
        time_bin_sse(iK) = sum(time_sumd);
    end

    figure(fig1);
    subplot(4,4,session_count)
    plot(1:no_k,time_bin_sse)
    xlabel('Number of Clusters')
    ylabel('Sum of Squared Error')
    title(['Session ',num2str(iS)])

    figure(fig2);
    subplot(4,4,session_count)

    [R,P] = corrcoef(norm_cluster_responses);
    imagesc(R)
    colorbar
    clim([0,0.75])
    colormap(flip(gray(256)))
    title(['Session ',num2str(iS),' All Laps'])
    color = hsv(2);
    label_index = time_bin_idx(:,2);
        for i = 1:length(label_index)
            text(0.5,i,num2str(label_index(i)),'Color',color(label_index(i),:),'FontSize',8)
        end
        for i = 1:length(start_time_all)/10
            xline(start_time_all((i-1)*10+1)/time_bin_size,'r')
        end
    
    time_all_idx{session_count} = time_bin_idx;

    all_clusters = sum(norm_cluster_responses,1);
    figure(fig3);
    subplot(4,4,session_count)
    plot(time_bin_size/2:time_bin_size:floor(max(tvec))-time_bin_size/2,all_clusters);

end
saveas(fig2,fullfile(base_folder,'correlation_map','time_all_clusters_correlation'),'png')
saveas(fig3,fullfile(base_folder,'correlation_map','sum_norm_responses_over_time'),'png')


%% manually setting selected_laps based on k-means early_lap_idx
base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
load(fullfile(base_folder, 'clusters_all_ks3.mat'));

clusters_all.selected_laps_idx = cell(19,1);
early_lap_idx = clusters_all.early_lap_index_all;

%session 1 M230310710 early part having problem
idx = early_lap_idx{1};
lap_ids = 1:length(clusters_all.lap_ID_all{1});
rewarded_laps = clusters_all.rewarded_lap_id{1};
selected_laps = rewarded_laps(idx(:,2) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{1} = selected_laps_idx;

%session 2 M230310711 ok
idx = early_lap_idx{2};
lap_ids = 1:length(clusters_all.lap_ID_all{2});
rewarded_laps = clusters_all.rewarded_lap_id{2};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{2} = selected_laps_idx;

%session 3 M230320712 later part having problem
idx = early_lap_idx{3};
lap_ids = 1:length(clusters_all.lap_ID_all{3});
rewarded_laps = clusters_all.rewarded_lap_id{3};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{3} = selected_laps_idx;

%session 4 M230320713 later part having problem
idx = early_lap_idx{4};
lap_ids = 1:length(clusters_all.lap_ID_all{4});
rewarded_laps = clusters_all.rewarded_lap_id{4};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{4} = selected_laps_idx;

%session 5 M230320718 later part having problem
idx = early_lap_idx{5};
lap_ids = 1:length(clusters_all.lap_ID_all{5});
rewarded_laps = clusters_all.rewarded_lap_id{5};
selected_laps = rewarded_laps(idx(:,2) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{5} = selected_laps_idx;

%session 6 M230320719 lots of early ones... half split
idx = early_lap_idx{6};
lap_ids = 1:length(clusters_all.lap_ID_all{6});
rewarded_laps = clusters_all.rewarded_lap_id{6};
selected_laps = rewarded_laps(idx(:,2) == 2);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{6} = selected_laps_idx;

%session 7 M230320720 V1 looks ok but MEC data is difficult to interpret so
%keep everything
idx = early_lap_idx{7};
lap_ids = 1:length(clusters_all.lap_ID_all{7});
rewarded_laps = clusters_all.rewarded_lap_id{7};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{7} = selected_laps_idx;

%session 8 M230320721 looks ok 
idx = early_lap_idx{8};
lap_ids = 1:length(clusters_all.lap_ID_all{8});
rewarded_laps = clusters_all.rewarded_lap_id{8};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{8} = selected_laps_idx;

%session 9 M230320722 looks ok 
idx = early_lap_idx{9};
lap_ids = 1:length(clusters_all.lap_ID_all{9});
rewarded_laps = clusters_all.rewarded_lap_id{9};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{9} = selected_laps_idx;


%session 10 M230340804 looks ok 
idx = early_lap_idx{10};
lap_ids = 1:length(clusters_all.lap_ID_all{10});
rewarded_laps = clusters_all.rewarded_lap_id{10};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{10} = selected_laps_idx;

%session 11 M230340805 looks ok 
idx = early_lap_idx{11};
lap_ids = 1:length(clusters_all.lap_ID_all{11});
rewarded_laps = clusters_all.rewarded_lap_id{11};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{11} = selected_laps_idx;

%session 12 M230340806 looks ok 
idx = early_lap_idx{12};
lap_ids = 1:length(clusters_all.lap_ID_all{12});
rewarded_laps = clusters_all.rewarded_lap_id{12};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{12} = selected_laps_idx;

%session 13 M230340807 looks alright - fewer neurons and laps so difficult
%to tell
idx = early_lap_idx{13};
lap_ids = 1:length(clusters_all.lap_ID_all{13});
rewarded_laps = clusters_all.rewarded_lap_id{13};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{13} = selected_laps_idx;

%session 14 M230370810 only the first part a bit so let it through
idx = early_lap_idx{14};
lap_ids = 1:length(clusters_all.lap_ID_all{14});
rewarded_laps = clusters_all.rewarded_lap_id{14};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{14} = selected_laps_idx;

%session 15 M230370811 first 46 laps a bit problematic
idx = early_lap_idx{15};
lap_ids = 1:length(clusters_all.lap_ID_all{15});
rewarded_laps = clusters_all.rewarded_lap_id{15};
selected_laps = rewarded_laps(idx(:,2) == 2);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{15} = selected_laps_idx;

%session 16 M230370812 alright but MEC a bit for the first 30 laps
idx = early_lap_idx{16};
lap_ids = 1:length(clusters_all.lap_ID_all{16});
rewarded_laps = clusters_all.rewarded_lap_id{16};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{16} = selected_laps_idx;

%session 17 M230370813 looks ok
idx = early_lap_idx{17};
lap_ids = 1:length(clusters_all.lap_ID_all{17});
rewarded_laps = clusters_all.rewarded_lap_id{17};
selected_laps = rewarded_laps(idx(:,1) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{17} = selected_laps_idx;

%session 18 M230380816 first 20 laps
idx = early_lap_idx{18};
lap_ids = 1:length(clusters_all.lap_ID_all{18});
rewarded_laps = clusters_all.rewarded_lap_id{18};
selected_laps = rewarded_laps(idx(:,2) == 2);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{18} = selected_laps_idx;

%session 19 M230380817 first 20 laps
idx = early_lap_idx{19};
lap_ids = 1:length(clusters_all.lap_ID_all{19});
rewarded_laps = clusters_all.rewarded_lap_id{19};
selected_laps = rewarded_laps(idx(:,2) == 1);
selected_laps_idx = ismember(lap_ids,selected_laps);
clusters_all.selected_laps_idx{19} = selected_laps_idx;

save(fullfile(base_folder,'clusters_all_ks3'),"clusters_all",'-v7.3')


%% run shuffles on selected laps

if ismac
    base_folder = '/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data';
    addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
else
    base_folder = 'C:\Users\adam.tong\OneDrive - University College London\data\';
    addpath(genpath('C:\Users\adam.tong\Documents\GitHub\V1_MEC_acute_MAT'))
end
%% setting metrics to screen good clusters

params.amplitude_cutoff = @(x) x<=0.01; %0.01 if strict and removes lots of units
params.sliding_rp_violation = @(x) x<=0.1; % chosen as 10% at IBL
params.num_negative_peaks = @(x) x<=1;
params.num_positive_peaks = @(x) x<=2;
params.peak_to_valley = @(x) x<=0.0008 & x>= 0.0002;
params.amplitude_cv_median = @(x) x<=0.7;
params.amplitude_cv_range = @(x) x<=0.7;
params.firing_rate = @(x) x>= 0.05;
params.presence_ratio = @(x) x>= 0.5;
params.sd_ratio = @(x) x<=4;
params.snr = @(x) x>= 2;
params.half_width = @(x) x<= 0.0003;
%% organize clusters from each session
sessions_to_run = [1,5,6,15,18,19];
SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};
option = 'V1-MEC';
session_count = 0;
for iSub = 1:5
    load(fullfile(base_folder,SUBJECTS{iSub},'analysis','experiment_info.mat'))
    Stimulus_type = 'Track';

    for nsession = 1:length(experiment_info)
        session_count = session_count + 1;
        if ismember(session_count,sessions_to_run)
        session_clusters = struct();
        iDate = nsession;
        session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        session_folder = fullfile(base_folder,SUBJECTS{iSub},'analysis',session_info.probe(1).SESSION);

        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        ANALYSIS_DATAPATH = session_info(1).probe(1).ANALYSIS_DATAPATH;
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH, 'Z:\ibn-vision\DATA\SUBJECTS', base_folder);
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH,'\', '/');
        session_info.probe(1).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
        session_info.probe(2).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
        load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));
        options = session_info.probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'));
        clusters = clusters_ks3;
        sorting_option = 'spikeinterface';
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        receptive_field_file =fullfile(session_folder,'SparseNoise_fullscreen','receptiveFields_ks3.mat');
        if exist(receptive_field_file,"file")
            load(receptive_field_file);
        else
            receptive_field = [];
        end
        for nprobe = 1:length(clusters)
            if isempty(receptive_field)
                clusters(nprobe).receptive_field = cell(size(clusters(nprobe).cluster_id));
            else
                clusters(nprobe).receptive_field = receptive_field{nprobe};
            end
        end
                task_info_fields = fieldnames(Task_info);
        for iField =1:length(task_info_fields)
            session_clusters.(task_info_fields{iField})= {Task_info.(task_info_fields{iField})};
        end
        % add behaviour info for each session
        behaviour_fields = fieldnames(Behaviour);
        for iF = 1:length(behaviour_fields)
            if ~(contains(behaviour_fields{iF},'eye_coordinates'))
                session_clusters.(behaviour_fields{iF}) = {Behaviour.(behaviour_fields{iF})};
            elseif contains(behaviour_fields{iF},'face_motion_SVD')
                session_clusters.eye_coordinates = {Behaviour.face_motion_SVD(:,10)};
            end
        end
        clusters = rmfield(clusters,'sorter');
        clusters = rmfield(clusters,'probe_id');

        for nprobe = 1:length(clusters)
            clusters(nprobe).region = strings(length(clusters(nprobe).cluster_id),1);
            clusters(nprobe).region(:) = 'n.a';
            options = session_info.probe(nprobe);
            if options.probe_id == options.probe_V1
                V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                HPC_channels = determine_region_channels(best_channels{nprobe},options,'region','HPC','group','by probe');
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HPC_channels))) = 'HPC';

            elseif options.probe_id == options.probe_MEC
                %             V1_channels = determine_region_channels(best_channels{nprobe},options,'region','V1','group','by probe');
                MEC_channels = determine_region_channels(best_channels{nprobe},options,'region','MEC_entry','group','by probe');
                HVA_channels = determine_region_channels(best_channels{nprobe},options,'region','HVA','group','by probe');
                %             clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,V1_channels))) = 'V1';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,MEC_channels))) = 'MEC';
                clusters(nprobe).region(find(ismember(clusters(nprobe).peak_channel,HVA_channels))) = 'HVA';
            end
        end
        field_names = fieldnames(clusters);

        for nprobe = 1:length(clusters)
            options = session_info.probe(nprobe);
            clusters_probe = select_clusters(clusters(nprobe),params); %only look at good clusters
            single_cluster_spike_id = cell(size(clusters_probe.cluster_id));
            single_cluster_spike_times = cell(size(clusters_probe.cluster_id));
            unique_spike_id = clusters_probe.spike_id + nprobe*10^4 + iDate* 10^6 + str2double(SUBJECTS{iSub}(2:end))*10^8;
            unique_cluster_id = clusters_probe.cluster_id + nprobe*10^4 + iDate* 10^6 + str2double(SUBJECTS{iSub}(2:end))*10^8;
            [unique_cluster_ids,first_index] = unique(unique_cluster_id);
            for ncell =1:length(unique_cluster_ids)
                single_cluster_spike_id{ncell,1} = unique_spike_id(unique_spike_id == unique_cluster_ids(ncell));
                single_cluster_spike_times{ncell,1} = clusters_probe.spike_times(unique_spike_id == unique_cluster_ids(ncell));
            end
            clusters_probe.cluster_id = unique_cluster_id(first_index);
            clusters_probe.spike_id = single_cluster_spike_id;
            clusters_probe.spike_times = single_cluster_spike_times;

            for iField =1:length(field_names)
                if nprobe == 1
                    session_clusters.(field_names{iField}) = [clusters_probe.(field_names{iField})];

                else
                    session_clusters.(field_names{iField}) = [session_clusters.(field_names{iField});clusters_probe.(field_names{iField})];

                end
            end
        end
        clear clusters clusters_ks3 clusters_probe
        selected_laps = clusters_all.selected_laps_idx{session_count};
        tvec = Behaviour.tvec;
        track_ID =Behaviour.track_ID;
        position = Behaviour.position;
        speed = Behaviour.speed;
        track_ID_all = Task_info.track_ID_all(selected_laps);
        start_time_all = Task_info.start_time_all(selected_laps);
        end_time_all = Task_info.end_time_all(selected_laps);
        no_clusters_this_session = length(session_clusters.cluster_id);
        spatial_response = cell(no_clusters_this_session,2);
        delay = 0;
        within_track_corr = cell(no_clusters_this_session,2);
        first_second_stability = zeros(no_clusters_this_session,2);
        odd_even_stability = zeros(no_clusters_this_session,2);
        peak_percentile = zeros(no_clusters_this_session,2);
        reliablity = zeros(no_clusters_this_session,2);
        odd_even_corr_shuffled = cell(no_clusters_this_session,2);
        first_second_corr_shuffled = cell(no_clusters_this_session,2);


        clear Behaviour Task_info

        all_spike_times = vertcat(session_clusters.spike_times{:});
        cluster_id = session_clusters.cluster_id;
        all_spike_id = vertcat(session_clusters.spike_id{:});
        spatial_response = fast_spatial_response(all_spike_id,cluster_id,all_spike_times,tvec,...
            position,speed,track_ID_all,start_time_all,end_time_all);

        nbatch = 10;
        % Calculate the indices to split the clusters into batches
        clusters_index_batch = round(linspace(0, no_clusters_this_session, nbatch + 1));

        position_shuffled = cell(1000,1);   %shuffled positions for shuffles

        parfor nshuffle = 1:1000
            %     tic
            disp(sprintf('Shuffle %i',nshuffle))
            position_shuffled{nshuffle} = position;

            for nlap = 1:length(start_time_all)
                s = RandStream('mrg32k3a','Seed',nlap+nshuffle*1000); % Set random seed for resampling

                position_shuffled{nshuffle}(tvec  >= start_time_all(nlap) & tvec <= end_time_all(nlap)) =...
                    circshift(position(tvec  >= start_time_all(nlap) & tvec <= end_time_all(nlap))...
                    ,randi(s,length(tvec  >= start_time_all(nlap) & tvec <= end_time_all(nlap))));
            end

        end
        for iBatch = 1:nbatch
            disp(['batch: ',num2str(iBatch)])
            all_spike_times_batch = vertcat(session_clusters.spike_times{clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1)});
            cluster_id_batch = session_clusters.cluster_id(clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1));
            all_spike_id_batch = vertcat(session_clusters.spike_id{clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1)});
            spatial_response_shuffled = cell(clusters_index_batch(iBatch+1)-clusters_index_batch(iBatch),2,1000); %pre-set shuffled response profiles
            parfor nshuffle = 1: 1000
                disp(sprintf('Shuffle %i',nshuffle))
                spatial_response_shuffled(:,:,nshuffle) = fast_spatial_response(all_spike_id_batch,cluster_id_batch,all_spike_times_batch,tvec...
                    ,position_shuffled{nshuffle},speed,track_ID_all,start_time_all,end_time_all);
            end
            clear all_spike_times_batch cluster_id_batch all_spike_id_batch

            for iCluster = clusters_index_batch(iBatch)+1:clusters_index_batch(iBatch+1)
                tic
                % Stability
                for track_id = 1:max(track_ID)
                    start_times = start_time_all(track_ID_all == track_id);

                    lap_correlation = corr(normalize(spatial_response{iCluster,track_id}','range'),...
                        normalize(spatial_response{iCluster,track_id}','range')); % lap by lap correlation

                    first_second_corr = mean(mean(lap_correlation(1:2:end,2:2:end),'omitnan'),'omitnan'); % First half vs Second half
                    odd_even_corr = mean(mean(lap_correlation(1:round(length(start_times)/2),round(length(start_times)/2)+1:end),'omitnan'),'omitnan'); % Odd laps vs even laps

                    within_track_corr{iCluster,track_id} = lap_correlation;

                    normalised_ratemap = normalize(spatial_response{iCluster,track_id}','range');

                    peak_shuffled = zeros(1000,1);
                    odd_even_corr_shuffled_temp = zeros(1000,1);
                    first_second_corr_shuffled_temp = zeros(1000,1);
                    sr_shuffled_temp = spatial_response_shuffled(iCluster-clusters_index_batch(iBatch),track_id,:);
                    parfor nshuffle = 1:1000
                        % Range between min and max FR
                        peak_shuffled(nshuffle,1) = max(mean(sr_shuffled_temp{1,1,nshuffle}))-min(mean(sr_shuffled_temp{1,1,nshuffle}));

                        lap_correlation_shuffled = corr(normalised_ratemap,...
                            normalize(sr_shuffled_temp{1,1,nshuffle}','range')); % original X shuffled correlation

                        odd_even_corr_shuffled_temp(nshuffle) = mean(mean(lap_correlation_shuffled(1:2:end,2:2:end),'omitnan'),'omitnan'); % odd lap from
                        first_second_corr_shuffled_temp(nshuffle) = mean(mean(lap_correlation_shuffled(1:round(length(start_times)/2),round(length(start_times)/2)+1:end),'omitnan'),'omitnan');
                    end
                    clear lap_correlation_shuffled
                    first_second_stability(iCluster,track_id) = sum(first_second_corr > first_second_corr_shuffled_temp)/length(first_second_corr_shuffled_temp);
                    first_second_corr_shuffled{iCluster,track_id} = first_second_corr_shuffled_temp;

                    odd_even_stability(iCluster,track_id) = sum(odd_even_corr > odd_even_corr_shuffled_temp)/length(odd_even_corr_shuffled_temp);
                    odd_even_corr_shuffled{iCluster,track_id} = odd_even_corr_shuffled_temp;



                    peak_percentile(iCluster,track_id) = sum(max(mean(spatial_response{iCluster,track_id}))-min(mean(spatial_response{iCluster,track_id}))>...
                        peak_shuffled)/length(peak_shuffled); % max-min FR relative to max-min of shuffled data


                end
                toc
            end
        end
        clear spatial_response_shuffled
        session_clusters.within_track_corr = within_track_corr;
        session_clusters.first_second_stability = first_second_stability;
        session_clusters.odd_even_stability = odd_even_stability;
        session_clusters.peak_percentile = peak_percentile;
        session_clusters.reliablity = reliablity;
        session_clusters.odd_even_corr_shuffled = odd_even_corr_shuffled;
        session_clusters.first_second_corr_shuffled = first_second_corr_shuffled;
        % add task info for each session


        save(fullfile(session_folder,'session_clusters.mat'),'session_clusters')
        end
    end

end

clusters_all_old = clusters_all;
clusters_all.reliability = clusters_all_old.reliability;
clusters_all.selected_laps_idx = clusters_all_old.selected_laps_idx;



