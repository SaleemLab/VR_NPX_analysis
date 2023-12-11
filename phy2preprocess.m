%%
% 
% spike_templates = readNPY('template_features.npy');
% 
% amplitudes = readNPY('amplitudes.npy');
% 
% 
% cluster_info = tdfread('cluster_info.tsv');
% 
% clear all

working_directory = 'D:\Neuropixel_recording\M22069\20221130\kilosort';
% cd('X:\BendorLab\Drobo\Lab Members\Ben\Neural data\R886\Recording data\23-05-22_kilosort_output\TT21')
folder_info = dir(working_directory);
channels_to_use = [7 16 21 22];
SUA_spike_time = [];
SUA_spike_id = [];
MUA_spike_time = [];
MUA_spike_id = [];

cd(working_directory)
spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');
cluster_group = tdfread('cluster_group.tsv'); % Information about cluster (good,mua or noise)
mean_waveforms = readNPY('mean_waveforms.npy'); % Information about mean waveform


temp = squeeze (mean_waveforms(1,1,:))

count = 1;
% Find MUA and SUA activity
unit_id = [];
unit_type = [];
for unit = 1:length(cluster_group.group)
    if strcmp(cluster_group.group(unit),'g')
        unit_type(count) = 1;
        unit_id(count,:) = [cluster_group.cluster_id(unit) 1000*nchannel+cluster_group.cluster_id(unit)];
        count = count + 1;
    elseif strcmp(cluster_group.group(unit),'m')
        unit_type(count) = 2;
        unit_id(count,:) = [cluster_group.cluster_id(unit) 1000*nchannel+cluster_group.cluster_id(unit)];
        %         elseif strcmp(cluster_group.group(unit),'n')
        count = count + 1;
    end

end

% Find MUA and SUA activity
index = find(unit_type==1);% Find SUA
if ~isempty(index)% If there are SUA
    data.extracted_clusters.(sprintf('TT%i',nchannel)).SUA_clusters = unit_id(index,:);
    SUA_index = [];
    for i = 1:length(index)
        SUA_index = [SUA_index; find(spike_clusters == unit_id(index(i),1))];
    end
    SUA_spike_time = [SUA_spike_time; spike_times(SUA_index)];
    SUA_spike_id = [SUA_spike_id; 1000*nchannel+spike_clusters(SUA_index)];
end

index = find(unit_type==2);% Find MUA
if ~isempty(index) %If there are MUA
    data.extracted_clusters.(sprintf('TT%i',nchannel)).MUA_clusters = unit_id(index,:);
    MUA_index = [];
    for i = 1:length(index)
        MUA_index = [MUA_index; find(spike_clusters == unit_id(index(i),1))];
    end
    %1 is spike time and 2 is cell id
    MUA_spike_time = [MUA_spike_time; spike_times(MUA_index)];
    MUA_spike_id = [MUA_spike_id; 1000*nchannel+spike_clusters(MUA_index)];
end

% Sort spike time from start to end
[~,index] = sort(SUA_spike_time);
data.spikes.SUA(:,1) = SUA_spike_id(index);
data.spikes.SUA(:,2) = SUA_spike_time(index);

[~,index] = sort(MUA_spike_time);
data.spikes.MUA(:,1) = MUA_spike_id(index);
data.spikes.MUA(:,2) = MUA_spike_time(index);

cluster_data = data;
cd('D:\Neuropixel_recording\M22069\20221130')
save cluster_data cluster_data -v7.3
%%


% Extract position data from video
cd('X:\BendorLab\Drobo\Lab Members\Ben\Neural data\R886\Recording data\23-05-22_raw_data');
extract_video;
% batch_analysis('EXTRACT_VIDEO')
process_positions_PRE;


% POST processing
rexposure_exception=[];  %if you have four tracks but they are not rexposures

cd('X:\BendorLab\Drobo\Lab Members\Ben\Neural data\R886\Recording data\23-05-22_raw_data')
load position_data;
    number_of_tracks=length(position.linear);
    if number_of_tracks<4 | ~isempty(rexposure_exception)
        rexposure=[];
    else
        rexposure=[1,3; 2,4];
    end
    extract_events;
    extract_dropped_samples;
    process_clusters;
    process_positions_POST(rexposure);  %input is empty unless you have track rexposure
    sleep_stager('auto');  %process sleep stages with default thresholds

% % % % Skip sleep and oscilation analysis
disp('processing CSC data')
% cd('X:\BendorLab\Drobo\Lab Members\Ben\Neural data\R885\Recording data\06-04-22_raw_data')
p = gcp; % Starting new parallel pool
if ~isempty(p)
    parallel_extract_PSD('sleep');
else
    disp('parallel processing not possible');
    extract_PSD([],'sleep');
end
% plot_PSDs;
best_channels = determine_best_channel('hpc');
extract_CSC('hpc');

%% Place field

cd('X:\BendorLab\Drobo\Lab Members\Ben\Neural data\R886\Recording data\23-05-22_place_field_calculation')
% EXTRACT PLACE FIELDS
disp('processing place_field data')
parameters=list_of_parameters;
calculate_place_fields(parameters.x_bins_width_bayesian);
calculate_place_fields(parameters.x_bins_width);
spike_count([],[],[],'Y');
bayesian_decoding([],[],'Y');
extract_replay_events; %finds onset and offset of replay events

load extracted_place_fields_BAYESIAN
fig = figure(1)
fig.Position = [720 140 1130 841]
for n = 1:16
    subplot(4,4,n)
    title(sprintf('Cell %i',n))
    plot(place_fields_BAYESIAN.track(1).raw{n})
%     plot(place_fields.track(1).smooth{n})
    hold on
%      plot(place_fields.track(2).smooth{n},'r')
     plot(place_fields_BAYESIAN.track(2).raw{n},'r')
%      legend('Track 1','Track 2')
end
saveas(gcf,'place cells 1_no legend.fig')
saveas(gcf,'place cells 1.pdf')
count = 1;

fig = figure(2)
fig.Position = [720 140 1130 841]
for n = 17:32
    subplot(4,4,count)
    title(sprintf('Cell %i',n))
    plot(place_fields_BAYESIAN.track(1).raw{n})
    hold on
    plot(place_fields_BAYESIAN.track(2).raw{n},'r')
%       legend('Track 1','Track 2')
      count = count + 1;
end
saveas(gcf,'place cells 2_no legend.fig')
saveas(gcf,'place cells 2.pdf')


fig = figure(2)
fig.Position = [720 140 1130 841]
for n = 1:16
    subplot(4,4,n)
    title(sprintf('Cell %i',n))
    plot(place_fields_BAYESIAN.track(3).raw{n})
    hold on
     plot(place_fields_BAYESIAN.track(4).raw{n},'r')
%      legend('Track 1','Track 2')
end

saveas(gcf,'Re place cells 1_no legend.fig')
saveas(gcf,'Re place cells 1.pdf')
count = 1;

fig = figure(2)
fig.Position = [720 140 1130 841]
for n = 17:32
    subplot(4,4,count)
    title(sprintf('Cell %i',n))
    plot(place_fields_BAYESIAN.track(3).raw{n})
    hold on
    plot(place_fields_BAYESIAN.track(4).raw{n},'r')
%       legend('Reexposure Track 1','Reexposure Track 2')
      count = count + 1;
end
saveas(gcf,'Re place cells 2_no legend.fig')
saveas(gcf,'Re place cells 2.pdf')

% Find cells 
% SU_id = unique(cluster_data.spikes.SUA(:,1));

% Use find() function to find


% 
% bayesian_decoding_error();

 % EXTRACT REPLAY EVENTS and BAYESIAN DECODING
    disp('processing replay events')
    replay_decoding; %extract and decodes replay events
    % SCORING METHODS: TEST SIGNIFICANCE ON REPLAY EVENTS
    disp('scoring replay events')
    scored_replay = replay_scoring([],[1 1 1 1]);
    save scored_replay scored_replay;
    
% RUN SHUFFLES
    disp('running shuffles')
    num_shuffles=1000;
%     analysis_type=[1 1 1 0];  %just linear fit, weighted correlation and pacman
    analysis_type=[0 1 1 0];  %just weighted correlation and pacman
    load decoded_replay_events
    p = []; %gcp; % Starting new parallel pool
    shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift'};
    
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events);
        end
    end
    save shuffled_tracks shuffle_type;
    
  % Evaluate significance
    load('scored_replay.mat');
    load('shuffled_tracks.mat');
    scored_replay=replay_significance(scored_replay, shuffle_type);
    save scored_replay scored_replay
    
    %%%%%%analyze segments%%%%%%%%%%
    %splitting replay events
    replay_decoding_split_events;
    load decoded_replay_events_segments;
    scored_replay1 = replay_scoring(decoded_replay_events1,[0 1 1 1]);
    scored_replay2 = replay_scoring(decoded_replay_events2,[0 1 1 1]);
    save scored_replay_segments scored_replay1 scored_replay2;
    num_shuffles=1000;
    analysis_type=[0 1 1 0];  %just weighted correlation and pacman
    load decoded_replay_events_segments;
    p = [];%gcp; % Starting new parallel pool
    shuffle_choice={'PRE spike_train_circular_shift','PRE place_field_circular_shift', 'POST place bin circular shift'};
    if ~isempty(p)
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
            shuffle_type2{shuffle_id}.shuffled_track = parallel_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
        end
    else
        disp('parallel processing not possible');
        for shuffle_id=1:length(shuffle_choice)
            shuffle_type1{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events1);
            shuffle_type2{shuffle_id}.shuffled_track = run_shuffles(shuffle_choice{shuffle_id},analysis_type,num_shuffles,decoded_replay_events2);
        end
    end
    save shuffled_tracks_segments shuffle_type1 shuffle_type2;
    load scored_replay_segments; load shuffled_tracks_segments;
    scored_replay1=replay_significance(scored_replay1, shuffle_type1);
    scored_replay2=replay_significance(scored_replay2, shuffle_type2);
    save scored_replay_segments scored_replay1 scored_replay2


  %EXTRACT SLEEP
    sleep_stager('manual');

  %ANALYSIS
    rexposure_exception=[];
    current_directory=pwd;
    load scored_replay;
    number_of_tracks=length(scored_replay);
    if number_of_tracks<4 | ~isempty(rexposure_exception)
        rexposure=[];
    else
        rexposure=[1,3; 2,4];
    end
    
%     scored_replay_first_exposure(1) = scored_replay(1)
%     scored_replay_first_exposure(2) = scored_replay(2)
%     
%     scored_replay_second_exposure(1) = scored_replay(3)
%     scored_replay_second_exposure(2) = scored_replay(4)
    rexposure_option = 2;
    
    significant_replay_events=number_of_significant_replays(0.05,3,'wcorr',rexposure_option);
    [sorted_replay_index,time_range]=sort_replay_events(rexposure_option,'wcorr');
    cumulative_replay([]);
    cd(current_directory);
    
  %COMPARE METHODS
    current_directory=pwd;
    load scored_replay;
    number_of_tracks=length(scored_replay);
    if number_of_tracks<4 | ~isempty(rexposure_exception)
        rexposure=[];
    else
        rexposure=[1,3; 2,4];
    end
    significant_replay_events=number_of_significant_replays(0.05,3,'linear',rexposure);
    plot_significant_events;
    significant_replay_events=number_of_significant_replays(0.05,3,'wcorr',rexposure);
    plot_significant_events;
    significant_replay_events=number_of_significant_replays(0.05,3,'path',rexposure);
    plot_significant_events;
    significant_replay_events=number_of_significant_replays(0.05,3,'spearman',rexposure);
    plot_significant_events;
    cd(current_directory);
    
%PLOT CLEANING STEPS
    plot_cleaning_steps;
    
%PLOT PLACE FIELDS
    plot_place_fields;
    
%PLOT SIGNIFICANT EVENTS
    plot_significant_events;
    
%PLOT CUMULATIVE REPLAY
    plot_cumulative_replay;
    
%PLOT RATE REMAP
    current_directory=pwd;
    %out=rate_remapping_ONE_TRACK([]);
    %plot_rate_remapping('ONE_TRACK');
    cd(current_directory)
    out=rate_remapping_TRACK_PAIRS([]);
    plot_rate_remapping('TRACK_PAIRS');
    cd(current_directory)
    
    
 extract_laps();   
 bayesian_decoding_error('method','cross_tracks');
 