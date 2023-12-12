%% Putative pipeline for processing NPX1 data recorded from hippocampus and V1
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)


%% Set the data folders and processing parameters
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% addpath('Z:\ibn-vision\USERS\Masa\code\Masa_utility')
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\NPXAnalysis\NPXAnalysis2022'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\visual_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'));



if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
%     ROOTPATH = 'X:\ibn-vision';
    ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive
%     ROOTPATH = '/research';
end

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'Masa2tracks';

%% import and align and store Bonsai data
for nsession = 1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            session_info(n).probe(nprobe).task_type = stimulus_name{n};

        end

        [peripherals,photodiodeData,eyeData] = import_and_align_Masa_VR_Bonsai(StimulusName,session_info(n).probe);
        [] = align_probes_NX1(session_info.probe);
        disp('Bonsai data loaded and alinged to both ephys sync pulse and photodiode signal')

        % figure
        % hold on
        % plot(MousePos.sglxTime,MousePos.pos)
        % scatter(MousePos.stimuli_onset(MousePos.stimuli_track == 1),1000*MousePos.stimuli_track(MousePos.stimuli_track == 1),'r')
        % scatter(MousePos.stimuli_onset(MousePos.stimuli_track == 2),100*MousePos.stimuli_track(MousePos.stimuli_track == 2),'b')
    end

end


%% Load LFP data
Stimulus_type = 'RUN'; % extract LFP during RUN
for nsession = 1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,'RUN'));
   

    for nprobe = 1:length(session_info.probe)
        options = session_info.probe(nprobe);

        options.importMode = 'LF';
        column = 1;
        [LF_FILE imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);
        [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column);
        [PSD power best_channels] = calculate_channel_PSD(raw_LFP,SR,sorted_config,options,'plot_option',1)
        
        
        % [gamma_coherence gamma_phase_coherence] = gamma_coherence_analysis(raw_LFP,tvec,SR,best_channels,sorted_config) % This function
        cd(options.ANALYSIS_DATAPATH)
        save extracted_PSD PSD power
        save best_channels best_channels
    end
end

%% Load Spike train data
options.importMode = 'KS';

if ~isempty(tvec)
    LFP_tvec = tvec;
else
    LFP_tvec = [];
end

% Load all spike data sorted according to the channel position
[SUA chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'group','by channel');
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
[L4_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L4_channel-5 best_channels.L4_channel+5],'group','by region');
% [L4_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L4_channel-10 best_channels.L4_channel],'group','by region');
% L5 spike data (roughly 200 micron or 20 channels))
[L5_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L5_channel-10 best_channels.L5_channel+10],'group','by region');

% all V1 spike data (roughly 300 micron or 30 channels))
[superficial_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.first_in_brain_channel-30 best_channels.first_in_brain_channel],'group','by region');
% [superficial_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.first_in_brain_channel-15 best_channels.first_in_brain_channel],'group','by region');

% all V1 spike data (roughly 300 micron or 30 channels))
[V1_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L5_channel-10 best_channels.first_in_brain_channel],'group','by region');


% CA1 spike data (roughly 200 micron or 20 channels))
[CA1_clusters chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.CA1_channel-10 best_channels.CA1_channel+10],'group','by region');



%% Detect Behavioural State

speed_interp = interp1(peripherals.sglxTime,peripherals.speed,tvec','linear');
speedTreshold = 1;

[freezing,quietWake,SWS,REM,movement] = detect_behavioural_states_masa(...
    [tvec' raw_LFP(find(sorted_config.Channel == best_channels.L5_channel),:)'],[tvec' raw_LFP(find(sorted_config.Channel == best_channels.CA1_channel),:)'],...
    [tvec' speed_interp],speedTreshold);

behavioural_state.freezing = freezing;
behavioural_state.quietWake = quietWake;
behavioural_state.SWS = SWS;
behavioural_state.REM = REM;
behavioural_state.movement = movement;

cd(options.EPHYS_DATAPATH)
save behavioural_state behavioural_state

%% Detect ripple and cortical slow wave oscillation (SO) and cortical spindles

% Detect CA1 populational bursting events (Candidate events)
zscore_min = 0;
zscore_max = 3;

cd(options.EPHYS_DATAPATH)
channel_to_use = find(sorted_config.Channel == best_channels.CA1_channel);
[replay,reactivations] = detect_candidate_events_masa(tvec,raw_LFP(channel_to_use,:),...
    CA1_clusters.MUA_zscore,[CA1_clusters.spike_id CA1_clusters.spike_times],peripherals,zscore_min,zscore_max,options) 
save extracted_candidate_events replay reactivations

channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);
[~,L5_reactivations] = detect_candidate_events_masa(tvec,raw_LFP(channel_to_use,:),...
    L5_clusters.MUA_zscore,[L5_clusters.spike_id L5_clusters.spike_times],peripherals,zscore_min,zscore_max,options) 
save extracted_L5_candidate_events L5_reactivations

% Detect CA1 ripple events
channel_to_use = find(sorted_config.Channel == best_channels.CA1_channel);
[ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,...
    'noise',raw_LFP(2,:)','passband',[125 300],'thresholds',[3 5])

figure
[ripples.SWS_offset,ripples.SWS_index] = RestrictInts(ripples.offset,behavioural_state.SWS);
ripples.SWS_onset = ripples.onset(ripples.SWS_index);
ripples.SWS_peaktimes = ripples.peaktimes(ripples.SWS_index);

histogram(abs(ripples.SWS_onset-ripples.SWS_offset)',0:0.005:0.2,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','Normalization','probability')

[ripples.awake_offset,ripples.awake_index] = RestrictInts(ripples.offset,behavioural_state.quietWake);
ripples.awake_onset = ripples.onset(ripples.awake_index);
ripples.awake_peaktimes = ripples.peaktimes(ripples.awake_index);
hold on
histogram(abs(ripples.awake_onset-ripples.awake_offset)',0:0.005:0.2,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','Normalization','probability')
legend('NREM Ripples','awake Ripples')
ylabel('Probability')
xlabel('Duration (sec)')

save extracted_ripples_events ripples

% Detect Cortical ripple events
channel_to_use = find(sorted_config.Channel == best_channels.first_in_brain_channel -12);% 240 micron
[V1_ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[125 300])
[V1_ripples.SWS_offset,V1_ripples.SWS_index] = RestrictInts(V1_ripples.offset,behavioural_state.SWS);
V1_ripples.SWS_onset = V1_ripples.onset(V1_ripples.SWS_index);
V1_ripples.SWS_peaktimes = V1_ripples.peaktimes(V1_ripples.SWS_index);

[V1_ripples.awake_offset,V1_ripples.awake_index] = RestrictInts(V1_ripples.offset,behavioural_state.quietWake);
V1_ripples.awake_onset = V1_ripples.onset(V1_ripples.awake_index);
V1_ripples.awake_peaktimes = V1_ripples.peaktimes(V1_ripples.awake_index);

V1_superficial_ripples = V1_ripples;

channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);% 240 micron
[V1_ripples] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[125 300])
[V1_ripples.SWS_offset,V1_ripples.SWS_index] = RestrictInts(V1_ripples.offset,behavioural_state.SWS);
V1_ripples.SWS_onset = V1_ripples.onset(V1_ripples.SWS_index);
V1_ripples.SWS_peaktimes = V1_ripples.peaktimes(V1_ripples.SWS_index);

[V1_ripples.awake_offset,V1_ripples.awake_index] = RestrictInts(V1_ripples.offset,behavioural_state.quietWake);
V1_ripples.awake_onset = V1_ripples.onset(V1_ripples.awake_index);
V1_ripples.awake_peaktimes = V1_ripples.peaktimes(V1_ripples.awake_index);

V1_deep_ripples = V1_ripples;
save extracted_V1_ripples_events V1_superficial_ripples V1_deep_ripples

figure
histogram(abs(V1_ripples.SWS_onset-V1_ripples.SWS_offset)',0:0.005:0.2,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','Normalization','probability')
hold on
histogram(abs(V1_ripples.awake_onset-V1_ripples.awake_offset)',0:0.005:0.2,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','Normalization','probability')
legend('NREM cortical Ripples','awake cortical Ripples')
ylabel('Probability')
xlabel('Duration (sec)')

% Detect cortical gamma events
channel_to_use = find(sorted_config.Channel == best_channels.first_in_brain_channel -12);% 240 micron
[gamma_events] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[60 100])
[gamma_events.SWS_offset,gamma_events.SWS_index] = RestrictInts(gamma_events.offset,behavioural_state.SWS);
gamma_events.SWS_onset = gamma_events.onset(gamma_events.SWS_index);
gamma_events.SWS_peaktimes = gamma_events.peaktimes(gamma_events.SWS_index);

[gamma_events.awake_offset,gamma_events.awake_index] = RestrictInts(gamma_events.offset,behavioural_state.quietWake);
gamma_events.awake_onset = gamma_events.onset(gamma_events.awake_index);
gamma_events.awake_peaktimes = gamma_events.peaktimes(gamma_events.awake_index);

V1_superficial_gamma_events = gamma_events;

channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);
[gamma_events] = FindRipples_masa(raw_LFP(channel_to_use,:)',tvec','minDuration',20,'durations',[30 200],'frequency',SR,'noise',raw_LFP(1,:)','passband',[60 100])
[gamma_events.SWS_offset,gamma_events.SWS_index] = RestrictInts(gamma_events.offset,behavioural_state.SWS);
gamma_events.SWS_onset = gamma_events.onset(gamma_events.SWS_index);
gamma_events.SWS_peaktimes = gamma_events.peaktimes(gamma_events.SWS_index);

[gamma_events.awake_offset,gamma_events.awake_index] = RestrictInts(gamma_events.offset,behavioural_state.quietWake);
gamma_events.awake_onset = gamma_events.onset(gamma_events.awake_index);
gamma_events.awake_peaktimes = gamma_events.peaktimes(gamma_events.awake_index);

V1_deep_gamma_events = gamma_events;

save extracted_V1_gamma_events V1_superficial_gamma_events V1_deep_gamma_events

% Detect Cortical spindle events
[spindles] = FindSpindles_masa(raw_LFP(channel_to_use,:)',tvec','durations',[400 3000],'frequency',SR,'noise',raw_LFP(1,:)','passband',[9 17],'thresholds',[1.5 3])
[spindles.SWS_offset,spindles.SWS_index] = RestrictInts(spindles.offset,SWS);
spindles.SWS_onset = spindles.onset(spindles.SWS_index);
spindles.SWS_peaktimes = spindles.peaktimes(spindles.SWS_index);

[spindles.awake_offset,spindles.awake_index] = RestrictInts(spindles.offset,behavioural_state.quietWake);
spindles.awake_onset = spindles.onset(spindles.awake_index);
spindles.awake_peaktimes = spindles.peaktimes(spindles.awake_index);
save extracted_spindles_events spindles

% Detect Slow wave Up and Down states (Using all layer 5)
options.importMode = 'KS';
[spikes chan_config_KS sorted_config_KS] = load_KS_NPX1(options,column,'LFP_tvec',LFP_tvec,'selected_channels',[best_channels.L5_channel-10 best_channels.L4_channel + 10] ,'group','Buz style');

channel_to_use = find(sorted_config.Channel == best_channels.L5_channel);
slow_waves= DetectSlowWaves_masa('time',tvec,'lfp',raw_LFP(channel_to_use,:)','NREMInts',behavioural_state.SWS,'spikes',spikes);

% save best_channels best_spindle_channel best_spindle_channel best_ripple_channel
save extracted_slow_waves slow_waves

%% Peri-event LFP amplitude and phase

% Stimulus
if strcmp(StimulusName,'replay_Masa2tracks') ;
    all_events = {MousePos.stimuli_onset(MousePos.stimuli_track == 1)',MousePos.stimuli_onset(MousePos.stimuli_track == 2)'};
    event_group = {'T1 stimuli','T2 stimuli'};


    lfpAvg = [];
    csd = [];
    lfpAvg.filter_type = {'SO','all'};
    lfpAvg.event_group = event_group;

    % Filtered at broad band (0.5Hz - 300Hz)
    [ csd.all, lfpAvg.all ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
        'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.2 0.8],'filter',[0.5 300]);

    % Filtered at slow wave oscilation band (0.5Hz - 4Hz)
    [ csd.SO, lfpAvg.SO ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
        'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.2 0.8],'filter',[0.5 4]);

    plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power,chan_config,sorted_config,best_channels);
    save visual_scene_LFP lfpAvg csd
end

% Brain events
all_events = {ripples.SWS_peaktimes',slow_waves.ints.UP(:,1)',spindles.SWS_peaktimes,V1_superficial_ripples.SWS_peaktimes',V1_deep_ripples.SWS_peaktimes'};
event_group = {'Ripple','UP','Spindle','V1 Superficial Ripple','V1 Deep Ripple'};

lfpAvg = [];
csd = [];
lfpAvg.event_group = event_group;
% lfpAvg.filter_type = {'SO','spindle','gamma','ripple','all'};
lfpAvg.filter_type = {'SO','spindle','gamma','ripple','all'};

% Filtered at slow wave oscilation band (0.5Hz - 4Hz)
[ csd.SO, lfpAvg.SO ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[0.5 4]);

% [ csd1, lfp1 ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
%     'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[1 1],'filter',[]);

% Filtered at spindle range (9Hz - 17Hz)
[ csd.spindle, lfpAvg.spindle ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[9 17]);

% Filtered at gamma range (30Hz - 100Hz)
[ csd.gamma, lfpAvg.gamma ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[30 100]);


% Filtered at ripple range (125Hz - 300Hz)
[ csd.ripple, lfpAvg.ripple ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[125 300]);

% Filtered at broad band (0.5Hz - 300Hz)
[ csd.all, lfpAvg.all ]  = perievent_CSD_LFP_amplitude_phase(raw_LFP',tvec',all_events,...
    'channels',1:1:size(raw_LFP,1),'samplingRate',SR,'twin',[0.5 0.5],'filter',[0.5 300]);

save peri_event_LFP lfpAvg csd

lfpAvg.filter_type = {'SO'}
plot_perievent_CSD_LFP_amplitude_phase(lfpAvg,csd,power,chan_config,sorted_config,best_channels);



%% Peri-event spike time histogram

% event = slow_waves.ints.DOWN(:,1);
% event = ripples.onset
% event = slow_waves.ints.UP(:,1);
% event = spindles.onset;
% event = V1_ripples.onset;
% event = reactivations.onset;
% event = MousePos.stimuli_onset(MousePos.stimuli_track == 1);
% events = slow_waves.ints.UP(:,1);

all_spike_data{1} = [superficial_clusters.spike_id superficial_clusters.spike_times];
all_spike_data{2} = [L4_clusters.spike_id L4_clusters.spike_times];
all_spike_data{3} = [L5_clusters.spike_id L5_clusters.spike_times];
all_spike_data{4} = [CA1_clusters.spike_id CA1_clusters.spike_times];

group_name = {'superficial','L4','L5','CA1'};

% V1 ripple
group_type = 'by cell zscore';

group_type = 'by region';
events = V1_superficial_ripples.SWS_peaktimes;
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 superficial layer NREM Ripple','twin',[-1.2 1.2])
events = V1_deep_ripples.SWS_peaktimes;
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 deep layer NREM Ripple','twin',[-1.2 1.2])

events = V1_superficial_ripples.awake_peaktimes;
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 superficial layer awake Ripple','twin',[-1.2 1.2])
events = V1_deep_ripples.awake_peaktimes;
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 deep layer awake Ripple','twin',[-1.2 1.2])


% Up state
events = slow_waves.ints.UP(:,1);
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by region','group_name',group_name,'event_name','Up state')
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by cell zscore','group_name',group_name,'event_name','Up state')

% Spindle
events = spindles.SWS_peaktimes;
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by region','group_name',group_name,'event_name','Spindle')
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by cell','group_name',group_name,'event_name','Spindle')

% Ripple
events = ripples.awake_onset
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by cell zscore','group_name',group_name,'event_name','CA1 awake Ripple','twin',[-0.6 0.6])
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by region','group_name',group_name,'event_name','CA1 awake Ripple','twin',[-1 1])
events = ripples.SWS_onset
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by cell zscore','group_name',group_name,'event_name','CA1 NREM Ripple','twin',[-0.6 0.6])
plot_perievent_spiketime_histogram(all_spike_data,events,'group','by region','group_name',group_name,'event_name','CA1 NREM Ripple','twin',[-1 1])

% Cortical gamma
group_type = 'by region';
events = V1_deep_gamma_events.awake_onset
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 deep layer awake Gamma','twin',[-1.1 1.1])
events = V1_deep_gamma_events.SWS_onset
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 deep layer NREM Gamma','twin',[-1.1 1.1])

events = V1_superficial_gamma_events.awake_onset
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 superficial layer awake Gamma','twin',[-1.1 1.1])
events = V1_superficial_gamma_events.SWS_onset
plot_perievent_spiketime_histogram(all_spike_data,events,'group',group_type,'group_name',group_name,'event_name','V1 superficial layer NREM Gamma','twin',[-1.1 1.1])


% Track stimuli
if strcmp(StimulusName,'replay_Masa2tracks') ;

    events = MousePos.stimuli_onset(MousePos.stimuli_track == 1);
    plot_perievent_spiketime_histogram(SUA,events,'group','by channel','event_name','Track 1 stimuli','twin',[0 0.8],'channel_map',chan_config_KS.Ks_ycoord')
    plot_perievent_spiketime_histogram(all_spike_data,events,'group','by region','group_name',group_name,'event_name','Track 1 stimuli','twin',[-0.2 0.8],'channel_map',chan_config_KS.Ks_ycoord')

    events = MousePos.stimuli_onset(MousePos.stimuli_track == 2);
    plot_perievent_spiketime_histogram(SUA,events,'group','by channel','event_name','Track 2 stimuli','twin',[0 0.8],'channel_map',chan_config_KS.Ks_ycoord')
    plot_perievent_spiketime_histogram(all_spike_data,events,'group','by region','group_name',group_name,'event_name','Track 2 stimuli','twin',[-0.2 0.8],'channel_map',chan_config_KS.Ks_ycoord')
end

% SUA chan_config_KS sorted_config_KS.Ks_ycoord
sorted_config_KS.Ks_ycoord

%% UP and ripple relationship

ripple_UP_relationship_masa


%% Peri event event histogram

[spindles.onset_with_UP,spindles.index_with_UP] = RestrictInts(spindles.SWS_onset,slow_waves.ints.UP);
spindles.offset_with_UP = spindles.onset(spindles.index_with_UP);
spindles.peaktimes_with_UP = spindles.peaktimes(spindles.index_with_UP);

event1 = {slow_waves.ints.UP(:,1)',spindles.peaktimes,V1_superficial_ripples.peaktimes,V1_deep_ripples.peaktimes};
event2 = {ripples.peaktimes',ripples.peaktimes',ripples.peaktimes',ripples.peaktimes'};
event1_name = {'UP state','Spindle','V1 superficial Ripple','V1 Deep Ripple','UP state','Spindle','Spindle + UP','Up State'};
event2_name = {'Ripple','Ripple','Ripple','Ripple','V1 Ripple','Ripple','Spindle'};

event1 = {slow_waves.ints.UP(:,1)',spindles.SWS_peaktimes,V1_superficial_ripples.SWS_peaktimes,V1_deep_ripples.SWS_peaktimes,V1_deep_ripples.SWS_peaktimes};
event2 = {ripples.SWS_peaktimes',ripples.SWS_peaktimes',ripples.SWS_peaktimes',ripples.SWS_peaktimes',spindles.SWS_peaktimes};
event1_name = {'UP state','Spindle','V1 superficial Ripple','V1 deep Ripple','V1 deep Ripple','Spindle','Spindle + UP','Up State'};
event2_name = {'Ripple','Ripple','Ripple','Ripple','Spindle','Spindle','Spindle'};

event1 = {slow_waves.ints.UP(:,1)',spindles.SWS_peaktimes,V1_superficial_gamma_events.SWS_peaktimes,V1_deep_gamma_events.SWS_peaktimes};
event2 = {ripples.SWS_peaktimes',ripples.SWS_peaktimes',ripples.SWS_peaktimes',ripples.SWS_peaktimes'};
event1_name = {'UP state','Spindle','V1 siperficial Gamma','V1 deep Gamma','UP state','Spindle','Spindle + UP','Up State'};
event2_name = {'Ripple','Ripple','Ripple','Ripple','V1 Ripple','Ripple','Spindle'};

% event1 = {slow_waves.ints.UP(:,1)',spindles.awake_peaktimes,V1_superficial_ripples.awake_peaktimes,V1_deep_ripples.awake_peaktimes};
% event2 = {ripples.awake_peaktimes',ripples.awake_peaktimes',ripples.awake_peaktimes',ripples.awake_peaktimes'};
% event1_name = {'UP state','Spindle','V1 siperficial Ripple','V1 deep Ripple','UP state','Spindle','Spindle + UP','Up State'};
% event2_name = {'Ripple','Ripple','Ripple','Ripple','V1 Ripple','Ripple','Spindle'};

figure
for n = 1:length(event1)
    nexttile
    plot_perievent_event_histogram(event1{n},event2{n},'twin',[-1.1 1.1])
title(sprintf('%s relative to %s',event1_name{n},event2_name{n}))
end
sgtitle('NREM Sleep')

event1 = {slow_waves.ints.UP(:,1)',slow_waves.timestamps',spindles.peaktimes,V1_ripples.peaktimes,slow_waves.ints.UP(:,1)',spindles.peaktimes,spindles.peaktimes_with_UP',slow_waves.ints.UP(:,1)'};
event2 = {ripples.peaktimes',ripples.peaktimes',ripples.peaktimes',ripples.peaktimes',V1_ripples.peaktimes',V1_ripples.peaktimes',ripples.peaktimes',spindles.peaktimes};
event1_name = {'UP state','DOWN state (Delta peak)','Spindle','V1 Ripple','UP state','Spindle','Spindle + UP','Up State'};
event2_name = {'Ripple','Ripple','Ripple','Ripple','V1 Ripple','V1 Ripple','Ripple','Spindle'};

figure
for n = 1:length(event1)
    nexttile
    plot_perievent_event_histogram(event1{n},event2{n},'twin',[-1.2 1.2])
title(sprintf('%s relative to %s',event1_name{n},event2_name{n}))
end
sgtitle('Awake')


%% Draw boundary based on gamma coherence analysis and Down-UP LFP/CSD
% Currently not used

%% Spatial representation
load extracted_position
load extracted_lap


place_fields_even = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'even laps')
place_fields_odd = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,'odd laps')
save extracted_place_fields_V1 place_fields

x_bins_width = 10;
cluster_name = 'V1';
clusters = V1_clusters;
place_fields_BAYESIAN = calculate_place_fields_masa_NPX(x_bins_width,position,clusters,[]);
V1_track_1_cells =  find(place_fields_BAYESIAN.track(1).mean_rate_track > 1 & place_fields_BAYESIAN.track(2).mean_rate_track < 0.5);
[C indexA indexB] = intersect(clusters.id_conversion(V1_track_1_cells,2),find(nominal_KSLabel=='good'))

figure
for cell = 1:length(indexB)
    imagesc(flip(initMap{indexB(cell)}(:,:,6,1)))
    colorbar
    nexttile
end


figure
subplot(2,2,1)

bar(place_fields_BAYESIAN.track(1).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate')
ylim([0 20])
title('Mean Firing Rate on Track 1 (when moving)')

subplot(2,2,2)
bar(place_fields_BAYESIAN.track(2).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate')
ylim([0 20])
title('Mean Firing Rate on Track 2 (when moving)')

subplot(2,2,3)
bar(place_fields_BAYESIAN.track(1).mean_rate_track(V1_track_1_cells) - place_fields_BAYESIAN.track(2).mean_rate_track(V1_track_1_cells))
ylabel('Mean Firing Rate Difference')
ylim([0 20])
title('Firing rate difference (when moving)')
sgtitle('Track 1 ''preferred'' V1 cells')

lap_place_fields_map = [];
for track = 1 : length(lap_times)
    potential_cells = find(place_fields_BAYESIAN.track(track).skaggs_info > 0.3);  
    for i = 1 : length(lap_times(track).lap)

        disp([num2str(i) ' out of ' num2str(length(lap_times(track).lap))])

        % Extract place field of each complete lap
        lap_start_time = lap_times(track).start(i);
        lap_end_time = lap_times(track).end(i);

        place_fields = get_lap_place_fields_masa(x_bins_width,place_fields_BAYESIAN,clusters,track,lap_start_time,lap_end_time);
        lap_place_fields(track).lap{i} = place_fields;
        clear place_fields
        for cell = 1:length(potential_cells)
        lap_place_fields_map{track}{cell}(i,:) = lap_place_fields(track).lap{i}.raw{potential_cells(cell)};
        end
    end

    figure
    for cell = 1:length(potential_cells)
        nexttile
        imagesc(lap_place_fields_map{track}{cell})
        colorbar
        xlabel('Position bin')
        ylabel('Lap')
        title(sprintf('unit %i',potential_cells(cell)))
    end
    sgtitle(sprintf('%s Track %i representation',cluster_name,track))
end

 %% Biasing track selective neruons in V1 by visual stimuli
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

%%
events = ripples;
events = reactivations;
CA1_track1_SWR_spikes = [];
V1_track1_SWR_spikes = [];
V1_track2_SWR_spikes = [];
active_cells_SWR = [];
UP_T1_ripple_number = [];
UP_ripple_time = [];
DOWN_ripple_time = [];
count = 1;
events = slow_waves.ints.DOWN;
 for event = 1:length(events)
    [ripple_this_UP,ripple_index_this_UP] = RestrictInts(ripples.SWS_peaktimes,[events(event,1) events(event,2)]);
    UP_ripple_number(event) = length(ripple_this_UP);
    UP_ripple_time = [UP_ripple_time; ripple_this_UP -  events(event,1)];
    DOWN_ripple_time = [DOWN_ripple_time; events(event,2) - ripple_this_UP];

    UP_T1_ripple_number(event) = 0;
    if ~isempty(ripple_this_UP)
        ripple_index_this_UP = find(ripple_index_this_UP == 1);
        active_cells_SWR = [];
        for n = 1:length(ripple_index_this_UP)
            onset = ripples.onset(ripple_index_this_UP(n));
            offset = ripples.offset(ripple_index_this_UP(n));
            spike_index = find(V1_clusters.spike_times > onset & V1_clusters.spike_times < offset);
            for cell = 1:length(V1_track_1_cells)
                V1_track1_SWR_spikes(cell,count) = sum(find(V1_clusters.spike_id(spike_index,1) == V1_track_1_cells(cell)));
            end
           active_cells_SWR (n) =  sum(V1_track1_SWR_spikes(:,count) > 0) >=3;
            count = count + 1;
            
        end
        UP_T1_ripple_number(n) = sum(active_cells_SWR);

    end
     
 end

 figure
     subplot(2,1,1)
 histogram(UP_ripple_time,[0:0.02:0.2],'Normalization','probability')
 xlabel('Seconds after each DOWN-UP transition (Ripple probability)')
 ylabel('Probability')
 hold on
    subplot(2,1,2)
 histogram(DOWN_ripple_time,[0:0.02:0.2],'Normalization','probability')
 xlabel('Seconds before each DOWN-UP transition (Ripple probability)')
 ylabel('Probability')

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


