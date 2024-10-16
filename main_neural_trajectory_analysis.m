
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))

%% Example 3: Single-trial neural trajectories
% The dataset corresponds to Section 3.2 in the DataHigh JNE paper.  
% The data includes two conditions with 15 trials each.  The space 
% was reduced from 61 neurons to 15 latent dimensions with 
% Gaussian-process factor analysis (GPFA) (Yu et al., 2009).
%
%  Detailed instructions can be found in the User Guide and the website.

cd ..
load('./data/ex2_rawspiketrains.mat');
% D(itrial).data : (num_neurons x num_1ms_bins)
DataHigh(D, 'DimReduce');
cd ./examples

cd ..
load('./data/ex2_singletrialtrajs.mat');
DataHigh(D);
% D(itrial).data : (num_latents x num_20ms_timebins)
cd ./examples


%% Load data for extracting ripple neural trajectory and then DLAG analysis
clear all
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\buzcode\externalPackages'))
% addpath(genpath('X:\ibn-vision\USERS\Masa\code\spikes'))
addpath(genpath('P:\corticohippocampal_replay\code\buzcode\externalPackages'))
addpath(genpath('P:\corticohippocampal_replay\code\spikes'))
addpath(genpath('Z:\ibn-vision\USERS\Masa\code\DLAG'))


SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';

Stimulus_types_all = {'RUN','POST'};
% Stimulus_type = 'RUN'; % extract LFP during RUN
% Stimulus_type = 'POST'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
c = 1;

for epoch = 1:length(Stimulus_types_all)
    Stimulus_type= Stimulus_types_all{epoch};

    for nsession =1:length(experiment_info)
        tic
        session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        if isempty(session_info)
            continue
        end
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

        for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
            cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
            load best_channels
            load extracted_PSD
            load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_laps')
            column = 1;

%             load(sprintf('extracted_CA1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_HPC_clusters_combined%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_V1_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            %             load(sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks')));
            load(sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load(sprintf('decoded_ripple_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')))
            load('extracted_V1_place_fields.mat')
            load('extracted_HPC_place_fields_combined.mat')

            tic
            % Reactivation log odds decoding
            for mprobe = 1:length(session_info(n).probe)
                % probe for ripples
                ripple_probe_no = session_info(n).probe(mprobe).probe_id + 1;
                ripple_probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                replay = ripples.probe(ripple_probe_no);
                %                 replay.offset = replay.onset(event) + 0.1;
                spike_count = [];

                for nprobe = 1:length(session_info(n).probe)
                    options = session_info(n).probe(nprobe);
                    options.importMode = 'KS';
                    options.gFileNum = gFileNum(n);
                    probe_no = session_info(n).probe(nprobe).probe_id + 1;
                    options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
                    options.ROOTPATH = ROOTPATH;
                    probe_hemisphere = options.probe_hemisphere;
                
                    %                     [place_fields_BAYESIAN,decoded_events,probability_ratio_original] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay,position,[],stimulus_name{n},timebin);
                    
                    % V1 cells
                    V1_spatial_cell_id = V1_clusters.probe(probe_no).id_conversion(V1_place_fields.probe(probe_no).good_place_cells_LIBERAL,2);
                    
                    for ncell = 1:length(V1_place_fields.probe(probe_no).good_place_cells_LIBERAL)
                        this_cell_spike_times = V1_clusters.probe(probe_no).spike_times(V1_clusters.probe(probe_no).spike_id== V1_spatial_cell_id(ncell));
                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(this_cell_spike_times,replay.onset, [-0.5 0.5], 0.02);
                        V1_spike_count(ncell,:,:) = binnedArray;
                    end
                    
                    % HPC cells
                    HPC_spatial_cell_id = HPC_clusters_combined.id_conversion(HPC_place_fields_combined.good_place_cells_LIBERAL,2);


                    for event = 1:length(replay.onset)

                        
                        replay1.offset = replay.onset(event) + 0.6;
                        replay1.onset = replay.onset(event);

                        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_data{ripple_probe_hemisphere}{probe_hemisphere}(spike_index,:),...
                            [zeros(1,event-1) event_index(event)], [-0.5 0.5], 0.001);

                        [place_fields_BAYESIAN,decoded_events,~,decoded_events_shuffled] = log_odds_reactivation_analysis(clusters,place_fields_BAYESIAN,replay1,position,[],stimulus_name{n},timebin);

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).replay = decoded_events(1).replay_events.replay;
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).replay = decoded_events(2).replay_events.replay;

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability = sum(decoded_events(1).replay_events.replay);
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability = sum(decoded_events(2).replay_events.replay);

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).spikes = decoded_events(1).replay_events.spikes;
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).spikes = decoded_events(2).replay_events.spikes;

                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(1).replay_events(event).timebins_edges = decoded_events(1).replay_events.timebins_edges;
                        decoded_ripple_events_V1(ripple_probe_no).probe(probe_no).track(2).replay_events(event).timebins_edges = decoded_events(2).replay_events.timebins_edges;

                        for nshuffle = 1:1000
                            decoded_ripple_events_V1_shuffled(ripple_probe_no).probe(probe_no).track(1).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(1).replay);
                            decoded_ripple_events_V1_shuffled(ripple_probe_no).probe(probe_no).track(2).replay_events(event).summed_probability(nshuffle,:) = sum(decoded_events_shuffled{nshuffle}(2).replay);
                        end

                    end

                end
            end
            toc


            save(sprintf('decoded_ripple_events_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_V1');
            %             save(sprintf('probability_ratio_global_remapped_V1%s.mat',erase(stimulus_name{n},'Masa2tracks')),'probability_ratio_global_remapped_V1');
            save(sprintf('decoded_ripple_events_V1_shuffled%s.mat',erase(stimulus_name{n},'Masa2tracks')),'decoded_ripple_events_V1_shuffled');
            clear decoded_ripple_events_V1_shuffled decoded_ripple_events_V1
        end
    end
end


%% Load data for extracting peri-ripple neural trajectory and then DLAP analysis


%%


%% Load data for extracting single lap neural trajectory and then DLAG analysis

clear all
SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision';
nsession = 8;

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels;
        load extracted_PSD;
        load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')));

        column = 1;
        
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            options = session_info(n).probe(nprobe);
            options.importMode = 'KS';
            options.gFileNum = gFileNum(n);
            probe_no = session_info(n).probe(nprobe).probe_id + 1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            options.probe_no = probe_no;
            options.ROOTPATH = ROOTPATH;
            % Load all spike data sorted according to the channel position
            



        end

    end

end





