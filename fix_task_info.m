%% Putative pipeline for processing NPX1 data recorded from hippocampus and V1
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)

% main_NPX_data_preprocessing

% Go to SUA_analysis_masa for kilosort + cell explorer cell classification
% Go to CellExplorerTest_masa_adapted for cell exploerer pipeline


%% Set the data folders and processing parameters
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))

addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
% addpath('Z:\ibn-vision\USERS\Masa\code\Masa_utility')
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\NPXAnalysis\NPXAnalysis2022'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\visual_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'));


%% import and align and store Bonsai and cluster spike data

addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
%%%%%% Option 1 use subject_session_stimuli_mapping for all animals you
%%%%%% want to process.
SUBJECTS = {'M23031','M23032','M23034','M23037','M23038'};
options = 'V1-MEC';
ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive

SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
SUBJECTS = {'M23028','M23087','M23153'};
options = 'bilateral';
ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive

 Stimulus_type = 'Track';
Stimulus_type = 'Masa2tracks';
% Stimulus_type = 'Checkerboard';

experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);

for nsession = 2%11:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    if isempty(stimulus_name) % If stimulus not existing for this session
        continue
    end

    for n = 1:length(session_info) % just in case there might be multiple recording for the same stimulus type (e.g. PRE RUN POST Masa2tracks)

        options = session_info(n).probe(1);

        if contains(Stimulus_type,'Masa2tracks')
            % If Masa2tracks, PRE, RUN and/or POST saved in one folder
            load(fullfile(options.ANALYSIS_DATAPATH,...
                sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},Stimulus_type))),'Behaviour')
            load(fullfile(options.ANALYSIS_DATAPATH,...
                sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},Stimulus_type))),'Task_info')
            %             load(fullfile(options.ANALYSIS_DATAPATH,...
            %                 sprintf('extracted_peripherals%s.mat',erase(stimulus_name{n},Stimulus_type))),'Peripherals')
        else
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'),'Task_info')
            %             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_peripherals.mat'),'Peripherals')
        end
        task_info = Task_info;
        behaviour = Behaviour;
        %% Extract laps info
        task_info.complete_laps_id = [];
        task_info.aborted_laps_id = [];

        if ~isempty(task_info.start_time_all)

            start_indices= [];

            x = behaviour.position;  % position data
            t= behaviour.sglxTime; % time

            for nlap = 1:length(task_info.track_ID_all)
                start_indices(nlap) = find(t == task_info.start_time_all(nlap));
            end

            for nlap = 1:length(start_indices)

                if nlap < length(start_indices)
                    current_lap_x = x(start_indices(nlap):start_indices(nlap+1));
                    current_lap_t = t(start_indices(nlap):start_indices(nlap+1));
                else
                    current_lap_x = x(start_indices(nlap):end);
                    current_lap_t = t(start_indices(nlap):end);
                end


                if length(current_lap_x) > 1 && length(current_lap_x(~isnan(current_lap_x))) >1% Only if getting more than 1 datapoint (maybe noise)
                    on_track_x = current_lap_x(~isnan(current_lap_x));
                    on_track_t = current_lap_t(~isnan(current_lap_x));

                    if sum(on_track_x==0)>0
                        start_position = find(on_track_x == 0);
                        if start_position < length(on_track_x)-10 % if start position happens around the end, ignore this
                            on_track_x = on_track_x(start_position(1):end);
                            on_track_t = on_track_t(start_position(1):end);
                        end
                    end

                    if on_track_x(end) ~= on_track_x(end-1)
                        on_track_x(end) = on_track_x(end-1);
                        on_track_t(end) = on_track_t(end-1);
                    end

                    [last_position last_position_index] = max(on_track_x);
                    %             task_info.end_time_all(nlap) = on_track_t(last_position_index);

                    if last_position_index*mean(diff(on_track_t)) < 0.1 % if 140 is reached in less than 0.1 second   (usually first point is still 140 )
                        on_track_x(1:last_position_index) = [];
                        on_track_t(1:last_position_index) = [];
                        [last_position last_position_index] = max(on_track_x);
                    end

                    if isfield(task_info,'reward_delivery_time')
                        if ~isempty(ismember(task_info.reward_delivery_time,on_track_t'))
                            task_info.rewarded_lap_id(ismember(task_info.reward_delivery_time,on_track_t')') = nlap;
                        end
                    end

                    if last_position >= 139 % sometimes last lap ends before 140cm
                        task_info.end_time_all(nlap) = on_track_t(last_position_index); % End time in terms of reaching end of track
                        task_info.complete_laps_id = [task_info.complete_laps_id nlap];% Lap id here is all laps (not track-specific id)
                    else
                        % If never reached the end, then end time is the end of the
                        % entire lap (jumped to next lap after reaching time threshold or lick threshold)
                        if last_position_index*mean(diff(on_track_t)) < 0.1 % if 140 is reached in less than 0.1 second   (usually first point is still 140 )
                            on_track_x(1:last_position_index) = [];
                            on_track_t(1:last_position_index) = [];
                            [last_position last_position_index] = max(on_track_x);
                        end

                        task_info.end_time_all(nlap) = on_track_t(end);
                        task_info.aborted_laps_id = [task_info.aborted_laps_id nlap]; % Lap id here is all laps (not track-specific id)
                    end
                    %                 lap_count = lap_count + 1;
                end
            end
        end

        Task_info = task_info;
        if contains(Stimulus_type,'Masa2tracks')
            % If Masa2tracks, PRE, RUN and/or POST saved in one folder
            %             save(fullfile(options.ANALYSIS_DATAPATH,...
            %                 sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},Stimulus_type))),'Behaviour')
            save(fullfile(options.ANALYSIS_DATAPATH,...
                sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},Stimulus_type))),'Task_info')
            %             save(fullfile(options.ANALYSIS_DATAPATH,...
            %                 sprintf('extracted_peripherals%s.mat',erase(stimulus_name{n},Stimulus_type))),'Peripherals')
        else
            %             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'),'Behaviour')
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'),'Task_info')
            %             save(fullfile(options.ANALYSIS_DATAPATH,'extracted_peripherals.mat'),'Peripherals')
        end
    end
end


