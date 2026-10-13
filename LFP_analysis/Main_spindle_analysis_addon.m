% MAIN_SPINDLE_ANALYSIS_ADDON
% Re-detect spindle events across all sessions using the improved
% DetectSpindles_masa algorithm (with relative sigma power, cycle count
% validation, and delta co-occurrence rejection).
%
% Loops through sessions, loads extracted LFP, loops through probes,
% runs spindle detection, and saves as detected_spindle_events.mat.
%
% Simplified from detect_behavioural_and_brain_states.m
%
% Masahiro Takigawa, 2026

clear all
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))

%% ==================== Configuration ====================
SUBJECTS = {'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS, option);
experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';

% Spindle detection parameters
spindle_params.passband      = [9 17];      % Narrower band to reduce delta harmonic leakage
spindle_params.thresholds    = [1.5 3];      % [onset/offset, peak] in SD
spindle_params.durations     = [400 3000];   % [min, max] in ms
spindle_params.min_rel_power = 0.15;         % Minimum relative sigma power (YASA-inspired)
spindle_params.min_cycles    = 3;            % Minimum oscillatory cycles
spindle_params.delta_reject  = false;         % Reject events with high delta co-occurrence
spindle_params.delta_thresh  = 3;            % Delta zscore threshold
spindle_params.show          = 'on';         % 'on' or 'off'

%% ==================== Main Loop ====================
for nsession = 1:length(experiment_info)

    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName, Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName, Stimulus_type));

    if isempty(stimulus_name)
        continue
    end

    %% Handle multiple recordings of same stimulus
    if length(stimulus_name) > 1
        if contains(Stimulus_type, 'PRE')
            n = find(contains(stimulus_name, '_2'));
        else
            session_info = session_info(~contains(stimulus_name, 'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name, 'PRE'));
            if length(stimulus_name) > 1
                n = find(contains(stimulus_name, '_2'));
            else
                n = 1;
            end
        end
    else
        n = 1;
    end

    options = session_info(n).probe(1);

    % Check that extracted data exists
    DIR = dir(fullfile(options.ANALYSIS_DATAPATH, 'extracted_clusters*.mat'));
    if isempty(DIR)
        disp(['Skipping session ' num2str(nsession) ': no extracted clusters.']);
        continue
    end

    fprintf('\n========================================\n');
    fprintf('Session %d/%d: %s %s\n', nsession, length(experiment_info), options.SUBJECT, options.SESSION);
    fprintf('========================================\n');

    %% Load LFP and Behaviour
    if contains(stimulus_name{n}, 'Masa2tracks')
        load(fullfile(options.ANALYSIS_DATAPATH, sprintf('extracted_behaviour%s.mat', erase(stimulus_name{n}, 'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH, sprintf('extracted_LFP%s.mat', erase(stimulus_name{n}, 'Masa2tracks'))), 'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters%s.mat',erase(stimulus_name,'Masa2tracks'))),'session_clusters');
    else
        load(fullfile(options.ANALYSIS_DATAPATH, 'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH, 'extracted_LFP.mat'), 'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters'); % Session clusters for SUA
    end

    %% Load best channels
    best_channels_file = fullfile(options.ANALYSIS_DATAPATH, '..', 'best_channels.mat');
    if exist(best_channels_file, 'file')
        load(best_channels_file);
    end

    %% Load behavioural state (for restricting to SWS if available)
    beh_state_file = fullfile(options.ANALYSIS_DATAPATH, 'behavioural_state_merged.mat');
    if exist(beh_state_file, 'file')
        load(beh_state_file, 'behavioural_state_merged');
    else
        behavioural_state_merged = struct();
    end

    %% Loop through probes
    clear spindles_new spindles_default
    if contains(stimulus_name,'Sleep')
        speed = double(session_clusters.mobility_thresholded{1});
        Behaviour.mobility_zscore=session_clusters.mobility_zscore{1};
    else
        speed = Behaviour.speed;
        speed(isnan(speed))=0;
        w = gausswin(9);
        w = w / sum(w);
        speed = filtfilt(w,1,speed')';
    end

    %     speed = interp1(Behaviour.tvec,speed,tvec,'linear');

    for nprobe = 1:length(session_info(n).probe)
        probe_no = session_info(n).probe(nprobe).probe_id + 1;

        fprintf('  Probe %d (hemisphere %d)...\n', nprobe, probe_no);

        tvec = LFP(nprobe).tvec;
        fs = mean(1 ./ diff(tvec));

        % Select best V1 channel for spindle detection
        % Priority: best_SO_V1_trough_peak_ratio > best_V1_high_freq_power
        if isfield(LFP(nprobe), 'best_SO_V1_trough_peak_ratio') && ~isempty(LFP(nprobe).best_SO_V1_trough_peak_ratio)
            [~, V1_ch] = max(LFP(nprobe).best_SO_V1_trough_peak_ratio);
        elseif isfield(LFP(nprobe), 'best_V1_high_freq_power') && ~isempty(LFP(nprobe).best_V1_high_freq_power)
            [~, V1_ch] = max(LFP(nprobe).best_V1_high_freq_power(:, 7));
        else
            V1_ch = 1;
            warning('No V1 channel selection criteria found. Using channel 1.');
        end

        %Get LFP data from best V1 channel
        if isfield(LFP(nprobe), 'best_V1_high_freq') && ~isempty(LFP(nprobe).best_V1_high_freq)
            lfp_data = LFP(nprobe).best_V1_high_freq(V1_ch, :)';
        elseif isfield(LFP(nprobe), 'best_SO_V1') && ~isempty(LFP(nprobe).best_SO_V1)
            lfp_data = LFP(nprobe).best_SO_V1(V1_ch, :)';
        else
            warning('  No V1 LFP data found for probe %d. Skipping.', nprobe);
            spindles_new(probe_no) = DetectSpindles_masa(zeros(10,1), (1:10)'/fs, ...
                'frequency', fs, 'best_channel', V1_ch);
            spindles_default(probe_no) = spindles_new(probe_no);
            continue
        end

        %Get actual channel number for metadata
        if isfield(LFP(nprobe), 'best_V1_high_freq_channel')
            actual_channel = LFP(nprobe).best_V1_high_freq_channel(V1_ch);
        else
            actual_channel = V1_ch;
        end

        % Loop over detection methods
        methods = {'new', 'default'};
        for m = 1:length(methods)
            current_method = methods{m};

            if strcmp(current_method, 'new')
                current_thresholds = [1.5 3];
                current_delta_reject = spindle_params.delta_reject;
            else
                current_thresholds = [2.5 5];
                current_delta_reject = false;
            end

            %Run improved spindle detection
            [temp_spindles] = DetectSpindles_masa(lfp_data, tvec', ...
                'method',          current_method, ...
                'frequency',       fs, ...
                'passband',        spindle_params.passband, ...
                'thresholds',      current_thresholds, ...
                'durations',       spindle_params.durations, ...
                'min_rel_power',   spindle_params.min_rel_power, ...
                'min_cycles',      spindle_params.min_cycles, ...
                'delta_reject',    current_delta_reject, ...
                'delta_thresh',    spindle_params.delta_thresh, ...
                'behaviour',       Behaviour, ...
                'noise',           [], ...
                'show',            spindle_params.show, ...
                'best_channel',    actual_channel);

            fprintf('  Probe %d (%s): %d spindles detected.\n', nprobe, current_method, length(temp_spindles.onset));

            %Add speed info per spindle
            if ~isempty(temp_spindles.onset)
                if isfield(Behaviour, 'speed')
                    for event = 1:length(temp_spindles.onset)
                        idx = Behaviour.sglxTime >= temp_spindles.onset(event) & ...
                            Behaviour.sglxTime <= temp_spindles.offset(event);
                        temp_spindles.speed(event) = mean(Behaviour.speed(idx));
                    end
                elseif isfield(Behaviour, 'mobility_zscore')
                    for event = 1:length(temp_spindles.onset)
                        idx = Behaviour.sglxTime >= temp_spindles.onset(event) & ...
                            Behaviour.sglxTime <= temp_spindles.offset(event);
                        temp_spindles.mobility_zscore(event) = mean(Behaviour.mobility_zscore(idx));
                    end
                end
            end

            % Restrict to SWS if behavioural states available
            if isfield(behavioural_state_merged, 'SWS') && ~isempty(behavioural_state_merged.SWS)
                if ~isempty(temp_spindles.onset)
                    [temp_spindles.SWS_offset, temp_spindles.SWS_index] = ...
                        RestrictInts(temp_spindles.offset, behavioural_state_merged.SWS);
                    temp_spindles.SWS_onset = temp_spindles.onset(temp_spindles.SWS_index);
                    temp_spindles.SWS_peaktimes = temp_spindles.peaktimes(temp_spindles.SWS_index);
                    fprintf('  Probe %d (%s): %d spindles during SWS.\n', nprobe, current_method, length(temp_spindles.SWS_onset));
                end
            end

            % Restrict to quiet wake if available
            if isfield(behavioural_state_merged, 'quietWake') && ~isempty(behavioural_state_merged.quietWake)
                if ~isempty(temp_spindles.onset)
                    [temp_spindles.awake_offset, temp_spindles.awake_index] = ...
                        RestrictInts(temp_spindles.offset, behavioural_state_merged.quietWake);
                    temp_spindles.awake_onset = temp_spindles.onset(temp_spindles.awake_index);
                    temp_spindles.awake_peaktimes = temp_spindles.peaktimes(temp_spindles.awake_index);
                end
            end

            if strcmp(current_method, 'new')
                spindles_new(probe_no) = temp_spindles;
            else
                spindles_default(probe_no) = temp_spindles;
            end

        end % end methods loop

        if  contains(stimulus_name{n},'RUN1')|contains(stimulus_name{n},'RUN2')
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(stimulus_name{n},'Masa2tracks_')),sprintf('Probe%i',probe_no)))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',erase(stimulus_name{n},'Masa2tracks_')),sprintf('Probe%i',probe_no)),[])
        else
            mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',stimulus_name{n}),sprintf('Probe%i',probe_no)))
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures',sprintf('LFP_events_%s',stimulus_name{n}),sprintf('Probe%i',probe_no)),[])
        end
    end  % probe loop

    %     plot(tvec',zscore(lfp_data));hold on;
    %     scatter(spindles(2).SWS_onset,1.1*ones(1,length(spindles(2).SWS_onset)),'k');hold on;scatter(spindles_new(2).SWS_onset,ones(1,length(spindles_new(2).SWS_onset)),'r')
    %     ylim([0.9 1.2])
    %% Save detected spindle events
    if contains(stimulus_name{n}, 'Masa2tracks')
        savefile_new = fullfile(options.ANALYSIS_DATAPATH, ...
            sprintf('detected_spindle_events_new%s.mat', erase(stimulus_name{n}, 'Masa2tracks')));
        savefile_default = fullfile(options.ANALYSIS_DATAPATH, ...
            sprintf('detected_spindle_events_default%s.mat', erase(stimulus_name{n}, 'Masa2tracks')));
    else
        savefile_new = fullfile(options.ANALYSIS_DATAPATH, 'detected_spindle_events_new.mat');
        savefile_default = fullfile(options.ANALYSIS_DATAPATH, 'detected_spindle_events_default.mat');
    end

    spindles = spindles_new;
    save(savefile_new, 'spindles');
    fprintf('  Saved: %s\n', savefile_new);

    spindles = spindles_default;
    save(savefile_default, 'spindles');
    fprintf('  Saved: %s\n', savefile_default);


    %% Save figures
    fig_dir = fullfile(options.ANALYSIS_DATAPATH, '..', 'figures', 'spindle_detection');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    if exist('save_all_figures', 'file')
        save_all_figures(fig_dir, []);
    end
    close all

end  % session loop

fprintf('\n=== Spindle detection addon complete. ===\n');
