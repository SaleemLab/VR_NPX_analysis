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
spindle_params.passband      = [10 16];      % Narrower band to reduce delta harmonic leakage
spindle_params.thresholds    = [1.5 3];      % [onset/offset, peak] in SD
spindle_params.durations     = [400 3000];   % [min, max] in ms
spindle_params.min_rel_power = 0.15;         % Minimum relative sigma power (YASA-inspired)
spindle_params.min_cycles    = 3;            % Minimum oscillatory cycles
spindle_params.delta_reject  = true;         % Reject events with high delta co-occurrence
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
    else
        load(fullfile(options.ANALYSIS_DATAPATH, 'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH, 'extracted_LFP.mat'), 'LFP');
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
    clear spindles
    
    for nprobe = 1:length(session_info(n).probe)
        probe_no = session_info(n).probe(nprobe).probe_id + 1;
        
        fprintf('  Probe %d (hemisphere %d)...\n', nprobe, probe_no);
        
        tvec = LFP(nprobe).tvec;
        fs = mean(1 ./ diff(tvec));
        
        %% Select best V1 channel for spindle detection
        % Priority: best_SO_V1_trough_peak_ratio > best_V1_high_freq_power
        if isfield(LFP(nprobe), 'best_SO_V1_trough_peak_ratio') && ~isempty(LFP(nprobe).best_SO_V1_trough_peak_ratio)
            [~, V1_ch] = max(LFP(nprobe).best_SO_V1_trough_peak_ratio);
        elseif isfield(LFP(nprobe), 'best_V1_high_freq_power') && ~isempty(LFP(nprobe).best_V1_high_freq_power)
            [~, V1_ch] = max(LFP(nprobe).best_V1_high_freq_power(:, 7));
        else
            V1_ch = 1;
            warning('No V1 channel selection criteria found. Using channel 1.');
        end
        
        %% Get LFP data from best V1 channel
        if isfield(LFP(nprobe), 'best_V1_high_freq') && ~isempty(LFP(nprobe).best_V1_high_freq)
            lfp_data = LFP(nprobe).best_V1_high_freq(V1_ch, :)';
        elseif isfield(LFP(nprobe), 'best_SO_V1') && ~isempty(LFP(nprobe).best_SO_V1)
            lfp_data = LFP(nprobe).best_SO_V1(V1_ch, :)';
        else
            warning('  No V1 LFP data found for probe %d. Skipping.', nprobe);
            spindles(probe_no) = DetectSpindles_masa(zeros(10,1), (1:10)'/fs, ...
                'frequency', fs, 'best_channel', V1_ch);
            continue
        end
        
        %% Get actual channel number for metadata
        if isfield(LFP(nprobe), 'best_V1_high_freq_channel')
            actual_channel = LFP(nprobe).best_V1_high_freq_channel(V1_ch);
        else
            actual_channel = V1_ch;
        end
        
        %% Run improved spindle detection
        [spindles(probe_no)] = DetectSpindles_masa(lfp_data, tvec', ...
            'frequency',       fs, ...
            'passband',        spindle_params.passband, ...
            'thresholds',      spindle_params.thresholds, ...
            'durations',       spindle_params.durations, ...
            'min_rel_power',   spindle_params.min_rel_power, ...
            'min_cycles',      spindle_params.min_cycles, ...
            'delta_reject',    spindle_params.delta_reject, ...
            'delta_thresh',    spindle_params.delta_thresh, ...
            'behaviour',       Behaviour, ...
            'noise',           [], ...
            'show',            spindle_params.show, ...
            'best_channel',    actual_channel);
        
        fprintf('  Probe %d: %d spindles detected.\n', nprobe, length(spindles(probe_no).onset));
        
        %% Add speed info per spindle
        if ~isempty(spindles(probe_no).onset)
            if isfield(Behaviour, 'speed')
                for event = 1:length(spindles(probe_no).onset)
                    idx = Behaviour.sglxTime >= spindles(probe_no).onset(event) & ...
                          Behaviour.sglxTime <= spindles(probe_no).offset(event);
                    spindles(probe_no).speed(event) = mean(Behaviour.speed(idx));
                end
            elseif isfield(Behaviour, 'mobility_zscore')
                for event = 1:length(spindles(probe_no).onset)
                    idx = Behaviour.sglxTime >= spindles(probe_no).onset(event) & ...
                          Behaviour.sglxTime <= spindles(probe_no).offset(event);
                    spindles(probe_no).mobility_zscore(event) = mean(Behaviour.mobility_zscore(idx));
                end
            end
        end
        
        %% Restrict to SWS if behavioural states available
        if isfield(behavioural_state_merged, 'SWS') && ~isempty(behavioural_state_merged.SWS)
            if ~isempty(spindles(probe_no).onset)
                [spindles(probe_no).SWS_offset, spindles(probe_no).SWS_index] = ...
                    RestrictInts(spindles(probe_no).offset, behavioural_state_merged.SWS);
                spindles(probe_no).SWS_onset = spindles(probe_no).onset(spindles(probe_no).SWS_index);
                spindles(probe_no).SWS_peaktimes = spindles(probe_no).peaktimes(spindles(probe_no).SWS_index);
                fprintf('  Probe %d: %d spindles during SWS.\n', nprobe, length(spindles(probe_no).SWS_onset));
            end
        end
        
        %% Restrict to quiet wake if available
        if isfield(behavioural_state_merged, 'quietWake') && ~isempty(behavioural_state_merged.quietWake)
            if ~isempty(spindles(probe_no).onset)
                [spindles(probe_no).awake_offset, spindles(probe_no).awake_index] = ...
                    RestrictInts(spindles(probe_no).offset, behavioural_state_merged.quietWake);
                spindles(probe_no).awake_onset = spindles(probe_no).onset(spindles(probe_no).awake_index);
                spindles(probe_no).awake_peaktimes = spindles(probe_no).peaktimes(spindles(probe_no).awake_index);
            end
        end
        
    end  % probe loop
    
    %% Save detected spindle events
    if contains(stimulus_name{n}, 'Masa2tracks')
        savefile = fullfile(options.ANALYSIS_DATAPATH, ...
            sprintf('detected_spindle_events%s.mat', erase(stimulus_name{n}, 'Masa2tracks')));
    else
        savefile = fullfile(options.ANALYSIS_DATAPATH, 'detected_spindle_events.mat');
    end
    
    save(savefile, 'spindles');
    fprintf('  Saved: %s\n', savefile);
    
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
