import re

file_path = r'c:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis\LFP_analysis\Main_spindle_analysis_addon.m'
with open(file_path, 'r') as f:
    code = f.read()

# 1. Update clear spindles
code = code.replace("clear spindles", "clear spindles_new spindles_default")

# 2. Replace the probe loop logic
old_logic = """        %% Run improved spindle detection
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
        
        fprintf('  Probe %d: %d spindles detected.\\n', nprobe, length(spindles(probe_no).onset));
        
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
                fprintf('  Probe %d: %d spindles during SWS.\\n', nprobe, length(spindles(probe_no).SWS_onset));
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
        
    end  % probe loop"""

new_logic = """        %% Loop over detection methods
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
            
            %% Run improved spindle detection
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
            
            fprintf('  Probe %d (%s): %d spindles detected.\\n', nprobe, current_method, length(temp_spindles.onset));
            
            %% Add speed info per spindle
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
            
            %% Restrict to SWS if behavioural states available
            if isfield(behavioural_state_merged, 'SWS') && ~isempty(behavioural_state_merged.SWS)
                if ~isempty(temp_spindles.onset)
                    [temp_spindles.SWS_offset, temp_spindles.SWS_index] = ...
                        RestrictInts(temp_spindles.offset, behavioural_state_merged.SWS);
                    temp_spindles.SWS_onset = temp_spindles.onset(temp_spindles.SWS_index);
                    temp_spindles.SWS_peaktimes = temp_spindles.peaktimes(temp_spindles.SWS_index);
                    fprintf('  Probe %d (%s): %d spindles during SWS.\\n', nprobe, current_method, length(temp_spindles.SWS_onset));
                end
            end
            
            %% Restrict to quiet wake if available
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
        
    end  % probe loop"""

code = code.replace(old_logic, new_logic)

# 3. Save logic update
save_old = """    %% Save detected spindle events
    if contains(stimulus_name{n}, 'Masa2tracks')
        savefile = fullfile(options.ANALYSIS_DATAPATH, ...
            sprintf('detected_spindle_events%s.mat', erase(stimulus_name{n}, 'Masa2tracks')));
    else
        savefile = fullfile(options.ANALYSIS_DATAPATH, 'detected_spindle_events.mat');
    end
    
    save(savefile, 'spindles');
    fprintf('  Saved: %s\\n', savefile);"""

save_new = """    %% Save detected spindle events
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
    fprintf('  Saved: %s\\n', savefile_new);
    
    spindles = spindles_default;
    save(savefile_default, 'spindles');
    fprintf('  Saved: %s\\n', savefile_default);"""

code = code.replace(save_old, save_new)

# Handle the early-exit default empty struct
code = code.replace("spindles(probe_no) = DetectSpindles_masa(zeros(10,1), (1:10)'/fs, ...\n                'frequency', fs, 'best_channel', V1_ch);", 
                    "spindles_new(probe_no) = DetectSpindles_masa(zeros(10,1), (1:10)'/fs, 'frequency', fs, 'best_channel', V1_ch);\n            spindles_default(probe_no) = spindles_new(probe_no);")

with open(file_path, 'w') as f:
    f.write(code)

print("Main script updated.")
