%%

addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
clear all

SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'Sleep';

for nsession = 1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    if isempty(stimulus_name)
        continue
    end
    
    if length(stimulus_name)>1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2'));
        else
            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));
            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2'));
            else
                n = 1;
            end
        end
    else
        n = 1;
    end
    
    options = session_info(n).probe(1);
    
    % Initialize periUP_DOWN_LFP_info_V1
    periUP_DOWN_LFP_info_V1 = struct();
    for nprobe = 1:2
        periUP_DOWN_LFP_info_V1(nprobe).UP_SO_amplitude{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_SO_amplitude{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude_smoothed{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude_smoothed{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_SO_phase{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_SO_phase{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_phase{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_phase{2} = [];

        periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_amplitude{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_amplitude{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude_smoothed{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude_smoothed{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_phase{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_phase{2} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_phase{1} = [];
        periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_phase{2} = [];

        periUP_DOWN_LFP_info_V1(nprobe).tvec = -1:0.02:1;
    end
    
    %% Load LFP and Slow Waves
    try
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat')); % loads slow_waves
    catch
        disp(['Missing files for session ' options.SESSION])
        continue
    end
    
    if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
        time_idx = 1:length(LFP(1).best_HPC);
    elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)
        time_idx = 1:length(LFP(2).best_HPC);
    else
        time_idx = 1:length(LFP(1).best_HPC);
    end
    
    probe_no = session_info(n).probe(1).probe_id+1;
    lfp = struct();
    lfp.samplingRate = single(round(1/mean(diff(LFP(probe_no).tvec))));
    lfp.timestamps = single(LFP(probe_no).tvec(time_idx));
    

    lfp_V1 = lfp;
    lfp_V1.data = [];
            
    for iprobe = 1:length(LFP)
        if isfield(LFP(iprobe), 'best_SO_V1_trough_peak_ratio') && ~isempty(LFP(iprobe).best_SO_V1_trough_peak_ratio)
            [~, V1_ch] = max(LFP(iprobe).best_SO_V1_trough_peak_ratio);
        elseif isfield(LFP(iprobe), 'best_V1_high_freq_power') && ~isempty(LFP(iprobe).best_V1_high_freq_power)
            [~, V1_ch] = max(LFP(iprobe).best_V1_high_freq_power(:, 7));
        else
            V1_ch = 1;
            warning(['No V1 channel selection criteria found for probe ' num2str(iprobe) '. Using channel 1.']);
        end
        
        if isfield(LFP(iprobe), 'best_SO_V1')
            % Take only the single best channel per probe
            lfp_V1.data = [lfp_V1.data, single(LFP(iprobe).best_SO_V1(V1_ch, time_idx)')];
        end
    end
    
    if isempty(lfp_V1.data)
        disp('lfp_V1.data is empty')
        continue
    end

    %% Filter for SO and Spindles
    filter_type = 'bandpass';
    filter_order = 3;
    
    % SO (1 - 4 Hz)
    freq_range = [1 4];
    norm_freq_range = freq_range / (lfp.samplingRate / 2);
    b_SO = fir1(filter_order, norm_freq_range, filter_type);
    
    % Spindles (9 - 17 Hz)
    freq_range = [9 17];
    norm_freq_range = freq_range / (lfp.samplingRate / 2);
    b_spindle = fir1(filter_order, norm_freq_range, filter_type);
    
    signal_SO = [];
    signal_spindle = [];
    for nShank = 1:size(lfp_V1.data, 2)
        signal_SO(:,nShank) = filtfilt(b_SO, 1, lfp_V1.data(:,nShank));
        signal_spindle(:,nShank) = filtfilt(b_spindle, 1, lfp_V1.data(:,nShank));
    end
    
    lfp_V1.SO_amplitude_LFP = zscore(abs(hilbert(signal_SO)));
    lfp_V1.SO_phase_LFP = angle(hilbert(signal_SO));
    
    lfp_V1.spindle_amplitude_LFP = zscore(abs(hilbert(signal_spindle)));
    lfp_V1.spindle_phase_LFP = angle(hilbert(signal_spindle));
    
    % --- Add smoothed spindle amplitude ---
    windowLength = round(lfp.samplingRate * 0.320);
    smooth_win = ones(windowLength, 1) / windowLength;
    sigma_envelope = abs(hilbert(signal_spindle));
    sigma_envelope_smoothed = zeros(size(sigma_envelope));
    for nShank = 1:size(sigma_envelope, 2)
        sigma_envelope_smoothed(:,nShank) = filtfilt(smooth_win, 1, sigma_envelope(:,nShank));
    end
    lfp_V1.spindle_amplitude_smoothed_LFP = zscore(sigma_envelope_smoothed);
    % ---------------------------------------------------
    
    clear signal_SO signal_spindle LFP
    
    %% Extract periUP_DOWN segments
    fs = single(lfp_V1.samplingRate);
    win_full_samp = round([-1:0.02:1] * fs);
    
    for nprobe = 1:2
        if ~isfield(slow_waves(nprobe), 'ints') || isempty(slow_waves(nprobe).ints)
            continue
        end
        
        % --- Extract for UP onset ---
        UP_onset_times = slow_waves(nprobe).ints.UP(:,1);
        event_sample_idx_UP = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), UP_onset_times, 'nearest', 'extrap'));
        
        for nevent = 1:length(event_sample_idx_UP)
            event_windows = event_sample_idx_UP(nevent) + win_full_samp; 
            if sum(event_windows > length(lfp_V1.SO_amplitude_LFP))>0
                event_windows(sum(event_windows <= length(lfp_V1.SO_amplitude_LFP)):end) = nan;
            elseif sum(event_windows < 1)>0
                event_windows(event_windows<1) = nan;
            end
            valid_bins = ~isnan(event_windows);

            for mprobe = 1:2
                if mprobe > size(lfp_V1.data, 2)
                    continue
                end
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_SO_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_SO_amplitude{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_smoothed_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude_smoothed{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_amplitude_smoothed{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_SO_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_SO_phase{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).UP_spindle_phase{mprobe} output];
            end
        end

        % --- Extract for DOWN onset (UP offset) ---
        DOWN_onset_times = slow_waves(nprobe).ints.DOWN(:,1);
        event_sample_idx_DOWN = round(interp1(lfp_V1.timestamps, 1:length(lfp_V1.timestamps), DOWN_onset_times, 'nearest', 'extrap'));
        
        for nevent = 1:length(event_sample_idx_DOWN)
            event_windows = event_sample_idx_DOWN(nevent) + win_full_samp; 
            if sum(event_windows > length(lfp_V1.SO_amplitude_LFP))>0
                event_windows(sum(event_windows <= length(lfp_V1.SO_amplitude_LFP)):end) = nan;
            elseif sum(event_windows < 1)>0
                event_windows(event_windows<1) = nan;
            end
            valid_bins = ~isnan(event_windows);

            for mprobe = 1:2
                if mprobe > size(lfp_V1.data, 2)
                    continue
                end
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_amplitude{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_amplitude_smoothed_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude_smoothed{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_amplitude_smoothed{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.SO_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_SO_phase{mprobe} output];
                
                output = nan(length(event_windows),1);
                output(valid_bins) = lfp_V1.spindle_phase_LFP(event_windows(valid_bins),mprobe);
                periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_phase{mprobe} = [periUP_DOWN_LFP_info_V1(nprobe).DOWN_spindle_phase{mprobe} output];
            end
        end
    end
    
    %% Save
    if contains(Stimulus_type,'Sleep')
        save(fullfile(options.ANALYSIS_DATAPATH,'periUP_DOWN_LFP_info_V1_best.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');
%     elseif contains(Stimulus_type,'PRE')
%         save(fullfile(analysis_folder,'periUP_DOWN_LFP_info_V1_PRE.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');
%     else
%         save(fullfile(options.ANALYSIS_DATAPATH,'periUP_DOWN_LFP_info_V1_RUN_best.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');
    end
%     if contains(Stimulus_type,'Sleep')
%         save(fullfile(options.ANALYSIS_DATAPATH,'periUP_DOWN_LFP_info_V1.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');
% %     elseif contains(Stimulus_type,'PRE')
% %         save(fullfile(analysis_folder,'periUP_DOWN_LFP_info_V1_PRE.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');
%     else
%         save(fullfile(options.ANALYSIS_DATAPATH,'periUP_DOWN_LFP_info_V1_RUN.mat'),'periUP_DOWN_LFP_info_V1','-v7.3');
%     end
    
    disp(['Finished session ' options.SESSION])
end
