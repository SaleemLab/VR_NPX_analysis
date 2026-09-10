%%% For analysis of trial-by-trial pupil size vs "omission responses". ERB 2026
% V1 L4 depth is identified via CSD analysis (for which, run CSD_GAVNIKstims_afterFILT_ellie)

clear all
addpath(genpath('C:\Users\eleanor.benoit\Documents\GitHub\VR_NPX_analysis'))

% 1/5 
shank_to_plot = 3; %for NPX1.0, set as 1.
stim_duration = 150;
subjects_dir = 'V:\Ellie\DATA\SUBJECTS';
SUBJECTS = {'M00069'};
option = 'V1-HPC';
experiment_info = subject_session_stimuli_mapping_Ellie(SUBJECTS, option);

%%% 2/5
Stimulus_type = 'GAVNIK_A_CD'; % 'GAVNIK_A_CD' 'GAVNIK_ABCD_1' 
if contains(Stimulus_type, '_A_CD')
    stimulus_color = [1.0,  0.0,  0.0];  %  red - GAVNIK_A_CD  
elseif contains(Stimulus_type, '_E_CD') 
    stimulus_color = [0.0,  1.0,  0.0];  % green - GAVNIK_E_CD  
elseif contains(Stimulus_type, '_ABCD')    
    stimulus_color = [0.0,  0.0,  1.0];  % blue - GAVNIK_ABCD
end  

% 3/5 files will be saved here in the cd
analysis_dir = fullfile(subjects_dir, SUBJECTS{1}, ...
    'analysis', '20251003', Stimulus_type);
cd(analysis_dir);

%% SET THIS 4/5***** For NPX2.0 you will use a different L4 channel for each shank. Use CSD to estimate the best channel to use in L4
probe_type = 1; % NPX1.0 is type 0, NPX2.0 is type 1.

% 5/5 
for nsession = 8 % row number of recording date in "experiment_info" 
    session_info = experiment_info(nsession).session(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(strcmp(experiment_info(nsession).StimulusName,Stimulus_type));
    
    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        subject_number = session_info(n).probe(1).SUBJECT;
        
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH, '..', 'earliest_V1sink_CSD.mat'))
        
        % Load the channel map here 
        files = dir(fullfile(options.EPHYS_DATAPATH, '*ChanMap*.mat')); %the channel map has the y coordinate of each channel. The dir function lists the contents of a folder or provides information about files and directories matching a specified pattern
        % Construct the full path to the channel map file
        channel_map_filename = fullfile(files(1).folder, files(1).name); % files(1) refers to the first file in the list of matches
        load(channel_map_filename);

        if probe_type == 0
            layerfour_channels = earliest_V1sink_CSD.overall_best_halfmax_channel;
            layerfour_channels = layerfour_channels(:);
            shank_ids = ones(size(layerfour_channels)); % dummy shank ID = 1
        
            assert(isfield(earliest_V1sink_CSD, 'overall_best_halfmax_depth'), ...
                'Missing overall_best_halfmax_depth in earliest_V1sink_CSD');
        
            earliest_V1sink_CSD = struct( ...
                'shank_id', 1, ...
                'best_channel_this_shank', layerfour_channels(1), ...
                'best_depth_this_shank', earliest_V1sink_CSD(1).overall_best_halfmax_depth ...
            );
        elseif probe_type == 1
            % NPX2.0 → multiple shanks
            % earliest_V1sink_CSD is a struct array, extract fields safely
            layerfour_channels = [earliest_V1sink_CSD.best_channel_this_shank]; 
            shank_ids = [earliest_V1sink_CSD.shank_id]; 
            % ensure column vectors
            layerfour_channels = layerfour_channels(:); 
            shank_ids = shank_ids(:);
        end

        
        % Select the L4 channel for the shank being analysed
        this_shank_idx = find(shank_ids == shank_to_plot, 1);      
        L4_channel = layerfour_channels(this_shank_idx);
        fprintf('Using shank %d, L4 channel %d\n', shank_to_plot, L4_channel);
        
        %% Load LFP data

        DIR = dir(fullfile(options.EPHYS_DATAPATH, '*lf.bin'));
        file_to_use = DIR(1).name;
        
        % Read metadata
        imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH, file_to_use));
        
        % LFP sampling rate
        lfp_sample_rate = SampRate(imecMeta);
        
        % Number of electrophysiology channels
        chanTypes = str2double(strsplit(imecMeta.acqApLfSy, ','));
        nEPhysChan = chanTypes(1);
        
        % Downsampling
        BinWidth = 1/1250;
        downSampleRate = round(lfp_sample_rate * BinWidth);
        
        % Only read the L4 channel
        selected_channels = L4_channel;

        if stim_duration == 150
            extractedTimeWindow = [-1 1.35];
        elseif stim_duration == 200
            extractedTimeWindow = [-1 1.55];
        elseif stim_duration == 250
            extractedTimeWindow = [-1 2.5];
        end
        
        nSamp = fix(lfp_sample_rate * range(extractedTimeWindow));
        % A is the first stimulus in every sequence
        ordered_oris = unique(Task_info.stim_orientation, 'stable');               
        ori = 1; % plot from onset of first stimulus in sequence
        stim_onsets = Task_info.stim_onset(Task_info.stim_orientation == ordered_oris(ori));
        
        nTrials = length(stim_onsets);
        
        
        % Work out how many samples remain after downsampling
        nTimepoints = ceil(nSamp / downSampleRate);
        
        % L4 LFP: time × trial
        resps = zeros(1, nTimepoints, nTrials, 'single');
        
        for thisTrial = 1:nTrials
        
            stim_on = stim_onsets(thisTrial);
        
            % Read LFP around stimulus onset
            tresps = ReadBin( ...
                fix(lfp_sample_rate * (stim_on + extractedTimeWindow(1))), ...
                nSamp, ...
                imecMeta, ...
                file_to_use, ...
                options.EPHYS_DATAPATH);
        
            % Keep only electrophysiology channels
            tresps(nEPhysChan+1:end,:) = [];
        
            % Gain correction
            if strcmp(imecMeta.typeThis, 'imec')
                tresps = GainCorrectIM(tresps, 1:nEPhysChan, imecMeta);
            else
                tresps = GainCorrectNI(tresps, 1:nEPhysChan, imecMeta);
            end
        
            % Select L4 channel
            tresps = tresps(L4_channel,:);
        
            % Zero mean
            tresps = tresps - mean(tresps);
        
            % Downsample
            tresps = downsample(tresps, downSampleRate);
        
            % Store
            resps(1,:,thisTrial) = single(tresps);
        
        end
        
        % Sampling rate after downsampling
        sampling_rate = SampRate(imecMeta) / downSampleRate;
        
        % Time vector corresponding to the downsampled LFP
        timebins = extractedTimeWindow(1) + 1/sampling_rate/2 : ...
                   1/sampling_rate : ...
                   extractedTimeWindow(2) - 1/sampling_rate/2;
        % Plot window
        if stim_duration == 250
            tidx = timebins >= -0.3 & timebins <= 2.5;
        else
            tidx = timebins >= -0.3 & timebins <= 1.35;
        end
        
        % Baseline: 1 second before stimulus onset
        baseline_duration = 1;
        baseline_samples = round(baseline_duration * sampling_rate);
        
        % Extract L4 LFP
        lfp_responses = squeeze(resps(1,:,:));   % time x trials
        
        % Baseline-correct each trial separately
        baseline_mean = mean(lfp_responses(1:baseline_samples,:), 1);
        lfp_responses = lfp_responses - baseline_mean;
        % Convert every trial to microvolts
        lfp_responses = lfp_responses * 1e6;

        % Mean across trials
        mean_lfp_trace = mean(lfp_responses, 2, 'omitnan');
        
        % SEM across trials
        sem_lfp_trace = std(lfp_responses, 0, 2, 'omitnan') ./ ...
                        sqrt(sum(isfinite(lfp_responses), 2));
               
        %% Plot
        
        figure;
        % Define grey intervals
        if (contains(Stimulus_type, 'GAVNIK250_A_CD'))
            grey_intervals = [-0.5 0; 1.0 2.0];
            element2_onset = 0.25;
        else
            grey_intervals = [-0.5 0; 0.6 2.0];
            element2_onset = 0.15;
        end    

        for i = 1:size(grey_intervals,1)
            h = xregion(grey_intervals(i,1), grey_intervals(i,2), ...
                FaceColor=[0.7 0.7 0.7], ...
                FaceAlpha=0.3);
            h.HandleVisibility = 'off';
        end     

        hold on;
        
        plot(timebins(tidx), mean_lfp_trace(tidx), ...
            'Color', stimulus_color, ...
            'LineWidth', 2);
        
        h = xline(element2_onset, 'k', 'Element 2 onset');
            h.LabelHorizontalAlignment = 'left';
            xline(element2_onset + 0.03, 'k--', '+30 ms');
            xline(element2_onset + 0.08, 'k--', '+80 ms');
        
        xlabel('Time from stimulus onset (s)');
        ylabel('LFP (\muV)');
        
        title(sprintf('%s Shank %d: Mean L4 LFP response for sense-checking: %s', ...
            Stimulus_type, shank_to_plot, subject_number), ...
            'Interpreter', 'none');
        
        set(gca, 'FontSize', 24);
        box off;
        
        hold off;        
        
     
        %% Extract trial-by-trial LFP "omission" measurements
                    
        % Preallocate
        LFP_at_element2onset = nan(nTrials, 1);
        Min_LFP_element2onset = nan(nTrials, 1);
        Min_LFP_time = nan(nTrials, 1);
        element2_response = nan(nTrials, 1);
        
        if stim_duration == 150
            
            % LFP value closest to 150 ms
            [~, idx_element2] = min(abs(timebins - 0.150));
            
            % Window for finding minimum LFP
            window_idx = timebins >= 0.180 & timebins <= 0.230;
            window_times = timebins(window_idx);
            
            target_time = 0.150;
            min_window = [0.180 0.230];
        
        elseif stim_duration == 250
            
            % LFP value closest to 250 ms
            [~, idx_element2] = min(abs(timebins - 0.250));
            
            % Window for finding minimum LFP
            window_idx = timebins >= 0.280 & timebins <= 0.330;
            window_times = timebins(window_idx);
            
            target_time = 0.250;
            min_window = [0.280 0.330];
        
        end
        
        
        %% Calculate omission response separately for each trial
        
        for trial = 1:nTrials
            
            % Extract this trial's LFP trace
            this_lfp = lfp_responses(:, trial);
            
            % Remove singleton dimensions
            this_lfp = squeeze(this_lfp);
            
            % LFP at element 2 onset
            LFP_at_element2onset(trial) = this_lfp(idx_element2);
            
            % Find minimum LFP in the specified window
            window_values = this_lfp(window_idx);
            
            [Min_LFP_element2onset(trial), idx_min] = min(window_values);
            
            % Time of minimum
            Min_LFP_time(trial) = window_times(idx_min);
            
            % Calculate omission response
            element2_response(trial) = ...
                LFP_at_element2onset(trial) - Min_LFP_element2onset(trial);
            
        end



        % Pupil data
        pupil_area = Behaviour.pupil_area_sixtyHz;
        pupil_time = Behaviour.sglxTime;
        
        disp(size(pupil_area))
        disp(size(pupil_time))
        disp([pupil_time(1), pupil_time(end)])
        disp([min(stim_onsets), max(stim_onsets)])
        

        % Pupil area at onset of A and element2 for each trial
        pupil_at_A_onset = nan(nTrials,1);
        pupil_at_element2_onset = nan(nTrials,1);
        % interpolate pupil_area from the two values at the two pupil_times that straddle the stim_onset time
        % pupil is only sampled at 60Hz (every 16.7ms)
        for trial = 1:nTrials
            % Pupil area at A onset
            pupil_at_A_onset(trial) = interp1( ...
                pupil_time, pupil_area, ...
                stim_onsets(trial), ...
                'linear');
            % Pupil area at element 2 onset
            pupil_at_element2_onset(trial) = interp1( ...
                pupil_time, pupil_area, ...
                stim_onsets(trial) + element2_onset, ...
                'linear');
        end
        
        % Change in pupil area from A onset to element 2 onset
        % Negative = pupil became smaller
        % Positive = pupil became larger
        pupil_change = pupil_at_element2_onset - pupil_at_A_onset;

        disp(size(pupil_at_A_onset))
        disp(size(element2_response))
        
        disp([pupil_at_A_onset(1:10), element2_response(1:10)])
        
        %% ---------------- PLOT 1: Pupil at A onset ----------------
        valid_trials = isfinite(element2_response) & isfinite(pupil_at_A_onset);

        x = element2_response(valid_trials);
        y = pupil_at_A_onset(valid_trials);
        
        fig = figure;
        scatter(x, y, 45, stimulus_color, 'filled');
        
        xlabel('Trough element2 LFP vs LFP at element2 onset');
        ylabel('Pupil area at A onset, arb. units');
        sgtitle(sprintf('%s Shank %d - %s: Pupil area at A onset vs. element2 L4 LFP', subject_number, shank_to_plot, Stimulus_type), 'Interpreter', 'none');
        %xlim([-0.5 5]);
        ylim([0 3000]);

        grid off;
        set(gca, 'FontSize', 24);  % Tick labels font size
        box off;
        
        %% Correlation and line of best fit

        [r, p_corr] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
        
        coeffs = polyfit(x, y, 1);
        x_fit = linspace(min(x), max(x), 100);
        y_fit = polyval(coeffs, x_fit);
        
        hold on;
        plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
        hold off;
        
        fprintf('Pearson correlation: r = %.3f, p = %.4g\n', r, p_corr);
        
        % Print correlation statistics on plot
        n_datapoints = numel(x);
        
        text(0.05, 0.95, ...
            sprintf('Pearson r = %.2f\np = %.3g\nn = %d', r, p_corr, n_datapoints), ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 12);
                        
        box off;
        
        % Calculate means and SEMs
        mean_omission = mean(element2_response, 'omitnan');
                
        mean_pupil = mean(pupil_at_A_onset, 'omitnan');
        
        hold on;
        plot(mean_omission, mean_pupil, ...
            'ko', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);

        % Also save as .fig
        fig_filename = sprintf('%s Shank %d Pupil at A onset vs element2 L4 LFP response %s.fig', subject_number, shank_to_plot, Stimulus_type);
        fig_save_path = fullfile(pwd, fig_filename);
        savefig(fig, fig_save_path);
        hold off;   
        %% ---------------- PLOT 2: Pupil at element 2 onset ----------------

        valid_trials = isfinite(element2_response) & ...
                       isfinite(pupil_at_element2_onset);
        
        x = element2_response(valid_trials);
        y = pupil_at_element2_onset(valid_trials);
        
        fig = figure;
        scatter(x, y, 45, stimulus_color, 'filled');
        
        xlabel('Downward V1 L4 LFP deflection following element 2 onset (µV)');
        ylabel('Pupil area at element 2 onset, arb. units');
        %xlim([-0.5 5]);
        ylim([0 3000]);

        grid off;
        set(gca, 'FontSize', 24);
        box off;
        
        %% Correlation and line of best fit
        
        [r, p_corr] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
        
        coeffs = polyfit(x, y, 1);
        x_fit = linspace(min(x), max(x), 100);
        y_fit = polyval(coeffs, x_fit);
        
        hold on;
        plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
        
        fprintf('Element 2 pupil: Pearson r = %.3f, p = %.4g\n', r, p_corr);
        
        n_datapoints = numel(x);

        text(0.05, 0.95, ...
            sprintf('Pearson r = %.2f\np = %.3g\nn = %d', r, p_corr, n_datapoints), ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 12);     


        %% Save mouse-level correlation result
        
        element2_pupil_vsLFP = struct();
        
        % --- Subject/session information ---
        element2_pupil_vsLFP.subject_number = subject_number;
        element2_pupil_vsLFP.Stimulus_type = Stimulus_type;
        element2_pupil_vsLFP.shank_to_plot = shank_to_plot;
        element2_pupil_vsLFP.L4_channel = L4_channel;
        
        % --- Correlation ---
        element2_pupil_vsLFP.pearson_r = r;
        element2_pupil_vsLFP.pearson_p = p_corr;
        element2_pupil_vsLFP.n_trials = n_datapoints;
        
        % --- Pupil measurement ---
        element2_pupil_vsLFP.pupil_measurement = 'Pupil area at element 2 onset';
        element2_pupil_vsLFP.pupil_time = element2_onset;
        
        % --- LFP omission response ---
        element2_pupil_vsLFP.LFP_measurement = ...
            'LFP at element 2 onset minus minimum LFP 30-80 ms after element 2 onset';
        
        element2_pupil_vsLFP.element2_onset = element2_onset;
        
        if stim_duration == 150
            element2_pupil_vsLFP.LFP_min_window = [0.180 0.230];
        elseif stim_duration == 250
            element2_pupil_vsLFP.LFP_min_window = [0.280 0.330];
        end
        
        element2_pupil_vsLFP.LFP_response_units = 'uV';
        element2_pupil_vsLFP.pupil_units = 'arbitrary units';
        
        % --- Trial-level data used for the correlation ---
        element2_pupil_vsLFP.element2_response = x;
        element2_pupil_vsLFP.pupil_at_element2_onset = y;
        
        % --- Analysis information ---
        element2_pupil_vsLFP.correlation_type = 'Pearson';
        element2_pupil_vsLFP.analysis_date = datetime("now");
        
        %% Save result struct

        result_filename = sprintf('element2_pupil_vsLFP_%s_%s.mat', ...
            subject_number, Stimulus_type);

        save(fullfile(subjects_dir, result_filename), 'element2_pupil_vsLFP');
               
        fprintf('Saved correlation result to:\n%s\n', ...
            fullfile(subjects_dir, result_filename));


        % Mean ± SEM
        mean_omission = mean(x, 'omitnan');
        %sem_omission = std(x, 'omitnan') / sqrt(sum(isfinite(x)));
        
        mean_pupil = mean(y, 'omitnan');
        %sem_pupil = std(y, 'omitnan') / sqrt(sum(isfinite(y)));
        
        plot(mean_omission, mean_pupil, ...
            'ko', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);
        
        title(sprintf('%s Shank %d - %s: Pupil area at element 2 onset vs. element2 LFP reponse', ...
            subject_number, shank_to_plot, Stimulus_type), ...
            'Interpreter', 'none');

        % Also save as .fig
        fig_filename = sprintf('%s Shank %d Pupil at element2 onset vs element2 LFP response %s.fig', subject_number, shank_to_plot, Stimulus_type);
        fig_save_path = fullfile(pwd, fig_filename);
        savefig(fig, fig_save_path);
        hold off; 
        
        
        
        %% ---------------- PLOT 3: Change in pupil area ----------------
        
        valid_trials = isfinite(element2_response) & ...
                       isfinite(pupil_change);
        
        x = element2_response(valid_trials);
        y = pupil_change(valid_trials);
        
        fig = figure;
        scatter(x, y, 45, stimulus_color, 'filled');
        
        xlabel('Trough element2 LFP vs LFP at element2 onset');
        ylabel('\Delta pupil area (element 2 onset - A onset)');
        %xlim([-0.5 5]);
        
        grid off;
        set(gca, 'FontSize', 24);
        box off;
        
        %% Correlation and line of best fit
        
        [r, p_corr] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
        
        coeffs = polyfit(x, y, 1);
        x_fit = linspace(min(x), max(x), 100);
        y_fit = polyval(coeffs, x_fit);
        
        hold on;
        plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
        
        fprintf('Pupil change: Pearson r = %.3f, p = %.4g\n', r, p_corr);
        
        n_datapoints = numel(x);

        text(0.05, 0.95, ...
            sprintf('Pearson r = %.2f\np = %.3g\nn = %d', r, p_corr, n_datapoints), ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 12);                
        % Mean ± SEM
        mean_omission = mean(x, 'omitnan');
        %sem_omission = std(x, 'omitnan') / sqrt(sum(isfinite(x)));
        
        mean_pupil_change = mean(y, 'omitnan');
        %sem_pupil_change = std(y, 'omitnan') / sqrt(sum(isfinite(y)));
        
        plot(mean_omission, mean_pupil_change, ...
            'ko', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);
        
        title(sprintf('%s Shank %d - %s: Pupil change vs. element2 LFP deflection', ...
            subject_number, shank_to_plot, Stimulus_type), ...
            'Interpreter', 'none');
        % Also save as .fig
        fig_filename = sprintf('%s Shank %d Pupil at element2 onset less at Aonset vs element2 LFP response %s.fig', subject_number, shank_to_plot, Stimulus_type);
        fig_save_path = fullfile(pwd, fig_filename);
        savefig(fig, fig_save_path);
        hold off;

        %% Sense check: pupil area across entire session

        fig = figure;
        plot(pupil_time, pupil_area, 'k-', 'LineWidth', 1);
        hold on;
        
        % Overlay A stimulus onsets
        for i = 1:length(stim_onsets)
            xline(stim_onsets(i), 'r-', 'HandleVisibility', 'off');
        end
        
        xlabel('Time (s)');
        ylabel('Pupil area (a.u.)');
        title(sprintf('%s %s: Pupil area across session', subject_number, Stimulus_type), ...
            'Interpreter', 'none');
        
        set(gca, 'FontSize', 18);
        box off;
        fig_filename = sprintf('%s Shank %d Pupil area across session %s.fig', subject_number, shank_to_plot, Stimulus_type);
        fig_save_path = fullfile(pwd, fig_filename);
        savefig(fig, fig_save_path);
        
    end
end

